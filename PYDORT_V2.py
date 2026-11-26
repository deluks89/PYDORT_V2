import numpy as np
from scipy.optimize import fsolve, curve_fit
from scipy.interpolate import interp1d, interpn
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve
import os
import subprocess
import shutil
import pandas as pd
import h5py
from soltrack import SolTrack
from pathlib import Path
import warnings
import copy
import re
from enum import Enum, IntEnum

# BASIC SETTINGS
################

# Directory of VLIDORT-PC-Library
MAIN_PATH = Path(__file__).resolve().parent.as_posix() + '/'

FORTRAN_PATH = MAIN_PATH + 'Fortran/'
QS_BACKEND_PATH = FORTRAN_PATH + 'vlidort_qs/'
STANDARD_BACKEND_PATH = FORTRAN_PATH + 'vlidort_standard/'
MIE_MODEL_PATH = FORTRAN_PATH + 'mie_model/'


class RadiativeTransferModel(str, Enum):
    """Available Fortran radiative-transfer backends."""
    STANDARD = 'vlidort_standard'
    AIRCRAFT = 'vlidort_ac'
    QS = 'vlidort_qs'


class SphericityCorrection(IntEnum):
    """Public spherical-correction accuracy levels."""
    EXPERIMENTAL = 0
    NONE = 1
    TWO_POINT = 2
    THREE_POINT = 3
    FULL = 4

# Sets the Fortran instance to use. See section "Run multiple Pydort instances in parallel" in the manual.
USE_FORTRAN_INSTANCE = 0

# Choice of solar irradiance literature spectrum. Either 'kurucz' (Chance and Kurucz, 2010) or 'tsis' (Coddington et al., 2023). 
SOLAR_IRRADIANCE_SPECTRUM = 'tsis' 


# If True, console output of Fortran executables will be printed in Python console. Gets messy very quickly though.
REPRINT_SUBPROCESS_OUTPUT = False

# ADVANCED SETTINGS: 
####################

DEFAULT_REF_WL = 400
DEFAULT_WAVELENGTHS = np.array([400])
DEFAULT_SPECTRAL_RESOLUTION = 0.5

PROFILE_HIGH_RESOLUTION = 1.
WAVELENGTH_HIGH_RESOLUTION = 0.002
SMOOTHING_FWHM_SAMPLING = 10

MAX_PF_MOMENTS = 500 # MUST NOT BE LARGER THAN MAXMOMENTS_INPUT IN Fortran/vlidort_qs/vqs_def/MSSeg_pars.f90
PF_MOMENT_CUTOFF = 1e-15 # PHASE FUNCTION MOMENTS WITH SMALLER VALUES WILL BE DISCARDED.

# Lower limits applied when writing RT input files (avoid numerical issues):
# Numerical background for active aerosol profiles; see aerosol_floor_1e8_validation.
LOWER_LIMIT_AEROSOL_EXTINCTION = 1e-8 # [1/km], at each type's profile wavelength
LOWER_LIMIT_GAS_CONCENTRATION  = 1.    # [molec/cm3]
LOWER_LIMIT_SURFACE_ALBEDO     = 1e-10

KB = 1.38065e-23 #Boltzmann constant
O2_MOLAR_FRACTION = 0.2095 #Molar fraction of oxygen in the atmosphere

####################################
# GENERAL UTILITY/HELPER FUNCTIONS #
####################################

PATH_PHYSICS = MAIN_PATH + 'physics_data/'

# Note that vertical resolution can be reduced at any point using the change_vertical_grids() method of the Pydort() class.
DEFAULT_VERTICAL_GRID = np.loadtxt(PATH_PHYSICS + 'default_vertical_grid.dat').astype(int)

def get_latest_tropomi_surface_albedo_file():
    """
    Return the newest available TROPOMI surface-albedo netCDF file in
    `physics_data/surface_albedo/`, based on the version string in the
    filename (e.g. v2.1, v2.2, ...).
    """
    surfalb_dir = Path(PATH_PHYSICS) / 'surface_albedo'
    files = list(surfalb_dir.glob('TROPOMI_Sentinel-5P_0125x0125_surface_DLER_v*.nc'))
    if not files:
        raise FileNotFoundError(
            f'No TROPOMI surface albedo netCDF file found in {surfalb_dir}'
        )

    def _version_key(path):
        m = re.search(r'_v(\d+(?:\.\d+)*)\.nc$', path.name)
        if m is None:
            return ()
        return tuple(int(part) for part in m.group(1).split('.'))

    return max(files, key=_version_key).as_posix()

def get_fortran_instance_suffix():
    return '_' + str(USE_FORTRAN_INSTANCE) if USE_FORTRAN_INSTANCE else ''


def get_qs_backend_path():
    """Return the active VLIDORT-QS backend path, including optional parallel instance suffix."""
    return QS_BACKEND_PATH[:-1] + get_fortran_instance_suffix() + '/'


def get_standard_backend_path():
    """Return the active Standard/VLIDORT-AC backend path, including optional parallel instance suffix."""
    return STANDARD_BACKEND_PATH[:-1] + get_fortran_instance_suffix() + '/'


def get_mie_model_path():
    """Return the active Mie-model path, including optional parallel instance suffix."""
    return MIE_MODEL_PATH[:-1] + get_fortran_instance_suffix() + '/'

def run_subprocess(executable):
    """ Helper function to call fortran executable and print progress to spyder console if desired """
    wd = os.getcwd()
    wd_local = '/'.join(executable.split('/')[:-1])
    try:
        os.chdir(wd_local)
        if not REPRINT_SUBPROCESS_OUTPUT:
            subprocess.run(executable, check=True)
        else:
            proc = subprocess.Popen([executable], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            while True:
                  line = proc.stdout.readline()
                  if not line:
                      break
                  print(line, end='')
            returncode = proc.wait()
            if returncode:
                raise subprocess.CalledProcessError(returncode, executable)
    finally:
        os.chdir(wd)


def gaussian_func(x, mean, fwhm, magnitude=1):
    sig = fwhm / 2 / np.sqrt(2*np.log(2))
    return magnitude*np.exp(-0.5*((x-mean)/sig)**2)


def angstrom_func(wavelength, reference_wavelength, ref_value, angstrom_coeffs):
    """ Implementation of the multi-order angstrom wavelength dependence as introduced in the paper section 'Wavelength-dependent parameters' """
    lnx = np.log(ref_value)
    for i, a in enumerate(angstrom_coeffs):
        lnx += a * np.log(reference_wavelength/wavelength)**(i+1)
    return np.exp(lnx)


def centers(array):
    return (array[1:] + array[:-1])/2


def logspace(mini, maxi, n, rounded = True):
    """ Similar to numpy.linspace but distributes values logarithmically within the desired range."""
    mi = np.log10(mini)
    ma = np.log10(maxi)
    vals = np.logspace(mi,ma,n)
    if rounded: vals = np.round(vals).astype(int)
    return vals


def expgrid(n_layers, alt_toa, surface_layer_height, round_to = 0, return_constants = False): #Create logarithmic altitude grid
    """
    Calculates and returns an exponential altitude grid.
    The boundariy altitudes are calculated using the equation
    h_n = a * (b**n - 1) with n being the boundary number starting from
    BOA (n = 0 is the surface, n = 1 is the upper boundary of the 1st layer
    and so on... ). Constants a and b are calculated such that the desired
    criteria (n_layers, alt_toa, surface_layer _height) are met.

    Proposed values for trade-off between computational performance
    and accuracy:
    n_layers = 50
    alt_toa = 60000
    surface_layer_height = 50

    Parameters
    ----------
    n_layers : int
        total number of layers of the new grid.
        Maximum is 70. Final grid has
        n_layers + 1 entries, since grid represents layer boundaries.
    alt_toa : float
        Upper altitude limit (TOA) in [m]
    surface_layer_height: float
        Thickness of the lowest layer in  [m]. Defines the exponential
        coefficient and thus the "bending" of the grid.
    return_constants : bool, default = False
        If True, besides the new grid, also the constants a and b of the
        above equation are returned.
    round_to : int or None
        If None: No rounding of altitude values.
        If int: Rounding rules of numpy.round() apply.

    Returns
    -------
    Array of layer boundaries of the exponential grid.
    If return_constants == True, returns tuple (boundaries, a, b)

    """
    N = n_layers
    H = alt_toa
    S = surface_layer_height

    n = np.arange(N+1) #enumerator of layer borders
    #b needs to be derived numerically:
    f_b = lambda x: sum(x**n[:-1]) - H/S
    b0 = (H/S)**(1/N)
    b = fsolve(f_b,b0)
    #calculation of a:
    a = S / (b - 1)
    #calculation of the acutal layer boundary heights
    h_n = a * (b**n - 1)
    if not round_to is None: h_n = np.round(h_n,round_to)
    if round_to == 0: h_n = h_n.astype(int)
    if return_constants: return h_n[::-1], a , b
    return h_n #turn  around to meet VLIDORT sorting convention


def interp(x, xp, yp, kind='cubic', 
           fill_value='extrapolate',
           **interp1d_kwargs):
    """
    Enhances scipy interp1d function to perform stable standard interpolation.
    Automatically reduces interpolation kind in case of few data points.
    For single value yp, this constant value will be used and fill_value
    will be ignored.
    
    Parameters
    ----------
    x : scalar or array, values to interpolate to.
    xp : array. x data.
    yp : array. y data.
    kind : string, optional. The default is 'cubic'.
        Defines the interpolation kind (see scipy.interp1d documentation online)
    fill_value : None, scalar, 2-tuple, 'constant' or 'extrapolate'
        What to do if x_new is outside interpolation range. Options are:
            None : no extrapolation, error is thrown if x_new is outisde range 
            Scalar : use this value to extrapolate both directions
            2-tuple, e.g. (10,0): use these values to interpolate lower and upper direction, repsectively.
            'constant': fill with nearest value.
            'extrapolate': values are extrapolated based on interpolation kind.
        The default is 'extrapolate'.
    **interp1d_kwargs : keyword arguments forwarded to scipy.interp1d.

    Returns
    -------
    1D array. Interpolated value(s)
    """  
    
    # Find right interp1d settings based on fill_value
    if fill_value == None: 
        fv = None
        be = True
    else: 
        be = False 
        fv = fill_value
        if fv == 'constant': fv = (yp[0],yp[-1])
    # Find interpolation kind based on number of values
    if len(xp) == 3 and kind=='cubic': kind = 'quadratic'
    if len(xp) == 2 and kind!= 'linear': kind = 'linear'
    if len(xp) == 1: kind = None
    # interpolate
    if kind == None:
        # warnings.warn('Interpolation on single data point not possible. Using constant value.')
        y = np.array([yp[0]]*len(np.atleast_1d(x)))
        if len(y) == 1: y = y[0]
    else:
        f = interp1d(xp, yp, kind = kind, fill_value=fv, bounds_error=be, **interp1d_kwargs)  
        y = f(x)
        if np.isscalar(x): y = y[()]
    return y


def log_interp(x, xp, yp, 
               xlog=True, ylog=True, 
               kind='linear', fill_value='extrapolate',
               **interp1d_kwargs):
    """ Uses interp() to perform a linear interpolation in logarithmic-space, e.g. for wavelength-dependent parameters with angstrom dependence"""
    xpt = np.log(xp) if xlog else xp
    xt = np.log(x) if xlog else x    
    ypt = np.log(yp) if ylog else yp
    yt = interp(xt, xpt, ypt, kind = kind, fill_value=fill_value, **interp1d_kwargs)  
    if ylog: yt = np.exp(yt)
    return yt

def henyey_greenstein_pf(angles, asymmetry_parameter):
    """ Angles in degree"""
    a = angles/180*np.pi
    g = asymmetry_parameter
    return 1/4/np.pi * (1-g**2) / (1 + g**2 - 2*g*np.cos(a))**(3/2)


def get_henyey_greenstein_moments(asymmetry_parameter, moment_order):
    """ Calculates the moment or order moment_order for a Henyey-Greenstein SPF """
    return (2*moment_order+1) * asymmetry_parameter**moment_order


def smooth_gaussian(array, fwhm):
    """
    Smooths data in an array by convolution with gaussian.
    FWHM is in number of elements --> data must be equidistant!
    """
    sigma = fwhm/2.355
    region = 5*sigma # considered wavelength interval in each direction for kernel
    kernel_x = np.arange(-region,region+1,1)
    kernel = np.exp(-kernel_x**2/sigma**2/2)
    kernel /= kernel.sum()
    return np.convolve(array, kernel, mode='same')


def array_to_rt_input(filename, data, header_comment):
    """ Writes a 1D or 2D array in the common RT-driver input format. """
    d = np.array(data.copy())
    if d.dtype == bool: d = np.array(['T' if i else 'F' for i in d])
    if d.ndim == 1: d = np.atleast_2d(d).T   
    if d.dtype == float: assert ~np.isnan(d).any(), "Nan encountered in RT input array! Not allowed."
    nrows, ncols = d.shape
    header = str(nrows) + ' ! n_rows\n'
    header += str(ncols) + ' ! n_cols\n'
    header += str(header_comment)
    fpath = get_qs_backend_path() + 'INPUT/' + filename
    if d.dtype == float: fmt = '%.15e'
    elif d.dtype == int: fmt = '%i'
    else: fmt = '%s' 
    np.savetxt(fpath, d, header=header, comments='', fmt=fmt)                                                           

def read_array_from_file(f):
    """ Read numpy arrays from csv (.inp, .out, .dat, ...) or binary (.npy) file """
    if not os.path.isfile(f): return None # Check if file exists
    allowed_ftypes = ['inp','out','dat','npy']
    ftype = f.split('.')[-1]
    if not ftype in allowed_ftypes: return None
    func = np.load if ftype == 'npy' else np.genfromtxt # pick right read function to read the file
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        arr = func(f)
    if not len(arr.flatten()): return None # return none if file is empty
    return np.atleast_1d(arr)

    
def save_array_to_file(f):
    """ Same as read array_from_file but for saving """
    location = f.split('/')[:-1]
    Path(location).mkdir(parents=True, exist_ok=True)
    ftype = f.split('.')[-1]
    func = np.save if ftype == 'npy' else np.savetxt # pick right read function to read the file
    return func(f)


def get_stepped_profile(boundaries,values):
    """
    Converts any profile into a step-shaped-profile, by doubling
    the values in boundaries and values. 
    Useful for interpolation and plotting of vertical profiles.
    """
    stepped_boundaries = np.repeat(boundaries,2)[1:-1]
    stepped_values = np.repeat(values,2)
    return stepped_boundaries, stepped_values


def merge_aerosol_bulk_properties(bulkprops, fractions):
    """
    Calculates effective aerosol bulk properties from bulk properties of several aerosol types.    
    Assumes that aerosols share a common profile shape and profile reference wavelength!
    Equations are in the Notes on Aerosol Retrievals, section 'Merging of Aerosol optical properties'.
    
    Parameters
    ----------
    bulkprop_dataframes : list of dataframes with bulk_properties of each aerosol type. 
        Each DataFrame must contain columns 'asymmetry_parameter', 'extinction_multiplier' and 'single_scattering_albedo'
    weights : list of floats 0 < f < 1.0
        Relative contribution of each aerosol to the total extinction at the aerosols' common profile-reference wavelength.

    Returns
    -------
    Effective aerosol bulk properties.
    """
    fractions = np.array(fractions)
    fractions = fractions/fractions.sum()
    exts = pd.concat([df.extinction_multiplier for df in bulkprops],axis=1)
    ext = (exts * fractions).sum(axis=1)
    ssas = pd.concat([df.extinction_multiplier*df.single_scattering_albedo/ext for df in bulkprops],axis=1)
    ssa = (ssas * fractions).sum(axis=1)
    asyms = pd.concat([df.extinction_multiplier*df.single_scattering_albedo*df.asymmetry_parameter/ext/ssa for df in bulkprops],axis=1)
    asym = (asyms * fractions).sum(axis=1)
    df = pd.concat([asym, ext, ssa], axis=1)
    df.columns = ['asymmetry_parameter', 'extinction_multiplier', 'single_scattering_albedo']
    return df


""" Following functions help to retrieve PSD parameters from a given measured PSD """

def psd_monomodal(r,n,rm,sig): # It is sig = log(S). Mie model requires S!!!
    return n / np.sqrt(2*np.pi) / sig / r * np.exp(-np.log(r/rm)**2 / 2 / sig**2)
def psd_bimodal(r, n1, n2, rm1, rm2, sig1, sig2):
    return psd_monomodal(r, n1, rm1, sig1) + psd_monomodal(r, n2, rm2, sig2)
def psd_trimodal(r, n1, n2, n3, rm1, rm2, rm3, sig1, sig2, sig3):
    return psd_monomodal(r, n1, rm1, sig1) + psd_monomodal(r, n2, rm2, sig2) + psd_monomodal(r, n3, rm3, sig3)

def sort_mode_params(p, s):
    """ 
    Sorts a set of retrieved PSD parameters by radius.
    p : list of parameters
    s : list of parameter standard deviations 
    """
    n = int(len(p)/3) # get number of modes from array length
    rs, ns, sigs, srs, sns, ssigs = zip(*sorted(zip(p[n:n*2],p[0:n],p[n*2:n*3],
                                                    s[n:n*2],s[0:n],s[n*2:n*3])))
    return np.hstack([ns,rs,sigs]),np.hstack([sns,srs,ssigs])

def get_monomodal_params_from_psd(r,dNdr):    
    scaler = np.trapz(dNdr,r)
    n_0 = 1
    rm_0 = r[dNdr.argmax()]
    sig_0 = 0.5
    params, cov = curve_fit(psd_monomodal,r,dNdr/scaler,p0=[n_0,rm_0,sig_0])
    sdevs = np.sqrt(np.diag(cov))
    params[0] *= scaler
    return params, sdevs

def get_bimodal_params_from_psd(r,dNdr):    
    scaler = np.trapz(dNdr,r)
    n1_0, n2_0 = 0.5, 0.5
    rm1_0, rm2_0 = r[dNdr.argmax()]/1.5, r[dNdr.argmax()]*1.5 
    sig1_0, sig2_0 = 0.5, 0.5
    params, cov = curve_fit(psd_bimodal,r,dNdr/scaler,p0=[n1_0, n2_0, rm1_0, rm2_0, sig1_0, sig2_0])
    sdevs = np.sqrt(np.diag(cov))
    params, sdevs = sort_mode_params(params, sdevs)
    params[0] *= scaler
    params[1] *= scaler
    return params, sdevs

def get_trimodal_params_from_psd(r,dNdr):    
    scaler = np.trapz(dNdr,r)
    n1_0, n2_0, n3_0 = 0.3, 0.3, 0.3
    rm1_0, rm2_0, rm3_0 = r[dNdr.argmax()]/1.5, r[dNdr.argmax()], r[dNdr.argmax()]*1.5 
    sig1_0, sig2_0, sig3_0 = 0.5, 0.5, 0.5
    params, cov = curve_fit(psd_trimodal, r, dNdr/scaler, p0=[n1_0, n2_0, n3_0, rm1_0, rm2_0, rm3_0, sig1_0, sig2_0, sig3_0])
    sdevs = np.sqrt(np.diag(cov))
    params[0] *= scaler
    params[1] *= scaler
    params[2] *= scaler
    return params, sdevs


def calculate_lapse_rate(surface_relative_humidity): 
    """
    Calculates effective lapse rate in [K/km] from a relative humidity surface measurement in [%],
    according to Beirle, 2022 (https://amt.copernicus.org/articles/15/987/2022/)
    """
    return 7.709 - 4.09 * surface_relative_humidity/100

def get_solar_geometry(lat, lon, datetime_utc):
    """
    Uses the Python soltrack library to infer the solar geometry for given geolocation and time.

    Parameters
    ----------
    lat : float, latitude
    lon : float, longitude
    datetime_utc : string.
        Date and time in string format, e.g. '2024-07-30 11:00:00'        
    """
    st = SolTrack(lon, lat, 
                  use_degrees=True,
                  use_north_equals_zero=True)
    st.set_date_time(datetime_utc)  
    st.compute_position()
    sza = 90 - st._altitude_uncorr # Use solar altitude without refractive correction. Refraction is done in VLIDORT.
    saa = st.azimuth
    return sza, saa

def get_rayleigh_extinction_coeff(wavelength, temp, pres, h2o):
    """
    Calculates the Rayleigh extinction coefficient (km^-1) for moist air.
    
    Parameters:
    wavelength (float): Wavelength in nanometers.
    temp (float): Temperature in Kelvin.
    pres (float): Total atmospheric pressure in hPa.
    h2o (float): Water vapor number density in molecules/cm^3.
    
    Returns:
    float: Extinction coefficient in km^-1.
    """
    # 1. Physical Constants
    kb = 1.380649e-23  # Boltzmann constant (J/K)
    wl_um = wavelength / 1000.0
    wl_cm = wavelength * 1e-7
    
    # 2. Calculate Total Molecular Number Density (N_total) in cm^-3
    # P(Pa) = P(hPa)*100. 1e-6 converts m^-3 to cm^-3.
    n_total = ((pres * 100) / (kb * temp)) * 1e-6
    
    # 3. Separate Dry Air and Water Vapor Densities
    n_dry = n_total - h2o
    if n_dry < 0:
        raise ValueError("Water vapor density cannot exceed total molecular density.")

    # 4. Refractive Index of Dry Air (n-1) at standard conditions
    # Peck and Reeder (1972) formula
    inv_wl_sq = 1 / (wl_um**2)
    n_minus_1_std = (5792105 / (238.0185 - inv_wl_sq) + 
                     167917 / (57.362 - inv_wl_sq)) * 1e-8
    
    # 5. King Factor for Dry Air (Depolarization correction)
    f_air = 1.034 + 3.17e-4 * inv_wl_sq
    
    # 6. Rayleigh Cross-section for Dry Air (sigma_dry) in cm^2/molec
    # N_std = 2.54689e19 molec/cm3 (at 1013.25 hPa, 288.15 K)
    n_std = 2.54689e19
    term_num = 8 * np.pi**3 * ( (n_minus_1_std**2 + 2*n_minus_1_std)**2 )
    term_den = 3 * n_std**2 * wl_cm**4
    sigma_dry = (term_num / term_den) * f_air
    
    # 7. Rayleigh Cross-section for Water Vapor (sigma_h2o)
    # H2O molecules are smaller/less polarizable, scattering ~85% of dry air
    sigma_h2o = 0.85 * sigma_dry 

    # 8. Calculate Total Extinction (beta)
    # beta = (sigma_dry * n_dry) + (sigma_h2o * n_h2o)
    beta_cm_inv = (sigma_dry * n_dry) + (sigma_h2o * h2o)
    
    # Convert cm^-1 to km^-1 (multiply by 100,000)
    beta_km_inv = beta_cm_inv * 1e5
    
    return beta_km_inv

##################
# PARENT CLASSES
##################
""" 
Define basic attributes and methods, inherited by classes actually used in the VLIDORT-PC model.
E.g.: Classes in the Pydort model like PressureProfile or TemperatureProfile inherit from VerticalGrid and GenericProfile.
""" 

class VerticalGrid():
    """
    Defines basic properties for vertical grid and profile type classes.
    Is initialized with the default vertical grid loaded from physics_data/default_vertical_grid.dat.
    Details in the docstrings of the inheriting classes.
    
    Important Attributes
    --------------------
    vertical_grid: 1D array, default defined via DEFAULT_VERTICAL_GRID. 
        Contains altitudes [m] above ground of the grid's layer boundaries. Monotonically increasing.
    """
    def __init__(self):
        self._vertical_grid = DEFAULT_VERTICAL_GRID.copy()
    
    def __repr__(self):
        return 'Vertical grid with layer boundaries:\nvertical_grid: ' + str(list(self.vertical_grid))

    def copy(self):
        return copy.deepcopy(self)
    
        
    @property
    def vertical_grid(self):
        return self._vertical_grid
        
    @vertical_grid.setter
    def vertical_grid(self, value):
        self._vertical_grid = np.array(value)
        
    def expgrid(self, n_layers, alt_toa, surface_layer_height, round_to = 0):
        self._vertical_grid = expgrid(n_layers, alt_toa, surface_layer_height, round_to, return_constants = False)
        
    def get_n_boundaries(self): return len(self.vertical_grid)

    def get_n_layers(self): return self.get_n_boundaries()-1
    
    def get_layer_thicknesses(self): 
        """ Calculates the thicknesses of each layer [m]"""
        return np.diff(self.vertical_grid)

    def get_layer_middles(self): 
        """ Calculates altitudes of the layer middles [m]"""
        return centers(self.vertical_grid)
    
    def get_layer_bottoms(self): 
        return self.vertical_grid[:-1]

    def get_layer_tops(self): 
        return self.vertical_grid[1:]
    
    def get_closest_layer_boundary(self, altitude):
        """ Returns index and altitude of the layer boundary closest to a given altitude [m]"""
        ix = np.argmin(np.abs(self.vertical_grid-altitude))
        alt = self.vertical_grid[ix]
        return ix, alt
    
    def get_closest_layer_middle(self, altitude):
        """ Returns index and altitude of the layer boundary closest to a given altitude [m]"""
        ix = np.argmin(np.abs(self.get_layer_middles()-altitude))
        alt = self.get_layer_middles()[ix]
        return ix, alt


class GenericProfile(VerticalGrid):
    """
    Extends the VerticalGrid base class to a general vertical profile.
    Class contains altitude grid, layer values and additional methods e.g. for interpolation.
    ExtinctionProfile, ConcentrationProfile, etc. inherit from this class.
    
    Important Attributes
    --------------------
    vertical_grid: 1D array, default defined via DEFAULT_VERTICAL_GRID. 
        Contains altitudes [m] above ground of the grid's layer boundaries. Monotonically increasing.
    values: 1d array, length is one less than vertical grid, default is zero array.
        Contains average values for each layer (e.g. temperature, pressure, gas concentration, ...).
    """
    
    def __init__(self):
        super().__init__()
        self.values = np.zeros(len(self.vertical_grid)-1)
        
    def __repr__(self):
        s = 'Vertical profile object with:\n'
        s += self.to_dataframe().__repr__()
        return s
        
    def copy(self):
        return copy.deepcopy(self)
      
    def fill_constant(self, value):
        """ Fills profile with cinstant value """
        self.values[:] = value
      
    def interpolate_from(self, altitudes, values, kind = 'cubic',
                         fill_value = None,
                         inplace=True):
        """
        Interpolates the provided profile information to the current grid 
        (vertical_grid) of the GenericProfile instance and inserts it by
        replacing values.
        Applies scipy.interpolate.interp1d to interpolate provided profile 
        data to high resolution first (defined by PROFILE_HIGH_RESOLUTION), 
        then calculates averages for each layer.

        Parameters
        ----------
        altitudes : 1d-array 
            Altitudes of the new profile (altitude of measurement, not layer boundaries)
        values : 1d-array
            Value at each altitude
        kind : str
            Interpolation kind to apply in the interp1d method (see scipy.interp1d documentation online)
    fill_value : None, scalar, 2-tuple or string.
        What to do if current grid exceeds the data interpolation range. Options are:
            None : no extrapolation, error is thrown if interpolation range is exceeded range 
            Scalar : use this value to extrapolate both directions
            2-tuple, e.g. (10,0): use these values to interpolate lower and upper direction, repsectively.
            'constant': fill with nearest value.
            'extrapolate': values are extrapolated based on interpolation kind.
        The default is 'extrapolate'.
    **interp1d_kwargs : keyword arguments forwarded to scipy.interp1d.
        inplace : bool, default = True
            If True, data will be interpolated into values.
            If False, interpolated data is returned, values remains unchanged.
        """
        xhr = np.arange(self.vertical_grid[0], self.vertical_grid[-1], PROFILE_HIGH_RESOLUTION)
        yhr = interp(xhr, altitudes, values, kind = kind, fill_value=fill_value)
        digitized = np.digitize(xhr, self.vertical_grid[1:])
        ynew = [yhr[digitized == i].mean() for i in range(0, self.get_n_layers())]
        ynew = np.array(ynew)
        if inplace: self.values = ynew
        else: return ynew
        
    def get_interpolated_value(self, altitude, kind = 'linear'):
        """ Get interpolated value at altitude [m]. Values (layer averages) are associated with layer middles for interpolation. """
        xp = self.get_layer_middles()
        yp = self.values
        y = interp(altitude, xp, yp, kind=kind, fill_value='extrapolate')
        return y

    def plot(self, new_fig = True, scaling=1, boundary_indexed = False, **plot_kwargs):
        """
        Plots the profile 

        Parameters
        ----------
        new_fig : bool
            Shall a new figure be created for the plot? If False, profile will be plotted into existing figure.
        scaling: scalar
            Scale data prior to plotting. Useful for plotting of multiple profiles into single plot.
        boundary_indexed: boolean,
            If True, vertical axis will show layer indexes instead of altitude.
        **plot_kwargs: kwargs
            Keyword arguments passed to matplotlib.pyplot.plot()
        """
        if boundary_indexed: x = np.arange(len(self.vertical_grid))
        else: x = self.vertical_grid/1e3
        if new_fig: plt.figure()
        x, y = get_stepped_profile(x, self.values)
        #plot layer borders
        for a in x: plt.axhline(a, color = 'k', alpha = 0.1, lw=0.5)
        plt.plot(y*scaling, x, **plot_kwargs)
        #global plot properties
        plt.ylabel('Altitude [km]')

    def change_vertical_grid(self, new_vertical_grid, kind = 'linear', upper_fill_value = 'extrapolate'):
        """
        Changes vertical_grid to a new grid and interpolates values correspondingly (inplace).
        The old profile is assumed to be step shaped and values are 'rebinned' into the new grid. 
        After interpolation, values below zero are set to zero!

        Parameters:
        -----------
        new_grid : Numpy Array
            New layer BOUNDARIES (hence, the new vertical_grid)
            that the profile will be interpolated to. Must be monotonic.
        kind : string, default is 'linear' (recommended since safe)
            Passed to the scipy interp1d function
        upper_fill_value: string or float,
            Defines how to fill high altitude layers, if former profile does not provide any data. Options are:
            'extrapolate' : values are extrapolated (can be dangerous!) Values below zero, will be set to zero.
            'last' : layers are filled with last value (constant)
            float : layers are filled with provided value
        """
        # vc = self.get_vertical_column()
        assert len(self.vertical_grid) == len(self.values)+1, 'Interpolation of profile not possible. Number of layers and number of values not consistent.'
        xnew = np.array(new_vertical_grid)
        if (np.diff(xnew) < 0).all(): #if its the wrong order turn it around (low to high altitudes)
            xnew = xnew[::-1]
        assert (np.diff(xnew) > 0).any(), 'Provide monotonically increasing/decreasing grid!'
        #Strange rebinning procedure based on the idea of integrating over the step profile:
        xhr = np.arange(xnew[0], xnew[-1], PROFILE_HIGH_RESOLUTION)
        x, y = self.get_layer_middles(), self.values
        fint = interp1d(x, y, kind = kind, fill_value = 'extrapolate', bounds_error = False)
        yhr = fint(xhr)
        if not upper_fill_value == 'extrapolate':
            if upper_fill_value == 'last': upper_fill_value = y[-1]
            yhr[xhr > x.max()] = upper_fill_value
        yhr[yhr<0] = 0.
        digitized = np.digitize(xhr, xnew[1:])
        ynew = np.array([yhr[digitized == i].mean() for i in range(0, len(xnew)-1)])
        self.vertical_grid = xnew
        self.values = ynew
    
    def to_dataframe(self):
        df =  pd.DataFrame({'lower_boundary': self.vertical_grid[:-1],
                            'upper_boundary': self.vertical_grid[1:],
                            'layer_middle': self.get_layer_middles(),
                            'value' : self.values})
        df.index.name = 'layer_index'
        return df
        


class Container():
    """
    Simple container for multiple variables, typically scalars, arrays or dataframes.
    Very similar to a python dictionary but uses the dot syntax and enables autocompletion in IDEs.
    Used e.g. for settings, result data, etc...
    """
    
    def __init__(self, attributes = []):
        for att in attributes:
            setattr(self, att, None)

    def __repr__(self):
        d = self.to_dictionary()
        s =  type(self).__name__ + ' object with members:'
        s += ''.join(['\n'+i+': '+str(d[i]) for i in d if i[0]!='_'])
        return s

    def copy(self):
        return copy.deepcopy(self)
            
    def to_dictionary(self):
        return self.__dict__
    


class FrozenContainer(Container):
    """
    Like container, but after initialisation, no further attributes can be added.
    This is useful for simulation settings: if not frozen, typos create new attributes 
    instead of setting existing ones, leading to unexpected behaviour.
    """
    def __init__(self, attributes = [], properties = []):
        super().__setattr__('_allowed_attributes', attributes+properties)
        for att in attributes:
            super().__setattr__(att, None)
            
    def copy(self):
        return copy.deepcopy(self)
    
    def __setattr__(self, name, value):
        if not name in self._allowed_attributes: raise Exception('Adding new attributes is not allowed!')
        else: super().__setattr__(name, value)



#################
# MODEL CLASSES #
#################


class ViewingGeometry(FrozenContainer):
    """
    Defines viewing geometry and observer height to apply in the simulation. 
    Atmospheric refraction is considered within VLIDORT.
    For simulation of actinic fluxes, only observer height and SZA are of relevance.
    
    Important Attributes:
    ---------------------
    
    observer_height_agl: integer, must be greater > 0.
        Height of observer above ground level in meters.
        Above TOA, the Fortran drivers trace the viewing ray to its TOA entry
        and transform the viewing and solar angles from the observer's local
        frame. Supported for radiance only, for rays entering the atmosphere
        on the dayside. Angles supplied here are always at the actual observer.
        
    sza: float, 0 <= sza <= 90 
        Solar zenith angle in degree. Zero is zenith. 90 is horizon.
    
    saa: float, arbitray value (cyclic).
        Solar azimuth angle in degree. Reference (zero) direction and rotation direction must be consistent with vaa.
    
    vza: float, 0 <= vza <= 180
        Viewing azimuth angle in degree. Zero is zenith, 90 is limb, 180 is nadir.
        
    vaa: float, any value possible, reference direction and rotation direction must be consistent with saa. 
        Viewing azimuth angle in degree.
    """
    def __init__(self):
        super().__init__(['sza',
                          'saa',
                          'vza',
                          'vaa',
                          'observer_height_agl'])
        self.sza = 30
        self.saa = 0
        self.vza = 60
        self.vaa = 90
        self.observer_height_agl = 0
    
    def copy(self):
        return copy.deepcopy(self)
        
    def get_solar_geometry(self, lat, lon, datetime_utc, inplace = True):
        """
        Uses the Python soltrack library to infer the solar geometry for given geolocation and time.

        Parameters
        ----------
        lat : float, latitude
        lon : float, longitude
        datetime_utc : string.
            Date and time in string format, e.g. '2024-07-30 11:00:00'
        inplace : bool, default is True
            If True, sza and saa are inserted into the current viewing_geometry object.
            If False, Mathod returns a (sza, saa) tuple. 
        """
        st = SolTrack(lon, lat, 
                      use_degrees=True,
                      use_north_equals_zero=True)
        st.set_date_time(datetime_utc)  
        st.compute_position()
        sza, saa = get_solar_geometry(lat, lon, datetime_utc)
        if inplace:
            self.sza = sza
            self.saa = saa
        else: 
            return sza, saa
        
    def get_rt_raa(self):
        raa = (self.vaa - self.saa)%360
        raa_mirror_flag = False
        if raa > 180:
            raa = 360 - raa
            raa_mirror_flag = True
        return raa, raa_mirror_flag

    def check_ranges(self):
        assert self.sza >= 0, "Ensure SZA >= 0!"
        assert self.sza <= 90, "Ensure SZA <= 90!"
        assert self.vza >= 0, "Ensure VZA >= 0!"
        assert self.vza <= 180, "Ensure VZA <= 180!"
        assert self.observer_height_agl >= 0, "Ensure observer_height_agl >= 0 m!"

#TEST:
# vg = ViewingGeometry()
# vg.get_solar_geometry(45, 180, '2000-03-21 00:00:00', False)

# vg.saa = 359
# vg.vaa = 181
# vg.get_rt_raa()


class SettingsRT(FrozenContainer):
    """
    A frozen container with settings controlling the radiative transfer
    backend routing and simulation options.
    
    Important Attributes
    --------------------
    
    sphericity_correction: int, default is 1.
        Selects the spherical-correction accuracy and the Fortran backend:
        0 uses experimental sparse, transmission-normalized aircraft MS sources.
        Native spherical FO and batched Standard MS share one driver call.
        Supports sunlit limb/turning rays and observers between layer boundaries;
        it remains an approximation, especially near the low-altitude horizon.
        1 uses Standard VLIDORT with spherical first-order scattering but
        no multiple-scattering spherical correction, 2 uses the
        2-point interpolated correction, 3 uses the 3-point interpolated
        correction at/above TOA and VLIDORT-QS below TOA, and 4 or higher uses
        VLIDORT-QS full spherical correction. Unsupported original 2-point geometries
        (including rays whose extension misses the surface, and observers
        between layer boundaries) fall back to QS with a warning.
    do_actinic_flux_at_observer: bool, default is False.
        Set to true if you require actinic fluxes at the exact observer altitude only.
    do_irradiance_at_observer: bool, default is False.
        Calculate cosine-weighted irradiance on a local horizontal plane,
        at the observer altitude (at or below TOA only).
    do_irradiance_jacobians: bool, default is False.
        Calculate layer-resolved irradiance derivatives for jacobian_list.
        Requires do_irradiance_at_observer=True.
    do_irradiance_scattering_components: bool, default is False.
        Separate first-order and multiply scattered irradiance. One surface
        reflection counts as first order. Requires a full calculation without
        delta-M scaling; adds four weak-scattering flux solves per wavelength.
    irradiance_scattering_step: float, default is 1e-3.
        Expansion step (1e-6 through 1e-2). Two extrapolations are compared
        internally; a failed convergence check raises a backend error.
    do_radiance: bool, default is True.
        Set to true if you require radiances at the observer and the specified viewing geometry.
    do_actinic_flux_jacobians: bool, default is False.
        Set to true 
    do_vectorized_simulation: bool, default is False.
        If True, radiative transfer considers full stokes vector (radiances and light polarisation state), instead of scalar radiances only.
        Much slower than scalar simulation and typically not required, unless polarization state shall be analyzed. 
    do_single_scattering_only: bool, default is False.
        If True, the contribution from multiply scattered light is omitted.
    n_full_sphericity_fine_layers: int, default is 4.
        Subdivision of atmospheric layers used by full-sphericity calculations.
    n_aircraft_source_points: int, default is 3.
        Number of local MS solutions for experimental correction 0 (2 to 33).
    jacobian_layer_flags: String ('all') or list of booleans with length equal to the number of profile layers.
        Limits the jacobian calculation to distinct layers to improve computational efficiency.
        if 'all', jacobian calculation will be performed inall layers.
        Otherwise, provide list of jacobian calculation flag for each layer starting from the bottom. 
        E.g. [True, ... , True, False, ... False].
        Note: calculation of jacobians can be                                         
    """
    def __init__(self):
        super().__init__(['sphericity_correction',
                          'n_aircraft_source_points',
                          'do_actinic_flux_at_observer',
                          'do_irradiance_at_observer',
                          'do_irradiance_jacobians',
                          'do_irradiance_scattering_components',
                          'irradiance_scattering_step',
                          'do_radiance',
                          'do_actinic_flux_jacobians',
                          'do_radiance_jacobians',
                          'jacobian_list',
                          'jacobian_layer_flags',
                          'n_halfspace_streams',
                          'n_max_pf_moments',
                          'do_o4_absorption',
                          'do_deltam_phase_function_scaling',
                          'planet_radius',
                          'do_vectorized_simulation',
                          'do_single_scattering_only',
                          'n_full_sphericity_fine_layers'])
        self.sphericity_correction = SphericityCorrection.NONE.value
        self.n_aircraft_source_points = 3
        self.do_actinic_flux_at_observer = False
        self.do_irradiance_at_observer = False
        self.do_irradiance_jacobians = False
        self.do_irradiance_scattering_components = False
        self.irradiance_scattering_step = 1e-3
        self.do_radiance = True
        self.do_actinic_flux_jacobians = False
        self.do_radiance_jacobians = False
        self.jacobian_list = ['absorption']
        self.jacobian_layer_flags = 'all'
        self.n_halfspace_streams = 4
        self.n_full_sphericity_fine_layers = 4
        self.n_max_pf_moments = MAX_PF_MOMENTS
        self.do_o4_absorption = True
        self.do_deltam_phase_function_scaling = True
        self.planet_radius = 6371.0
        self.do_vectorized_simulation = False
        self.do_single_scattering_only = False

    def copy(self):
        return copy.deepcopy(self)

    def validate_backend_settings(self):
        if self.do_irradiance_scattering_components:
            if not self.do_irradiance_at_observer or self.do_single_scattering_only:
                raise ValueError('Scattering components require irradiance and a full-scattering calculation.')
            if not 1e-6 <= self.irradiance_scattering_step <= 1e-2:
                raise ValueError('irradiance_scattering_step must be between 1e-6 and 1e-2.')
        if self.do_irradiance_jacobians and not self.do_irradiance_at_observer:
            raise ValueError('Enable do_irradiance_at_observer when requesting irradiance Jacobians.')
        raw_correction = self.sphericity_correction
        if isinstance(raw_correction, bool):
            raise AssertionError(
                'settings_rt.sphericity_correction must be an integer: '
                '0=experimental sparse-source correction, '
                '1=no correction, 2=2-point interpolation, '
                '3=3-point interpolation, >=4=full spherical correction.'
            )
        try:
            correction = int(raw_correction)
        except (TypeError, ValueError):
            raise AssertionError(
                'settings_rt.sphericity_correction must be an integer: '
                '0=experimental sparse-source correction, '
                '1=no correction, 2=2-point interpolation, '
                '3=3-point interpolation, >=4=full spherical correction.'
            )
        if isinstance(raw_correction, str):
            assert raw_correction.strip() == str(correction), (
                'settings_rt.sphericity_correction must be an integer: '
                '0=experimental sparse-source correction, '
                '1=no correction, 2=2-point interpolation, '
                '3=3-point interpolation, >=4=full spherical correction.'
            )
        else:
            assert float(raw_correction) == float(correction), (
                'settings_rt.sphericity_correction must be an integer: '
                '0=experimental sparse-source correction, '
                '1=no correction, 2=2-point interpolation, '
                '3=3-point interpolation, >=4=full spherical correction.'
            )
        assert correction >= 0, (
            'settings_rt.sphericity_correction must be >= 0: '
            '0=experimental sparse-source aircraft correction, '
            '1=no correction, 2=2-point interpolation, '
            '3=3-point interpolation, >=4=full spherical correction.'
        )
        self.sphericity_correction = correction
        points = self.n_aircraft_source_points
        assert not isinstance(points, bool) and isinstance(points, (int, np.integer)) and 2 <= points <= 33, (
            'settings_rt.n_aircraft_source_points must be an integer from 2 to 33.'
        )


class PressureProfile(GenericProfile):
    """
    Contains information on the atmospheric pressure profile. Inherits from the GenericProfile class.
    
    Important Attributes
    --------------------
    
    vertical_grid: 1D array, default defined via DEFAULT_VERTICAL_GRID. 
        Contains altitudes [m] above ground of the grid's layer boundaries. Monotonically increasing.
    values: 1d array, length is one less than vertical grid, default is zero array.
        Contains average pressures [hPa] in each layer.
    
    """
    def __init__(self):
        super().__init__()
        d = np.loadtxt(PATH_PHYSICS + 'profiles/pressure_profile.dat', comments = '#')
        self.interpolate_from(d[:,0],d[:,1])

    def copy(self):
        return copy.deepcopy(self)
    
    def fill_exponential_profile(self, scale_height = 8000, surface_value = 1013):
        """
        Fills values with an exponentially decaying profile assuming the current vertical grid.
        
        Parameters
        ----------
        scale_height : scalar
            Scale height of exponential profile [m]
        surface_value : scalar
            Surface pressure [hPa]
        """
        #Create high resolution exponential profile
        xhr = np.arange(self.vertical_grid.min(), self.vertical_grid.max() + PROFILE_HIGH_RESOLUTION, PROFILE_HIGH_RESOLUTION)
        yhr = surface_value*np.exp(-xhr/scale_height)
        self.interpolate_from(xhr, yhr)
    
    def scale_from_point_measurement(self, meas_value, meas_altitude):
        """
        Scales the whole pressure profile based on a pressure measurement meas_value [hPa] at altitude meas_altitude [m].
        """
        sf = meas_value/self.get_interpolated_value(meas_altitude,'cubic')
        self.values *= sf

class TemperatureProfile(GenericProfile):
    """
    Contains information on the atmospheric temperature profile. Inherits from the GenericProfile class.
    
    Important Attributes
    --------------------
    
    vertical_grid: 1D array, default defined via DEFAULT_VERTICAL_GRID. 
        Contains altitudes [m] above ground of the grid's layer boundaries. Monotonically increasing.
    values: 1d array, length is one less than vertical grid, default is zero array.
        Contains average temperatures [K] in each layer.
    
    """
    def __init__(self):
        super().__init__()
        d = np.loadtxt(PATH_PHYSICS + 'profiles/temperature_profile.dat', comments = '#')
        self.interpolate_from(d[:,0],d[:,1])
        
    def copy(self):
        return copy.deepcopy(self)
    
    def fill_constant(self, temperature = 300):
        self.values = np.array([temperature]*self.get_n_layers())
                
    def scale_troposphere_from_point_measurement(self, meas_value, meas_altitude, lapse_rate = 6.5, tropopause_temperature = 216.5, transition = (13000,18000)):
        """
        transition : tuple, defines altitude interval in which the new tropospheric profile  
        """
        meas_altitude = round(meas_altitude)
        lapse_rate /= 1e3 # to [K/m]
        surface_temp = meas_value + lapse_rate * meas_altitude
        alts = self.get_layer_middles()
        # Create high resolution linear tropospheric profile with low limit:
        tp_lin = surface_temp - lapse_rate*alts
        tp_lin[tp_lin<tropopause_temperature] = tropopause_temperature
        # Interpolate existing profile to hr:
        tp_old = self.values
        # Create transition weights:
        a = 1/(transition[0]-transition[1])
        b = - a * transition[1]
        w = a*alts + b
        w[w>1] = 1
        w[w<0] = 0
        # Create new hr profile:
        tp = w*tp_lin + (1-w)*tp_old
        self.values = tp
        
        
        
class ExtinctionProfile(GenericProfile):
    """
    Describes an aerosol vertical distribution in terms of extinction 
    coefficient [1/km] at a specific reference wavelength. The reference wavelength 
    defined later within the parent OpticalAerosol or MicrophysicalAerosol object. 
    Inherits from the GenericProfile class.
    
    Important Attributes
    --------------------
    
    vertical_grid: 1D array, default defined via DEFAULT_VERTICAL_GRID. 
        Contains altitudes [m] above ground of the grid's layer boundaries. Monotonically increasing.
    values: 1d array, length is one less than vertical grid, default is zero array.
        Contains average extinction coeffcients [1/km] for each layer.
    
    """

    def copy(self):
        return copy.deepcopy(self)
         
    def get_vertical_column(self):
        """ Calculates the optical depth (AOD) of the current profile """
        y = self.values
        dx = self.get_layer_thicknesses()/1e3
        return (y*dx).sum()

    def set_vertical_column(self, optical_depth = 1.0):
        """ Rescales the profile to the provided optical depth (AOD) """
        if not self.get_vertical_column(): self.fill_constant(1)
        self.values = self.values / self.get_vertical_column() * optical_depth
    
    def fill_exponential_profile(self, scale_height = 1000, optical_depth = 0.2):
        """
        Fills values with an exponentially decaying profile.
        
        Parameters
        ----------
        scale_height : scalar
            Scale height of exponential profile [m].
        optical_depth : scalar
            Optical depth (AOD) of total profile.
        """
        #Create high resolution exponential profile
        xhr = np.arange(self.vertical_grid.min(), self.vertical_grid.max() + PROFILE_HIGH_RESOLUTION, PROFILE_HIGH_RESOLUTION)
        yhr = np.exp(-xhr/scale_height)
        self.interpolate_from(xhr, yhr)
        self.set_vertical_column(optical_depth)

    def fill_gaussian_profile(self, center_altitude = 1000, fwhm = 1000, optical_depth = 0.2):
        """
        Fills values with a gaussian shaped profile (in place).

        Parameters
        ----------
        center_altitude : float
            Center altitude [m] of the gaussian profile.
        fwhm : float
            Full width half maximum of the gaussian profile [m]
        optical_depth : float
            Optical depth (AOD) of total profile.
        """
        xhr = np.arange(self.vertical_grid.min(), self.vertical_grid.max() + PROFILE_HIGH_RESOLUTION, PROFILE_HIGH_RESOLUTION)
        yhr = gaussian_func(xhr, center_altitude,fwhm)
        self.interpolate_from(xhr, yhr)
        self.set_vertical_column(optical_depth)

    def fill_fermi_profile(self, height = 1000, blur = 100, optical_depth = 0.2):
        """
        Fills values with a fermi shaped profile (in place).
        Equation: 1/(exp((h - height)/blur) + 1), with h being altitude.

        Parameters
        ----------
        height : float
            Height of the layer.
        blur : float
            Blur of the layer edge.
        optical_depth : float
            Optical depth (AOD) of total profile.
        """
        #Create high resolution exponential profile
        xhr = np.arange(self.vertical_grid.min(), self.vertical_grid.max() + PROFILE_HIGH_RESOLUTION, PROFILE_HIGH_RESOLUTION)
        #Precalulate exponent and calculate different regions separately to avoid overflow:
        if not blur: blur = 1e-3
        exponents = (xhr - height)/blur
        yhr = np.empty(len(exponents))
        yhr[exponents <= -20] = 1.
        yhr[exponents >= 20]  = 0.
        m = (exponents > -20) & (exponents < 20)
        yhr[m] = 1/(np.exp(exponents[m]) + 1)  
        self.interpolate_from(xhr, yhr)
        self.set_vertical_column(optical_depth)
    

class GasProfile(GenericProfile):
    """
    Describes a trace gas vertical distribution in terms of number concentrations 
    [molec/cm3]. Inherits from the GenericProfile class.
    
    Important Attributes
    --------------------
    
    vertical_grid: 1D array, default defined via DEFAULT_VERTICAL_GRID. 
        Contains altitudes [m] above ground of the grid's layer boundaries. Monotonically increasing.
    values: 1d array, length is one less than vertical grid, default is zero array.
        Contains average number concentrations [molec/cm3] for each layer.
    
    """

    def copy(self):
        return copy.deepcopy(self)
 
    def get_vertical_column(self):
        y = self.values
        dx = self.get_layer_thicknesses()*1e2
        return (y*dx).sum()

    def set_vertical_column(self, vertical_column = 1e19):
        """
        Normalizes the profile to the provided vertical column (in place).
        """
        if not self.get_vertical_column(): self.fill_constant(1)
        self.values = self.values / self.get_vertical_column() * vertical_column
    
    def fill_constant(self, number_concentration):
        self.values = np.array([number_concentration]*self.get_n_layers())

    def fill_exponential_profile(self, scale_height = 1000, vertical_column = 1e19):
        """
        Parameters
        ----------
        scale_height : scalar
            Scale height of exponential profile [m]
        vertical_column : scalar
            Vertical column density of total profile [1/cm2].
        """
        #Create high resolution exponential profile
        xhr = np.arange(self.vertical_grid.min(), self.vertical_grid.max() + PROFILE_HIGH_RESOLUTION, PROFILE_HIGH_RESOLUTION)
        yhr = np.exp(-xhr/scale_height)
        self.interpolate_from(xhr, yhr)
        self.set_vertical_column(vertical_column)

    def fill_gaussian_profile(self, center_altitude = 1000, fwhm = 1000, vertical_column = 1e19):
        """
        Fills values with a gaussian shaped profile (in place).

        Parameters
        ----------
        center_altitude : float
            Center altitude [m] of the gaussian profile.
        fwhm : float
            Full width half maximum of the gaussian profile [m]
        vertical_column : float
            Vertical column density of total profile [1/cm2].
        """
        #Create high resolution exponential profile
        xhr = np.arange(self.vertical_grid.min(), self.vertical_grid.max() + PROFILE_HIGH_RESOLUTION, PROFILE_HIGH_RESOLUTION)
        yhr = gaussian_func(xhr, center_altitude,fwhm)
        self.interpolate_from(xhr, yhr)
        self.set_vertical_column(vertical_column)

    def fill_fermi_profile(self, height = 1000, blur = 100, vertical_column = 1e19):
        """
        Fills values with a fermi shaped profile (in place).
        Equation: 1/(exp((h - height)/blur) + 1), with h being altitude.

        Parameters
        ----------
        height : float
            Height of the layer.
        blur : float
            Blur of the layer edge.
        vertical_column : float
            Vertical column density of total profile [1/cm2].
        """
        #Create high resolution exponential profile
        xhr = np.arange(self.vertical_grid.min(), self.vertical_grid.max() + PROFILE_HIGH_RESOLUTION, PROFILE_HIGH_RESOLUTION)
        #Precalulate exponent and calculate different regions separately to avoid overflow:
        if not blur: blur = 1e-3
        exponents = (xhr - height)/blur
        yhr = np.empty(len(exponents))
        yhr[exponents <= -20] = 1.
        yhr[exponents >= 20]  = 0.
        m = (exponents > -20) & (exponents < 20)
        yhr[m] = 1/(np.exp(exponents[m]) + 1)  
        self.interpolate_from(xhr, yhr)
        self.set_vertical_column(vertical_column)
        
    def get_partial_pressures(self, temperature):
        """ Converts the number concentration profile to a partial pressure profile, for a given temperature [K] or a temperature profile. """
        return self.values*1e6 * KB * temperature / 1e2
        

# Tests:
# p = PressureProfile()
# p.fill_exponential_profile(8000,1013)
# p.plot()
# p.interpolate_from([0,60000], [1013,0], kind='linear')
# p.plot(new_fig=False)
# p.change_vertical_grid(expgrid(10,60000,200))
# p.plot(new_fig=False) 

# p = GasProfile()
# p.fill_fermi_profile()
# p.plot()

class DiscreteSpectralParameter():

    """
    Describes a wavelength-dependent parameter, by providing values wavelength
    by wavelength. When calling pyd.run_rt(), this parameter must contain values
    for all simulation wavelenths (as defined by pyd.wavelengths). Interpolation
    to these wavelengths is facilitated by the change_wavelengths method.

    Important attributes:
    ---------------------
    
    wavelengths : 1D-array
        Wavelengths of each datapoint.
    values : 1D array
        Values associated with the wavelengths.
    """
    def __init__(self,
                 wavelengths = DEFAULT_WAVELENGTHS.copy(),
                 values = None): 
        self._wavelengths = wavelengths
        if values is None: self._values = np.zeros(len(wavelengths))
        else: self._values = np.array(values)
        
    def __repr__(self):
        s = ''.join(['DiscreteSpectralParameter with:',
                     '\n\twavelengths: ' + str(self.wavelengths),
                     '\n\tvalues: ' + str(self.values)])
        return s

    def copy(self):
        return copy.deepcopy(self)
     
    @property
    def wavelengths(self):
        return self._wavelengths
        
    @wavelengths.setter
    def wavelengths(self, value):
        self._wavelengths = np.array(value)
        
    @property
    def values(self):
        return self._values
        
    @values.setter
    def values(self, value):
        self._values = np.array(value)  
    
    def _get_data_from_default_file(self, file): #Import of profile from a file
        """
        Imports data from a file. Assumes following file format:
            First column: wavelength [nm]
            Second column: parameter value
            Whitespaces as delimiter.
            Comments at beginning of file (header) start with '#'
            
        file: relative path (with respect to the copied physics_data folder) and filename.
        """
        d = np.loadtxt(PATH_PHYSICS + file, comments = '#')
        self.wavelengths = d[:,0]
        self.values = d[:,1]
    
    def change_wavelengths(self, new_wavelengths, inplace = True):
        """ Interpolates data to a new wavelength grid. Uses log-log-interpolation """
        ynew = log_interp(new_wavelengths, self.wavelengths, self.values)
        if inplace:
            self.wavelengths = new_wavelengths.copy()
            self.values = ynew
        else:
            return ynew
    
    def get_value(self, wavelength):
        """ Return value for specific wavelength(s). 
        wavelength : scalar or 1D-array.
            Wavelength(s) must be included in self.wavelengths """
        wl = np.atleast_1d(wavelength)
        assert all([w in self.wavelengths for w in wl]), 'DiscreteSpectralParameter does not contain requested wavelength. Use change_wavelengths() to interpolate?'
        m = [w in wl for w in self.wavelengths]
        vals = self.values[m]
        if len(vals) == 1: vals = vals[0] #Convert to scalar if single value was requested.
        return vals
    
    def interpolate_from(self, wavelengths, values, kind='cubic', 
                         fill_value='extrapolate', **interp1d_kwargs):
        """
        Interpolates the provided data to the current wavelengths 
        of the DiscreteSpectralParameter object replacing its values.
        Uses log_log_interpolation.

        Parameters
        ----------
        wavelengths : 1d-array 
            Wavelengths of the dataset to interpolate from.
        values : 1d-array
            Values at each wavelength.
        kind : string, optional. The default is 'cubic'.
            Defines the interpolation kind (see scipy.interp1d documentation online).
            Note: is applied in log-log-space
        fill_value : None, scalar, 2-tuple, 'constant' or 'extrapolate'
            What to do if x_new is outside interpolation range. Options are:
                None : no extrapolation, error is thrown if x_new is outisde range 
                Scalar : use this value to extrapolate both directions
                2-tuple, e.g. (10,0): use these values to interpolate lower and upper direction, repsectively.
                'constant': fill with nearest value.
                'extrapolate': values are extrapolated based on interpolation kind.
            The default is 'extrapolate'.
        **interp1d_kwargs : keyword arguments forwarded to scipy.interp1d.
        """
        self.values = log_interp(self.wavelengths, wavelengths, values, 
                                 kind=kind, fill_value=fill_value, **interp1d_kwargs)

    def set_constant(self, value):
        """ Fills values with a constant for all wavelengths """
        self.values = np.ones(len(self.wavelengths))*value
        
    def plot(self, new_fig = True):
        if new_fig: plt.figure()
        plt.plot(self.wavelengths,self.values,'o')
        plt.xlabel('Wavelength [nm]')
        
# Tests:
# p = DiscreteSpectralParameter([0.8,0.5,0.3,0.2],[300,400,500,600])
# p.plot()
# p.change_wavelengths(np.linspace(200,700))
# p.plot(new_fig=False)
# plt.close('all')
# p = SpectralParameter(0.1,DEFAULT_REF_WL,[5,1])
# p.plot()


class SpectralParameter():
    
    """
    Describes a wavelength-dependent parameter, with a smooth (log-polynomial) 
    wavelength-dependence, according to.
    ln(x) = ln(x_0) + sum [ alpha_i * ln(lambda_0/lambda)**i ]
    For more details, see manual and paper sections 
    "Wavelength-dependent parameters".

    Important attributes:
    ---------------------
    
    value : scalar
        Value of the parameter x_0 at the reference wavelength lmabda0.
    reference_wavelength : scalar
        Reference wavelength lambda_0[nm] of the value.
    angstrom_parameters : list
        List of angstrom parameters alpha_i (increasing order), describing the 
        wavelength dependence of the parameter. 
        Provide a single value (e.g. angstrom_parameters = [1.7]),
        to obtain a classic angstrom dependence.
        Provide zero (angstrom_parameters = [0]), for spectrally constant parameters.
    """
    def __init__(self,
                 value = 1,
                 reference_wavelength = DEFAULT_REF_WL,
                 angstrom_parameters = [0]):
        self.value = value
        self.reference_wavelength = reference_wavelength
        self.angstrom_parameters = np.atleast_1d(angstrom_parameters)
        
    def __repr__(self):
        s = ''.join(['SpectralParameter with:',
                     '\n\tvalue: ' + str(self.value),
                     '\n\treference_wavelength: ' + str(self.reference_wavelength),
                     '\n\tangstrom_parameters: ' + str(self.angstrom_parameters)])
        return s
        

    def copy(self):
        return copy.deepcopy(self)
     
    def get_value(self, wavelength):
        """ Calclate the parameter's value for an arbitrary wavelength(s) """
        return angstrom_func(np.array(wavelength), self.reference_wavelength, self.value, self.angstrom_parameters)
    
    def plot(self, wl_min = 300, wl_max = 640, new_fig = True, color='k', **plot_kwargs):
        if new_fig: plt.figure()
        plt.plot([self.reference_wavelength],[self.value],'o', markersize=10, markerfacecolor='none', lw = 0, color=color)
        wls = np.linspace(wl_min, wl_max, 100)
        plt.plot(wls, self.get_value(wls), color=color, **plot_kwargs)
        plt.xlabel('Wavelength [nm]')
        
    def set_constant(self, value):
        """ Resets the value and set wavelength dependence to zero """
        self.value = value
        self.angstrom_parameters = np.array([0])
        
    def get_properties_from_data_fit(self, data_wavelengths, data_values, 
                                     deg = 1, reference_wavelength = None, 
                                     plot=False, relative_weighting = False,
                                     inplace = True):
        """
        Fits the equation above to a set of datapoints provided by the user, 
        to infer value and angstrom_parameters (by default, this happens in place!).
        
        Parameters
        ----------
        data_wavelengths : 1D-array.
            Wavelengths of the datapoints.
        data_values : 1D-array.
            Values of the datapoints.
        deg : integer, optional, the default is 1.
            Degree of the wavelength dependences (i.e. number of angstrom prameters).
            One corresponds to classic angstrom dependence.
        reference_wavelength : scalar, optional, default is None.
            Defines the reference wavelenth.
            If None, current reference_wavelength of the SpectralParameter object will be used.
        plot : bool, optional.
            If True, a figure with datapoints and fit results will be produced for inspection. 
            The default is False.
        relative_weighting : bool, optional.
            If True, data_values during the polynomial fit are weighted by 1/data_values.
            This results in kind of relative weighting, i.e. large values allow 
            larger absolute residuals.
        inplace : bool, optional.
            If True, new value, angstrom_parameters and reference_wavelength
            will be put into the current SpectralParameter object.
            If False, method returns new value and angstrom_parameters.
            Default is True.

        """
        if reference_wavelength == None: reference_wavelength = self.reference_wavelength
        x = np.log(reference_wavelength/data_wavelengths)
        y = np.log(data_values)
        w = None if relative_weighting else data_values
        params = np.polyfit(x, y, deg,w=w)[::-1]
        ref_value = np.exp(params[0])
        angstrom_coeffs = params[1:]
        angstrom_coeffs[np.abs(angstrom_coeffs) < 1e-8] = 0.
        if plot:
            plt.figure()
            plt.plot(data_wavelengths, data_values, '.')
            wls = np.linspace(data_wavelengths.min(),data_wavelengths.max())
            plt.plot(wls, angstrom_func(wls,reference_wavelength,ref_value,angstrom_coeffs))
        if inplace:
            self.reference_wavelength = reference_wavelength
            self.value = ref_value
            self.angstrom_parameters = angstrom_coeffs.copy()
        if not inplace: return ref_value, angstrom_coeffs
        
    def shift_reference_wavelength(self, new_reference_wavelength):
        """
        Shifts the reference_wavelength adapting the value and angstrom parameter 
        attributes such that the parameter's spectral properties are not affected.
        """
        deg = len(self.angstrom_parameters)
        x = np.linspace(new_reference_wavelength,new_reference_wavelength/2,deg*2)
        y = self.get_value(x)
        self.get_properties_from_data_fit(x, y, deg, new_reference_wavelength)
        
    def to_discrete_parameter(self, wavelengths, inplace=False):
        """
        Convertes the SpectralParameter object to an DiscreteSpectralParameter object.
        Datapoints in the DiscreteSpectralParameter object will be calculated for the provided 
        wavelengths, based on the SpectralParameter object's properties.
        Not in place!
        """
        if inplace:
            self = DiscreteSpectralParameter(wavelengths, self.get_value(wavelengths))
            return None
        else:
            return DiscreteSpectralParameter(wavelengths, self.get_value(wavelengths))
    
    def get_angstrom(self, wavelength):
        """ Calculates the angstrom parameter dln(x)/dln(wl) at a given wavelength """
        ln = np.log(self.reference_wavelength/wavelength)
        angstrom = 0
        for i,a in enumerate(self.angstrom_parameters):
            angstrom += (i+1) * a * ln**i
        return angstrom
        
# sp = SpectralParameter(1, 400, [1,2,3])
# sp.plot()
# sp.shift_reference_wavelength(532)
# sp.plot()





class ParticleSizeDistribution(FrozenContainer):
    
    """
    Describes the particle size distribution of MicrophysicalAerosol objects.
    Uses a trimodal size distribution according to:
        dN/dr = Normalization * (
                 f1 * exp( -0.5 * [ (ln(r) - ln(r1)) / sig1 ]**2)
        +        f2 * exp( -0.5 * [ (ln(r) - ln(r2)) / sig2 ]**2)
        + (1-f1-f2) * exp( -0.5 * [ (ln(r) - ln(r3)) / sig3 ]**2)
    Ensure f1+f2 < 1.
    
    For large radii and width, the Mie model fails. Presumably because hard
    coded Fortran parameters have to be increased. Check in the future!
    
    Attributes
    ----------

    1.) Mie aerosol parameters:

        r1: float
            number median radius of the first aerosol mode [um]
        r2: float
            number median radius of the second aerosol mode [um]
        r3: float
            number median radius of the third aerosol mode [um]
        sig1: float
            width of the first aerosol mode.
        sig2: float
            width of the second aerosol mode.
        sig3: float
            width of the third aerosol mode.
        f1: float
            fraction of particles within the first mode
        f2: float
            fraction of particles within the second mode
        modal_fraction_mode, string:
            Determines if the modal fraction parameters f1 and f2 above, 
            refer to particle number, surface area, or particle volume.
            Allowed values are 'number', 'area', and 'volume'.
            When changing this parameter, the current fractions in the psd 
            object are automatically converted (!). So if you want to provide 
            area or volume modal fractions to the model, change 
            modal_fraction_mode before (!) setting f1 and f2.
                Example code:
                    aer = MicrophysicalAerosol()
                    aer.psd.modal_fraction_mode = 'volume'
                    aer.psd.f1, aer.psd.f2 = 0.4, 0.5
        calcrange_exponent: float, default = 3
            Defines the radius integration window used for the Mie calculations.
            Trade-off between accuracy (higher values) and preformance (lower values). 
            For each mode, lower limit rX - (sigX * calcrange_exponent)
            and upper limit rX + (sigX * calcrange_exponent) are determined. 
            Minimum and maximum of obtained values are used as integration limits in the Mie model.
            This calculation is made based on the aerosol area concentration dA/dr. 
    """
    
    def __init__(self):
        super().__init__(['r1','r2','r3',
                          'sig1','sig2','sig3',
                          'f1','f2',
                          'calcrange_exponent','_modal_fraction_mode'],
                         ['modal_fraction_mode'])
        self.f1, self.f2 = 1, 0
        self.r1, self.r2, self.r3 = 0.1, 0.5, 0.05
        self.sig1, self.sig2, self.sig3 = 0.4, 0.4, 0.1
        self.calcrange_exponent = 3
        self._modal_fraction_mode = 'number'

    @property
    def modal_fraction_mode(self):
        return self._modal_fraction_mode
            
    @modal_fraction_mode.setter
    def modal_fraction_mode(self, value):
        assert value in ['area','number','volume'], "ParticleSizeDistribution().modal_fraction_mode must be 'number', 'area' or 'volume'"
        # warnings.warn('Setting ParticelSizeDistribution().modal_fraction_mode converts modal fraction values based on current PSD parameters!')
        if value == 'number': self.f1, self.f2 = self.get_number_modal_fractions()
        elif value == 'area': self.f1, self.f2 = self.get_area_modal_fractions()
        elif value == 'volume': self.f1, self.f2 = self.get_volume_modal_fractions()
        self._modal_fraction_mode = value

    def copy(self):
        return copy.deepcopy(self)
         
    def __repr__(self):
        return 'Particle size distribution:' + \
                '\n\tModal fractions (f1, f2):\t' + str([self.f1, self.f2])[1:-1] + \
                '\n\tNumber median radii (r1, r2, r3):\t' + str([self.r1, self.r2, self.r3])[1:-1] + \
                '\n\tModal widths (sig1, sig2, sig3):\t' + str([self.sig1,self.sig2, self.sig3])[1:-1] + \
                '\n\tCalculation range exponent (calcrange_exponent):\t' + str(self.calcrange_exponent) + \
                '\n\tModal fraction refers to (modal_fraction_mode):\t' + str(self.modal_fraction_mode)

    def _number_to_area_modal_fraction(self,f1,f2):
        a1 = f1 * self.r1**2 * np.exp(2*self.sig1**2)
        a2 = f2 * self.r2**2 * np.exp(2*self.sig2**2)
        a3 = (1 - f1 - f2) * self.r3**2 * np.exp(2*self.sig3**2)
        tot = a1+a2+a3
        return a1/tot, a2/tot
    
    def _area_to_number_modal_fraction(self,f1,f2):
        n1 = f1 / self.r1**2 / np.exp(2*self.sig1**2)
        n2 = f2 / self.r2**2 / np.exp(2*self.sig2**2)
        n3 = (1 - f1 - f2) / self.r3**2 / np.exp(2*self.sig3**2)
        tot = n1+n2+n3
        return n1/tot, n2/tot

    def _number_to_volume_modal_fraction(self,f1,f2):
        a1 = f1 * self.r1**3 * np.exp(4.5*self.sig1**2)
        a2 = f2 * self.r2**3 * np.exp(4.5*self.sig2**2)
        a3 = (1 - f1 - f2) * self.r3**3 * np.exp(4.5*self.sig3**2)
        tot = a1+a2+a3
        return a1/tot, a2/tot

    def _volume_to_number_modal_fraction(self,f1,f2):
        n1 = f1 / self.r1**3 / np.exp(4.5*self.sig1**2)
        n2 = f2 / self.r2**3 / np.exp(4.5*self.sig2**2)
        n3 = (1 - f1 - f2) / self.r3**3 / np.exp(4.5*self.sig3**2)
        tot = n1+n2+n3
        return n1/tot, n2/tot

    def get_area_modal_fractions(self):
        if self.modal_fraction_mode == 'area': 
            return self.f1, self.f2
        elif self.modal_fraction_mode == 'number': 
            return self._number_to_area_modal_fraction(self.f1, self.f2)
        elif self.modal_fraction_mode == 'volume':
            fn1, fn2 = self._volume_to_number_modal_fraction(self.f1, self.f2)
            return self._number_to_area_modal_fraction(fn1, fn2)

    def get_volume_modal_fractions(self):
        if self.modal_fraction_mode == 'volume': 
            return self.f1, self.f2
        elif self.modal_fraction_mode == 'number': 
            return self._number_to_volume_modal_fraction(self.f1, self.f2)
        elif self.modal_fraction_mode == 'area':
            fn1, fn2 = self._area_to_number_modal_fraction(self.f1, self.f2)
            return self._number_to_volume_modal_fraction(fn1, fn2)
    
    def get_number_modal_fractions(self):
        if self.modal_fraction_mode == 'number': return self.f1, self.f2
        elif self.modal_fraction_mode == 'area':
            return self._area_to_number_modal_fraction(self.f1, self.f2)
        elif self.modal_fraction_mode == 'volume':
            return self._volume_to_number_modal_fraction(self.f1, self.f2)


    def get_rmin(self, calcrange_exponent = None):
        """
        Provides Mie model lower integration limits for each mode of the
        current size distribution (see calcrange_parameter - docstring above).
        """
        if calcrange_exponent is None: calcrange_exponent = self.calcrange_exponent
        sig = np.array([self.sig1, self.sig2, self.sig3])
        ra = np.array([self.r1, self.r2, self.r3]) * np.exp(sig**2)
        lims = (ra / np.exp(sig)**calcrange_exponent).round(4)
        lims[lims < 0.0001] = 0.0001 #avoid very small numbers
        return list(lims)


    def get_rmax(self, calcrange_exponent = None):
        """
        Provides Mie model upper integration limits for each mode of the
        current size distribution (see calcrange_parameter - docstring above).
        """
        if calcrange_exponent is None: calcrange_exponent = self.calcrange_exponent
        sig = np.array([self.sig1, self.sig2, self.sig3])
        ra = np.array([self.r1, self.r2, self.r3]) * np.exp(sig**2)
        return list( (ra * np.exp(sig)**calcrange_exponent).round(4) )


    def get_size_distribution(self,
                              N0 = 1.0,
                              ordinate = 'number',
                              dlnr = False,
                              rescale_modes = [1,1,1],
                              calcrange_exponent = None,
                              renormalise = False,
                              r_resolution = 1000):
        """
        Calculates the PSD, for a high resolution radius grid.
        Each Mode is calculated within the limits defined by calcrange_exponent.

        Parameters:
        -----------
        N0 : float, default = 1.0
            total number of particles
        rescale_coarse_mode : float, default = 1.0
            Scaling factor applied on the coarse mode. Useful for visibility in
            during plotting.
        ordinate : String default = 'number'
            Whether to plot aerosol number, surface area or volume distribution:
            'number' --> dN/dr,
            'area' --> dA/dr,
            'volume' --> dV/dr
            Note that this changes amplitude and position of the modes.
        dlnr: bool, default = False
            if True plot dX/dlnr instead of dX/dr
        calcrange_exponent : float, default = None
            Provide a number here to use another calcrange_exponent value than
            AerosolType().calcrange_exponent
        renormalise : bool, default = False
            If False: Any size distribution (number, area and volume) will be normalised to N0 particles.
            If True: Size distribution integral will be normalised to unity: N --> 1 particle/cm^3, A --> 1 um^2/cm^3, V --> 1 um^3/cm^3
        r_resolution : int, default = 1000
            Number of calculated points for each mode (logarithmic space).

        Returns:
        --------
        Pandas Dataframe
            Index: radius r
            Columns: the three modes
        """

        rmins = self.get_rmin(calcrange_exponent)
        rmaxs = self.get_rmax(calcrange_exponent)
        r = np.logspace(np.log10(np.min(rmins)),np.log10(np.max(rmaxs)), r_resolution)
        match ordinate: 
            case 'number': weight = np.ones(r_resolution)
            case 'area': weight = 4*np.pi*r**2
            case 'volume': weight = 4/3*np.pi*r**3
            case _: Exception("Ordinate must be 'number', 'area' or 'volume'.")
        if dlnr: weight *= r
        buff = np.zeros((len(r),3)) * np.nan
        f1,f2 = self.get_number_modal_fractions()
        if f1:
            mask = (r >= rmins[0]) & (r <= rmaxs[0])
            vals = psd_monomodal(r[mask], f1 * N0, self.r1, self.sig1) * weight[mask] * rescale_modes[0]
            buff[mask,0] = vals
        if f2:
            mask = (r >= rmins[1]) & (r <= rmaxs[1])
            vals = psd_monomodal(r[mask], f2 * N0, self.r2, self.sig2) * weight[mask] * rescale_modes[1]
            buff[mask,1] = vals
        f3 = 1.0 - self.f1 - self.f2
        if f3:
            mask = (r >= rmins[2]) & (r <= rmaxs[2])
            vals = psd_monomodal(r[mask], f3 * N0, self.r3, self.sig3) * weight[mask]  * rescale_modes[2]
            buff[mask,2] = vals
        #Renormalisation:
        if renormalise:
            total_sd = np.nan_to_num(buff).sum(axis=1)
            drs = np.diff(r)
            if dlnr: drs /= centers(r)
            sd_vals = centers(total_sd)
            norm = (sd_vals*drs).sum()
            buff /= norm
        return pd.DataFrame(buff, index = r, columns=['mode1', 'mode2', 'mode3'])

    
    def plot(self,
             N0 = 1.,
             ordinate = 'number',
             xscale = 'log',
             dlnr = False,
             rescale_modes = [1,1,1],
             color = 'k',
             calcrange_exponent = None,
             ax = None,
             legend = True,
             tight_layout = True,
             show_single_modes = True,
             renormalise=False,
             r_resolution = 1000,
             **total_dist_plot_kwargs):
        """
        Plot the aerosol size distribution

        Parameters:
        -----------
        xscale : str, default = 'log'
            X-scale of the plot, either 'linear' or 'log'
        N0 : float, default = 1.0
            total number of particles
        rescale_modes : list of numbers, default = [1,1,1]
            Scaling factors applied to each mode. Useful for better visibility.
        ordinate : String default = 'number'
            Whether to plot aerosol number, surface area or volume distribution:
            'number' --> dN/dr,
            'area' --> dA/dr,
            'volume' --> dV/dr
            Note that this changes amplitude and position of the modes.
        color: str, default = 'k'
            Line color. Color convention of Matplotlib applies.
        calcrange_exponent : float, default = None
            Provide a number here to use another calcrange_exponent value than
            AerosolType().calcrange_exponent
        dlnr: bool, default = False
            if True plot dX/dlnr instead of dX/dr
        calcrange_exponent : float, default = None
            Provide a number here to use another calcrange_exponent value than
            AerosolType().calcrange_exponent
        ax : Matplotlib axis object
            Axis to plot into.
        legend : bool, default = True
            Show legend or not.
        tight_layout :  bool
            Apply matplotlib.pyplot.tight_layout() or not.
        show_single_modes : bool, default = True:
            Show fine and coarse mode separately or not.
        renormalise : bool, default = False
            If False: Any size distribution (number, area and volume) will be normalised to N0 particles.
            If True: Size distribution integral will be normalised to unity: N --> 1 particle/cm^3, A --> 1 um^2/cm^3, V --> 1 um^3/cm^3
        r_resolution : int, default = 1000
            Number of calculated points for each mode (logarithmic space).
        **total_dist_plot_kwargs:
            keyword arguments passed to matplotlib.pyplot.plot() for
            the total (bimodal) size distribution.

        Returns:
        --------
        Pandas Dataframe
            Index: radius r
            Columns: 'fine' and 'coarse' mode distributions
        """

        df = self.get_size_distribution(N0, ordinate, dlnr, rescale_modes, calcrange_exponent, renormalise, r_resolution)
        r, sd = df.index, df.values
        tot = np.nan_to_num(sd).sum(axis=1)
        if ax is None: _, ax = plt.subplots(1,1)
        if 'label' not in total_dist_plot_kwargs: total_dist_plot_kwargs['label'] = 'Bimodal'
        ax.plot(r, tot, color = color, lw = 2, **total_dist_plot_kwargs)
        if show_single_modes:
            for i in [0,1,2]:
                lbl = ['1st','2nd','3rd'][i] + ' mode'
                if rescale_modes[i] != 1: lbl += ' x' + str(rescale_modes[i])
                ax.plot(r, sd[:,i], color = color, ls=':', label = lbl )
        ax.set_xscale(xscale)
        ax.set_xlabel('r [µm]')
        #ylabel:
        denom = 'dlnr' if dlnr else 'dr'
        if ordinate == 'number':
            numerator = r'dN'
            unit = r' [µm$^{-1}\,$cm$^{-3}$]'
        elif ordinate == 'area':
            numerator = r'dA'
            unit = r' [µm$^2\,$cm$^{-3}\,$µm$^{-1}$]'
        elif ordinate == 'volume':
            numerator = r'dV'
            unit = r' [µm$^3\,$cm$^{-3}\,$µm$^{-1}$]'
        ylabel =  numerator + '/' + denom + unit
        ax.set_ylabel(ylabel)
        if legend: ax.legend()
        if tight_layout: plt.tight_layout()


    def get_effective_radius(self, mode='all'):
        """
        Numerical calculation of the effective radius of the size
        distribution according to Hansen and Travis 1974.
        mode : either 'all', 1, 2 or 3 
            Choose for which mode to calculate effective radius.
        """
        df = self.get_size_distribution(calcrange_exponent=5)
        r, sd = df.index, df.values
        sd = np.nan_to_num(sd)
        if mode=='all': sd = sd.sum(axis=1)
        else: sd = sd[:,mode-1]
        rc = centers(r)
        sdc = centers(sd)
        pir2 = np.pi * rc**2
        dr = np.diff(r)
        return sum(rc*pir2*sdc*dr)/sum(pir2*sdc*dr)

    def get_average_radius(self, mode='all'):
        """
        Numerical calculation of the mean particle radius of the size
        distribution.
        mode : either 'all', 1, 2 or 3 
            Choose for which mode to calculate effective radius.
        """
        df = self.get_size_distribution(calcrange_exponent=5)
        r, sd = df.index, df.values
        sd = np.nan_to_num(sd)
        if mode=='all': sd = sd.sum(axis=1)
        else: sd = sd[:,mode-1]
        rc = centers(r)
        sdc = centers(sd)
        dr = np.diff(r)
        return sum(rc*sdc*dr)/sum(sdc*dr)

    def get_mean_geometric_xs(self, mode='all'):
        """
        Numerical calculation of the mean geometric cross-sectional area of
        the bimodal size distribution.
        mode : either 'all', 1, 2 or 3 
            Choose for which mode to calculate effective geometric cross section.
        """
        df = self.get_size_distribution(calcrange_exponent=5)
        r, sd = df.index, df.values
        sd = np.nan_to_num(sd)
        if mode=='all': sd = sd.sum(axis=1)
        else: sd = sd[:,mode-1]
        rc = centers(r)
        nd = centers(sd)
        dr = np.diff(r)
        return sum(np.pi * rc**2 * nd * dr)/sum(nd * dr)
    
    def get_mean_particle_volume(self, mode='all'):
        """
        Numerical calculation of the mean geometric cross-sectional area of
        the bimodal size distribution.
        mode : either 'all', 1, 2 or 3 
            Choose for which mode to calculate effective geometric cross section.
        """
        df = self.get_size_distribution(calcrange_exponent=5)
        r, sd = df.index, df.values
        sd = np.nan_to_num(sd)
        if mode=='all': sd = sd.sum(axis=1)
        else: sd = sd[:,mode-1]
        rc = centers(r)
        nd = centers(sd)
        dr = np.diff(r)
        return sum(4/3*np.pi*rc**3 * nd*dr)/sum(nd*dr)
    
    def get_params_from_data_fit(self, data_radii, data_dNdr, n_modes = 2, plot=False):
        """
        Get PSD parameters by fitting a multimode log-normal distribution
        to some measured data.
        
        Parameters
        ----------
        data_radii : particle radii in um
        data_dNdr : particle number concentration at each radius
        n_modes : Number of modes to fit.
        plot : if True, a plot with the fit result will be created.
        """
        assert self.modal_fraction_mode == 'number', "ParticleSizeDistribution().get_params_from_data_fit() currently only works if modal_fraction_mode is 'number'"
            
        match n_modes:
            case 1: 
                p, s = get_monomodal_params_from_psd(data_radii, data_dNdr)
                self.f1, self.f2 = 1.0, 0.0
                self.r1, self.r2, self.r3 = p[1], 0.1, 0.1
                self.sig1, self.sig2, self.sig3 = p[2], 0.1, 0.1
            case 2:
                p, s = get_bimodal_params_from_psd(data_radii, data_dNdr)
                f1 = p[0]/(p[0]+p[1])
                self.f1, self.f2 = f1, 1-f1
                self.r1, self.r2, self.r3 = p[2], p[3], 0.1
                self.sig1, self.sig2, self.sig3 = p[4], p[5], 0.1
            case 3:
                p, s = get_trimodal_params_from_psd(data_radii, data_dNdr)
                n = p[0] + p[1] + p[2]
                self.f1, self.f2 = p[0]/n, p[1]/n
                self.r1, self.r2, self.r3 = p[3], p[4], p[5]
                self.sig1, self.sig2, self.sig3 = p[6], p[7], p[8]
            case '_': raise Exception('PSD data fit failed. Number of modes must be 1,2 or 3.')
        if plot:
            self.plot(calcrange_exponent=5)
            norm = sum(centers(data_dNdr)* np.diff(data_radii))
            plt.plot(data_radii, data_dNdr/norm)


# TESTS:
# psd = ParticleSizeDistribution(f1=0.3, f2=0.3, calcrange_exponent=5)
# psd.get_rmin()
# psd.get_rmax()
# psd.plot(dlnr=True)
# psd.get_effective_geometric_xs()

        
class MicrophysicalAerosol():
    """
    First way to describe an aerosol species. See also manual section
    'Microphysical and optical aerosols'.
    Aerosol properties are parameterised by a particle size distribution 
    and the material refracive indices. Bulk properties are calculated using 
    the Fortran Mie model.
    
    Important attributes:
    ---------------------
    
    profile : ExtinctionProfile object. 
        Aerosol veritcal distribution (see docstring of ExtinctionProfile class).
    profile_wavelength : scalar. 
        Reference wavelength for the extinction profile.
    psd : ParticleSizeDistribution object.
        Aerosol PSD (see docstring of ParticleSizeDistribution class)
    real_refractive_index_reduced : either SpectralParameter or DiscreteSpectralParameter object.
        Reduced real refractive index of the aerosol material (n-1).
        Not real refractive index (n) !!!
        See docstring of SpectralParameter or DiscreteSpectralParameter class,
        and manual section 'Wavelength-dependent parameters'.
    imag_refractive_index : either SpectralParameter or DiscreteSpectralParameter object.
        Imaginary refractive index of the aerosol material.
        See docstring of SpectralParameter or DiscreteSpectralParameter class,
        and manual section 'Wavelength-dependent parameters'.
    results : Container object, read only.
        Filled with Mie model results when self.run_model() is called.
        See docstring of self.run_model() method.
    """
    
    def __init__(self):
        self._results_container_key_list = ['bulk_properties',
                                            'phasefunction_moments',
                                            'phasefunction',
                                            'bulk_cross_sections']
        self.profile = ExtinctionProfile()
        self.profile_wavelength = DEFAULT_REF_WL
        self.psd = ParticleSizeDistribution()
        self.real_refractive_index_reduced = SpectralParameter(0.53, DEFAULT_REF_WL,0)
        self.imag_refractive_index = SpectralParameter(0.007,DEFAULT_REF_WL,0)
        self.results = FrozenContainer(self._results_container_key_list)
        self.jacobians = Container()

    def copy(self):
        return copy.deepcopy(self)
     
    def run_model(self, wavelengths):
        """
        Calls the Fortran/mie_model executable with the current settings and reads the 
        results into the 'result'-attribute.
        
        Parameters
        ----------
        wavelengths : list or 1D-array.
            Wavelengths to run the Mie model for. 
        
        Results
        -------
        Once self.run_model() was called, results are accessible via self.results.XXX
        with XXX being the following attributes:
        
            bulk_properties : pandas DataFrame.
                Index is wavelengths.
                Columns are:
                    asymmetry_parameter : Asymmetry parameter of the scattering phase function.
                    extinction_multiplier : describes wavelength [nm] dependence of the aerosol extinction coefficient.
                        It is unity at the profile_wavelength.
                    single_scattering_albeo : Aerosol single scattering albedo
            
            phasefunction_moments : pandas DataFrame.
                Legendre moments of the scattering phase matrix, as required by the RT backends.
                Index is two-level: wavelength and moment order.
                Columns are different matrix entries. Labels follow the nomenclature in the Mie model manual.

            phasefunction : pandas DataFrame.
                Scattering phase function and other non-zero phase matrix entries.
                Index is two-level: wavelength and scattering angle.
                Columns represent different matrix entries.
            
            bulk_cross_sections : pandas DataFrame.
                Index is wavelength.
                Columns are effective particle cross-sections [um2] for extinction, scattering and absorption.
        """
        
        f1, f2 = self.psd.get_number_modal_fractions()
        f3 = 1.0 - f1 - f2
        
        assert f3>-1e-10, 'Modal fractions add up to more than one!'
        for p in ['r1','r2','r3','sig1','sig2','sig3']:
            assert getattr(self.psd,p) > 0, 'Invalid PSD parameter ecnountered (less equal zero)'
        
        s1, s2, s3 = np.exp(self.psd.sig1), np.exp(self.psd.sig2), np.exp(self.psd.sig3) # Convert modal widths to mie model convention

        wls = np.array(wavelengths.copy())
        
        path_mie_model = get_mie_model_path()
         
        wstring = 'Mie model error. Microphysical aerosol profile_wavelength outside simulated wavelength range! ' \
                   + 'Use aerosol\'s shift_profile_wavelength() method or the Pydort object\'s shift_reference_wavelengths() method to adjust.'
        assert (self.profile_wavelength >= wls.min()) and (self.profile_wavelength <= wls.max()), wstring     

        
        # Delete all former input and output files. Otherwise it's hard to recognize if Fortran code failed.
        
        for f in Path(path_mie_model + 'INPUT').glob("**/*.inp"): f.unlink() 
        for f in Path(path_mie_model + 'OUTPUT').glob("**/*.out"): f.unlink() 
    
        # Write configuration file:
        # -------------------------
         
        config_string = '\n'.join([
            '! This configuration file is used and modified by the python wrapper!',
            '! Changes here will not be permanent.',
            '! See python wrapper for detailed parameter description',
            '',
            ', '.join([str(self.psd.r1),str(self.psd.r2),str(self.psd.r3)]) + '  ! modal radii',
            ', '.join([str(s1),str(s2),str(s3)]) + '   ! modal widths expressed as exp(sigma)',
            ', '.join([str(f1),str(f2),str(f3)]) + '   ! modal fractions',
            ', '.join(str(float(value)) for value in self.psd.get_rmin()) + '   ! rmin',
            ', '.join(str(float(value)) for value in self.psd.get_rmax()) + '   ! rmax'])
         
        with open(path_mie_model + 'INPUT/configuration.inp','w+') as f:
            f.write(config_string)
             
        # Write wavelengths and refractive indices to file:
        #--------------------------------------------------

        dstr = '' # String with data
        rris = self.real_refractive_index_reduced.get_value(wls)
        iris = self.imag_refractive_index.get_value(wls)

        
        assert (rris>0).all(), 'Real refractive index is less equal zero.'
        #assert (iris>0).all(), 'Imaginary refractive index is less equal zero.'
        if (rris>1).any():
            wstring = 'In Pydort, the reduced real refractive index must be provided (n-1). Did you accidentally provide the real refractive index (n) instead?'
            warnings.warn(wstring)

        for iwl,wl in enumerate(wls):
             dstr += str(wl) + ', '
             dstr += str(rris[iwl]+1) + ', '
             dstr += str(iris[iwl]) + '\n'
         
        refr_string = '\n'.join([
            str(len(wls)) + ' ! number of wavelengths\n',
            'wavelength, n_real, n_imag',
            dstr])
        
        with open(path_mie_model + 'INPUT/refractive_indexes.inp','w+') as f:
             f.write(refr_string)
             
         
        # Call Mie model executable:
        #---------------------------
        run_subprocess(path_mie_model + 'MIE_MODEL.exe')
         
         
        # Gather results:
        #----------------   
        # Read bulk properties and replace extinction coefficient by extinction multiplier based on profile reference wavelength:
        bp = pd.read_csv(path_mie_model + 'OUTPUT/bulk_properties.out', index_col=0, sep='\\s+')
        ec = bp.extinction_xs.copy()
        rec = log_interp(self.profile_wavelength, ec.index, ec.values, fill_value=None)
        bp['extinction_multiplier'] = ec/rec
        bc = pd.concat([bp.extinction_xs,
                        bp.extinction_xs * bp.single_scattering_albedo,
                        bp.extinction_xs * (1 - bp.single_scattering_albedo)],
                       axis=1)
        bc.columns = ['extinction_xs','scattering_xs','absorption_xs']
        self.results.bulk_cross_sections = bc.copy()
        self.results.bulk_properties = bp.drop('extinction_xs', axis=1).sort_index(axis=1)
        self.results.phasefunction_moments = pd.read_csv(path_mie_model + 'OUTPUT/phasefunction_moments.out', index_col=[0,1], sep='\\s+')
        self.results.phasefunction = pd.read_csv(path_mie_model + 'OUTPUT/phasefunction.out', index_col=[0,1], sep='\\s+') / 4 / np.pi
     
    def shift_profile_wavelength(self, new_profile_wavelength):
        """
        Shifts the profile reference wavelength (profile_wavelength) and adjusts 
        the extinction profile such that the aerosol's properties are maintained.
        Involves a Mie model run!
        """
        self.run_model([self.profile_wavelength, new_profile_wavelength])
        ems = self.results.bulk_properties.extinction_multiplier.values
        self.profile.values *= ems[1]/ems[0] # no ratio require
        self.profile_wavelength = new_profile_wavelength   
     
    def shift_reference_wavelengths(self, new_reference_wavelength):
        """
        Shifts all reference wavelengths of the aerosol object (extinction profile and spectral parameter),
        adjusting all related parameters such that the aerosol's properties are maintained.
        """
        self.shift_profile_wavelength(new_reference_wavelength)
        self.real_refractive_index_reduced.shift_reference_wavelength(new_reference_wavelength)
        self.imag_refractive_index.shift_reference_wavelength(new_reference_wavelength)
        
    def get_jacobians(self, parameter, wavelengths, dx=None, rerun_fd_base_simulation=True):
        if rerun_fd_base_simulation: self.run_model(wavelengths)
        aerd = self.copy()
        rdx = 1e-5
        match parameter:
            case 'f1':
                if dx is None: dx = - aerd.psd.f1 * rdx
                assert self.psd.f1 + self.psd.f2 == 1, "Mie model jacobians for modal fraction only available if f1 + f2 are one." 
                aerd.psd.f1 += dx 
                aerd.psd.f2 -= dx
            case 'r1':
                if dx is None: dx = aerd.psd.r1 * rdx
                aerd.psd.r1 += dx 
            case 'r2':
                if dx is None: dx = aerd.psd.r2 * rdx
                aerd.psd.r2 += dx 
            case 'sig1':
                if dx is None: dx = aerd.psd.sig1 * rdx
                aerd.psd.sig1 += dx 
            case 'sig2':
                if dx is None: dx = aerd.psd.sig2 * rdx 
                aerd.psd.sig2 += dx 
            case 'n0':
                if dx is None: dx = aerd.real_refractive_index_reduced.value * rdx
                aerd.real_refractive_index_reduced.value += dx 
            case 'k0':
                if dx is None: dx = aerd.imag_refractive_index.value * rdx
                aerd.imag_refractive_index.value += dx 
            case _:
                raise Exception('Parameter ' + parameter + ' not available or not implemented.')
        aerd.run_model(wavelengths)
        setattr(self.jacobians,parameter,FrozenContainer(self._results_container_key_list))
        for k in self._results_container_key_list:
            y0 = getattr(self.results,k)
            yd = getattr(aerd.results,k)
            dydx = (yd-y0)/dx
            dydx = dydx.reindex(y0.index).fillna(0)
            setattr(getattr(self.jacobians, parameter), k, dydx)
            
# a = MicrophysicalAerosol()
# a.profile.fill_exponential_profile(2000,0.5)
# a.profile.plot()
# a.imag_refractive_index.plot()
# a.shift_reference_wavelengths(400)
# a.profile.plot()
# a.imag_refractive_index.plot()
# a.run_model()
# a.bulk_properties


class OpticalAerosol():
    """
    Second way to describe an aerosol species. See also manual section
    'Microphysical and optical aerosols'.
    Aerosol properties are parameterised by an extinction wavelength dependence,
    single scattering absorption (= 1 - SSA) and asymmetry parameter 
    (assuming a Henyey-Greenstein scattering phase function).
    
    Important attributes:
    ---------------------
    
    profile : ExtinctionProfile object. 
        Aerosol veritcal distribution (see docstring of ExtinctionProfile class).
    profile_wavelength : scalar. 
        Reference wavelength for the extinction profile. In OpticalAerosols,
        this is linked to reference_wavelength of the extinction_multiplier attribute.
    extinction_multiplier : Either SpectralParameter or DiscreteSpectralParameter object.
        Describes the aerosol's extinction wavelength dependence.  
        See docstring of SpectralParameter or DiscreteSpectralParameter class,
        and manual section 'Wavelength-dependent parameters'.
        Be careful:
            extinction_multiplier.reference_wavelength is linked to profile_wavelength.
            extinction_multiplier.value should be unity to not accidently scale the aerosol extinction.
    single_scattering_absortion : either SpectralParameter or DiscreteSpectralParameter object.
        (1 - single scattering albedo) of the aerosol.
        See docstring of SpectralParameter or DiscreteSpectralParameter class,
        and manual section 'Wavelength-dependent parameters'.
    asymmetry_parameter : either SpectralParameter or DiscreteSpectralParameter object.
        Asymmetry parameter of the assumed Henyey-Greenstein scattering phase function.
        See docstring of SpectralParameter or DiscreteSpectralParameter class,
        and manual section 'Wavelength-dependent parameters'.
    results : Container object, read only.
        Filled with model results when self.run_model() is called.
        See docstring of self.run_model() method.
    """
    def __init__(self): 
        self._results_container_key_list = ['bulk_properties',
                                            'phasefunction_moments',
                                            'phasefunction']
        self.profile = ExtinctionProfile()
        self.extinction_multiplier = SpectralParameter(1.0, DEFAULT_REF_WL, 1.7)
        self.single_scattering_absorption = SpectralParameter(0.05, DEFAULT_REF_WL, -0.5)
        self.asymmetry_parameter = SpectralParameter(0.6, DEFAULT_REF_WL, 0.4)
        self.hg_moment_cutoff = 1e-8
        self.results = FrozenContainer(self._results_container_key_list)
        self.jacobians = Container()
        self._profile_wavelength = DEFAULT_REF_WL        

    def copy(self):
        return copy.deepcopy(self)
 
    # In contrast to Microphysical aerosol, profile_wavelength is a property here, since defined by extinction_multiplier.
    @property
    def profile_wavelength(self):
        if isinstance(self.extinction_multiplier, SpectralParameter): # if extinction is a SpectralParameter, extinction_multiplier.reference_wavelength and profile_wavelength must be linked
            return self.extinction_multiplier.reference_wavelength
        else: # DiscreteSpectralParameter has no reference wavelength. No linking required.
            return self._profile_wavelength
        
    @profile_wavelength.setter
    def profile_wavelength(self, value):
        if isinstance(self.extinction_multiplier, SpectralParameter):
            wstring = 'For optical aerosol, profile_wavelength and extinction_multiplier.reference_wavelength are linked! ' \
                        + 'Setting one of them sets the other as well.'
            warnings.warn(wstring)
            self.extinction_multiplier.reference_wavelength = value
        else:
            self._profile_wavelength = value
            
    def run_model(self, wavelengths):
        """
        Equivalent to the run_model() method of the MicrophysicalAerosol class.
        Here, the method does not really run a model but fills the  
        but creates results of the same format. 
        
        Parameters
        ----------
        wavelengths : list or 1D-array.
            Wavelengths to run the Mie model for. 
        
        Results
        -------
        Once self.run_model() was called, results are accessible via self.results.XXX
        with XXX being the following attributes:
        
            bulk_properties : pandas DataFrame.
                Index is wavelengths.
                Columns are:
                    asymmetry_parameter : Asymmetry parameter of the scattering phase function.
                    extinction_multiplier : describes wavelength [nm] dependence of the aerosol extinction coefficient.
                        It is unity at the profile_wavelength.
                    single_scattering_albedo : Aerosol single scattering albedo
            
            phasefunction_moments : pandas DataFrame.
                Legendre moments of the scattering phase matrix, as required by the RT backends.
                Index is two-level: wavelength and moment order.
                Columns are different matrix entries. Labels follow the nomenclature in the Mie model manual.

            phasefunction : pandas Series.
                Scattering phase function and other non-zero phase matrix entries.
                Index is two-level: wavelength and scattering angle.
        """
        idx = np.array(wavelengths.copy())

        # Aerosol bulk properties:
        em = np.atleast_1d(self.extinction_multiplier.get_value(idx))
        ssa = np.atleast_1d(1-self.single_scattering_absorption.get_value(idx))
        ap = np.atleast_1d(self.asymmetry_parameter.get_value(idx))
        self.results.bulk_properties = pd.DataFrame(index = idx,
                                            data={'extinction_multiplier':em,
                                                  'single_scattering_albedo':ssa,
                                                  'asymmetry_parameter':ap}).sort_index(axis=1)

        # Scattering phasefunction moments:
        dl = []
        for iwl, wl in enumerate(idx):
            moms = get_henyey_greenstein_moments(ap[iwl], np.arange(MAX_PF_MOMENTS))
            moms = moms[moms > self.hg_moment_cutoff]
            dl += [pd.Series(moms)]
        d = pd.concat(dl,keys=idx)
        d.index.set_names(['wavelength','moment_order'], inplace=True)
        d = d.rename('alpha1').to_frame()
        for c in ['alpha2', 'alpha3', 'alpha4', 'beta1', 'beta2']:
            d[c] = 0.0
        self.results.phasefunction_moments = d.copy()
        
        # Scattering phasefunction:
        dl = []
        a = np.arange(0,180.1,0.2)
        for iwl, wl in enumerate(idx):
            pf = henyey_greenstein_pf(a, self.asymmetry_parameter.get_value(wl))
            dl += [pd.Series(pf, index = a)]
        d = pd.concat(dl,keys=idx)
        d.index.set_names(['wavelength','scattering_angle'], inplace=True)
        self.results.phasefunction = d.copy()

        
    def shift_profile_wavelength(self, new_profile_wavelength):
        """
        Shifts the profile reference wavelength (extinction_multiplier.reference_wavelength) 
        and adjusts both, extinction profile and extinction multiplier such that
        properties are maintained and extinction_multiplier.value stays unity.
        """
        ems = self.extinction_multiplier.get_value([self.extinction_multiplier.reference_wavelength, new_profile_wavelength])        
        ratio = ems[1]/ems[0]
        self.extinction_multiplier.shift_reference_wavelength(new_profile_wavelength) # Shift wavelength without changing the wavelength dependence.
        self.extinction_multiplier.value /= ratio # reduce new value to one
        self.profile.values *= ratio # Scale profile to compensate for change in extinction multiplier
        
    def shift_reference_wavelengths(self, new_reference_wavelength):
        """
        Shifts all reference wavelengths of the aerosol object (extinction profile and spectral parameter),
        adjusting all related parameters such that the aerosol's properties are maintained.
        """
        self.shift_profile_wavelength(new_reference_wavelength)
        if isinstance(self.single_scattering_absorption, SpectralParameter):
            self.single_scattering_absorption.shift_reference_wavelength(new_reference_wavelength)
        if isinstance(self.asymmetry_parameter, SpectralParameter):
            self.asymmetry_parameter.shift_reference_wavelength(new_reference_wavelength)

    def get_jacobians(self, parameter, wavelengths, dx=None, rerun_fd_base_simulation=True):
        if rerun_fd_base_simulation: self.run_model(wavelengths)
        aerd = self.copy()
        rdx = 1e-5
        match parameter:
            case 'ang':
                if dx is None: dx = aerd.extinction_multiplier.angstrom_parameters[0] * rdx
                aerd.extinction_multiplier.angstrom_parameters[0] += dx
            case 'ssca0':
                if dx is None: dx = aerd.single_scattering_absorption.value * rdx
                aerd.single_scattering_absorption.value += dx 
            case 'asym0':
                if dx is None: dx = aerd.asymmetry_parameter.value * rdx
                aerd.asymmetry_parameter.value += dx 
            case _:
                raise Exception('Parameter ' + parameter + ' not available or not implemented.')
        aerd.run_model(wavelengths)
        setattr(self.jacobians,parameter,FrozenContainer(self._results_container_key_list))
        for k in self._results_container_key_list:
            y0 = getattr(self.results,k)
            yd = getattr(aerd.results,k)
            dydx = (yd-y0)/dx
            dydx = dydx.reindex(y0.index).fillna(0)
            setattr(getattr(self.jacobians, parameter), k, dydx)

class LiteratureSpectrum():
    """
    Stores literature spectrum information, e.g. solar irradiance spectrum or trace gas absorption cross-sections.
    Facilitates:
        Import of the spetral data from textfile.
        Gaussian smoothing of high resolution spectra (see paper section 'Spectral averaging and interpolation').
        Temperature dependent cross-sections.
    Object initialisation may take a moment, since a high resolution literature 
    spectrum is read from file and immediately smoothed using fast-fourier-transformation.
        
    Important attributes:
    ---------------------
    filename : String. Required during initialisation. E.g. 'data.dat'
        Name of the file containing the raw data.
        File must be located in physics_data/literature_spectra/
        The assumed format is:
            Comments/header lines start with '#'
            First column is wavelength in nm.
            Second column is values (cross section [cm2] or irradiance [photons/s/cm2/nm])
            Further columns (optional) contain temperature dependence coefficients (cross sections only)
            in units of cm2/K, cm2/K2, cm2/K3, ...
        See already existing files, e.g. physics_data/literature_spectra/o3.dat
    unit_string : String 
        The units of the data (e.g. 'cm2'). Only used for labeling during plotting.
    raw_data : 2D-array. Do not modify after initialisation!
        Contains the raw data as array, in the same format as in the file (see filename attribute above).
    smoothing_fwhm : scalar. Default is 1.0 
        Full-width half-maximum of Gaussian smoothing kernel [nm]. 
        Also defines the sampling resolution of the smoothed data, which is smoothing_fwhm/5.
        Is a Python property: everytime smoothing_fwhm is changed, smoothed_data will be recalculated immediately. 
    smoothed_data : 2D-array.
        Contains the smoothed data as array, in the same format as in the file (see filename attribute above).
    reference_temperature : scalar.
        Only of relevance for absorption cross-sections with temperature dependence (e.g. O3).
        Reference temperature for the temperature dependence coefficients.
    """
    
    
    def __init__(self,
                 filename,
                 unit_string = '',
                 smoothing_fwhm = DEFAULT_SPECTRAL_RESOLUTION,
                 spectral_smoothing_on_initialisation = True,
                 reference_temperature = 273.15):
        self.raw_data = np.genfromtxt(PATH_PHYSICS + 'literature_spectra/' +  filename)
        self._smoothing_fwhm = smoothing_fwhm
        self.smoothed_data = None
        if spectral_smoothing_on_initialisation: self.apply_smoothing()
        self.unit_string = unit_string
        self.reference_temperature = reference_temperature

    def copy(self):
        return copy.deepcopy(self)
            
    def apply_smoothing(self):
        """ Recalculates smoothed_data, based on current raw_data and smoothed_fwhm."""
        xs = self.raw_data
        sigma = self.smoothing_fwhm/2.355
        smoothing_resolution = self.smoothing_fwhm/SMOOTHING_FWHM_SAMPLING
        ncols = xs.shape[1] # some xs files also contain temperature coeffcients -> smooth all
        # Create equidistant wavelength grid
        wls_hr = np.arange(xs[:,0].min(),xs[:,0].max(),WAVELENGTH_HIGH_RESOLUTION)
        # Create kernel. Apply rounding to ensure gauss curve is centered. 
        region = int(5*self.smoothing_fwhm/WAVELENGTH_HIGH_RESOLUTION) * WAVELENGTH_HIGH_RESOLUTION
        kernel_wls = np.arange(-region,region+WAVELENGTH_HIGH_RESOLUTION/2,WAVELENGTH_HIGH_RESOLUTION)
        kernel = np.exp(-kernel_wls**2/sigma**2/2)
        kernel /= kernel.sum()
        # iterate over columns:
        cl = []
        for c in range(1,ncols):
            interp_kind = 'cubic' if c == 1 else 'linear' # some columns feature sudden jumps -> cubic spline introduces oscillation.
            fint = interp1d(xs[:,0],xs[:,c],kind=interp_kind,bounds_error=False,fill_value=0)
            col_hr = fint(wls_hr)
            col_conv = fftconvolve(col_hr, kernel,mode='valid')
            wls_conv = wls_hr[:len(col_conv)] + region
            # Interpolate to smoothed resolution:
            wls_sres = np.arange(wls_conv.min(),wls_conv.max(), smoothing_resolution)
            fint = interp1d(wls_conv, col_conv, kind=interp_kind, bounds_error=False, fill_value=0)
            col_sres = fint(wls_sres)
            if c == 1: col_sres[col_sres<0] = 0 # remove negative values in cross sections.
            cl += [col_sres]
        sd = np.vstack([wls_sres, *cl])
        self.smoothed_data = sd.T.copy()
        
    def get_values(self, wavelengths = None, temperatures = None, use_raw_data = False, as_dataframe = False):
        """
        Calculates cross-section values for specific wavelengths and optionally for specific temperatures.
        
        Parameters
        ----------
        wavelengths : 1D-array or None.
            Wavelengths to provide xs values for. 
            If None, use wavelengths from the literature data. The default is None.
        temperatures : float, 1D-array or None.
            Only of relevance for absorption cross sections with temperature dependence.
            Temperature(s) to provide cross section values for.
            If float or 1D-array, columns of returned object represent the different temperatures
            If None, columns will be the same as in the raw data, i.e. cross section values and temperature coefficients.
        use_raw_data : bool, optional
            If True, use raw_data instead of smoothed_data. The default is False.
        as_dataframe : bool.
            If True, data is returned as labeled dataframe instead of array.

        Returns
        -------
        1D-array, 2D-array or Pandas dataframe. Depends on temperatures and as_dataframe (see above). 
        """
        xs = self.raw_data.copy() if use_raw_data else self.smoothed_data.copy()
        if wavelengths is None: wls = xs[:,0]
        else:
            wls = np.atleast_1d(wavelengths)
            xs = self._interpolate_wavelengths(xs, wls)
        if temperatures is None:
            xs = xs[:,1:]
            cols = ['xs'] + ['temp_c' + str(i) for i in range(self.n_temperature_coeffs)]
            col_name = None
        if not temperatures is None: # If temperatures are provided, process data to cross-sections for specific temperatures.
            xs = xs.T # Transpose only for easier indexing.
            t = np.atleast_1d(temperatures)
            t = t - self.reference_temperature
            xsn = np.repeat(xs[1][:,np.newaxis],len(t),axis=1)
            for c in range(2,len(xs)): # iterate over coefficient columns
                xsn += xs[c][:,np.newaxis] * t**(c-1)
            xs = xsn
            cols = t+self.reference_temperature
            col_name = 'temperature'
        if not as_dataframe: 
            if xs.shape[1] == 1: xs = xs.flatten()
            return xs
        else:
            df = pd.DataFrame(xs, index=wls, columns=cols)
            df.columns.name = col_name
            df.index.name = 'wavelength'
            return df
        
    def plot(self, new_fig = True, color='k', alpha = 1.0, lw = 3, show_unsmoothed=True, **plot_kwargs):
        if new_fig: plt.figure()
        if show_unsmoothed:
            x,y = self.raw_data.T[0:2]
            plt.plot(x, y, lw = 1, alpha = alpha*0.5, color=color, **plot_kwargs)
        x,y = self.smoothed_data.T[0:2]
        plt.plot(x, y, lw = lw, color=color, alpha=alpha, **plot_kwargs)
        plt.grid()
        plt.xlabel('Wavelength [nm]')
        plt.ylabel(self.unit_string)
      
    @ property
    def n_temperature_coeffs(self): 
        return self.raw_data.shape[1]-2
      
    @ property
    def smoothing_fwhm(self):
        return self._smoothing_fwhm
    
    @ smoothing_fwhm.setter
    def smoothing_fwhm(self, value):
        self._smoothing_fwhm = value
        self.apply_smoothing()
        
    def _interpolate_wavelengths(self, xs_data, new_wavelengths):
        """
        Helper function to interpolate cross-section data to new wavelength grid.
        
        Parameters
        ----------
        xs_data : 2D array, columns are:
            1.  Wavelength
            2.  Cross section
            3+. Temperature or pressure dependence coefficients.
        
        Returns
        -------
        2D-array. Wavelength interpolated data, same column format as xs_data (wavelength column included).
        """
        xs = xs_data
        ncols = xs.shape[1] # some xs files also contain temperature coeffcients -> smooth all
        # iterate over columns:
        cl = []
        for c in range(1,ncols):
            interp_kind = 'cubic' if c == 1 else 'linear' # some columns feature sudden jumps -> cubic spline introduces oscillation.
            fint = interp1d(xs[:,0], xs[:,c], kind=interp_kind, bounds_error=False, fill_value=0)
            nc = fint(new_wavelengths)
            if c == 1: nc[nc<0] = 0 # remove negative values in cross sections.
            cl += [nc]
        return np.vstack([new_wavelengths, *cl]).T



class TraceGas():
    """
    Contains all information required to consider a trace gas species in the simulation.
    
    Important attributes:
    ---------------------
    name : String. 
        Name of the trace gas. During initialisation, profile and cross section 
            information will be read from <name>.dat files in:
            physics_data/profiles/
            physics_data/literature_spectra/
        Make sure to provide these files!
    xs : LiteratureSpectrum object.
        Gas absorption cross section information
    profile : GasProfile objet.
        Vertical gas concentration profile.
    """
  
    def __init__(self, name, spectral_smoothing_on_initialisation=True):
        self.xs = LiteratureSpectrum(name + '.dat','cm2',
                                     spectral_smoothing_on_initialisation=spectral_smoothing_on_initialisation) 
        self.profile = GasProfile()
        d = np.loadtxt(PATH_PHYSICS + 'profiles/' + name + '.dat', comments = '#')
        self.profile.interpolate_from(d[:,0],d[:,1])

    def copy(self):
        return copy.deepcopy(self)       

class TraceGasContainer():
    """
    Container for multiple TraceGas objects. By default, it is initialized containing
    TraceGas objects for O4, O3, NO2 and H2O.
    
    Important Attributes:
    ---------------------
    XXX: TraceGas objects for various gases. Can be removed or added.
    water_vapor_label: Identifier for water vapour. 
        Number concentration profile of this trace gas is used 
        for the (minor) water vapour correction of Rayleigh cross section and 
        scattering matrix.
    o4_label: Identifier for O4. The trace gas with this label is assumed to be O4. 
        Its profile will be recalculated from pressure and temperature profile before RT simulation.
    """    

    def __init__(self, 
                 names = ['o4','o3','no2','h2o','hcho','bro','io','chocho','hono'],
                 spectral_smoothing_on_initialisation = True):
        for g in names: setattr(self, g, TraceGas(g, spectral_smoothing_on_initialisation=spectral_smoothing_on_initialisation))
        self.water_vapor_label = 'h2o'
        self.o4_label = 'o4'

    def copy(self):
        return copy.deepcopy(self)
         
    def __repr__(self):
        s = 'Trace gas container with gases:\n'
        s += '\t\n'.join(self.get_name_list())
        return s
    
    def get_name_list(self):
        l = []
        for a in dir(self):
            att = getattr(self,a)
            if isinstance(att,TraceGas): l += [a]
        return l
        
    def to_dictionary(self):
        return {g:getattr(self, g) for g in self.get_name_list()}
        
    def remove_gas(self, name):
        delattr(self, name)
        
    def get_n_species(self):
        return len(self.get_name_list())
    
    
class AerosolContainer():
    """
    Container for multiple Aerosol species. By default, it is initialized containing
    two OpticalAerosol objects (typical tropospheric and stratospheric background aerosol),
    and two MicrophysicalAerosol objects (typical biomass burning black and brown carbon).
    
    Important Attributes:
    ---------------------
    XXX: Aerosol objects for various species. Can be removed or added.
    """    

    def __init__(self):        

        # Create typical plume aerosols:
        bc = MicrophysicalAerosol()
        bc.profile.fill_gaussian_profile(3000,1000,0.4)
        bc.psd.f1, bc.psd.f2 = 1., 0.
        bc.psd.r1 = 0.104 
        bc.psd.sig1 = 0.4
        bc.real_refractive_index_reduced.set_constant(0.7)
        bc.imag_refractive_index.set_constant(0.45)
        self.black_carbon = bc        

        brc = MicrophysicalAerosol()
        brc.profile.fill_gaussian_profile(3000,1000,2.6)
        brc.psd.f1, brc.psd.f2 = 1., 0.
        brc.psd.r1 = 0.104 
        brc.psd.sig1 = 0.4
        brc.real_refractive_index_reduced.set_constant(0.53)
        brc.imag_refractive_index.reference_wavelength = 532
        brc.imag_refractive_index.value = 1.613e-2
        brc.imag_refractive_index.angstrom_parameters = [5.25, -2.78]
        self.brown_carbon = brc
        
        # Create typical tropospheric aerosol:
        tropo = OpticalAerosol()
        tropo.profile.fill_fermi_profile(1000, 100, 0.1)
        tropo.asymmetry_parameter.reference_wavelength = 400
        tropo.asymmetry_parameter.value = 0.6
        tropo.asymmetry_parameter.angstrom_parameters = [0.4]
        tropo.single_scattering_absorption.reference_wavelength = 400
        tropo.single_scattering_absorption.value = 0.04
        tropo.single_scattering_absorption.angstrom_parameters = [0.5]
        tropo.extinction_multiplier.angstrom_parameters = [1.7]
        self.tropospheric_background = tropo
        
        # Create typical stratospheric aerosol:
        strato = OpticalAerosol()
        d = np.loadtxt(PATH_PHYSICS + 'profiles/stratospheric_aerosol.dat', comments = '#')
        strato.profile.interpolate_from(d[:,0],d[:,1])
        strato.asymmetry_parameter.reference_wavelength = 400
        strato.asymmetry_parameter.value = 0.6
        strato.asymmetry_parameter.angstrom_parameters = [0.4]
        strato.single_scattering_absorption.reference_wavelength = 400
        strato.single_scattering_absorption.set_constant(1e-5)
        strato.extinction_multiplier.angstrom_parameters = [1.5]
        self.stratospheric_background = strato

    def copy(self):
        return copy.deepcopy(self)
         
    def __repr__(self):
        s = 'Aerosol container with aerosol types:'
        adic = self.to_dictionary()
        for a in adic:
            at = adic[a]
            if isinstance(at, MicrophysicalAerosol): t = 'MicrophysicalAerosol' # Extract model type
            if isinstance(at, OpticalAerosol): t = 'OpticalAerosol'
            s += '\t\n' + a + ' (' + t +  ')'
        return s
        
    def get_name_list(self):
        l = []
        for a in dir(self):
            att = getattr(self,a)
            if isinstance(att, MicrophysicalAerosol) or isinstance(att, OpticalAerosol): l += [a]
        return l
            
    def to_dictionary(self):
        return {g:getattr(self, g) for g in self.get_name_list()}           
        
    def remove_aerosol(self, name):
        delattr(self, name)
        
    def remove_all(self):
        """ Removes all aerosol species """
        for a in self.get_name_list():
            self.remove_aerosol(a)
        
    def run_models(self, wavelengths):
        """ Runs Mie model and optical models, respectively, for all aerosol species """
        for a in self.get_name_list():
            aer = getattr(self, a)
            aer.run_model(wavelengths)
            
    def get_n_species(self):
        return len(self.get_name_list())        



class SurfaceAlbedo(DiscreteSpectralParameter):
    """
    Describes the spectrally dependent surface albedo (lambertian equivalent reflectivity).
    By default initialized from file, with TROPOMI-data for shady fire forest.
    But can easily be set for any location and month using get_tropomi_albedo().
    """
    
    def __init__(self):
        super(SurfaceAlbedo,self).__init__()
        self._get_data_from_default_file('surface_albedo/tropomi_shady_fire_forest_july.dat')

    def change_wavelengths(self, new_wavelengths, inplace = True):
        """ Overloads change_wavelengths from DiscreteSpectralParameter to change from log_interp to standard interp """
        ynew = interp(new_wavelengths, self.wavelengths, self.values)
        if inplace:
            self.wavelengths = new_wavelengths.copy()
            self.values = ynew
        else:
            return ynew

    def copy(self):
        return copy.deepcopy(self)
     
    def get_tropomi_albedo(self, lat, lon, month, inplace=True):
        """
        Reads and linearly interpolates TROPOMI-2 Lambertian equivalent reflectance (LER).
        Data gets interpolated onto the current wavelengths in self.wavelengths.
        Not working close to Earths poles and lon = 180.
        Surface albedos outside valid data range raises error.
        Values below 328 nm will be extrapolated (see code below for details).
        ----------
        lat : float, Latitude, -89 < lat < 89
        lon : float, Longitude, -179 < lon < 179
        month : int, Month of the year (starting at 1 = January), 1 <= month <= 12
        """
        wavelengths = self.wavelengths
        file = get_latest_tropomi_surface_albedo_file()
        d = h5py.File(file,'r')
        lats = np.array(d.get('latitude'))
        lons = np.array(d.get('longitude'))
        wls = np.array(d.get('wavelength'))
        # Prepare slices for efficient data reading:
        sl_lat = slice(np.argmax(lats > lat-1), np.argmax(lats > lat+1))
        sl_lon = slice(np.argmax(lons > lon-1), np.argmax(lons > lon+1))
        lats = lats[sl_lat]
        lons = lons[sl_lon]
        albs = np.array(d['minimum_LER_clear'][month-1,:,sl_lon,sl_lat]) # dims: month, wavelength, lon, lat
        # 2D-interpolation to exact location for all wavelengths:
        newx = [(wl,lon,lat) for wl in wls]
        albs = interpn((wls, lons, lats), albs, newx)
        wls_reg = np.arange(150,800,2)
        albs_reg = np.interp(wls_reg, wls, albs)
        albs_reg = smooth_gaussian(albs_reg,10)
        albs = np.interp(wavelengths, wls_reg, albs_reg)
        if inplace: self.values = albs
        else: return pd.Series(albs, index = self.wavelengths.copy())
        
class ResultsRT():
    """
    Contains and post-processes the results from the VLIDORT-QS RT simulation.
    The raw results (as read from the VLIDORT-QS Fortran executable output) are 
    contained in the _raw attribute. However, access to the results is typically
    done via the methods provided by this class (see below), that also post-process 
    the data.
    """
    
    def __init__(self,
                 parent_vpc):
        self._parent_vpc = parent_vpc
        self._raw = Container()        

    def _read_raw_results(self, directory):
        # Read raw RT results:
        p = directory
        for f in os.listdir(p):
            attr_name = f.split('.')[0]
            if not bool(attr_name): continue # Ignores .gitkeep and things.
            setattr(self._raw, attr_name, read_array_from_file(p + f))

    def copy(self):
        return copy.deepcopy(self)      

    def get_actinic_flux(self,
                         up_down_both = 'both', 
                         diffuse_direct_both = 'both',
                         normalized = False):
        """
        Returns the actinic flux [photons/s/cm2/nm]
        at the observer altitude, for the simulated wavelengths.

        Parameters
        ----------
        up_down_both : String, optional. The default is 'both'.
            Either:
                'up' : consider upwelling photons (coming from lower halfspace) only.
                'down' : consider downwelling photons (coming from upper halfspace) only.
                'both' : consider all directions.
        diffuse_direct_both : String, optional. The default is 'both'.
            Either :                 
                'direct' : consider photons from direct solar beam (unscattered) only.
                'diffuse' : consider scattered photons only.
                'both' : consider all photons. 
        normalized : bool, optional. The default is False.
            If False, actinic fluxes are returned in units of [photons/s/cm2/nm]
            If True, actinic fluxes are normalized to the TOA solar irradiance.

        Returns
        -------
        Pandas series with actinic fluxes
        """
        du = self._raw.flux_diffuse_upwelling
        dd = self._raw.flux_diffuse_downwelling
        di = self._raw.flux_direct
        assert not du is None, 'Fluxes at observer not simulated. Set settings_rt.do_actinic_flux_at_observer flag to True' 
        d = self._merge_actinic_flux_files(du,dd,di,
                                           up_down_both, 
                                           diffuse_direct_both)
        wls = self._parent_vpc.wavelengths.copy()
        sol = self._parent_vpc.solar_spectrum.get_values(wls).flatten()
        if not normalized: d *= sol
        ser = pd.Series(d, wls)
        ser.index.name = 'wavelength'
        return ser
    
    def _get_irradiance_array(self, direction, component, jac_label=None):
        if direction not in ('down', 'up', 'net', 'sum'):
            raise ValueError('direction must be down, up, net or sum')
        if component not in ('direct', 'diffuse', 'total', 'single', 'multiple'):
            raise ValueError('component must be direct, diffuse, total, single or multiple')
        prefix = 'irradiance_' if jac_label is None else 'jacs_irradiance_'
        suffix = '' if jac_label is None else '_' + jac_label
        arrays = []
        names = (component + '_upwelling', component + '_downwelling') if component in ('single', 'multiple') else (
            'diffuse_upwelling', 'diffuse_downwelling', 'direct')
        for name in names:
            value = getattr(self._raw, prefix + name + suffix, None)
            if value is None:
                raise ValueError('Requested irradiance or irradiance Jacobians were not calculated.')
            arrays.append(np.asarray(value, dtype=float).copy())
        if component in ('single', 'multiple'):
            up, down = arrays
            return {'up': up, 'down': down, 'net': down-up, 'sum': down+up}[direction]
        up, down, direct = arrays
        if component == 'direct':
            up, down = np.zeros_like(up), direct
        elif component == 'total':
            down = down + direct
        return {'up': up, 'down': down, 'net': down-up, 'sum': down+up}[direction]

    def _irradiance_scale(self, normalized, units):
        if units not in ('photons', 'energy'):
            raise ValueError('units must be photons or energy')
        wls = self._parent_vpc.wavelengths.copy()
        if normalized:
            return np.ones(len(wls))
        scale = self._parent_vpc.solar_spectrum.get_values(wls).flatten()
        if units == 'energy':
            # photons / (s cm2 nm) -> W / (m2 nm), using exact SI constants.
            scale = scale * 1e4 * 6.62607015e-34 * 299792458.0 / (wls * 1e-9)
        return scale

    def get_irradiance(self, direction='down', component='total', normalized=False, units='photons'):
        """Spectral irradiance on a local horizontal plane at the observer.

        direction: down/up (positive magnitudes), net (down minus up), or sum.
        component: direct/diffuse/total/single/multiple. Upward direct solar
        irradiance is zero. Single includes one surface reflection or one
        atmospheric scattering. Single/multiple require the scattering-
        components setting; direct + single + multiple equals total.
        units: photons [photons/s/cm2/nm] or energy [W/m2/nm]. With normalized=True,
        return the ratio to TOA solar beam irradiance, independent of units.
        Observer must be at or below TOA. Returns a wavelength-indexed Series.
        """
        values = self._get_irradiance_array(direction, component).reshape(-1)
        result = pd.Series(values * self._irradiance_scale(normalized, units),
                           index=self._parent_vpc.wavelengths.copy(), name='irradiance')
        result.index.name = 'wavelength'
        return result

    def get_irradiance_jacobians(self, jac_label='absorption', direction='down', component='total',
                                normalized=False, units='photons'):
        """Layer-resolved derivatives of get_irradiance, with the same scaling.

        Uses settings_rt.jacobian_list labels and jacobian_layer_flags. Rows
        are layer-middle altitudes [m]; columns are wavelengths [nm]. Physical
        derivatives use the requested irradiance units per parameter unit;
        normalized=True returns derivatives of solar-normalized irradiance.
        """
        values = self._get_irradiance_array(direction, component, jac_label)
        wls = self._parent_vpc.wavelengths.copy()
        values = values.reshape(len(wls), self._parent_vpc.get_n_layers()).T[::-1]
        result = pd.DataFrame(values * self._irradiance_scale(normalized, units),
                              index=self._parent_vpc.get_layer_middles(), columns=wls)
        result.index.name = 'altitude'
        result.columns.name = 'wavelength'
        return result

    def get_radiance(self, normalized = False):
        """
        Returns the radiances [photons/s/cm2/sr/nm] at the simulated wavelengths
        as Pandas dataframe. If normalized is True, radiances are normalized
        against the TOA solar irradiance spectrum.        
        """
        assert not self._raw.radiances is None, 'Radiances not simulated. Set settings_rt.do_radiance flag to True'
        wls = self._parent_vpc.wavelengths.copy()
        d = pd.Series(self._raw.radiances.copy(), wls)
        if not normalized: d *= self._parent_vpc.solar_spectrum.get_values(wls).flatten()
        return d 
    
    def get_bamfs(self):
        """
        Calculate box-airmass factors (bAMFs) from Jacobians w.r.t. to each layers absorption optical depth.
        Currently only works for radiances (radiance_or_flux = 'radiance').
        Returns a pandas dataframe with bAMFs for various altitudes [m] (index) and wavelengths [nm] (columns). 
        """
        fail = (self._raw.radiances is None) \
                or not hasattr(self._raw, 'jacs_radiance_absorption')
        assert not fail, 'Radiances or radiance Jacobians w.r.t. absorption optical depth not simulated. Check flags in pyd.settings_rt.'
        rads = self._raw.radiances.copy()
        d_rad_d_abs = np.atleast_2d(self._raw.jacs_radiance_absorption).T
        bamfs = - d_rad_d_abs / rads
        df = pd.DataFrame(bamfs[::-1], columns = self._parent_vpc.wavelengths.copy(), index=self._parent_vpc.get_layer_middles())
        df.index.name = 'altitude'
        df.columns.name = 'wavelength'
        return df
    
    def get_jacobians(self, jac_label):
        _f = 'jacs_radiance_' + jac_label
        fail = (self._raw.radiances is None) \
                or not hasattr(self._raw, _f ) 
        assert not fail, 'Radiances or radiance Jacobians w.r.t. the aerosol types extinction profile not simulated. Check flags in pyd.settings_rt.'
        _ = np.atleast_2d(getattr(self._raw, _f)).T
        df = pd.DataFrame(_[::-1], columns = self._parent_vpc.wavelengths.copy(), index=self._parent_vpc.get_layer_middles())
        df.index.name = 'altitude'
        df.columns.name = 'wavelength'
        return df

    def get_flux_bamfs(self,
                       up_down_both = 'both', 
                       diffuse_direct_both = 'both'):
        fail = (self._raw.flux_diffuse_downwelling is None) \
                or not hasattr(self._raw, 'jacs_flux_diffuse_downwelling_absorption')
        assert not fail, 'Fluxes or flux Jacobians w.r.t. gas absorption not simulated. Check flags in settings_rt.'
        af = self.get_actinic_flux(up_down_both, diffuse_direct_both, normalized=True).values
        jdu = self._raw.jacs_flux_diffuse_upwelling_absorption
        jdd = self._raw.jacs_flux_diffuse_downwelling_absorption
        jdi = self._raw.jacs_flux_direct_absorption
        d_af_d_abs = self._merge_actinic_flux_files(jdu,jdd,jdi,
                                                    up_down_both, 
                                                    diffuse_direct_both)
        d_af_d_abs = np.atleast_2d(d_af_d_abs).T
        bamfs = - d_af_d_abs / af
        df = pd.DataFrame(bamfs[::-1], columns = self._parent_vpc.wavelengths.copy(), index=self._parent_vpc.get_layer_middles())
        df.index.name = 'altitude'
        df.columns.name = 'wavelength'
        return df
    
    def plot_bamfs(self, wavelength, new_fig = True, boundary_indexed = False, **plot_kwargs):
        if boundary_indexed: x = np.arange(len(self._parent_vpc.vertical_grid))
        else: x = self._parent_vpc.vertical_grid/1e3
        if new_fig: plt.figure()
        bamfs = self.get_bamfs()[wavelength].values
        x,y = get_stepped_profile(x, bamfs)
        #plot layer borders
        for a in x: plt.axhline(a, color = 'k', alpha = 0.1)
        plt.plot(y, x, **plot_kwargs)
        #global plot properties
        plt.ylabel('Altitude [km]')

    def plot_flux_bamfs(self, wavelength, 
                        up_down_both = 'both',
                        diffuse_direct_both = 'both',
                        new_fig = True, boundary_indexed = False, **plot_kwargs):
        if boundary_indexed: x = np.arange(len(self._parent_vpc.vertical_grid))
        else: x = self._parent_vpc.vertical_grid/1e3
        if new_fig: plt.figure()
        bamfs = self.get_flux_bamfs(up_down_both,diffuse_direct_both)[wavelength].values
        x,y = get_stepped_profile(x, bamfs)
        #plot layer borders
        for a in x: plt.axhline(a, color = 'k', alpha = 0.1)
        plt.plot(y, x, **plot_kwargs)
        #global plot properties
        plt.ylabel('Altitude [km]')

    
    def get_scds(self, gas_label):
        """ 
        Calculates trace gas slant column densities (SCDs) [1/cm2], based on 
        box-airmass factors and current gas vertical profile. 
            gas_label : gas to get SCDs for. Use the same labels as used in
            pyd.trace_gases container. By default, these are 
            'o4', 'no2', 'h2o' and 'o3'. 
        """
        dh = self._parent_vpc.get_layer_thicknesses() * 100 #layer thicknesses in cm
        bamfs = self.get_bamfs()
        if (gas_label == self._parent_vpc.trace_gases.o4_label): 
            prof = self._parent_vpc.get_o4_profile()
        else:
            prof = getattr(self._parent_vpc.trace_gases, gas_label).profile.values
        return (prof*dh)@bamfs
    
    def get_amfs(self, gas_label):
        """ 
        Calculates trace gas slant airmass factors (AMFs), based on 
        box-airmass factors and current gas vertical profile. 
            gas_label : gas to get AMFs for. Use the same labels as used in
            pyd.trace_gases container. By default, these are 
            'o4', 'no2', 'h2o' and 'o3'. 
        """
        if (gas_label == self._parent_vpc.trace_gases.o4_label): 
            vcd = self._parent_vpc.get_o4_vcd()
        else:
            vcd = getattr(self._parent_vpc.trace_gases, gas_label).profile.get_vertical_column()
        return self.get_scds(gas_label) / vcd
    
    def _merge_actinic_flux_files(self,
                                  du,dd,di,
                                  up_down_both, 
                                  diffuse_direct_both):
        """ Helper function to merge actinic flux related output files. """
        match (up_down_both, diffuse_direct_both):
            case ('up','diffuse'):   d = du.copy()
            case ('up','direct'):    d = np.zeros(dd.shape)
            case ('up','both'):      d = du.copy()
            case ('down','diffuse'): d = dd.copy()
            case ('down','direct'):  d = di.copy()
            case ('down','both'):    d = dd + di
            case ('both','diffuse'): d = dd + du
            case ('both','direct'):  d = di.copy()
            case ('both','both'):    d = du + dd + di
            case _: raise Exception('Invalid option for up_down_both or diffuse_direct_both. Check docstring for valid values!')
        return d
    
class Pydort(VerticalGrid):
    
    """
    Pydort main class.
    
    Important attributes:
    ---------------------
    
    wavelengths : 1D-array of simulated wavelengths [nm].
        Note: profile_wavelengths and reference_wavelength attributes of 
        MicrophysicalAerosol, OpticalAerosol and SpectralParameter objects must
        be covered by wavelengths. This can always be achieved using the 
        shift_reference_wavelengths() method (see below).
    vertical_grid : 1D-array.
        Vertical grid (altitudes of layer boundaries [m], monotonically increasing).
        Must be consistent with vertical grids of all profile objects before calling run_rt().
        It is reommended to modify the vertical grid exclusively using the 
        change_vertical_grids() method (see below). 
        Use get_inconsistent_vertical_grids() method to reveal inconsistencies.
    pressure_profile : PressureProfile object.
    temperature_profile : TemperatureProfile object.
    trace_gases : TraceGasContainer object.
        Contains the trace gases considered in the simulation.
    aerosols : AerosolContainer object.
        Contains the aerosol species considered in the simulation.
    surface_albedo : SurfaceAlbedo object.
    viewing_geometry : ViewingGeometry object.
    solar_spectrum : LiteratureSpectrum object.
        Contains the spectral solar irradiance at TOA.
    settings_rt : SettingsRT object.
        Settings for the radiative transport calculation.
    results_rt : ResultsRT object.
        Will be filled with results once run_rt() is called.
        Contains methods for post-processing, including box-airmass factor and
        slant column density calculations (see docstrings of the ResultsRT class).
    """
    
    def __init__(self, 
                 trace_gases = ['o4','o3','no2','h2o','hcho','bro','io','chocho','hono'],
                 spectral_smoothing_on_initialisation = True):     
        self._wavelengths = DEFAULT_WAVELENGTHS.copy()
        self.trace_gases = TraceGasContainer(trace_gases, 
                                             spectral_smoothing_on_initialisation=spectral_smoothing_on_initialisation)
        self.aerosols = AerosolContainer()        
        self.pressure_profile = PressureProfile()
        self.temperature_profile = TemperatureProfile() 
        self._vertical_grid = DEFAULT_VERTICAL_GRID.copy()
        self.surface_albedo = SurfaceAlbedo()
        self.viewing_geometry = ViewingGeometry()
        self.solar_spectrum = LiteratureSpectrum(SOLAR_IRRADIANCE_SPECTRUM + '.dat', unit_string='photons/s/cm2/nm',
                                                 spectral_smoothing_on_initialisation = spectral_smoothing_on_initialisation)
        self.results_rt = None
        self.settings_rt = SettingsRT()
        self.results_rt = ResultsRT(self)

    def __repr__(self):
        s = 'Pydort object. For details on the object structure refer to the PYDORT_V2 manual.'
        # Add list of attributes:
        atts = list(vars(self).keys()).copy()
        for ia,a in enumerate(atts):
            if a[0]=='_': atts[ia] = a[1:]
        s += '\nAttributes:\n\t'
        s += '\n\t'.join(atts)
        # Add list of methods:
        meths = [m for m in dir(self) if callable(getattr(self, m)) and m[0] != '_']
        s += '\nMethods:\n\t'
        s += '\n\t'.join(meths)
        return s

    def copy(self):
        return copy.deepcopy(self)
     
    @property
    def wavelengths(self):
        return self._wavelengths
        
    @wavelengths.setter
    def wavelengths(self, value):
        self._wavelengths = np.array(value)
    
    @property
    def vertical_grid(self):
        return self._vertical_grid
        
    @vertical_grid.setter
    def vertical_grid(self, value):
        wstring = 'Changing the vertical_grid parameter of the Pydort class alone leads to inconsistency with grids in child profile classes (pressure, temperature, gases, ...). ' \
                    + 'Use the change_vertical_grids() method to consistently interpolate all profiles to a new grid.'
        warnings.warn(wstring)
        self._vertical_grid = np.array(value)
    
    def get_profile_attributes(self):
        """ Returns all profile objects of the current Pydort object as a dictionary. 
        Contains references to the objects, i.e. modifications in the dictionary
        are readily applied inside the Pydort object.
        """
        aerosol_dic = self.aerosols.to_dictionary()
        gas_dic = self.trace_gases.to_dictionary()
        plist = [self.pressure_profile,
                 self.temperature_profile]
        keys = ['pressure','temperature'] 
        for a in aerosol_dic: 
            plist += [aerosol_dic[a].profile]
            keys += [a]
        for g in gas_dic: 
            if hasattr(gas_dic[g], 'profile'): 
                plist += [gas_dic[g].profile]
                keys += [g]
        return dict(zip(keys, plist))

    def get_spectra_attributes(self):
        """ Returns all literature spectrum objects of the current Pydort object 
        as a dictionary. Contains references to the objects, i.e. modifications in the dictionary
        are readily applied inside the Pydort object.
        """
        gas_dic = self.trace_gases.to_dictionary()
        plist = [self.solar_spectrum]
        keys = ['solar_spectrum']
        for g in gas_dic: 
            plist += [gas_dic[g].xs]
            keys += ['xs_' + g]
        return dict(zip(keys, plist))

    def get_inconsistent_vertical_grids(self):
        """ Provides a list of inconsistent vertical grids between various profile objects within the Pydort object. """
        profs = self.get_profile_attributes() 
        refgrid = self.vertical_grid
        inconsistent_grids = []
        for ip,p in enumerate(profs):
            grid = profs[p].vertical_grid
            if len(grid) != len(refgrid):
                inconsistent_grids += [p]
                continue
            if (grid != refgrid).any():
                inconsistent_grids += [p]
                continue
        return inconsistent_grids
    
    def disable_jacobians_above(self, altitude, inplace=True):
        """
        Disables the Jacobians above a certain altitude by setting pyd.settings_rt.jacobian_layer_flags. 
        Any layer extending into the altitude range (between zero and altitude) will stay activated.
        altitude : altitude limit in meters
        """
        flags = [a < altitude for a in self.get_layer_bottoms()] 
        if inplace: 
            self.settings_rt.jacobian_layer_flags = flags
            return sum(flags) # return number of jacobian active layers
        else:
            return flags
    
    def disable_jacobians_below(self, altitude, inplace=True):
        """
        Disables the Jacobians below a certain altitude by setting pyd.settings_rt.jacobian_layer_flags. 
        Any layer extending into the altitude range (between TOA and altitude) will stay activated.
        altitude : altitude limit in meters
        """
        flags = [a > altitude for a in self.get_layer_tops()] 
        if inplace: 
            self.settings_rt.jacobian_layer_flags = flags
            return sum(flags) # return number of jacobian active layers
        else:
            return flags
    
    def enable_jacobians_between(self, min_altitude, max_altitude, inplace=True):
        f1 = self.disable_jacobians_above(max_altitude, False) 
        f2 = self.disable_jacobians_below(min_altitude, False) 
        flags = [a and b for a, b in zip(f1, f2)]
        if inplace: 
            self.settings_rt.jacobian_layer_flags = flags
            return sum(flags) # return number of jacobian active layers
        else:
            return flags
        
    def change_spectral_resolution(self, new_smoothing_fwhm):
        sas = self.get_spectra_attributes()
        for s in sas:
            sas[s].smoothing_fwhm = new_smoothing_fwhm

    def _observer_at_or_above_toa(self):
        observer_height = float(self.viewing_geometry.observer_height_agl)
        toa_height = float(np.asarray(self.vertical_grid)[-1])
        return observer_height >= toa_height - 1e-6

    def _resolve_rt_backend(self, warn=True):
        correction = int(self.settings_rt.sphericity_correction)
        observer_at_or_above_toa = self._observer_at_or_above_toa()

        def qs_fallback(reason):
            if warn:
                warnings.warn(
                    f'sphericity_correction={correction}: {reason}; '
                    'using VLIDORT-QS full spherical correction instead.',
                    RuntimeWarning, stacklevel=3,
                )
            return {
                'model': RadiativeTransferModel.QS.value,
                'backend_mode': RadiativeTransferModel.QS.value,
                'sphericity_interpolation_points': 4,
                'fallback_reason': reason,
            }

        height = float(self.viewing_geometry.observer_height_agl)
        grid = np.asarray(self.vertical_grid, dtype=float)
        radius = float(self.settings_rt.planet_radius) * 1000.0
        vza = float(self.viewing_geometry.vza)
        ray_radius = (radius + height) * np.sin(np.deg2rad(min(vza, 180.0-vza)))
        misses_surface = ray_radius >= radius + grid[0]

        if correction == SphericityCorrection.EXPERIMENTAL.value:
            return {
                'model': RadiativeTransferModel.AIRCRAFT.value,
                'backend_mode': 'aircraft_combined',
                'sphericity_interpolation_points': self.settings_rt.n_aircraft_source_points,
            }

        if correction == SphericityCorrection.NONE.value:
            if self.settings_rt.do_radiance and vza > 90 and misses_surface:
                return qs_fallback('the spherical first-order limb ray does not intersect the surface')
            return {
                'model': RadiativeTransferModel.STANDARD.value,
                'backend_mode': 'standard',
                'sphericity_interpolation_points': 2,
            }

        if correction == SphericityCorrection.TWO_POINT.value:
            if observer_at_or_above_toa:
                return {
                    'model': RadiativeTransferModel.STANDARD.value,
                    'backend_mode': 'standard_spherical',
                    'sphericity_interpolation_points': 2,
                }
            if self.settings_rt.do_radiance:
                if misses_surface:
                    return qs_fallback('the aircraft 2-point ray cannot be extended to the surface')
                if not np.any(np.isclose(grid[1:-1], height, rtol=0, atol=1e-6)):
                    return qs_fallback('the aircraft 2-point observer must be on an interior layer boundary')
            return {
                'model': RadiativeTransferModel.AIRCRAFT.value,
                'backend_mode': 'aircraft_spherical',
                'sphericity_interpolation_points': 2,
            }

        if correction == SphericityCorrection.THREE_POINT.value:
            if observer_at_or_above_toa:
                return {
                    'model': RadiativeTransferModel.STANDARD.value,
                    'backend_mode': 'standard_spherical',
                    'sphericity_interpolation_points': 3,
                }
            if warn:
                warnings.warn(
                    'sphericity_correction=3 is not available for aircraft-height '
                    'observers; using VLIDORT-QS full spherical correction instead.'
                )
            return {
                'model': RadiativeTransferModel.QS.value,
                'backend_mode': RadiativeTransferModel.QS.value,
                'sphericity_interpolation_points': 3,
            }

        return {
            'model': RadiativeTransferModel.QS.value,
            'backend_mode': RadiativeTransferModel.QS.value,
            'sphericity_interpolation_points': correction,
        }
        
    def _generate_rt_input(self, run_aerosol_models = True):
        
        #PREPARATIONS AND CHECKS:
        
        wls = np.array(self.wavelengths)
        path_qs_backend = get_qs_backend_path()

        # Check vertical grid consistency:
        incon = self.get_inconsistent_vertical_grids()
        s = 'Inconsistent vertical grids for:\n' + str(incon) \
            + '\nApply change_vertical_grid functions to ensure consistency!'
        assert len(incon)==0, s 

        gas_dic = self.trace_gases.to_dictionary()


        # Run aerosol models if necesary:
        aerosol_dic = self.aerosols.to_dictionary()        
        if run_aerosol_models: self.aerosols.run_models(self.wavelengths)
        
        # Calculate necessary aerosol-jacobians:
        if self.settings_rt.do_radiance_jacobians or self.settings_rt.do_actinic_flux_jacobians or self.settings_rt.do_irradiance_jacobians:
            for jlbl in self.settings_rt.jacobian_list:
                if jlbl in aerosol_dic or not '_' in jlbl: continue # skip non-aerosol-parameter jacobians
                if not run_aerosol_models: warnings.warn('Mie model jacobians are calculated based on results from old Mie model run!')
                albl,plbl = jlbl.rsplit('_', 1)
                assert albl in aerosol_dic, jlbl + ' does not refer to a known aerosol label.'
                a = aerosol_dic[albl]
                if isinstance(a,MicrophysicalAerosol):
                    a.get_jacobians(plbl, wls, rerun_fd_base_simulation=False)
                elif isinstance(a,OpticalAerosol):
                    a.get_jacobians(plbl, wls)
        
        # If no aerosol provided, create dummy aerosol (with zero profile), so the RT drivers are happy:
        if not len(aerosol_dic):
            oa = OpticalAerosol()
            oa.profile.change_vertical_grid(self.vertical_grid)
            oa.profile.fill_constant(0)
            oa.shift_reference_wavelengths(self.wavelengths[-1])
            oa.run_model(self.wavelengths)
            aerosol_dic['_zero_dummy_'] = oa
        
        # Check viewing geometries:
        self.viewing_geometry.check_ranges()

        active_backend = getattr(self, '_active_rt_backend_settings', None)
        if active_backend is None:
            active_backend = self._resolve_rt_backend(warn=False)
        
        # Check consistency of Jacobian settings:
        # !!! REVISE ONCE JACOBIAN LABELS ARE FINALIZED
        # _anames = self.aerosols.get_name_list()
        # _valid = ['absorption'] + _anames + [a+'_ssa' for a in _anames]
        # for q in self.settings_rt.jacobian_list:
        #     assert q in _valid, q + ' is not a valid jacobian identifier. Check pyd.settings_rt.jacobian_list. Valid entries are\n' + str(_valid)
        if not len(self.settings_rt.jacobian_list):
            if self.settings_rt.do_radiance_jacobians or self.settings_rt.do_actinic_flux_jacobians or self.settings_rt.do_irradiance_jacobians:
                warnings.warn('Jacobian list empty. All Jacobian calculation flags will be set to False.')
                self.settings_rt.do_actinic_flux_jacobians = self.settings_rt.do_radiance_jacobians = False
                self.settings_rt.do_irradiance_jacobians = False
        
        # Delete old data specific gas species and aerosol types:
        for f in Path(path_qs_backend + 'INPUT').glob("**/*.inp"): f.unlink()                 

        # WRITE DATA TO COMMON RT INPUT FILES:
        
        # Vertial grid:
        grid = self.vertical_grid[::-1]
        header = str(len(grid)-1) + ' ! n_layers (n_vertical_grid-1)\n'
        header += 'Layer boundary altitudes [m]'
        assert ~np.isnan(grid).any(), 'Vertical grid contains NAN values.'
        np.savetxt(path_qs_backend + 'INPUT/vertical_grid.inp', grid, 
                   header=header, comments='')
            
        # Wavelength grid:
        array_to_rt_input('wavelengths.inp', wls,
                           'Wavelength [nm]')  
        
        # Aerosol type names:
        with open(path_qs_backend + 'INPUT/aertype_names.inp', 'w') as f:
            f.write(str(len(aerosol_dic)) + ' ! n_aertypes\n')
            for a in aerosol_dic:
                f.write(f"{a}\n")
          
        # Aerosol profiles:
        d_list = []
        for a in aerosol_dic:
             aer = aerosol_dic[a]
             d_list += [aer.profile.values]
        d = np.vstack(d_list).T[::-1]
        assert bool(len(d_list)), 'Provide at least one (dummy) aerosol type!'
        is_aerosol_zero = not bool(d.sum())  # needed later for configuration file
        # Raise the background only for aerosol types that are actually present.
        # Preserve the legacy tiny guard for inactive types / the Rayleigh dummy.
        floors = np.where(np.any(d > 0, axis=0), LOWER_LIMIT_AEROSOL_EXTINCTION,
                          min(LOWER_LIMIT_AEROSOL_EXTINCTION, 1e-10))
        d = np.maximum(d, floors)
        header_comment = ' '.join(list(aerosol_dic.keys()))
        array_to_rt_input('aerosol_extinction_profiles.inp', d, header_comment)
          
        # Aerosol single scattering albedos:
        d_list = []
        for a in aerosol_dic:
              d = aerosol_dic[a].results.bulk_properties['single_scattering_albedo']
              d_list += [d.copy()]
        d = np.vstack(d_list).T
        header_comment = ' '.join(list(aerosol_dic.keys()))
        array_to_rt_input('aerosol_single_scattering_albedos.inp', d, header_comment)
   
        # Aerosol extinction multipliers:
        d_list = []
        for a in aerosol_dic:
              d = aerosol_dic[a].results.bulk_properties['extinction_multiplier']
              d_list += [d.copy()]
        d = np.vstack(d_list).T
        header_comment = ' '.join(list(aerosol_dic.keys()))
        array_to_rt_input('aerosol_extinction_multipliers.inp', d, header_comment)
        
        # Aerosol moments:
        for a in aerosol_dic:
            d = aerosol_dic[a].results.phasefunction_moments.copy()
            # d = d.fillna(0)
            m1 = d.index.get_level_values(1) <= self.settings_rt.n_max_pf_moments # mask to remove unwanted additional moments.
            m2 = d.iloc(1)[0].abs().values >= PF_MOMENT_CUTOFF # mask to remove moments with negligible impact
            d = d.loc[m1 & m2]
            d = d.reset_index()
            # Replace wavelengths in multiindex by wavelength index (1,2,3,...):
            d.iloc(1)[0] = ((d.iloc(1)[0].diff()>0).cumsum()+1).astype(int)
            d = d.values
            assert ~np.isnan(d).any(), 'Aerosol moments contain NAN values.'
            header = str(len(d)) + ' ! n_all_moments\n'
            header += str('wavelength_index moment_index moment_value')
            f = path_qs_backend + 'INPUT/aerosol_moments/' + a + '.inp'
            np.savetxt(f,d, header=header, comments='', 
                       fmt='%i %i %.15e %.15e %.15e %.15e %.15e %.15e')            
            
        # Pressure profile:
        d = self.pressure_profile.values[::-1] # upside down to meet RT-driver convention
        array_to_rt_input('pressure_profile.inp', d, 'Pressure [hPa]')

        # Temperature profile:
        d = self.temperature_profile.values[::-1] # upside down to meet RT-driver convention
        array_to_rt_input('temperature_profile.inp', d, 'Temperature [K]')
        
        # Water vapor pressure profile:    
        l = self.trace_gases.water_vapor_label
        lbls = self.trace_gases.get_name_list()
        if not l in lbls: 
            warnings.warn('trace_gases.water_vapor_label not in the list of gases! Water vapor pressure is set to pressure x 0.01!')
            d = self.pressure_profile.values[::-1] * 1e-2
        else:
            d = gas_dic[l].profile.get_partial_pressures(self.temperature_profile.values)[::-1]
        array_to_rt_input('vapor_pressure_profile.inp', d, 
                           'Water vapor partial pressure [hPa], required for Rayleigh scattering properties.')
        
        # Gas species list:
        with open(path_qs_backend + 'INPUT/gas_names.inp', 'w') as f:
            f.write(str(len(gas_dic)) + ' ! n_gases\n')
            for g in gas_dic:
                f.write(f"{g}\n")
                
        # Gas absorption cross sections:
        for g in gas_dic:
            xs = gas_dic[g].xs
            d = xs.get_values(self.wavelengths)
            header_comment = 'CrossSection [cm2]    TempCoeffs'
            array_to_rt_input('xs/' + g + '.inp', d, header_comment)         
        
        # Gas vertical profiles:
        o4_lbl = self.trace_gases.o4_label
        d_list = []
        for g in gas_dic:
            tg = gas_dic[g]
            if g == o4_lbl: 
                _ = self.get_o4_profile()
                if not self.settings_rt.do_o4_absorption: _ *= 0.
                d_list += [_.copy()]
            else: d_list += [tg.profile.values]
        assert bool(len(d_list)), 'Provide at least one (dummy) trace gas!'
        d = np.vstack(d_list).T[::-1]
        d[d<LOWER_LIMIT_GAS_CONCENTRATION] = LOWER_LIMIT_GAS_CONCENTRATION
        header_comment = ' '.join(list(gas_dic.keys())) + ', all in number concentrations [cm-3]'
        array_to_rt_input('gas_concentration_profiles.inp', d, header_comment)
        
        # Surface albedo:
        d = np.atleast_1d(self.surface_albedo.get_value(wls))
        d[d<LOWER_LIMIT_SURFACE_ALBEDO] = LOWER_LIMIT_SURFACE_ALBEDO
        array_to_rt_input('surface_albedos.inp', d,
                           'Surface albedo for each wavelength')
        
        # Jacobian list:
        with open(path_qs_backend + 'INPUT/jacobian_list.inp', 'w') as f:
            f.write(str(len(self.settings_rt.jacobian_list)) + ' ! n_jacobians\n')
            for s in self.settings_rt.jacobian_list:
                f.write(f"{s}\n")
        
        # Aerosol model Jacobians:
        if self.settings_rt.do_radiance_jacobians or self.settings_rt.do_actinic_flux_jacobians or self.settings_rt.do_irradiance_jacobians:
            for jlbl in self.settings_rt.jacobian_list:
                if jlbl in aerosol_dic or not '_' in jlbl:continue
                albl,plbl = jlbl.rsplit('_', 1)
                assert albl in aerosol_dic, jlbl + ' does not refer to a known aerosol label.'
                jacs = getattr(aerosol_dic[albl].jacobians, plbl)
                # Export d_extmul_d_param:
                d = jacs.bulk_properties.extinction_multiplier.values.T
                array_to_rt_input('aerosol_model_jacobians/' + jlbl + '_extmul.inp', d, 
                                   'd_extmul_d_' + plbl + ' at each wavelength')
                # Export d_ssa_d_param:
                d = jacs.bulk_properties.single_scattering_albedo.values.T
                array_to_rt_input('aerosol_model_jacobians/' + jlbl + '_ssa.inp', d, 
                                   'd_ssa_d_' + plbl + ' at each wavelength')
                # Export d_moms_d_param:
                # Create masks based on moment values:
                d = aerosol_dic[albl].results.phasefunction_moments.copy()
                d = d.fillna(0)
                m1 = d.index.get_level_values(1) <= self.settings_rt.n_max_pf_moments # mask to remove unwanted additional moments.
                m2 = d.iloc(1)[0].values >= PF_MOMENT_CUTOFF # mask to remove moments with negligible impact
                # Read in actual moment jacobians:
                d = jacs.phasefunction_moments.copy()
                d = d.loc[m1 & m2]
                d = d.reset_index()
                # Replace wavelengths in multiindex by wavelength index (1,2,3,...):
                d.iloc(1)[0] = ((d.iloc(1)[0].diff()>0).cumsum()+1).astype(int)
                d = d.values
                assert ~np.isnan(d).any(), 'Aerosol moments contain NAN values.'
                header = str(len(d)) + ' ! n_all_moments\n'
                header += str('wavelength_index moment_index moment_value')
                f = path_qs_backend + 'INPUT/aerosol_model_jacobians/' + jlbl + '_moms.inp'
                np.savetxt(f,d, header=header, comments='', 
                       fmt='%i %i %.12e %.12e %.12e %.12e %.12e %.12e')
        
        # # Jacobian aertype index list:
        # d = []
        # if 'absorption'in self.settings_rt.jacobian_list: d+=[0]
        # for i, a in enumerate(aerosol_dic):
        #     if a in self.settings_rt.jacobian_list: d+=[i]
        # array_to_rt_input('jacobian_list_aertype_index.inp', d, 'Aerosol type index for requested Jacobians. For absorption it is zero.')
        
        # # Jacobian aertype flags:
        # jaf = self.settings_rt.jacobian_aertype_flags        
        # if isinstance(jaf, str):
        #     if jaf == 'all': d = np.array([True]*self.aerosols.get_n_species())
        # elif isinstance(jaf, list):
        #     d = np.array(jaf.copy())
        #     assert len(d) == self.aerosols.get_n_species(), 'Length of settings_rt.jacobian_aertype_flags must match number of aerosols!'
        # array_to_rt_input('jacobian_aertype_flags.inp', d, 
        #                    'Jacobian calculation flag for each aerosol type.')

        # Jacobian layer flags:
        #!!! Are these considered in the result read-in-routine??? Might not work...
        jlf = self.settings_rt.jacobian_layer_flags        
        if isinstance(jlf, str):
            assert jlf == 'all', 'settings_rt.jacobian_layer_flags must be "all" or iterable of booleans.'   
            d = np.array([True]*self.get_n_layers())
        else:
            d = np.array(jlf.copy())
            assert len(d) == self.get_n_layers(), 'Length of settings_rt.jacobian_layer_flags must match number of layers!'
        d = d[::-1]
        array_to_rt_input('jacobian_layer_flags.inp', d, 
                           'Jacobian calculation flag for each layer, starting from bottom.')
        
        
        #CREATE MAIN CONFIGURATION FILE:
            
        do_deltam = self.settings_rt.do_deltam_phase_function_scaling
        do_rayleigh_only = False
        if is_aerosol_zero: 
            do_deltam = False
            do_rayleigh_only = True

        if self.settings_rt.do_irradiance_scattering_components and do_deltam:
            raise ValueError('Scattering-order irradiances require do_deltam_phase_function_scaling=False.')

        values = [str(self.settings_rt.do_single_scattering_only),
                  str(self.settings_rt.do_vectorized_simulation),
                  str(self.settings_rt.do_radiance),
                  str(self.settings_rt.do_actinic_flux_at_observer),
                  str(False),
                  str(do_rayleigh_only),
                  str(do_deltam),
                  str(self.settings_rt.n_halfspace_streams),
                  str(self.settings_rt.n_full_sphericity_fine_layers),
                  str(float(self.viewing_geometry.observer_height_agl)),
                  str(self.settings_rt.planet_radius),
                  str(self.viewing_geometry.sza),
                  str(self.viewing_geometry.vza),
                  str(self.viewing_geometry.get_rt_raa()[0]),
                  str(self.settings_rt.do_radiance_jacobians),
                  str(self.settings_rt.do_actinic_flux_jacobians),
                  active_backend['backend_mode'],
                  str(int(active_backend['sphericity_interpolation_points'])),
                  str(self.settings_rt.do_irradiance_at_observer),
                  str(self.settings_rt.do_irradiance_jacobians),
                  str(self.settings_rt.do_irradiance_scattering_components),
                  str(self.settings_rt.irradiance_scattering_step)]
        comments = ['Flag for performing only single scatter (FO) radiative transfer, no multiple scattering (MS)',
                    'Flag for performing a vector calculation',
                    'Flag for calculating the radiance at the desired viewing geometry',
                    'Flag for calculating actinic fluxes at observer altitude additionally.',
                    'Flag for calculating actinic flux profiles (disabled by the Python wrapper)',
                    'Flag for Rayleigh only. If true, only moment order 0 to 2 are considered for scattering phasefunction.',
                    'Flag for using deltaM Scaling (normally true with MS contributions and aerosols)',
                    'Number of discrete ordinate streams',
                    'Subdivision of layers for full-sphericity calculations',
                    'Observer height [m]',
                    'Planet radius [km]',
                    'Solar zenith angle',
                    'Viewing zenith angle',
                    'Relative azimuth angle',
                    'Flag for radiance jacobian calculations',
                    'Flag for actinic flux jacobian calculation',
                    'Resolved radiative-transfer backend mode',
                    'Resolved sphericity interpolation points for Standard/AC corrections',
                    'Flag for horizontal irradiance at observer altitude',
                    'Flag for irradiance Jacobians',
                    'Flag for irradiance scattering-order decomposition',
                    'Scattering-strength expansion step']
        with open(path_qs_backend + 'INPUT/configuration.inp', 'w') as f:
            for v,c in zip(values,comments):
                f.write(v + '\t! ' + c + '\n')
    
    def _run_qs_backend(self, run_aerosol_models=True):
        """Run the VLIDORT-QS executable with the common RT input files."""
        self.settings_rt.validate_backend_settings()
        path_qs_backend = get_qs_backend_path()
        for f in Path(path_qs_backend + 'OUTPUT').glob("**/*.out"): f.unlink() # Delete all former output files.
        self._generate_rt_input(run_aerosol_models=run_aerosol_models)
        run_subprocess(path_qs_backend + 'DRIVER.exe')
        self.results_rt = ResultsRT(self) # Create new clean result instance to prevent mixing old and new results
        self.results_rt._read_raw_results(path_qs_backend + 'OUTPUT/')

    def _stage_standard_backend_input(self):
        source_input = Path(get_qs_backend_path()) / 'INPUT'
        standard_path = Path(get_standard_backend_path())
        target_input = standard_path / 'INPUT'
        target_output = standard_path / 'OUTPUT'
        if not standard_path.exists():
            raise FileNotFoundError('Standard/AC backend path does not exist: ' + standard_path.as_posix())
        if not (standard_path / 'DRIVER.exe').exists():
            raise FileNotFoundError(
                'DRIVER.exe not found in ' + standard_path.as_posix()
                + '. Build it with "mingw32-make -C VLIDORT_RT\\Fortran\\vlidort_standard all".'
            )
        target_input.mkdir(parents=True, exist_ok=True)
        target_output.mkdir(parents=True, exist_ok=True)
        for f in target_input.glob("**/*.inp"):
            f.unlink()
        shutil.copytree(source_input, target_input, dirs_exist_ok=True)

    def _run_standard_backend(self, run_aerosol_models=True):
        self._generate_rt_input(run_aerosol_models=run_aerosol_models)
        self._stage_standard_backend_input()
        path_standard = get_standard_backend_path()
        for f in Path(path_standard + 'OUTPUT').glob("**/*.out"):
            f.unlink()
        run_subprocess(path_standard + 'DRIVER.exe')
        self.results_rt = ResultsRT(self)
        self.results_rt._read_raw_results(path_standard + 'OUTPUT/')

    def run_rt(self, run_aerosol_models=True):
        """
        Run the selected radiative-transfer backend.
        """
        self.settings_rt.validate_backend_settings()
        active_backend = self._resolve_rt_backend(warn=True)
        self._active_rt_backend_settings = active_backend
        try:
            if active_backend['model'] == RadiativeTransferModel.QS.value:
                return self._run_qs_backend(run_aerosol_models=run_aerosol_models)
            if active_backend['model'] in [
                    RadiativeTransferModel.STANDARD.value,
                    RadiativeTransferModel.AIRCRAFT.value]:
                return self._run_standard_backend(run_aerosol_models=run_aerosol_models)
            raise RuntimeError('Unexpected radiative transfer model: ' + str(active_backend['model']))
        finally:
            self._active_rt_backend_settings = None

  
    def shift_reference_wavelengths(self, new_reference_wavelength):
        """
        Shifts all reference wavelengths within the Pydort object (i.e. profile_wavelength attributes 
        in OpticalAerosol and MicrophysicalAerosol objects, and reference_wavelength attributes in 
        SpectralParameter objects) to new_reference_wavelength. 
        All related parameters are adjusted to maintain the aerosol/parameter properties.
        """
        aer_dic = self.aerosols.to_dictionary()
        for a in aer_dic:
            aer = aer_dic[a]
            aer.shift_reference_wavelengths(new_reference_wavelength)
        if isinstance(self.surface_albedo, SpectralParameter):
            self.surface_albedo.shift_reference_wavelength(new_reference_wavelength)
            
    def change_vertical_grids(self, new_vertical_grid):
        """
        Interpolates vertical_grid attributes of all profiles in the 
        current Pydort object (PressureProfile, TemperatureProfile, GasProfile and 
        ExtinctionProfile objects) to new_vertical_grid.
        """
        self._vertical_grid = np.array(new_vertical_grid).copy()
        prof_attr = self.get_profile_attributes()
        for a in prof_attr:
            prof_attr[a].change_vertical_grid(new_vertical_grid)
            
    def elevate_ground_to_layer_boundary(self,layer_boundary):
        """
        Changes the ground elevation by removing the lower layer_boundary layers from all profiles.
        Observer altitude (viewing_geometry.observer_height_agl) is shifted to maintain barometric altitude.
        See also manual section 'Ground elevation'.
        
        Parameters
        ----------
        layer_boundary : integer
            The number of layers that will be removed.
        """
        assert len(self.get_inconsistent_vertical_grids()) == 0,\
            'Ground cannot be elevated as profile grids are insconsistent.\
            Use change_vertical_grids() to make grids consistent.' 
        alt = self.vertical_grid[layer_boundary]
        profs = self.get_profile_attributes() 
        for ip,p in enumerate(profs):
            profs[p].vertical_grid = profs[p].vertical_grid[layer_boundary:]
            profs[p].vertical_grid -= alt
            profs[p].values = profs[p].values[layer_boundary:]
        self.viewing_geometry.observer_height_agl -= alt
        self.vertical_grid = self.vertical_grid[layer_boundary:]
        self.vertical_grid -= alt
  
    def elevate_ground(self, ground_elevation, grid_mode='compress'):
        """
        Changes the ground elevation by removing the ground_elevation lower metres 
        of the atmosphere and squeezing the vertical grid correspondingly.
        The procedure is skectched in the manual section 'Ground elevation'.
        Notes:
            The number of layers is maintained.
            The grid's TOA (AGL) changes.
            Observer altitude (viewing_geometry.observer_height_agl) is shifted 
            to maintain barometric altitude (!!!).
        
        Parameters
        ----------
        ground_elevation : integer
            Height by which the ground shall be elevated [m].
        grid_mode : string, how shall the vertical grid be treated?
            'compress' : grid is compressed against TOA. Number of layers is conserved. Levels change.
            'cutoff' : grid will be cut off at the bottom. Levels are conserved. Number of layers changes.
        
        returns:
            if grid_mode == 'cutoff': number of surface layers that where entirely removed.
            else: None
        """
        ground_elevation = round(ground_elevation)

        if ground_elevation == 0: return 0 if grid_mode == 'cutoff' else None # Do nothing if ground elevation is zero
        assert ground_elevation > 0, 'Negative ground elevations are not allowed.'
        assert len(self.get_inconsistent_vertical_grids()) == 0,\
            'Ground cannot be elevated as profile grids are insconsistent. Use change_vertical_grids() to make grids consistent.' 
        assert grid_mode in ['compress','cutoff'],\
            "Allowed values for grid_mode are 'compress' and 'cutoff'."
        profs = self.get_profile_attributes()
        
        if grid_mode == 'compress':
            span  = max(self.vertical_grid) - min(self.vertical_grid)
            new_grid = np.round(self.vertical_grid * (1 - ground_elevation/span) + ground_elevation).astype(int) # Create squeezed grid
            self.change_vertical_grids(new_grid) # Interpolate all profiles to squeezed grid. Includes altitude offset.
            # Remove altitude offest from all profile grids:       
            self._vertical_grid = self._vertical_grid - ground_elevation
            for p in profs:
                profs[p].vertical_grid = profs[p].vertical_grid - ground_elevation
        
        elif grid_mode == 'cutoff':
            ilc = np.where(self.vertical_grid > ground_elevation)[0][0] # layer number where to cut
            for p in profs.values():
                p.vertical_grid = np.hstack([0, p.vertical_grid[ilc:]-ground_elevation]) 
                p.values = p.values[ilc-1:]
            self._vertical_grid = np.hstack([0, self.vertical_grid[ilc:]-ground_elevation]) 
        
        self.viewing_geometry.observer_height_agl = self.viewing_geometry.observer_height_agl - ground_elevation  
        return ilc-1 if grid_mode == 'cutoff' else None

    def get_rayleigh_extinction_profile(self, wavelength):
        return np.array([get_rayleigh_extinction_coeff(wavelength, t, p, h2o) 
                  for t, p, h2o in zip(self.temperature_profile.values,
                                       self.pressure_profile.values, 
                                       self.trace_gases.h2o.profile.values)])
   
    def get_o4_profile(self):
        """ Calculates the O4 number concentration profile [molec/cm3] from 
        pressure and temperature profile. Returns it as 1D array."""
        t = self.temperature_profile.values
        p = self.pressure_profile.values
        assert len(t) == len(p), 'O4 profile calculation requires pressure and temperature profiles on same grid!' 
        o4_conc = (p*1e2/KB/t * 1e-6 * O2_MOLAR_FRACTION)**2 #molecule number concentration [molecs/cm^3]
        return o4_conc

    def get_o4_vcd(self):
        """ Calculates the O4 vertical column density [molec^2/cm^5] from temperature and pressure profiles. """
        h = self.pressure_profile.get_layer_thicknesses()*1e2
        o4 = self.get_o4_profile()
        return (h*o4).sum()
    
    # def save_raw_results_rt(self, location):
    #     """ Writes the raw results (self.results_rt._raw) from VLIDORT-QS RT simulation to a user defined folder. 
    #     Not crazily helpful, though. Better store the results you actually need writing your own functions. """
    #     if not location[-1] == '/': location += '/' 
    #     for a in dir(self.results_rt):
    #         if (not 'flux' in a) or (not 'radiance' in a): continue
    #         arr = getattr(self.results_rt._raw,a)
    #         save_array_to_file(location + a + '.npy', arr)
    
    # def load_raw_results_rt(self, location):
    #     """ Reads raw results from VLIDORT-QS RT simulation from a user defined folder into self.results_rt._raw. 
    #     Only useful if the pyd object's settings are the same as during simulation, otherwise post-processing fails. """
    #     if not location[-1] == '/': location += '/'
    #     self.results_rt = ResultsRT(parent_vpc = self, raw_data_path = location, ftype = '.npy')
#%% TEST

# DEFAULT_VERTICAL_GRID is loaded from physics_data/default_vertical_grid.dat.

# pyd = Pydort()

# pyd.aerosols.remove_all()

# pyd.aerosols.aer = MicrophysicalAerosol()
# pyd.aerosols.aer.profile.fill_exponential_profile(1000,0.1)

# ilay = 0

# pyd.wavelengths = [300,532]
# pyd.settings_rt.do_radiance = True
# pyd.settings_rt.do_actinic_flux_at_observer = False
# pyd.settings_rt.do_radiance_jacobians = True
# pyd.settings_rt.jacobian_list = ['aer','aer_f','aer_r1']
# pyd.viewing_geometry.vza = 89

# pyd.run_rt()
