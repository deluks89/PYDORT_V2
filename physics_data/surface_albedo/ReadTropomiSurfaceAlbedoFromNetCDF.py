import numpy as np
import h5py
from scipy.interpolate import interpn

path_to_netcdf = 'C:/VLIDORT-PC/physics_data/''surface_albedo/TROPOMI_Sentinel-5P_0125x0125_surface_DLER_v2.0.nc'

def get_tropomi_albedo(lat, lon, month):
    """
    Reads and linearly interpolates TROPOMI-2 Lambertian equivalent reflectance (LER).
    Not working close to Earths poles and lon = 180.
    Surface albedos outside valid data range raises error.
    ----------
    lat : float, Latitude, -89 < lat < 89
    lon : float, Longitude, -179 < lon < 179
    month : integer, e.g. 8 for August
    """
    d = h5py.File(path_to_netcdf,'r')
    lats = np.array(d.get('latitude'))
    lons = np.array(d.get('longitude'))
    wls = np.array(d.get('wavelength'))
    # Prepare slices for efficient data reading:
    sl_lat = slice(np.argmax(lats > lat-1), np.argmax(lats > lat+1))
    sl_lon = slice(np.argmax(lons > lon-1), np.argmax(lons > lon+1))
    lats = lats[sl_lat]
    lons = lons[sl_lon]
    albs = np.array(d['minimum_LER_clear'][month, :, sl_lon, sl_lat]) # dims: month, wavelength, lon, lat
    # 2D-interpolation to exact location for all wavelengths:
    newx = [(wl,lon,lat) for wl in wls]
    albs = interpn((wls, lons, lats), albs, newx)
    return wls, albs

ler = get_tropomi_albedo(40.8079,-112.088,8) 
