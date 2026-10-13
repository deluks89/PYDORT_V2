import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
sys.path.append('C:/PYDORT_V2/')
import PYDORT_V2 as pydort
from scipy.interpolate import interp1d

xnew = np.arange(0,70001,10)

#%% Temperature profile:

d = pd.read_csv('temperature_profile.dat', delimiter = '\t', index_col=0).iloc(1)[0]
d.index = d.index*1e3

x = (d.index[:-1]+d.index[1:])/2
y = d.values[1:]
x = x[::-1]
y = y[::-1]

fint = interp1d(x, y, kind = 'linear', fill_value = 'extrapolate', bounds_error = False)
ynew = fint(xnew)

data = np.vstack([xnew,ynew])
header='Altitude[m] Temperature[K]'
np.savetxt('../temperature_profile.dat', data.T, header=header)


#%%

d = pd.read_csv('ozone_profile.dat', sep='\s+', index_col=0, comment='#', header=None)[1]
d.index = d.index*1e3
x,y = d.index, d.values

fint = interp1d(x, y, kind = 'cubic', fill_value = 'extrapolate', bounds_error = False)
ynew = fint(xnew)

data = np.vstack([xnew,ynew])

header='Altitude[m] Number_concentration [cm-3]'
np.savetxt('../o3.dat', data.T, header=header)

#%% H2O profile:

d = pd.read_csv('atmmod.afglmw.100.dat', sep='\s+', index_col=0, comment='#')
d.index = d.index*1e3
d = d.sort_index()

x,y = d.index.copy(), d['h2o(cm-3)'].values.copy()

fint = interp1d(x, y, kind = 'cubic', fill_value = 'extrapolate', bounds_error = False)
ynew = fint(xnew)

data = np.vstack([xnew,ynew])

header='Layer_boundaries_altitudes[m] Layer_center_number_concentration [cm-3]'
np.savetxt('../h2o.dat', data.T, header=header)


#%% Pressure

d = pd.read_csv('pressure_profile.dat', index_col=0, comment='#')
d.index = d.index*1e3
d /= 1e2
d = d.sort_index()

x,y = d.index.copy(), d.pres.values.copy()

fint = interp1d(x, y, kind = 'cubic', fill_value = 'extrapolate', bounds_error = False)
ynew = fint(xnew)

data = np.vstack([xnew,ynew])

header='Altitudes[m] Pressure [hPa]'
np.savetxt('../pressure_profile.dat', data.T, header=header)


#%%

# Add stratospheric layer: 
# Shape based on Hendrick, 2004
# VCD based on Dirksen, 2011
gp = pydort.GasProfile()
gp.change_vertical_grid(xnew)
gp.fill_gaussian_profile(center_altitude=25000,
                          fwhm = 15000,
                          vertical_column = 2.5e15)
y = gp.values.copy()

# Add troposheric layer. 
gp.fill_exponential_profile(scale_height=1000,
                          vertical_column = 2e15)
y = y + gp.values.copy()

# Interpolate
x = gp.get_layer_middles()
fint = interp1d(x, y, kind = 'cubic', fill_value = 'extrapolate', bounds_error = False)
ynew = fint(xnew)

data = np.vstack([xnew,ynew])

header='Altitudes[m] Number_concentration [cm-3]'
np.savetxt('../no2.dat', data.T, header=header)

#%% Stratospheric aerosol profile:

d = pd.read_csv('stratospheric_aerosol.dat', skiprows = 1, index_col=0).iloc(1)[0]

d = d.sort_index()

x,y = d.index.copy(), d.values.copy()

fint = interp1d(x, y, kind = 'linear', fill_value = 0, bounds_error = False)
ynew = fint(xnew)

ynew = ynew * (532/400)**1.5

data = np.vstack([xnew,ynew])

header='Altitudes[m] Extinction@400nm[km-1]'
np.savetxt('../stratospheric_aerosol.dat', data.T, header=header)

