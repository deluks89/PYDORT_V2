"""
Example script:

    SIMULATE OBSERVER-ALTITUDE ACTINIC FLUXES WITH PYDORT_V2

This example shows how to configure a small atmosphere, run radiative transfer
at selected wavelengths, and post-process actinic fluxes at the observer
altitude.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import PYDORT_V2 as pydort


plt.close("all")
pydort.DEFAULT_VERTICAL_RESOLUTION = "HR"
pydort.REPRINT_SUBPROCESS_OUTPUT = True

EXAMPLE_DATA = Path(__file__).resolve().parent / "example_data"


# INITIALISE A PYDORT OBJECT WITH DEFAULT VALUES
pyd = pydort.Pydort()
pyd.wavelengths = np.array([310.0, 340.0, 400.0])
pyd.shift_reference_wavelengths(400.0)


# REMOVE AND ADD AEROSOL SPECIES
print("\nDefault aerosols:")
print(pyd.aerosols)

pyd.aerosols.remove_all()
pyd.aerosols.plume = pydort.OpticalAerosol()
pyd.aerosols.dust_like = pydort.MicrophysicalAerosol()


# ADAPT VIEWING GEOMETRY
pyd.viewing_geometry.sza = 30.0
pyd.viewing_geometry.saa = 0.0
pyd.viewing_geometry.vza = 0.0
pyd.viewing_geometry.vaa = 0.0
pyd.viewing_geometry.observer_height_agl = 3000.0


# MODIFY VERTICAL PROFILES
pyd.pressure_profile.values *= 1.01

pyd.aerosols.dust_like.profile_wavelength = 400.0
pyd.aerosols.dust_like.profile.fill_gaussian_profile(
    center_altitude=3000.0,
    fwhm=1000.0,
    optical_depth=0.15,
)

alts, vals = np.genfromtxt(EXAMPLE_DATA / "plume_profile_example_data.dat").T
pyd.aerosols.plume.profile.interpolate_from(alts, vals, fill_value=0.0)

pyd.aerosols.dust_like.profile.plot(scaling=10.0, label="dust_like x10")
pyd.aerosols.plume.profile.plot(new_fig=False, label="plume")
plt.legend()


# MODIFY MICROPHYSICAL AND SPECTRAL AEROSOL PARAMETERS
psd = pyd.aerosols.dust_like.psd
psd.f1 = 0.9
psd.f2 = 0.1
psd.r2 = 0.3
psd.sig2 = 0.5

pyd.aerosols.plume.single_scattering_absorption.reference_wavelength = 400.0
pyd.aerosols.plume.single_scattering_absorption.value = 0.1
pyd.aerosols.plume.single_scattering_absorption.angstrom_parameters = [2.0, 1.0]

wls, vals = np.genfromtxt(EXAMPLE_DATA / "refractive_index_example_data.dat").T
pyd.aerosols.dust_like.imag_refractive_index.get_properties_from_data_fit(
    wls,
    vals,
    deg=3,
    plot=True,
)


# MODIFY TRACE GASES
print("\nTrace gases:")
print(pyd.trace_gases)
pyd.trace_gases.no2.profile.fill_fermi_profile(1500.0, 200.0, 5e16)
pyd.trace_gases.o3.profile.set_vertical_column(1e19)


# OPTIONAL: ELEVATE GROUND AND REDUCE VERTICAL GRID RESOLUTION
pyd.elevate_ground(1000.0)
pyd.change_vertical_grids(pydort.expgrid(40, 60000.0, 100.0))


# RADIATIVE TRANSFER SETTINGS
pyd.settings_rt.sphericity_correction = 1
pyd.settings_rt.do_actinic_flux_at_observer = True
pyd.settings_rt.do_radiance = False
pyd.settings_rt.n_halfspace_streams = 2


# RUN RADIATIVE TRANSFER SIMULATION
print("\nRunning the simulation. Here comes the Fortran console output:")
pyd.run_rt(run_aerosol_models=True)


# ACCESS OBSERVER-ALTITUDE ACTINIC FLUXES
actinic_flux = pyd.results_rt.get_actinic_flux()
actinic_flux_normalized = pyd.results_rt.get_actinic_flux(normalized=True)

print("\nActinic flux [photons s-1 cm-2 nm-1]:")
print(actinic_flux)

plt.figure()
actinic_flux.plot(marker="o")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Actinic flux [photons s-1 cm-2 nm-1]")
plt.tight_layout()

plt.figure()
actinic_flux_normalized.plot(marker="o")
plt.xlabel("Wavelength [nm]")
plt.ylabel("Normalized actinic flux")
plt.tight_layout()
