# """
# Example script:

#     SIMULATE RADIANCES, (BOX-)AIRMASS FACTORS, AND GAS SLANT COLUMNS

# If you are new to PYDORT_V2, start with Examples/actinic_fluxes.py. This
# script focuses on radiances and radiance Jacobians.
# """

# import matplotlib.pyplot as plt
# import numpy as np

# import PYDORT_V2 as pydort


# plt.close("all")
# pydort.REPRINT_SUBPROCESS_OUTPUT = False


# # INITIALISE A PYDORT OBJECT WITH DEFAULT VALUES
# pyd = pydort.Pydort()


# # CONFIGURE SIMULATION WAVELENGTHS
# pyd.wavelengths = np.array([360.0, 477.0])
# pyd.shift_reference_wavelengths(400.0)


# # ADAPT VIEWING GEOMETRY
# pyd.viewing_geometry.sza = 0.0
# pyd.viewing_geometry.saa = 0.0
# pyd.viewing_geometry.vza = 0.0
# pyd.viewing_geometry.vaa = 0.0
# pyd.viewing_geometry.observer_height_agl = 0.0

# #%%
# # RADIATIVE TRANSFER SETTINGS
# pyd.settings_rt.sphericity_correction = 1
# pyd.settings_rt.do_actinic_flux_at_observer = False
# pyd.settings_rt.do_radiance = True
# pyd.settings_rt.n_halfspace_streams = 4


# # REQUEST RADIANCE JACOBIANS
# pyd.settings_rt.do_radiance_jacobians = True
# pyd.settings_rt.do_actinic_flux_jacobians = False
# pyd.settings_rt.jacobian_list = ["absorption"]


# # RUN RADIATIVE TRANSFER SIMULATION
# pyd.run_rt()


# # ACCESS RADIATIVE TRANSFER RESULTS
# print("\nNormalized radiances:")
# print(pyd.results_rt.get_radiance(normalized=True))

# print("\nO4 AMF:")
# print(pyd.results_rt.get_amfs("o4"))

# print("\nNO2 SCD [molec cm-2]:")
# print(pyd.results_rt.get_scds("no2"))

# print("\nBox-AMFs:")
# print(pyd.results_rt.get_bamfs())


# # PLOT BOX-AIRMASS FACTORS
# pyd.results_rt.plot_bamfs(360.0, label="Default aerosol", color="C1")

# pyd_rayleigh = pyd.copy()
# pyd_rayleigh.aerosols.remove_all()
# pyd_rayleigh.run_rt()
# pyd_rayleigh.results_rt.plot_bamfs(360.0, new_fig=False, label="Rayleigh", color="C0")

# plt.legend()
# plt.tight_layout()
