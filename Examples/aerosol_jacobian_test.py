"""
Example script:

    COMPARE AEROSOL EXTINCTION JACOBIANS WITH FINITE DIFFERENCES

Radiance Jacobians are available through all radiance backends.
"""

import numpy as np

import PYDORT_V2 as pydort


pydort.DEFAULT_VERTICAL_RESOLUTION = "LR"

DX = 1e-7
LAYER_INDEX = 0
WAVELENGTHS = np.array([310.0, 400.0, 532.0])
AODS = [0.05, 0.1, 0.5, 1.0]


pyd = pydort.Pydort()
pyd.wavelengths = WAVELENGTHS
pyd.shift_reference_wavelengths(400.0)

pyd.settings_rt.sphericity_correction = 4
pyd.settings_rt.do_actinic_flux_at_observer = False
pyd.settings_rt.do_radiance = True
pyd.settings_rt.do_radiance_jacobians = True
pyd.settings_rt.jacobian_list = ["aer"]

pyd.aerosols.remove_all()
pyd.aerosols.aer = pydort.MicrophysicalAerosol()
pyd.aerosols.aer.imag_refractive_index.set_constant(0.5)
pyd.aerosols.run_models(WAVELENGTHS)



for aod in AODS:
    print("\nAOD:", aod)
    pyd.aerosols.aer.profile.fill_exponential_profile(1000.0, aod)

    # Disturbed finite-difference simulation.
    pyd_disturbed = pyd.copy()
    pyd_disturbed.settings_rt.do_radiance_jacobians = False
    pyd_disturbed.aerosols.aer.profile.values[LAYER_INDEX] += DX
    pyd_disturbed.run_rt(run_aerosol_models=False)
    radiance_disturbed = pyd_disturbed.results_rt.get_radiance(normalized=True)

    # Baseline simulation with analytic Jacobians.
    pyd.settings_rt.do_radiance_jacobians = True
    pyd.run_rt(run_aerosol_models=False)
    radiance_base = pyd.results_rt.get_radiance(normalized=True)
    jacobian = pyd.results_rt.get_jacobians("aer").iloc[LAYER_INDEX]

    finite_difference = (radiance_disturbed - radiance_base) / DX
    print(jacobian / finite_difference)
