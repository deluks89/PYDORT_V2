"""Public-API regression checks in isolated Fortran working directories.

Run as a Spyder script or with python Examples/testing/validate_rt_workflows.py.
Existing simulation inputs and comparison results are not changed.
"""

#%% Imports
import contextlib
import io
import sys
import warnings
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from validate_native_jacobians import isolated_backends
import PYDORT_V2 as rt

#%% Geometry, QS fallback, and spherical first-order checks
with isolated_backends():
    rt.REPRINT_SUBPROCESS_OUTPUT = False
    pyd = rt.Pydort()
    pyd.wavelengths = [477.0]
    pyd.aerosols.remove_all()
    for gas in pyd.trace_gases.to_dictionary().values():
        gas.xs.raw_data[:, 1:] = 0
        if gas.xs.smoothed_data is not None:
            gas.xs.smoothed_data[:, 1:] = 0
    pyd.surface_albedo.set_constant(0)
    pyd.settings_rt.do_radiance_jacobians = True
    pyd.settings_rt.jacobian_list = ['absorption']
    pyd.settings_rt.n_halfspace_streams = 4
    pyd.viewing_geometry.sza = 60
    pyd.viewing_geometry.saa = 0
    pyd.viewing_geometry.vaa = 90
    pyd.viewing_geometry.vza = 89
    results = {}
    for altitude, correction, ss in [(100,1,False), (100,2,False), (100,4,False),
                                    (15000,1,False), (15000,2,False), (15000,4,False),
                                    (15000,1,True), (15000,4,True)]:
        pyd.viewing_geometry.observer_height_agl = altitude
        pyd.settings_rt.sphericity_correction = correction
        pyd.settings_rt.do_single_scattering_only = ss
        with warnings.catch_warnings(record=True) as caught, contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter('always')
            pyd.run_rt(run_aerosol_models=False)
        radiance = float(pyd.results_rt.get_radiance(normalized=True).iloc[0])
        scd = float(pyd.results_rt.get_scds('o4').iloc[0])
        assert np.isfinite(radiance) and np.isfinite(scd)
        results[(altitude,correction,ss)] = (radiance, scd)
        print(altitude, correction, ss, radiance, scd, [str(w.message) for w in caught], flush=True)
    np.testing.assert_equal(results[(15000,2,False)], results[(15000,4,False)])
    np.testing.assert_allclose(results[(15000,1,True)], results[(15000,4,True)], rtol=5e-3)
    # Repeated wavelengths exercise the radiance/flux configuration resets.
    pyd.settings_rt.sphericity_correction = 1
    pyd.settings_rt.do_single_scattering_only = False
    pyd.wavelengths = [477., 500.]
    pyd.settings_rt.do_actinic_flux_at_observer = True
    pyd.settings_rt.do_actinic_flux_jacobians = True
    with contextlib.redirect_stdout(io.StringIO()):
        pyd.run_rt(run_aerosol_models=False)
    values = np.asarray(pyd.results_rt.get_radiance(normalized=True))
    assert values.shape == (2,) and np.all(np.isfinite(values))
    assert np.all(np.isfinite(pyd.results_rt.get_actinic_flux(normalized=True)))
    assert np.all(np.isfinite(pyd.results_rt.get_flux_bamfs()))
    print('multiple wavelengths and observer flux/Jacobian run', values, flush=True)
    pyd.settings_rt.sphericity_correction = 2
    pyd.viewing_geometry.observer_height_agl = 100
    with contextlib.redirect_stdout(io.StringIO()):
        pyd.run_rt(run_aerosol_models=False)
    assert np.all(np.isfinite(pyd.results_rt.get_radiance(normalized=True)))
    assert np.all(np.isfinite(pyd.results_rt.get_actinic_flux(normalized=True)))
    assert np.all(np.isfinite(pyd.results_rt.get_flux_bamfs()))
    print('PASS: public radiance, flux, native-Jacobian and fallback workflows.', flush=True)
