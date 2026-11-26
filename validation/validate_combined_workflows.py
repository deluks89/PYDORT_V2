"""Spyder-style experimental routing, aerosol-Jacobian and mixed-output tests.

Finite differences here validate native derivatives; they are never used by
the library. All generated input/output files live in temporary directories.
"""
#%% Imports
from contextlib import redirect_stdout
import io
from pathlib import Path
import sys
import warnings
import numpy as np

import hashlib
import validate_limb_hybrid as h
rt=h.rt

def isolated_backends():
    hashes={name:hashlib.sha256((h.ROOT/'Fortran'/name/'DRIVER.exe').read_bytes()).hexdigest()
            for name in ['vlidort_qs','vlidort_standard']}
    return h.isolated(hashes)

def run_executable(folder):
    return float(h.native(folder,h.ROOT/'tmp/combined_workflow_native.log')['radiance'][0])

write_table=h.write_csv


#%% Routing and validation
def check_routing():
    pyd = rt.Pydort()
    pyd.settings_rt.sphericity_correction = 0
    pyd.settings_rt.validate_backend_settings()
    for height, vza in [(15000,89),(550,89),(100,0),(100,180),(15000,91),(100,90)]:
        pyd.viewing_geometry.observer_height_agl=height
        pyd.viewing_geometry.vza=vza
        assert pyd._resolve_rt_backend()['backend_mode']=='aircraft_combined'
    for invalid in (True, 1, 34, 2.5):
        pyd.settings_rt.n_aircraft_source_points = invalid
        try:
            pyd.settings_rt.validate_backend_settings()
        except AssertionError:
            continue
        raise AssertionError(f'Accepted invalid point count: {invalid}')


#%% Aerosol extinction derivatives with nonzero delta-M truncation
def check_aerosol():
    rows = []
    with isolated_backends() as work:
        pyd = rt.Pydort()
        pyd.aerosols.remove_all()
        pyd.aerosols.haze = rt.OpticalAerosol()
        pyd.change_vertical_grids(np.array([0, 100, 1000, 3000, 6000, 10000, 20000, 70000]))
        pyd.aerosols.haze.profile.values[:] = 0.02*np.exp(-pyd.aerosols.haze.profile.get_layer_middles()/3000)
        pyd.aerosols.haze.asymmetry_parameter.set_constant(0.8)
        pyd.aerosols.haze.single_scattering_absorption.set_constant(0.05)
        pyd.surface_albedo.set_constant(0.2)
        pyd.wavelengths = [477.0]
        pyd.settings_rt.sphericity_correction = 0
        pyd.settings_rt.do_radiance_jacobians = True
        pyd.settings_rt.jacobian_list = ['absorption', 'haze']
        pyd.viewing_geometry.observer_height_agl = 550
        pyd.viewing_geometry.sza = 60
        pyd.viewing_geometry.saa = 0
        pyd.viewing_geometry.vaa = 90
        pyd.viewing_geometry.vza = 89
        with redirect_stdout(io.StringIO()):
            pyd.run_rt()
        folder = work / 'vlidort_standard'
        native = np.loadtxt(folder / 'OUTPUT/jacs_radiance_haze.out').reshape(-1)
        base = float(np.loadtxt(folder / 'OUTPUT/radiances.out'))
        config = folder / 'INPUT/configuration.inp'
        lines = config.read_text().splitlines()
        lines[14] = 'False ! independent forward check'
        config.write_text('\n'.join(lines)+'\n')
        np.testing.assert_allclose(run_executable(folder), base, rtol=1e-9)
        path = folder / 'INPUT/aerosol_extinction_profiles.inp'
        header = '\n'.join(path.read_text().splitlines()[:3])
        values = np.loadtxt(path, skiprows=3, ndmin=2)
        for layer in (2, 4, 6):
            for delta in (1e-6, 1e-7):
                intensities = []
                for multiplier in (1, 2, 3):
                    changed = values.copy()
                    changed[layer, 0] += multiplier*delta
                    np.savetxt(path, changed, header=header, comments='')
                    intensities.append(run_executable(folder))
                fd = (-5*intensities[0]+8*intensities[1]-3*intensities[2])/(2*delta)
                passed = bool(np.isclose(native[layer], fd, rtol=1e-4, atol=1e-7))
                rows.append(dict(layer=layer, step=delta, native=native[layer], finite_difference=fd, passed=passed))
                print(rows[-1], flush=True)
        # Public runs regenerate files and exercise state resets after an MS run.
        pyd.wavelengths = [477.0, 500.0]
        pyd.settings_rt.do_actinic_flux_at_observer = True
        pyd.settings_rt.do_actinic_flux_jacobians = True
        pyd.settings_rt.do_irradiance_at_observer = True
        pyd.settings_rt.do_irradiance_jacobians = True
        with redirect_stdout(io.StringIO()):
            pyd.run_rt()
        assert np.all(np.isfinite(pyd.results_rt.get_radiance(normalized=True)))
        assert np.all(np.isfinite(pyd.results_rt.get_actinic_flux(normalized=True)))
        assert np.all(np.isfinite(pyd.results_rt.get_flux_bamfs()))
        assert np.all(np.isfinite(pyd.results_rt.get_irradiance(normalized=True)))
        assert np.all(np.isfinite(pyd.results_rt.get_irradiance_jacobians('haze',normalized=True)))
        pyd.settings_rt.do_single_scattering_only = True
        pyd.settings_rt.do_actinic_flux_at_observer = False
        pyd.settings_rt.do_actinic_flux_jacobians = False
        pyd.settings_rt.do_irradiance_at_observer = False
        pyd.settings_rt.do_irradiance_jacobians = False
        with redirect_stdout(io.StringIO()):
            pyd.run_rt(run_aerosol_models=False)
        assert np.all(np.isfinite(pyd.results_rt.get_radiance(normalized=True)))
        # Direct FO also covers upward and surface-intersecting paths.
        for vza in [0,180]:
            pyd.viewing_geometry.observer_height_agl=100
            pyd.viewing_geometry.vza=vza
            pyd.settings_rt.sphericity_correction=0
            with redirect_stdout(io.StringIO()):
                pyd.run_rt(run_aerosol_models=False)
            ac=np.loadtxt(work/'vlidort_standard/OUTPUT/radiances.out')
            ac_jac=np.loadtxt(work/'vlidort_standard/OUTPUT/jacs_radiance_haze.out')
            pyd.settings_rt.sphericity_correction=4
            with redirect_stdout(io.StringIO()):
                pyd.run_rt(run_aerosol_models=False)
            qs=np.loadtxt(work/'vlidort_qs/OUTPUT/radiances.out')
            qs_jac=np.loadtxt(work/'vlidort_qs/OUTPUT/jacs_radiance_haze.out')
            np.testing.assert_allclose(ac,qs,rtol=1e-9,atol=1e-12)
            np.testing.assert_allclose(ac_jac,qs_jac,rtol=2e-7,atol=1e-9)
    output = WORKFLOW_OUTPUT
    output.mkdir(parents=True,exist_ok=True)
    write_table(output / 'aerosol.txt', rows)
    assert all(row['passed'] for row in rows), 'Aerosol derivative check failed'


#%% Run
WORKFLOW_OUTPUT = h.ROOT / 'Examples/testing/limb_hybrid/combined_workflows'

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=WORKFLOW_OUTPUT)
    WORKFLOW_OUTPUT = parser.parse_args().output
    check_routing()
    check_aerosol()
    print('PASS: combined AC routing, aerosol derivatives, wavelengths, actinic flux, irradiance and single scattering.')
