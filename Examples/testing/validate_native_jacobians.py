"""Validate native absorption Jacobians without replacing them at runtime.

Spyder: run the file for two short checks, or set FULL_SUITE below.
CLI: python Examples/testing/validate_native_jacobians.py --suite
All Fortran inputs/outputs are isolated; existing comparison results are kept.
"""

#%% Imports and configuration
from pathlib import Path
from contextlib import contextmanager, redirect_stdout
import argparse
import csv
import io
import shutil
import subprocess
import sys
import tempfile
import time
import warnings

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
import PYDORT_V2 as rt

FULL_SUITE = False
OUTPUT_DIR = Path(__file__).resolve().parent / 'native_jacobian_validation'
FD_STEPS = (1e-5, 1e-6)
RTOL = 1e-4
ATOL = 1e-7


#%% Isolated executable runs
@contextmanager
def isolated_backends():
    original = rt.QS_BACKEND_PATH, rt.STANDARD_BACKEND_PATH
    with tempfile.TemporaryDirectory(prefix='pydort_native_') as tmp:
        work = Path(tmp)
        for backend in ('vlidort_qs', 'vlidort_standard'):
            folder = work / backend
            for subdir in ('INPUT', 'INPUT/aerosol_moments', 'INPUT/xs',
                           'INPUT/aerosol_model_jacobians', 'OUTPUT'):
                (folder / subdir).mkdir(parents=True, exist_ok=True)
            shutil.copy2(PROJECT_ROOT / 'Fortran' / backend / 'DRIVER.exe', folder)
        rt.QS_BACKEND_PATH = (work / 'vlidort_qs').as_posix() + '/'
        rt.STANDARD_BACKEND_PATH = (work / 'vlidort_standard').as_posix() + '/'
        try:
            yield work
        finally:
            rt.QS_BACKEND_PATH, rt.STANDARD_BACKEND_PATH = original


def run_executable(folder):
    proc = subprocess.run([str(folder / 'DRIVER.exe')], cwd=folder,
                          capture_output=True, text=True, timeout=180)
    if proc.returncode:
        raise RuntimeError(proc.stdout + proc.stderr)
    radiance = float(np.loadtxt(folder / 'OUTPUT/radiances.out'))
    if not np.isfinite(radiance) or radiance <= 0:
        raise AssertionError(f'Nonpositive/nonfinite radiance: {radiance}')
    return radiance


#%% Native versus finite-difference verification
def check_case(correction, sza, altitude, vza=89, full_grid=False, vector=False,
               background=0.0, verify_derivatives=True):
    with isolated_backends() as work:
        pyd = rt.Pydort()
        if not full_grid:
            pyd.change_vertical_grids(np.array([0, 100, 1000, 3000, 6000, 10000, 20000, 70000]))
        pyd.wavelengths = [477.0]
        pyd.aerosols.remove_all()
        for gas in pyd.trace_gases.to_dictionary().values():
            gas.xs.raw_data[:, 1:] = 0
            if gas.xs.smoothed_data is not None:
                gas.xs.smoothed_data[:, 1:] = 0
        pyd.surface_albedo.set_constant(0)
        pyd.viewing_geometry.sza = sza
        pyd.viewing_geometry.saa = 0
        pyd.viewing_geometry.vaa = 90
        pyd.viewing_geometry.vza = vza
        pyd.viewing_geometry.observer_height_agl = altitude
        pyd.settings_rt.sphericity_correction = correction
        pyd.settings_rt.n_halfspace_streams = 4
        pyd.settings_rt.n_full_sphericity_fine_layers = 4
        pyd.settings_rt.do_vectorized_simulation = vector
        pyd.settings_rt.do_radiance = True
        pyd.settings_rt.do_radiance_jacobians = True
        pyd.settings_rt.do_actinic_flux_at_observer = False
        pyd.settings_rt.do_actinic_flux_jacobians = False
        pyd.settings_rt.jacobian_list = ['absorption']
        resolved = pyd._resolve_rt_backend()
        pyd._active_rt_backend_settings = resolved
        with redirect_stdout(io.StringIO()):
            pyd._generate_rt_input(run_aerosol_models=False)
        is_qs = resolved['model'] == rt.RadiativeTransferModel.QS.value
        folder = work / ('vlidort_qs' if is_qs else 'vlidort_standard')
        if not is_qs:
            pyd._stage_standard_backend_input()
        inp, out = folder / 'INPUT', folder / 'OUTPUT'

        # An otherwise transparent O3 tracer imposes known absorption OD.
        gas_file = inp / 'gas_concentration_profiles.inp'
        gas_header = '\n'.join(gas_file.read_text().splitlines()[:3])
        gases = np.loadtxt(gas_file, skiprows=3, ndmin=2)
        names = (inp / 'gas_names.inp').read_text().splitlines()[1:]
        tracer = names.index('o3')
        heights = np.loadtxt(inp / 'vertical_grid.inp', skiprows=2)
        thickness_cm = -np.diff(heights) * 100
        xs_file = inp / 'xs/o3.inp'
        xs_header = '\n'.join(xs_file.read_text().splitlines()[:3])
        xs = np.loadtxt(xs_file, skiprows=3, ndmin=2)
        xs[:, 0] = 1e-20
        np.savetxt(xs_file, xs, header=xs_header, comments='')

        def set_absorption(tau):
            gases[:, tracer] = np.asarray(tau) / (1e-20 * thickness_cm)
            np.savetxt(gas_file, gases, header=gas_header, comments='')

        set_absorption(background)
        start = time.perf_counter()
        baseline = run_executable(folder)
        seconds = time.perf_counter() - start
        jac = np.loadtxt(out / 'jacs_radiance_absorption.out').reshape(-1)
        assert len(jac) == len(thickness_cm) and np.all(np.isfinite(jac))
        config_file = inp / 'configuration.inp'
        config = config_file.read_text().splitlines()
        config[14] = 'False ! radiance Jacobians disabled for consistency check'
        config_file.write_text('\n'.join(config) + '\n')
        regular = run_executable(folder)
        np.testing.assert_allclose(regular, baseline, rtol=1e-7, atol=1e-11)
        # Keep Jacobians disabled for the independent forward perturbations;
        # the native derivative has already been saved above.

        columns = pyd.get_o4_profile()[::-1] * thickness_cm
        scd = float(-jac @ columns / baseline)
        summary = dict(correction=correction, backend=resolved['backend_mode'],
                       fallback_reason=resolved.get('fallback_reason', ''), sza_deg=sza,
                       altitude_m=altitude, vza_deg=vza, vector=vector, layers=len(jac),
                       background_layer_od=background, radiance=baseline,
                       radiance_without_jacobians=regular, o4_scd=scd, seconds=seconds)
        print(summary, flush=True)
        rows = []
        if not verify_derivatives:
            return summary, rows
        if is_qs:
            raise ValueError('Derivative regression must exercise a native Standard/AC backend, not a QS fallback.')

        layers = list(range(len(jac))) if len(jac) < 10 else [0, 10, 30, 50, len(jac)-2, len(jac)-1]
        weights = [(f'layer_{layer}', np.eye(1, len(jac), layer).reshape(-1)) for layer in layers]
        weights.append(('o4_weighted', columns / columns.max()))
        for name, weight in weights:
            native = float(jac @ weight)
            for step in FD_STEPS:
                values = []
                for multiplier in (1, 2, 3):
                    set_absorption(background + multiplier * step * weight)
                    values.append(run_executable(folder))
                # Extrapolate from positive ODs, avoiding the SSA limiter's
                # clipping kink at zero absorption. This is a TEST ONLY.
                fd = (-5*values[0] + 8*values[1] - 3*values[2]) / (2*step)
                passed = bool(np.isclose(native, fd, rtol=RTOL, atol=ATOL))
                row = dict(correction=correction, backend=resolved['backend_mode'],
                           sza_deg=sza, altitude_m=altitude, vza_deg=vza, vector=vector,
                           background_layer_od=background, quantity=name, step=step,
                           native=native, finite_difference=fd, passed=passed)
                rows.append(row)
                print(row, flush=True)
            set_absorption(background)
        summary['derivatives_passed'] = all(row['passed'] for row in rows)
        return summary, rows


#%% Routing regression (no Fortran execution)
def check_routing():
    pyd = rt.Pydort()
    pyd.viewing_geometry.vza = 89
    pyd.settings_rt.sphericity_correction = 2
    for altitude in (100, 500):
        pyd.viewing_geometry.observer_height_agl = altitude
        assert pyd._resolve_rt_backend(warn=False)['backend_mode'] == 'aircraft_spherical'
    for altitude in (1000, 8100, 15000, 150):
        pyd.viewing_geometry.observer_height_agl = altitude
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            resolved = pyd._resolve_rt_backend()
        assert resolved['model'] == rt.RadiativeTransferModel.QS.value
        assert any('VLIDORT-QS' in str(item.message) for item in caught)


#%% Run and store separate text tables
def write_table(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite', action='store_true', default=FULL_SUITE)
    parser.add_argument('--correction', type=int, choices=(1, 2))
    parser.add_argument('--sza', type=float, default=60)
    parser.add_argument('--altitude', type=float, default=100)
    parser.add_argument('--vza', type=float, default=89)
    parser.add_argument('--full-grid', action='store_true')
    parser.add_argument('--vector', action='store_true')
    parser.add_argument('--background', type=float, default=0)
    parser.add_argument('--output', type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    check_routing()
    if args.suite:
        cases = [(1, 0, 100, 89), (1, 60, 8100, 89), (1, 85, 15000, 89),
                 (2, 0, 100, 89), (2, 60, 500, 89), (2, 85, 15000, 85),
                 (1, 60, 550, 89)]
    else:
        cases = [(c, args.sza, args.altitude, args.vza) for c in
                 ([args.correction] if args.correction else [1, 2])]
    args.output.mkdir(parents=True, exist_ok=True)
    summaries, derivatives = [], []
    for correction, sza, altitude, vza in cases:
        summary, rows = check_case(correction, sza, altitude, vza,
                                   full_grid=args.suite or args.full_grid,
                                   vector=args.vector, background=args.background)
        summaries.append(summary)
        derivatives.extend(rows)
        write_table(args.output / 'summary.txt', summaries)
        write_table(args.output / 'derivatives.txt', derivatives)
    if not all(row['passed'] for row in derivatives):
        raise AssertionError('Native Jacobian regression failed; inspect derivatives.txt.')
    print(f'PASS: {len(summaries)} cases, {len(derivatives)} derivative checks. Results: {args.output}')


if __name__ == '__main__':
    main()
