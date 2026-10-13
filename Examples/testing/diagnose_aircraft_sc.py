"""Compare valid native AC geometries with Standard and QS, separating FO/MS.

Spyder: edit CASES and run this file. CLI: add --suite for the geometry sweep.
The optional Fortran diagnostics require PYDORT_RT_DIAGNOSTICS=1 (set here).
"""

#%% Imports and cases
import argparse
from contextlib import redirect_stdout
import io
import os
from pathlib import Path
import shutil
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import PYDORT_V2 as rt
from validate_native_jacobians import isolated_backends, write_table

CASES = [(100, 60, 89)]
OUTPUT = Path(__file__).resolve().parent / 'aircraft_diagnosis'


#%% Simulations with isolated inputs and outputs
def compare_case(altitude, sza, vza, output, streams=4, fine_layers=4, vector=False):
    rows = []
    with isolated_backends() as work:
        pyd = rt.Pydort()
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
        pyd.settings_rt.n_halfspace_streams = streams
        pyd.settings_rt.n_full_sphericity_fine_layers = fine_layers
        pyd.settings_rt.do_vectorized_simulation = vector
        pyd.settings_rt.do_radiance_jacobians = True
        pyd.settings_rt.jacobian_list = ['absorption']
        pyd.settings_rt.do_actinic_flux_at_observer = False
        pyd.settings_rt.do_actinic_flux_jacobians = False
        pyd.settings_rt.sphericity_correction = 2
        assert pyd._resolve_rt_backend()['backend_mode'] == 'aircraft_spherical'
        for correction, ss in [(1, False), (2, False), (4, False), (1, True), (4, True)]:
            pyd.settings_rt.sphericity_correction = correction
            pyd.settings_rt.do_single_scattering_only = ss
            with redirect_stdout(io.StringIO()):
                pyd.run_rt(run_aerosol_models=False)
            radiance = float(pyd.results_rt.get_radiance(normalized=True).iloc[0])
            scd = float(pyd.results_rt.get_scds('o4').iloc[0])
            row = dict(altitude_m=altitude, sza_deg=sza, vza_deg=vza,
                       correction=correction, single_scattering_only=ss,
                       streams=streams, fine_layers=fine_layers, vector=vector,
                       radiance=radiance, o4_scd=scd)
            case_dir = output / f'h{altitude:g}_sza{sza:g}_vza{vza:g}_c{correction}_ss{int(ss)}'
            case_dir.mkdir(parents=True, exist_ok=True)
            backend = work / ('vlidort_qs' if correction == 4 else 'vlidort_standard')
            if correction != 4 and not ss:
                components = np.loadtxt(backend / 'OUTPUT/rt_components.out')
                for key, value in zip(['original_fo', 'original_ms', 'original_total',
                                       'corrected_fo', 'corrected_ms', 'corrected_total'], components):
                    row[key] = float(value)
                for diagnostic in (backend / 'OUTPUT').glob('rt_*.out'):
                    shutil.copy2(diagnostic, case_dir)
            for name in ['radiances.out', 'jacs_radiance_absorption.out']:
                shutil.copy2(backend / 'OUTPUT' / name, case_dir)
            rows.append(row)
            print(row, flush=True)
    return rows


#%% Run and save text tables
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite', action='store_true')
    parser.add_argument('--output', type=Path, default=OUTPUT)
    parser.add_argument('--streams', type=int, default=4)
    parser.add_argument('--fine-layers', type=int, default=4)
    parser.add_argument('--vector', action='store_true')
    parser.add_argument('--samples', action='store_true', help='Evaluate native MS sources at each layer midpoint.')
    args = parser.parse_args()
    cases = CASES if not args.suite else [(h, s, 89) for h in (100, 500) for s in (0, 60, 85)]
    if args.suite:
        cases += [(15000, 60, 85), (100, 60, 70), (100, 60, 160)]
    args.output.mkdir(parents=True, exist_ok=True)
    saved = os.environ.get('PYDORT_RT_DIAGNOSTICS')
    saved_samples = os.environ.get('PYDORT_RT_SAMPLE_SOURCES')
    os.environ['PYDORT_RT_DIAGNOSTICS'] = '1'
    os.environ['PYDORT_RT_SAMPLE_SOURCES'] = '1' if args.samples else '0'
    try:
        rows = []
        for altitude, sza, vza in cases:
            rows.extend(compare_case(altitude, sza, vza, args.output,
                                     args.streams, args.fine_layers, args.vector))
            write_table(args.output / 'components.txt', rows)
    finally:
        if saved is None:
            os.environ.pop('PYDORT_RT_DIAGNOSTICS', None)
        else:
            os.environ['PYDORT_RT_DIAGNOSTICS'] = saved
        if saved_samples is None:
            os.environ.pop('PYDORT_RT_SAMPLE_SOURCES', None)
        else:
            os.environ['PYDORT_RT_SAMPLE_SOURCES'] = saved_samples


if __name__ == '__main__':
    main()
