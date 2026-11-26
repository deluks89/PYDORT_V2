"""Research harness: QS first order + Standard sparse spherical MS.

Public backend routing is intentionally unchanged. Native mode
``aircraft_limb_ms`` returns ONLY multiple scattering. Run from any directory:
    python validation/validate_limb_hybrid.py --quick
Outputs go to the ignored Examples/testing/limb_hybrid directory.
"""
from contextlib import contextmanager, redirect_stdout
from pathlib import Path
import argparse
import csv
import io
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import PYDORT_V2 as rt


@contextmanager
def isolated(expected_hashes):
    original = rt.QS_BACKEND_PATH, rt.STANDARD_BACKEND_PATH
    with tempfile.TemporaryDirectory(prefix='limb_', dir=ROOT / 'tmp') as temp:
        work = Path(temp)
        for name in ('vlidort_qs', 'vlidort_standard'):
            path = work / name
            for sub in ('INPUT/aerosol_moments', 'INPUT/aerosol_model_jacobians', 'INPUT/xs', 'OUTPUT'):
                (path / sub).mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / 'Fortran' / name / 'DRIVER.exe', path)
            if hashlib.sha256((path / 'DRIVER.exe').read_bytes()).hexdigest() != expected_hashes[name]:
                raise RuntimeError('A native executable changed during validation; rerun after the build finishes.')
        rt.QS_BACKEND_PATH = (work / 'vlidort_qs').as_posix() + '/'
        rt.STANDARD_BACKEND_PATH = (work / 'vlidort_standard').as_posix() + '/'
        try:
            yield work
        finally:
            rt.QS_BACKEND_PATH, rt.STANDARD_BACKEND_PATH = original


def native(folder, log):
    start = time.perf_counter()
    proc = subprocess.run([str(folder / 'DRIVER.exe')], cwd=folder,
                          capture_output=True, text=True, timeout=900)
    log.write_text(proc.stdout + proc.stderr)
    if proc.returncode:
        raise RuntimeError(f'Native calculation failed; see {log}:\n{proc.stdout[-2000:]}{proc.stderr}')
    output = folder / 'OUTPUT'
    result = dict(radiance=np.loadtxt(output / 'radiances.out').reshape(-1),
                  seconds=time.perf_counter()-start)
    for label in ('absorption', 'haze'):
        path = output / f'jacs_radiance_{label}.out'
        if path.exists():
            result[label] = np.loadtxt(path).reshape(-1)
    for key, value in result.items():
        if not np.all(np.isfinite(value)):
            raise AssertionError(f'Nonfinite {key}; see {log}')
    return result


def scene(height, vza, raa, vector, fine=False):
    p = rt.Pydort()
    p.aerosols.remove_all()
    p.aerosols.haze = rt.OpticalAerosol()
    grid = np.unique(np.r_[np.arange(0, 20001, 1000 if fine else 2000), 25000, 30000, 40000, 55000, 70000])
    p.change_vertical_grids(grid)
    aerosol = p.aerosols.haze
    aerosol.profile.values[:] = np.exp(-aerosol.profile.get_layer_middles()/2000)
    aerosol.profile.values[:] *= 0.1 / np.sum(aerosol.profile.values*np.diff(grid)/1000)
    aerosol.asymmetry_parameter.set_constant(0.7)
    aerosol.single_scattering_absorption.set_constant(0.05)
    for gas in p.trace_gases.to_dictionary().values():
        gas.xs.raw_data[:, 1:] = 0
        if gas.xs.smoothed_data is not None:
            gas.xs.smoothed_data[:, 1:] = 0
    p.wavelengths = [477.0]
    p.surface_albedo.set_constant(0.2)
    p.viewing_geometry.observer_height_agl = height
    p.viewing_geometry.sza = 60
    p.viewing_geometry.saa = 0
    p.viewing_geometry.vaa = raa
    p.viewing_geometry.vza = vza
    p.settings_rt.sphericity_correction = 4
    p.settings_rt.n_halfspace_streams = 4
    # QS SS-only disables delta-M truncation. Use the same unscaled optical
    # convention in both MS calculations so FO+MS is a valid decomposition.
    p.settings_rt.do_deltam_phase_function_scaling = False
    p.settings_rt.n_full_sphericity_fine_layers = 4
    p.settings_rt.do_vectorized_simulation = vector
    p.settings_rt.do_radiance_jacobians = True
    p.settings_rt.do_actinic_flux_at_observer = False
    p.settings_rt.do_actinic_flux_jacobians = False
    p.settings_rt.jacobian_list = ['absorption', 'haze']
    p.aerosols.run_models(p.wavelengths)
    return p


def generate(p, work, mode, points=9, ss=False):
    p.settings_rt.do_single_scattering_only = ss
    p._active_rt_backend_settings = dict(model='qs', backend_mode=mode,
                                        sphericity_interpolation_points=points)
    with redirect_stdout(io.StringIO()):
        p._generate_rt_input(run_aerosol_models=False)
    p._active_rt_backend_settings = None
    if mode in ('aircraft_limb_ms', 'aircraft_limb_raw'):
        p._stage_standard_backend_input()
        return work / 'vlidort_standard'
    return work / 'vlidort_qs'


def finite_difference(folder, base, output):
    output.mkdir(parents=True, exist_ok=True)
    config = folder / 'INPUT/configuration.inp'
    lines = config.read_text().splitlines()
    lines[14] = 'False ! independent forward finite difference'
    config.write_text('\n'.join(lines)+'\n')
    path = folder / 'INPUT/aerosol_extinction_profiles.inp'
    header = '\n'.join(path.read_text().splitlines()[:3])
    values = np.loadtxt(path, skiprows=3, ndmin=2)
    rows = []
    # Include aerosol-rich and nearly Rayleigh layers, above and below tangent.
    # Also perturb the observer's parent layer to test sublayer-Jacobian mapping.
    grid = np.loadtxt(folder / 'INPUT/vertical_grid.inp', skiprows=2).reshape(-1)
    height = float(lines[9].split()[0])
    containing = np.flatnonzero((grid[:-1] >= height) & (grid[1:] < height))
    layers = sorted(set([3, len(values)-3, len(values)-1] + containing.tolist()))
    for layer in layers:
        for step in (1e-6, 1e-7):
            intensities = []
            for multiplier in (1, 2, 3):
                changed = values.copy()
                changed[layer, 0] += multiplier*step
                np.savetxt(path, changed, header=header, comments='')
                intensities.append(native(folder, output / 'fd_latest.log')['radiance'][0])
            fd = (-5*intensities[0]+8*intensities[1]-3*intensities[2])/(2*step)
            analytic = base['haze'][layer]
            rows.append(dict(layer_top_down=layer, step=step, native=analytic, finite_difference=fd,
                             relative_error=(analytic-fd)/max(abs(fd), 1e-12),
                             passed=bool(np.isclose(analytic, fd, rtol=2e-3, atol=1e-7))))
    write_csv(output / 'finite_difference.csv', rows)
    assert all(row['passed'] for row in rows), rows


def check_guards(folder, output):
    path = folder / 'INPUT/configuration.inp'
    original = path.read_text().splitlines()
    checks = [('flux', 3, 'True', 'radiance-only'),
              ('single_scattering', 0, 'True', 'returns MS only'),
              ('above_toa', 9, '80000', 'at or below TOA'),
              ('surface_ray', 12, '180', 'misses BOA'),
              ('twilight', 11, '95', 'twilight')]
    results = []
    try:
        for name, index, value, expected in checks:
            lines = original.copy()
            lines[index] = value
            path.write_text('\n'.join(lines)+'\n')
            proc = subprocess.run([str(folder / 'DRIVER.exe')], cwd=folder,
                                  capture_output=True, text=True, timeout=180)
            message = proc.stdout+proc.stderr
            (output / f'guard_{name}.log').write_text(message)
            passed = proc.returncode != 0 and expected in message
            results.append(dict(case=name, passed=passed))
        write_csv(output / 'guards.csv', results)
        assert all(r['passed'] for r in results), results
    finally:
        path.write_text('\n'.join(original)+'\n')


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_results(output, rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cases = list(dict.fromkeys((r['height_m'], r['vza'], r['raa']) for r in rows))
    fig, axes = plt.subplots(len(cases), 2, figsize=(11, 3.2*len(cases)), squeeze=False)
    for index, (height, vza, raa) in enumerate(cases):
        selected = [r for r in rows if (r['height_m'], r['vza'], r['raa']) == (height, vza, raa)]
        for row in selected:
            points = row['source_points_per_leg']
            data = np.load(output / f'h{height:g}_v{vza:g}_a{raa:g}_p{points}.npz')
            grid = data['grid']
            z = (grid[1:]+grid[:-1])[::-1]/2000
            for column, label in enumerate(('absorption', 'haze')):
                ax = axes[index, column]
                if row is selected[0]:
                    ax.plot(data['qs_'+label], z, 'k-', lw=2, label='QS total')
                ax.plot(data['fo_'+label]+data['ms_'+label], z, '--', marker='.',
                        label=f'QS FO + Standard MS ({points} points/leg)')
                ax.set_xlabel(f'dI / d{label} (native units)')
                ax.set_ylabel('Altitude (km)')
                ax.set_title(f'h={height:g} m, VZA={vza:g}°, RAA={raa:g}°')
                ax.grid(alpha=0.25)
                ax.legend(handlelength=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(output / 'jacobian_profiles.png', dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 4))
    for points in sorted(set(r['source_points_per_leg'] for r in rows)):
        subset = [r for r in rows if r['source_points_per_leg'] == points]
        labels = [f"{r['height_m']:g} m / {r['vza']:g}° / {r['raa']:g}°" for r in subset]
        ax.plot(labels, [r['error_pct'] for r in subset], 'o-', label=f'{points} source points/leg')
    ax.axhline(0, color='k', lw=0.7)
    ax.set_ylabel('Hybrid radiance − QS (%)')
    ax.set_xlabel('Observer height / VZA / relative azimuth')
    ax.tick_params(axis='x', labelrotation=30)
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / 'comparison.png', dpi=160)
    plt.close(fig)


def main():
    (ROOT / 'tmp').mkdir(exist_ok=True)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true')
    parser.add_argument('--fd', action='store_true')
    parser.add_argument('--qs-fd', action='store_true', help='Also validate the QS reference Jacobians')
    parser.add_argument('--guards', action='store_true')
    parser.add_argument('--vector', action='store_true')
    parser.add_argument('--fine', action='store_true')
    parser.add_argument('--unsplit', action='store_true', help='Disable internal observer/tangent layer splitting')
    parser.add_argument('--case', type=float, nargs=3, action='append', metavar=('HEIGHT', 'VZA', 'RAA'))
    parser.add_argument('--points', type=int, nargs='+')
    parser.add_argument('--qs-fine', type=int, default=4, help='QS path quadrature points (Standard stays at 4)')
    parser.add_argument('--output', type=Path, default=ROOT / 'Examples/testing/limb_hybrid')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    hashes = {name: hashlib.sha256((ROOT / 'Fortran' / name / 'DRIVER.exe').read_bytes()).hexdigest()
              for name in ('vlidort_qs', 'vlidort_standard')}
    (args.output / 'manifest.json').write_text(json.dumps(dict(
        executable_sha256=hashes, arguments=sys.argv[1:],
        atmosphere='477 nm, SZA 60, surface albedo 0.2, aerosol AOD 0.1, scale height 2 km, SSA 0.95, g 0.7',
        delta_m=False, streams=4, scalar=not args.vector), indent=2))
    cases = [(15000, 91, 90), (550, 90, 90)] if args.quick else [
        (15000, 89, 90), (15000, 90, 90), (15000, 91, 90), (15000, 93, 90),
        (550, 90, 90), (550, 90.5, 90), (15000, 91, 0), (15000, 91, 180)]
    if args.case:
        cases = args.case
    rows = []
    for height, vza, raa in cases:
        name = f'h{height:g}_v{vza:g}_a{raa:g}'
        with isolated(hashes) as work:
            p = scene(height, vza, raa, args.vector, args.fine)
            p.settings_rt.n_full_sphericity_fine_layers = args.qs_fine
            qs = native(generate(p, work, 'qs'), args.output / f'{name}_qs.log')
            if args.qs_fd and height == 15000 and vza == 91 and raa == 90:
                finite_difference(work / 'vlidort_qs', qs, args.output / 'qs_fd')
            fo = native(generate(p, work, 'qs', ss=True), args.output / f'{name}_fo.log')
            printed_fo = re.search(r'Radiance first order:\s*([-+\d.EeDd]+)',
                                   (args.output / f'{name}_qs.log').read_text())
            assert printed_fo, 'QS did not report its first-order component'
            np.testing.assert_allclose(fo['radiance'][0], float(printed_fo[1].replace('D', 'E')), rtol=1e-11)
            p.settings_rt.n_full_sphericity_fine_layers = 4
            for points in (args.points or ((5, 9) if args.quick else (2, 5, 9, 17))):
                folder = generate(p, work, 'aircraft_limb_raw' if args.unsplit else 'aircraft_limb_ms', points)
                ms = native(folder, args.output / f'{name}_ms{points}.log')
                total = fo['radiance'][0] + ms['radiance'][0]
                row = dict(height_m=height, vza=vza, raa=raa, vector=args.vector, fine=args.fine, split=not args.unsplit,
                           source_points_per_leg=points, qs=qs['radiance'][0], hybrid=total,
                           qs_fo=fo['radiance'][0], hybrid_ms=ms['radiance'][0],
                           error_pct=100*(total/qs['radiance'][0]-1),
                           qs_seconds=qs['seconds'], fo_seconds=fo['seconds'], ms_seconds=ms['seconds'],
                           hybrid_seconds=fo['seconds']+ms['seconds'])
                for label in ('absorption', 'haze'):
                    reference = qs[label]
                    combined = fo[label]+ms[label]
                    row[label+'_jac_peak_error_pct'] = 100*np.max(np.abs(combined-reference))/np.max(np.abs(reference))
                rows.append(row)
                np.savez(args.output / f'{name}_p{points}.npz', grid=p.vertical_grid,
                         **{'qs_'+k: v for k, v in qs.items()}, **{'fo_'+k: v for k, v in fo.items()},
                         **{'ms_'+k: v for k, v in ms.items()})
                print(json.dumps(row), flush=True)
                write_csv(args.output / 'comparison.csv', rows)
            if args.fd and height == 15000 and vza == 91 and raa == 90:
                finite_difference(folder, ms, args.output)
            if args.guards and height == 15000 and vza == 91 and raa == 90:
                check_guards(folder, args.output)
    plot_results(args.output, rows)
    print('Validation complete.', flush=True)


if __name__ == '__main__':
    main()
