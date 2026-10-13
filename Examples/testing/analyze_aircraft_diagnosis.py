"""Analyze the component and layer dumps from diagnose_aircraft_sc.py (Spyder).

The midpoint calculation below is a diagnostic reconstruction, not an RT backend.
"""

#%% Imports and paths
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent / 'aircraft_diagnosis'


def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write_rows(path, rows):
    if not rows:
        return
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


#%% Component comparisons
def analyze(folder):
    data = read_rows(folder / 'components.txt')
    keys = sorted({(float(r['altitude_m']), float(r['sza_deg']), float(r['vza_deg'])) for r in data})
    comparisons, checks = [], []
    for altitude, sza, vza in keys:
        case = {(int(r['correction']), r['single_scattering_only'] == 'True'): r for r in data
                if (float(r['altitude_m']), float(r['sza_deg']), float(r['vza_deg'])) == (altitude, sza, vza)}
        standard, ac, qs = case[(1, False)], case[(2, False)], case[(4, False)]
        iqs = float(qs['radiance'])
        foqs = float(case[(4, True)]['radiance'])
        msqs = iqs - foqs
        comparison = dict(altitude_m=altitude, sza_deg=sza, vza_deg=vza,
                          standard_radiance=float(standard['radiance']), ac_radiance=float(ac['radiance']),
                          qs_radiance=iqs, standard_fo=float(standard['corrected_fo']),
                          ac_fo=float(ac['corrected_fo']), qs_fo=foqs,
                          standard_ms=float(standard['corrected_ms']), ac_ms=float(ac['corrected_ms']),
                          qs_ms=msqs,
                          standard_radiance_error_pct=100*(float(standard['radiance'])/iqs-1),
                          ac_radiance_error_pct=100*(float(ac['radiance'])/iqs-1),
                          ac_ms_error_pct=100*(float(ac['corrected_ms'])/msqs-1),
                          standard_o4_error_pct=100*(float(standard['o4_scd'])/float(qs['o4_scd'])-1),
                          ac_o4_error_pct=100*(float(ac['o4_scd'])/float(qs['o4_scd'])-1))
        comparisons.append(comparison)
        print(comparison)
        case_dir = folder / f'h{altitude:g}_sza{sza:g}_vza{vza:g}_c2_ss0'
        layers = np.loadtxt(case_dir / 'rt_layers.out', ndmin=2)
        if vza > 90:
            continue
        # All selected tests use the default Earth radius and no delta-M truncation.
        radius = 6371000.0
        b = (radius + altitude)*np.sin(np.deg2rad(vza))
        mu_ac = np.cos(np.deg2rad(vza))
        mu_toa = np.sqrt(1-(b/(radius+layers[0, 1]))**2)
        layers = layers[layers[:, 2] >= altitude-1e-6]
        mu_top = np.sqrt(1-(b/(radius+layers[:, 1]))**2)
        mu_bottom = np.sqrt(1-(b/(radius+layers[:, 2]))**2)
        lengths = np.sqrt((radius+layers[:, 1])**2-b*b)-np.sqrt((radius+layers[:, 2])**2-b*b)
        trans = np.exp(-layers[:, 3]*lengths/(layers[:, 1]-layers[:, 2]))
        weights = (mu_toa-layers[:, 4:6].mean(axis=1))/(mu_toa-mu_ac)
        interpolated = weights*layers[:, 8]+(1-weights)*layers[:, 9]
        downstream = np.cumprod(np.r_[1.0, layers[:0:-1, 6]])[::-1]
        reconstructed = float(interpolated @ downstream)
        check = dict(altitude_m=altitude, sza_deg=sza, vza_deg=vza,
                     max_path_mu_error=float(np.max(np.abs(layers[:, 4:6]-np.c_[mu_top, mu_bottom]))),
                     max_transmittance_error=float(np.max(np.abs(layers[:, 6]-trans))),
                     min_interpolation_weight=float(weights.min()), max_interpolation_weight=float(weights.max()),
                     reconstructed_ac_ms=reconstructed,
                     native_ac_ms=float(ac['corrected_ms']),
                     reconstructed_minus_native=reconstructed-float(ac['corrected_ms']))
        np.testing.assert_allclose(reconstructed, float(ac['corrected_ms']), rtol=1e-10)
        checks.append(check)
        print(check)

        sampled_path = case_dir / 'rt_sampled_sources.out'
        if not sampled_path.exists():
            continue
        sampled = np.loadtxt(sampled_path, ndmin=2)
        np.testing.assert_array_equal(sampled[:, 0], layers[:, 0])
        local = sampled[:, 4]
        midpoint_ms = float(local @ downstream)
        print('DIAGNOSTIC midpoint sources: MS=', midpoint_ms,
              'total error vs QS (%)=', 100*((midpoint_ms+float(ac['corrected_fo']))/iqs-1))
        contributions = []
        for i in range(len(layers)):
            contributions.append(dict(layer=int(layers[i, 0]), altitude_mid_m=float(layers[i, 1:3].mean()),
                                      local_vza_deg=sampled[i, 2], endpoint_weight=weights[i],
                                      ms_interpolated=interpolated[i], ms_direct_midpoint=local[i],
                                      downstream_transmittance=downstream[i],
                                      interpolation_excess_at_observer=(interpolated[i]-local[i])*downstream[i]))
        write_rows(case_dir / 'source_analysis.txt', contributions)
        height = layers[:, 1:3].mean(axis=1)/1000
        fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True, constrained_layout=True)
        axes[0].plot(interpolated*downstream, height, label='AC 2-point interpolation')
        axes[0].plot(local*downstream, height, label='Direct local-geometry evaluation')
        axes[0].set(xlabel='Layer MS contribution at observer', ylabel='Altitude (km)', ylim=(0, 30))
        axes[0].ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
        axes[0].locator_params(axis='x', nbins=5)
        axes[0].legend(fontsize=8)
        axes[1].plot(100*(interpolated/local-1), height, color='crimson')
        axes[1].set(xlabel='Interpolation excess over direct source (%)')
        for ax in axes:
            ax.grid(alpha=0.25)
        fig.suptitle(f'Valid aircraft geometry: {altitude:g} m, SZA {sza:g}, VZA {vza:g}')
        fig.savefig(case_dir / 'source_interpolation.png', dpi=160)
        plt.close(fig)
    write_rows(folder / 'comparison.txt', comparisons)
    write_rows(folder / 'geometry_checks.txt', checks)


#%% Run
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder', nargs='?', type=Path, default=ROOT / 'sampled')
    args = parser.parse_args()
    analyze(args.folder)
