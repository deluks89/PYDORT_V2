"""Check saved three-point curved-source AC accuracy and continuity across 90 degrees.
Run validate_shared_ac.py with three horizon cases and --fine-points 6 first.
"""
from pathlib import Path
import argparse
import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--input', type=Path, default=Path(__file__).resolve().parents[1] / 'Examples/testing/limb_hybrid/chord_horizon6')
a = p.parse_args()
values = []
for angle in ['89.99999', '90', '90.00001']:
    case = f'h550_v{angle}_a90'
    with np.load(a.input / (case + '_qs.npz')) as q, np.load(a.input / (case + '_p3_shared.npz')) as s:
        values.append({key: s[key].copy() for key in ['radiance', 'absorption', 'haze']})
        for key, limit in [('radiance', 1e-4), ('absorption', 1e-4), ('haze', 1e-3)]:
            assert np.all(np.isfinite(s[key])), (angle, key)
            error = np.max(np.abs(s[key] - q[key])) / np.max(np.abs(q[key]))
            assert error < limit, (angle, key, error, limit)
for key in values[0]:
    scale = np.max(np.abs(values[1][key]))
    jump = max(np.max(np.abs(values[i][key] - values[1][key])) for i in [0, 2]) / scale
    assert jump < 1e-5, (key, jump)
    print(f'{key}: maximum adjacent-angle change / peak = {jump:.8g}')
print('Curved-source horizon accuracy and continuity checks passed.')
