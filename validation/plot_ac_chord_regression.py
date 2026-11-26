"""Plot the saved before/after horizon regression against QS."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parents[1]
folder = root / 'Examples/testing/limb_hybrid/chord_horizon6'
with (folder / 'comparison.csv').open() as f:
    rows = list(csv.DictReader(f))
x = np.arange(len(rows))
fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
for ax, key, label in zip(axes, ['radiance', 'haze'], ['Radiance', 'Aerosol extinction Jacobian']):
    for prefix, name, style in [('previous', 'Previous AC', 'o--'), ('shared', 'Corrected AC', 's-')]:
        values = [100*float(r[prefix+'_'+key+'_qs_error']) for r in rows]
        ax.semilogy(x, values, style, label=name, linewidth=1.8)
    ax.set_xticks(x, ['89.99999', '90', '90.00001'])
    ax.set_xlabel('Viewing zenith angle (degrees)')
    ax.set_ylabel('Maximum difference / QS peak (%)')
    ax.set_title(label)
    ax.grid(alpha=.25, which='both')
    ax.legend(handlelength=3.5)
fig.suptitle('550 m observer, 477 nm: three angular points; six FO/QS fine points')
fig.savefig(folder / 'comparison.png', dpi=180)
print(folder / 'comparison.png')
