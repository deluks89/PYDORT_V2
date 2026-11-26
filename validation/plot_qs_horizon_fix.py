"""Plot saved before/after QS horizon checks (both use 16-point quadrature)."""
from pathlib import Path
import argparse
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1] / 'Examples/testing/limb_hybrid'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', type=Path, default=ROOT / 'horizon_fine16/comparison.csv')
    parser.add_argument('--after', type=Path, default=ROOT / 'corrected_horizon/comparison.csv')
    parser.add_argument('--output', type=Path, default=ROOT / 'qs_horizon_fix.png')
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for path, label, style in ((args.before, 'Original QS', 'r--^'),
                               (args.after, 'Corrected QS', 'k-o')):
        with path.open() as stream:
            rows = list(csv.DictReader(stream))
        x = [1000*(float(row['vza'])-90) for row in rows]
        for ax, key, title in zip(axes, ('qs', 'qs_fo'), ('Total radiance', 'First-order radiance')):
            ax.plot(x, [float(row[key]) for row in rows], style, lw=2, label=label)
            ax.set_title(title)
            ax.set_xlabel('Viewing angle − 90° (millidegrees)')
            ax.set_ylabel('Normalized radiance')
            ax.grid(alpha=0.25)
            ax.legend(handlelength=4)
    fig.suptitle('QS horizon geometry fix: observer 550 m, SZA 60°, RAA 90°, 477 nm')
    fig.tight_layout()
    fig.savefig(args.output, dpi=180)
    print(args.output)


if __name__ == '__main__':
    main()
