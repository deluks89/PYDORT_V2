"""
Visualize spherical-correction comparison results.

This script reads the text files written by sphericity_comparisons.py and
creates profile-style figures for limb radiances and O4 slant columns.
"""

#%% Imports

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


#%% Configuration

RESULT_DIR = Path(__file__).resolve().parent / "sphericity_comparison_results"

RADIANCE_FILE = RESULT_DIR / "limb_rayleigh_radiances.txt"
O4_SCD_FILE = RESULT_DIR / "limb_rayleigh_o4_scds.txt"

FIG_ABSOLUTE = RESULT_DIR / "sphericity_absolute_values.png"
FIG_RELATIVE = RESULT_DIR / "sphericity_relative_to_full_qs.png"

CASE_ORDER = ["none", "aircraft_2point", "full_qs"]
CASE_LABELS = {
    "none": "No correction",
    "aircraft_2point": "Aircraft 2-point",
    "full_qs": "Full QS",
}
CASE_COLORS = {
    "none": "C3",
    "aircraft_2point": "C0",
    "full_qs": "k",
}
CASE_MARKERS = {
    "none": "o",
    "aircraft_2point": "s",
    "full_qs": "^",
}

BASELINE_CASE = "full_qs"


#%% Load results

def require_file(path):
    if not path.exists():
        raise FileNotFoundError(
            f"Missing result file: {path}\n"
            "Run Examples/testing/sphericity_comparisons.py first."
        )


def load_results():
    require_file(RADIANCE_FILE)
    require_file(O4_SCD_FILE)

    radiances = pd.read_csv(RADIANCE_FILE, sep="\t")
    o4_scds = pd.read_csv(O4_SCD_FILE, sep="\t")

    radiance_columns = [
        "radiance_photons_s_cm2_sr_nm",
        "radiance_normalized",
    ]
    o4_scd_columns = ["o4_scd_molec2_cm5"]
    for column in radiance_columns:
        radiances[column] = pd.to_numeric(radiances[column], errors="coerce")
    for column in o4_scd_columns:
        o4_scds[column] = pd.to_numeric(o4_scds[column], errors="coerce")

    if not np.isfinite(radiances["radiance_normalized"]).any():
        raise RuntimeError("Radiance result file contains no finite simulations.")
    if not np.isfinite(o4_scds["o4_scd_molec2_cm5"]).any():
        raise RuntimeError("O4 SCD result file contains no finite simulations.")

    return radiances, o4_scds


radiances, o4_scds = load_results()


#%% Plot helpers

def sorted_unique(values):
    return sorted(pd.unique(values))


def case_sort_key(case):
    if case in CASE_ORDER:
        return CASE_ORDER.index(case)
    return len(CASE_ORDER)


def iter_cases(df):
    return sorted(pd.unique(df["case"]), key=case_sort_key)


def plot_quantity_profile(ax, df, value_column, sza, wavelength, xlabel, xscale=None):
    sub = df[
        (df["sza_deg"] == sza)
        & (df["wavelength_nm"] == wavelength)
    ]

    for case in iter_cases(sub):
        case_data = sub[sub["case"] == case].sort_values("observer_altitude_m")
        x = case_data[value_column].to_numpy()
        y = case_data["observer_altitude_m"].to_numpy() / 1000.0
        finite = np.isfinite(x) & np.isfinite(y)
        if not finite.any():
            ax.plot(
                [],
                [],
                marker=CASE_MARKERS.get(case, "o"),
                color=CASE_COLORS.get(case, None),
                linestyle=":",
                label=CASE_LABELS.get(case, case) + " (no finite data)",
            )
            continue
        ax.plot(
            x[finite],
            y[finite],
            marker=CASE_MARKERS.get(case, "o"),
            color=CASE_COLORS.get(case, None),
            label=CASE_LABELS.get(case, case),
        )

    if xscale is not None:
        ax.set_xscale(xscale)
    ax.set_title(f"SZA {sza:g} deg")
    ax.set_xlabel(xlabel)
    ax.grid(True, alpha=0.3)


def ratio_to_baseline(df, value_column):
    key_cols = ["sza_deg", "wavelength_nm", "observer_altitude_m"]
    baseline = df[df["case"] == BASELINE_CASE][key_cols + [value_column]].copy()
    baseline = baseline.rename(columns={value_column: "baseline_value"})

    merged = df.merge(baseline, on=key_cols, how="inner")
    merged = merged[merged["case"] != BASELINE_CASE].copy()
    merged["ratio"] = merged[value_column] / merged["baseline_value"]
    merged["percent_difference"] = (merged["ratio"] - 1.0) * 100.0
    return merged


def plot_ratio_profile(ax, df, sza, wavelength, xlabel):
    sub = df[
        (df["sza_deg"] == sza)
        & (df["wavelength_nm"] == wavelength)
    ]

    for case in iter_cases(sub):
        case_data = sub[sub["case"] == case].sort_values("observer_altitude_m")
        x = case_data["percent_difference"].to_numpy()
        y = case_data["observer_altitude_m"].to_numpy() / 1000.0
        finite = np.isfinite(x) & np.isfinite(y)
        if not finite.any():
            ax.plot(
                [],
                [],
                marker=CASE_MARKERS.get(case, "o"),
                color=CASE_COLORS.get(case, None),
                linestyle=":",
                label=CASE_LABELS.get(case, case) + " (no finite data)",
            )
            continue
        ax.plot(
            x[finite],
            y[finite],
            marker=CASE_MARKERS.get(case, "o"),
            color=CASE_COLORS.get(case, None),
            label=CASE_LABELS.get(case, case),
        )

    ax.axvline(0.0, color="0.4", lw=1)
    ax.set_title(f"SZA {sza:g} deg")
    ax.set_xlabel(xlabel)
    ax.grid(True, alpha=0.3)


def add_shared_legend(fig, axes):
    handles, labels = [], []
    for ax in np.ravel(axes):
        h, l = ax.get_legend_handles_labels()
        for handle, label in zip(h, l):
            if label not in labels:
                handles.append(handle)
                labels.append(label)
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=len(handles), frameon=False)


#%% Absolute-value figure

def make_absolute_figure(wavelength):
    szas = sorted_unique(radiances["sza_deg"])
    fig, axes = plt.subplots(
        2,
        len(szas),
        figsize=(4.5 * len(szas), 8),
        sharey=True,
        constrained_layout=True,
    )

    for col, sza in enumerate(szas):
        plot_quantity_profile(
            axes[0, col],
            radiances,
            "radiance_normalized",
            sza,
            wavelength,
            "Normalized radiance",
            xscale="log",
        )
        plot_quantity_profile(
            axes[1, col],
            o4_scds,
            "o4_scd_molec2_cm5",
            sza,
            wavelength,
            "O4 SCD [molec2 cm-5]",
            xscale="log",
        )

    axes[0, 0].set_ylabel("Observer altitude [km]")
    axes[1, 0].set_ylabel("Observer altitude [km]")
    fig.suptitle(f"Sphericity comparison at {wavelength:g} nm", y=1.04)
    add_shared_legend(fig, axes)
    return fig


#%% Relative-to-QS figure

def make_relative_figure(wavelength):
    radiance_ratios = ratio_to_baseline(radiances, "radiance_normalized")
    o4_scd_ratios = ratio_to_baseline(o4_scds, "o4_scd_molec2_cm5")
    szas = sorted_unique(radiances["sza_deg"])

    fig, axes = plt.subplots(
        2,
        len(szas),
        figsize=(4.5 * len(szas), 8),
        sharey=True,
        constrained_layout=True,
    )

    for col, sza in enumerate(szas):
        plot_ratio_profile(
            axes[0, col],
            radiance_ratios,
            sza,
            wavelength,
            "Radiance difference vs full QS [%]",
        )
        plot_ratio_profile(
            axes[1, col],
            o4_scd_ratios,
            sza,
            wavelength,
            "O4 SCD difference vs full QS [%]",
        )

    axes[0, 0].set_ylabel("Observer altitude [km]")
    axes[1, 0].set_ylabel("Observer altitude [km]")
    fig.suptitle(f"Sphericity differences relative to full QS at {wavelength:g} nm", y=1.04)
    add_shared_legend(fig, axes)
    return fig


#%% Create and save figures

wavelengths = sorted_unique(radiances["wavelength_nm"])

for wavelength in wavelengths:
    fig_abs = make_absolute_figure(wavelength)
    fig_rel = make_relative_figure(wavelength)

    if len(wavelengths) == 1:
        fig_abs.savefig(FIG_ABSOLUTE, dpi=200, bbox_inches="tight")
        fig_rel.savefig(FIG_RELATIVE, dpi=200, bbox_inches="tight")
        print("Wrote:", FIG_ABSOLUTE)
        print("Wrote:", FIG_RELATIVE)
    else:
        fig_abs_path = RESULT_DIR / f"sphericity_absolute_values_{wavelength:g}nm.png"
        fig_rel_path = RESULT_DIR / f"sphericity_relative_to_full_qs_{wavelength:g}nm.png"
        fig_abs.savefig(fig_abs_path, dpi=200, bbox_inches="tight")
        fig_rel.savefig(fig_rel_path, dpi=200, bbox_inches="tight")
        print("Wrote:", fig_abs_path)
        print("Wrote:", fig_rel_path)

plt.show()
