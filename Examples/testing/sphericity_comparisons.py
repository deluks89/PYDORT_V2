"""
Compare spherical-correction options for Rayleigh limb observations.

The script simulates limb radiances and O4 slant columns for ten observer
altitudes, three solar zenith angles, and three spherical-correction settings:

    1. no spherical correction
    2. automatic aircraft 2-point correction
    3. full VLIDORT-QS spherical correction

Results are written as tab-delimited text files in OUTPUT_DIR.
"""

#%% Imports

from pathlib import Path
import sys
import traceback

import numpy as np
import pandas as pd


#%% Project import

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import PYDORT_V2


#%% Configuration

OUTPUT_DIR = Path(__file__).resolve().parent / "sphericity_comparison_results"

WAVELENGTHS_NM = np.array([477.0])
OBSERVER_ALTITUDES_M = np.array([100,   500,  1000,  2000, 4000, 6000, 8100, 
                              10000, 11000, 12500, 15000])
SOLAR_ZENITH_ANGLES_DEG = [0.0, 60.0, 85.0]

SAA_DEG = 0.0
VAA_DEG = 90.0
VZA_DEG = 89

SURFACE_ALBEDO = 0.0
N_STREAMS = 4
N_FULL_SPHERICITY_FINE_LAYERS = 4

SPHERICITY_CASES = [
    {
        "label": "none",
        "sphericity_correction": 1,
        "description": "Standard VLIDORT with spherical FO but no MS spherical correction",
    },
    {
        "label": "aircraft_2point",
        "sphericity_correction": 2,
        "description": "VLIDORT-AC 2-point correction for below-TOA observers",
    },
    {
        "label": "full_qs",
        "sphericity_correction": 4,
        "description": "VLIDORT-QS full spherical correction",
    },
]


#%% Setup helpers

def zero_trace_gas_absorption(pyd):
    """Keep gas profiles available but remove their absorption from the RT run."""
    for gas in pyd.trace_gases.to_dictionary().values():
        gas.xs.raw_data[:, 1:] = 0.0
        if gas.xs.smoothed_data is not None:
            gas.xs.smoothed_data[:, 1:] = 0.0


def make_base_pydort():
    pyd = PYDORT_V2.Pydort()
    pyd.wavelengths = WAVELENGTHS_NM.copy()
    pyd.aerosols.remove_all()
    zero_trace_gas_absorption(pyd)
    pyd.surface_albedo.set_constant(SURFACE_ALBEDO)

    pyd.viewing_geometry.saa = SAA_DEG
    pyd.viewing_geometry.vaa = VAA_DEG
    pyd.viewing_geometry.vza = VZA_DEG

    pyd.settings_rt.do_radiance = True
    pyd.settings_rt.do_actinic_flux_at_observer = False
    pyd.settings_rt.do_radiance_jacobians = True
    pyd.settings_rt.do_actinic_flux_jacobians = False
    pyd.settings_rt.jacobian_list = ["absorption"]
    pyd.settings_rt.n_halfspace_streams = N_STREAMS
    pyd.settings_rt.n_full_sphericity_fine_layers = N_FULL_SPHERICITY_FINE_LAYERS

    return pyd


#%% Result helpers

def get_resolved_backend(pyd):
    pyd.settings_rt.validate_backend_settings()
    return pyd._resolve_rt_backend(warn=False)


def append_failure_rows(radiance_rows, o4_scd_rows, case, sza, altitude, error):
    common = {
        "case": case["label"],
        "sphericity_correction": case["sphericity_correction"],
        "description": case["description"],
        "sza_deg": sza,
        "saa_deg": SAA_DEG,
        "vza_deg": VZA_DEG,
        "vaa_deg": VAA_DEG,
        "observer_altitude_m": altitude,
        "status": "failed",
        "error": error,
    }
    for wavelength in WAVELENGTHS_NM:
        radiance_rows.append({
            **common,
            "resolved_model": "",
            "resolved_backend_mode": "",
            "sphericity_interpolation_points": np.nan,
            "wavelength_nm": wavelength,
            "radiance_photons_s_cm2_sr_nm": np.nan,
            "radiance_normalized": np.nan,
        })
        o4_scd_rows.append({
            **common,
            "resolved_model": "",
            "resolved_backend_mode": "",
            "sphericity_interpolation_points": np.nan,
            "wavelength_nm": wavelength,
            "o4_scd_molec2_cm5": np.nan,
        })


def append_success_rows(radiance_rows, o4_scd_rows, pyd, case, sza, altitude, backend):
    radiance = pyd.results_rt.get_radiance(normalized=False)
    normalized_radiance = pyd.results_rt.get_radiance(normalized=True)
    o4_scd = pyd.results_rt.get_scds("o4")

    common = {
        "case": case["label"],
        "sphericity_correction": case["sphericity_correction"],
        "description": case["description"],
        "resolved_model": backend["model"],
        "resolved_backend_mode": backend["backend_mode"],
        "fallback_reason": backend.get("fallback_reason", ""),
        "sphericity_interpolation_points": backend["sphericity_interpolation_points"],
        "sza_deg": sza,
        "saa_deg": SAA_DEG,
        "vza_deg": VZA_DEG,
        "vaa_deg": VAA_DEG,
        "observer_altitude_m": altitude,
        "status": "ok",
        "error": "",
    }

    new_radiance_rows = []
    new_o4_scd_rows = []
    for wavelength in WAVELENGTHS_NM:
        radiance_value = float(radiance.loc[wavelength])
        normalized_radiance_value = float(normalized_radiance.loc[wavelength])
        o4_scd_value = float(o4_scd.loc[wavelength])
        values = [radiance_value, normalized_radiance_value, o4_scd_value]
        if not np.all(np.isfinite(values)):
            raise ValueError(
                "Non-finite RT result encountered for "
                f"{case['label']}, SZA={sza:g} deg, "
                f"altitude={altitude:.1f} m, wavelength={wavelength:g} nm: "
                f"radiance={radiance_value}, "
                f"normalized_radiance={normalized_radiance_value}, "
                f"o4_scd={o4_scd_value}"
            )

        new_radiance_rows.append({
            **common,
            "wavelength_nm": wavelength,
            "radiance_photons_s_cm2_sr_nm": radiance_value,
            "radiance_normalized": normalized_radiance_value,
        })
        new_o4_scd_rows.append({
            **common,
            "wavelength_nm": wavelength,
            "o4_scd_molec2_cm5": o4_scd_value,
        })

    radiance_rows.extend(new_radiance_rows)
    o4_scd_rows.extend(new_o4_scd_rows)


#%% Main comparison

PYDORT_V2.REPRINT_SUBPROCESS_OUTPUT = False
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

base_pyd = make_base_pydort()
radiance_rows = []
o4_scd_rows = []

for case in SPHERICITY_CASES:
    for sza in SOLAR_ZENITH_ANGLES_DEG:
        for altitude in OBSERVER_ALTITUDES_M:
            pyd = base_pyd.copy()
            pyd.settings_rt.sphericity_correction = case["sphericity_correction"]
            pyd.viewing_geometry.sza = sza
            pyd.viewing_geometry.observer_height_agl = float(altitude)
            backend = get_resolved_backend(pyd)

            print(
                f"Running {case['label']}: "
                f"SZA={sza:g} deg, altitude={altitude:.1f} m, "
                f"backend={backend['backend_mode']}"
            )

            try:
                pyd.run_rt(run_aerosol_models=False)
                append_success_rows(
                    radiance_rows,
                    o4_scd_rows,
                    pyd,
                    case,
                    sza,
                    float(altitude),
                    backend,
                )
            except Exception as exc:
                error = repr(exc)
                print("  FAILED:", error)
                traceback.print_exc()
                append_failure_rows(
                    radiance_rows,
                    o4_scd_rows,
                    case,
                    sza,
                    float(altitude),
                    error,
                )

radiance_table = pd.DataFrame(radiance_rows)
o4_scd_table = pd.DataFrame(o4_scd_rows)

radiance_table.to_csv(
    OUTPUT_DIR / "limb_rayleigh_radiances.txt",
    sep="\t",
    index=False,
    float_format="%.10e",
)
o4_scd_table.to_csv(
    OUTPUT_DIR / "limb_rayleigh_o4_scds.txt",
    sep="\t",
    index=False,
    float_format="%.10e",
)

print("\nWrote:")
print(" ", OUTPUT_DIR / "limb_rayleigh_radiances.txt")
print(" ", OUTPUT_DIR / "limb_rayleigh_o4_scds.txt")

