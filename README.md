# VLIDORT_RT

This directory is the working integration tree for the multi-backend Python
wrapper. The original source folders are intentionally left untouched:

- `../VLIDORT-PC`
- `../VLIDORT_SCAirPlane_2p8p3_31Jan2025`

## Layout

- `PYDORT_V2.py` contains the Python wrapper and automatic backend routing.
- `physics_data/` contains the copied wrapper physics data.
- `Fortran/mie_model/` contains the copied Mie helper model.
- `Fortran/vlidort_qs/` contains the copied VLIDORT-QS backend.
- `Fortran/vlidort_standard/` contains the copied Standard VLIDORT 2.8.3 and
  VLIDORT-AC source tree.
- `Fortran/driver_support/` contains support routines shared by the copied
  executable drivers.

## Backend Routing

The Python wrapper exposes one public spherical-correction setting:

- `pyd.settings_rt.sphericity_correction = 1`: Standard VLIDORT without
  spherical correction.
- `pyd.settings_rt.sphericity_correction = 2`: 2-point interpolated correction;
  Standard VLIDORT at or above the top of the altitude grid, VLIDORT-AC below it.
- `pyd.settings_rt.sphericity_correction = 3`: 3-point interpolated correction
  with Standard VLIDORT at or above the top of the altitude grid; below it the
  wrapper falls back to VLIDORT-QS and issues a warning.
- `pyd.settings_rt.sphericity_correction >= 4`: VLIDORT-QS full spherical
  correction.

Use `pyd.run_rt()` as the backend-dispatching entry point.

## Build

The copied Fortran sources can be built from the shared Fortran entry point:

```sh
cd Fortran
mingw32-make
```

The common makefile builds the QS driver, the Standard/AC driver, and the Mie
model. Override the compiler with `FC=...` if needed.

## Actinic flux output

The wrapper supports actinic fluxes at the observer altitude. Full vertical
actinic-flux profiles have been removed from the Python API and are no longer
requested from the copied QS backend.

The copied Standard VLIDORT `vlidort_pars.f90` now sets `MAXLAYERS = 200` for
the wrapper's default vertical grids and `MAX_USER_LEVELS = 1` for single-level
output.
