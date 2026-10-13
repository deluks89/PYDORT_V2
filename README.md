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

- `pyd.settings_rt.sphericity_correction = 1`: Standard VLIDORT with spherical
  first-order scattering, but without the multiple-scattering spherical correction.
- `pyd.settings_rt.sphericity_correction = 2`: 2-point interpolated correction;
  Standard VLIDORT at or above the top of the altitude grid, VLIDORT-AC below it.
- `pyd.settings_rt.sphericity_correction = 3`: 3-point interpolated correction
  with Standard VLIDORT at or above the top of the altitude grid; below it the
  wrapper falls back to VLIDORT-QS and issues a warning.
- `pyd.settings_rt.sphericity_correction >= 4`: VLIDORT-QS full spherical
  correction.

Use `pyd.run_rt()` as the backend-dispatching entry point.

Aircraft 2-point calculations require an observer on an interior layer boundary
and a ray that can be extended to the surface. Unsupported geometries fall back
to QS with a warning, rather than silently rounding the observer height or
returning NaNs. At VZA 89 degrees on Earth, the surface-intersection limit is
approximately 970 m. The comparison tables record the resolved backend and
`fallback_reason`; these rows are QS results, not aircraft interpolation results.
Surface-missing, downward-looking limb rays also require QS for spherical FO.

## Build

The copied Fortran sources can be built from the shared Fortran entry point:

```sh
cd Fortran
mingw32-make
```

The common makefile builds the QS driver, the Standard/AC driver, and the Mie
model. The validated Standard/AC build uses gfortran and binary128 arithmetic
(`-freal-8-real-16`) throughout the driver, core, and bundled LAPACK routines.
This prevents catastrophic cancellation in native LPS Jacobians near conservative
scattering. The adapter also uses an ASYMTX tolerance of `1e-28`: `1e-20`
still left errors in conservative polarized calculations. The existing native
SSA limiter is unchanged; no background absorption is added, and the Jacobians
are not replaced by finite differences. QS and Mie are unchanged.

The makefile rebuilds all Standard/AC objects after changes to dimensions,
type definitions, or the makefile. When overriding compiler flags on the command
line, also use `-B`; mixing precision settings between objects is invalid.
Other compiler configurations have not been validated for conservative Jacobians.

## Jacobian regression

```sh
python Examples/testing/validate_native_jacobians.py
python Examples/testing/validate_native_jacobians.py --suite
python Examples/testing/validate_native_jacobians.py --suite --vector --output Examples/testing/native_jacobian_validation/vector_full_grid
python Examples/testing/validate_rt_workflows.py
```

The Spyder-style validation script compares native absorption Jacobians with
independent layer-OD perturbations and an O4-weighted perturbation, and checks
that enabling Jacobians does not change the radiance. It uses temporary Fortran
run directories and writes separate text tables under
`Examples/testing/native_jacobian_validation/`. Existing comparison results are
not overwritten. Finite differences are used only in this validation script.
The public-API workflow check also exercises QS fallback, spherical first-order
agreement with QS, and repeated wavelengths with both radiances and observer
actinic fluxes enabled. These checks cover absorption Jacobians of Stokes I
(including polarized RT), not every aerosol or surface Jacobian.

## Actinic flux output

The wrapper supports actinic fluxes at the observer altitude. Full vertical
actinic-flux profiles have been removed from the Python API and are no longer
requested from the copied QS backend.

The copied Standard VLIDORT `vlidort_pars.f90` sets `MAXLAYERS = 101` for
the wrapper's default vertical grids and `MAX_USER_LEVELS = 1` for single-level
output.
