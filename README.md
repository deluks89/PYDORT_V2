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

- `pyd.settings_rt.sphericity_correction = 0`: experimental sparse-source
  aircraft approximation (see below). Existing values retain their meanings.
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

## Experimental Aircraft Approximation

```python
pyd.settings_rt.sphericity_correction = 0
pyd.settings_rt.n_aircraft_source_points = 5  # integer, 2 through 33
pyd.run_rt()
```

Implemented in `Fortran/vlidort_standard/sparse_aircraft.f90`, this method
samples local pseudo-spherical VLIDORT MS solutions at a small number of
geometries along the observer ray. For each layer it divides the integrated
source by `1-exp(-tau/mu)` before interpolating in `log(mu)`. Here `tau` is
the delta-M-scaled vertical optical depth and `mu` the viewing cosine.
This removes much of the near-horizon path-length dependence responsible
for the original AC interpolation bias. The normalized sources are then
integrated using spherical shell intersections and exact segment transmission,
with a piecewise-constant effective source within each intersected layer.
Nodes are geometry-based, not adaptively chosen from a radiance error estimate.

Spherical first-order scattering remains enabled and is evaluated separately
along the observer ray. Native LPS source derivatives are propagated through
normalization, interpolation and transmission analytically. No finite-difference
Jacobian replacement is used. For aerosol delta-M calculations, the FO helper
uses a full local solve to retain consistent scaled extinction; its MS output
is discarded. This adds some overhead.

Supported radiance geometries include upward-looking rays (`0 <= VZA < 90`)
even when their backward extension misses the surface, and surface-intersecting
downward-looking rays (`90 < VZA <= 180`). Observers may lie between layer
boundaries. Exactly horizontal and downward-looking turning rays currently
fall back to QS with a warning. Twilight local-source geometries are rejected.
The atmosphere is horizontally uniform and refraction is not included.
Observer fluxes retain the native pseudo-spherical flux calculation; this
line-of-sight correction is not an angularly integrated flux correction.

The five-point Rayleigh test at 477 nm, 100 m, SZA=60, VZA=89 reduces radiance
error relative to QS from about 11.1% (old AC) to 0.0008%, and O4 SCD error
from 8.4% to 1.0% with four half-space streams. These are selected test results,
not accuracy guarantees. At 15 km, SZA=60, VZA=85 the radiance is extremely
close to QS but the O4 SCD differs by about -5.4%, worse than old AC in that
case. Angular-stream and vertical-grid convergence remain important.
Increasing source points alone does not guarantee closer agreement with QS.
Runtime grows with point count and is not always lower than QS, particularly
at high altitude or for short downward-looking paths. Treat this as an
experimental intermediate approximation, not a generally validated QS replacement.

Reproducible Spyder-style comparisons and native-derivative checks:

```sh
python Examples/testing/compare_sparse_aircraft.py --suite
python Examples/testing/validate_native_jacobians.py --correction 0
python Examples/testing/validate_native_jacobians.py --correction 0 --altitude 8100 --full-grid --output Examples/testing/sparse_native_checks/high_altitude
python Examples/testing/validate_native_jacobians.py --correction 0 --altitude 550 --vector --output Examples/testing/sparse_native_checks/vector_partial
python Examples/testing/validate_sparse_workflows.py
```

Comparisons run three, five and nine source points, recording the actual
backend, O4 SCDs, radiances and timings in `Examples/testing/sparse_aircraft_results/`.
They omit unsupported old-AC cases rather than mislabeling QS fallbacks.
Tests use temporary backend folders and leave normal simulation files untouched.

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

## Above-TOA radiance observers

Set `viewing_geometry.observer_height_agl` to the actual observer altitude,
with SZA, VZA and azimuths defined in the local frame at that observer.
The QS and Standard Fortran drivers automatically trace an above-TOA viewing
ray to the near atmospheric boundary and transform all three geometry angles
(SZA, VZA and relative azimuth) into the local frame at that entry point.
The conversion assumes straight propagation in vacuum and a distant Sun.
Radiance and its atmospheric Jacobians need no distance scaling.

Observers at or below TOA retain their existing behavior. QS still places
its numerical camera 1 mm inside TOA. Rays pointing away from the atmosphere,
missing or merely tangent to TOA, and entries with local SZA greater than
90 degrees are rejected with an error. Existing backend restrictions on
surface-missing limb rays still apply. Above-TOA actinic flux is rejected:
mapping one viewing ray does not transform an angularly integrated flux.
Disable `do_actinic_flux_at_observer` for satellite radiance simulations.

Regression: `python Examples/testing/validate_satellite_geometry.py` compares
external observers with independently constructed equivalent TOA geometries,
including native absorption Jacobians, using isolated backend directories.

## Irradiance output

Both drivers expose native cosine-weighted irradiance on a local horizontal
plane at the observer. This is distinct from actinic flux, which has no cosine
weighting. Both are restricted to observer heights at or below TOA.

```python
pyd.settings_rt.do_irradiance_at_observer = True
pyd.settings_rt.do_irradiance_jacobians = True  # optional
pyd.settings_rt.jacobian_list = ['absorption']
pyd.run_rt()

down = pyd.results_rt.get_irradiance()  # photons / (s cm2 nm)
net = pyd.results_rt.get_irradiance(direction='net', units='energy')  # W / (m2 nm)
direct = pyd.results_rt.get_irradiance(component='direct', normalized=True)
jac = pyd.results_rt.get_irradiance_jacobians('absorption', normalized=True)
```

`direction` accepts `down`, `up`, `net` (down minus up), or `sum` (down plus up).
`component` accepts `direct`, `diffuse`, or `total`. Upward direct solar
irradiance is zero; reflected sunlight belongs to the upward diffuse field.
`units` accepts `photons` (default) or `energy`. `normalized=True` returns
irradiance divided by the TOA solar beam spectrum, independent of the unit
selection; direct normalized irradiance is therefore cosine(SZA) times the
solar transmission. Energy output is spectral, not wavelength integrated.
Jacobians use the same options and scaling, with altitude rows and wavelength
columns, and support the existing Jacobian labels and layer flags.

Flux solves calculate both angular moments together; requesting irradiance
alongside actinic flux does not require additional native RT solves. Irradiance
alone also makes the companion actinic outputs available. Standard/AC flux
calculations retain their native pseudo-spherical treatment; the line-of-sight
spherical correction is not applied to these angular integrals. QS uses its
native regular-flux outputs. Viewing zenith and azimuth do not select the
orientation of the horizontal irradiance plane. Tilted planes and vertical
output profiles are not exposed by this API.

Run `python Examples/testing/validate_irradiance.py` for component, unit,
direct-beam normalization and finite-difference Jacobian checks in isolated
backend directories. Older Fortran configuration files remain valid: the two
new trailing irradiance flags default to false if absent.

The QS profile-linearized diffuse regular-flux quadrature now includes the
same `2*pi` azimuth-integration factor as the forward regular flux. This fixes
a native normalization omission exposed by the irradiance finite-difference
checks; direct flux, actinic flux and radiance formulas are unchanged.

### Single and multiple scattering irradiance

Enable the optional order split before running:

```python
pyd.settings_rt.do_irradiance_at_observer = True
pyd.settings_rt.do_deltam_phase_function_scaling = False
pyd.settings_rt.do_irradiance_scattering_components = True
pyd.run_rt()
single = pyd.results_rt.get_irradiance(component='single')
multiple = pyd.results_rt.get_irradiance(component='multiple')
```

One atmospheric scattering **or one surface reflection** counts as first order.
The unscattered solar beam remains `direct`; `diffuse = single + multiple` and
`total = direct + single + multiple`. The same components and directions are
available through `get_irradiance_jacobians` when irradiance Jacobians are enabled.
This option requires a full-scattering run (`do_single_scattering_only=False`)
and delta-M scaling disabled (Rayleigh-only runs already disable delta-M).
Delta-M redistributes forward scattering and would change the order convention.

The drivers extract the coefficient of scattering strength from their native
hemispheric flux solution: total extinction and phase matrices stay fixed while
all single-scattering albedos and the surface albedo are scaled together.
Four additional weak-scattering solves per wavelength give a quadratic
extrapolation of diffuse flux divided by scattering strength to zero strength.
Two step sizes must agree within `1e-9 + 5e-6*abs(single)` in normalized units;
otherwise the run fails. `irradiance_scattering_step` defaults to `1e-3` (allowed
range `1e-6` to `1e-2`). This is a numerical order extraction, not a separate
native first-order flux output. Native parameter Jacobians undergo the same
extrapolation. Multiple scattering is the original diffuse flux minus first
order; residuals are not clipped. Ordinary total, direct, actinic and radiance
outputs retain the baseline calculation.

QS uses scalar probe solves even for polarized totals: with the unpolarized
solar source and Lambertian surface exposed here, first-order Stokes I depends
only on P11. This avoids weak-scattering vector eigenmode degeneracy. The
baseline total and the resulting multiple-scattering residual retain the
requested polarization treatment.

QS now also honors disabled delta-M scaling in its first-order extinction and
linearized extinction setup, keeping the direct beam fixed during extraction.
At exactly ground level its upward flux is reconstructed from the Lambertian
boundary condition (albedo times total downward irradiance); previously the
zero-length upward path left that flux empty. This also corrects upward actinic
flux and its atmospheric Jacobians at that boundary. Standard input now pads
zero phase moments to the required quadrature order for unscaled isotropic
aerosols.
Run `python Examples/testing/validate_irradiance_orders.py` for an independent
single-layer first-order solution, surface reflection and Jacobian checks.

## Standard spherical-correction geometry

The Standard spherical-correction radiance driver converts observer-referenced
TOA angles to the ground-referenced angles required by the native SC master.
This conversion applies to both forward and profile-linearized calculations.
Exactly vertical SC rays are evaluated at an internal off-vertical angle of
`1e-5` degrees to avoid divisions by zero in the native three-point geometry
and interpolation. Public geometry and angularly integrated flux geometry
remain unchanged. AC uses its separate geometry handling.
Regression and aerosol-Jacobian diagnostics are in
`Examples/testing/investigate_toa_jacobians.py`; saved comparisons distinguish
the original SC geometry from the corrected geometry.
Standard uses a Fourier convergence tolerance of `1e-9`; the active QS Mark3
forward/profile solvers retain `1e-12`. Both require two successive converged
terms for aerosol calculations. These settings remove a reproduced SC AOD-sweep
outlier and improve aerosol Jacobian convergence. They can increase runtime,
particularly for Standard's binary128 calculations. This is a Fourier stopping
criterion, not a guarantee of `1e-12` accuracy of the complete RT solution.
The matched tolerance benchmark found `1e-9` adequate for Standard SC with an
8–9% runtime saving versus `1e-12`. QS optical-aerosol Jacobians changed by up
to 3.7 ppm at `1e-9`; `1e-10` also narrowly missed the predefined 1 ppm check
and saved only about 2% for the tested Jacobian cases, so QS retains `1e-12`.
See `Examples/testing/fourier_tolerance_comparison/decision.txt` for the tests.
The separate native polarized-Jacobian issue in nearly aerosol-free layers
remains: the tighter eigenvalue-degeneracy threshold is still experimental.
See `Examples/testing/toa_jacobian_investigation/report.txt` for the original
investigation and numerical evidence.

## Solar-beam geometry

Plain Standard (mode 1) and Standard SC both use pseudo-spherical solar-beam
attenuation (`TS_DO_PLANE_PARALLEL = .false.`). The SC option separately enables
the radiance line-of-sight correction; disabling SC no longer forces flat-layer
solar attenuation. The diffuse multiple-scattering solver remains based on a
plane-parallel formulation. See `Examples/testing/standard_profile_investigation/`
for the finite-difference and isolated geometry-control evidence.

## Aerosol extinction floor

Active aerosol profiles have a native-input extinction floor of `1e-8 km^-1`
per type, at that type's profile wavelength. This reduces the observed
near-Rayleigh Jacobian jump while adding very little optical depth. Inactive
types and the aerosol-free dummy retain the previous `1e-10 km^-1` guard.
The Python profile objects are not modified; native Jacobians describe the
floored atmosphere, so finite-difference checks must perturb from the effective
floor rather than from a smaller, clipped value. This is a numerical workaround,
not a general repair of eigenvalue degeneracy. The targeted two-aerosol tests
are saved in `Examples/testing/aerosol_floor_1e8_validation/`.

## Actinic flux output

The wrapper supports actinic fluxes at the observer altitude. Full vertical
actinic-flux profiles have been removed from the Python API and are no longer
requested from the copied QS backend.

The copied Standard VLIDORT `vlidort_pars.f90` sets `MAXLAYERS = 101` for
the wrapper's default vertical grids and `MAX_USER_LEVELS = 1` for single-level
output.
