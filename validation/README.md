> Current default (2026-11-25): Standard/AC uses binary64 and AC shares separated
> MS kernels across both ray legs, with curved within-layer source integration
> (see the latest section at the end). Earlier timings
> below describe the former binary128 build unless explicitly labeled otherwise.
> The production switch follows the precision benchmark near the end of this file.

# AC limb implementation and validation

## Combined native implementation

Experimental correction `sphericity_correction=0` now selects
`aircraft_combined`. A public `run_rt()` call launches the Standard executable
once for its entire wavelength grid. The executable evaluates spherical FO
and profile derivatives through the bundled geometry/FO kernels, shares the
QS separated MS kernels across sparse samples, and integrates the source
along the spherical ray. Turning rays retain both legs; partial observer and
tangent layers are integrated directly on the original atmosphere grid.

The default is three source points per leg (configurable from 2 to 33).
Both legs now share the homogeneous solution and BVP factorization once per
Fourier order; beam-dependent work remains per sample. No QS master or fine
LOS quadrature is called. The older Standard source-only master mode
(`SPHERCORR_TYPE=S`) remains available for research comparisons. Legacy
sphericity modes retain their existing input checks and FO calculations.

FO extinction and its derivatives use the same delta-M scaling as MS. SS-only
requests omit MS and use unscaled FO, matching the existing QS convention.
Actinic flux and irradiance still use the existing native Standard flux path;
this change is a radiance correction, not a spherical flux model. Twilight
MS paths and an exactly surface-tangent horizontal observer remain unsupported
and fail explicitly. Horizontal limb support requires an observer above BOA.

`validate_combined_ac.py` compares with full QS and a clearly defined baseline:
**the identical spherical FO plus ordinary local Standard MS**, with a
pseudo-spherical solar beam but no spherical LOS MS integration. This extends
the FO-only curvature baseline to turning rays where the original Standard
FO geometry cannot run. Exactly horizontal baseline MS uses the limiting
angle 89.99999 degrees, since the plane-parallel LOS is singular there.
The validation reports both process wall/CPU time and native AC FO/MS time,
so process startup is not mistaken for RT cost.

```sh
python validation/validate_combined_ac.py --points 3 5 --public-api
python validation/validate_combined_ac.py --quick --vector --delta-m --property-jacobian
python validation/validate_combined_ac.py --case 15000 91 90 --points 3 --fd
python validation/validate_combined_ac.py --case 15000 91 90 --wavelengths 8 --points 3 5
```

## Combined-driver measurements (earlier binary128 implementation)

The final scalar, four-stream suite (477 nm, optical aerosol AOD 0.1, two
profile-Jacobian families) gives these three-point results. Jacobian error is
max(abs(AC-QS))/max(abs(QS)) across the profile, not a pointwise ratio near zeros.

| Observer / VZA / RAA | AC radiance error | FO-only baseline radiance error | AC aerosol Jacobian error | Baseline aerosol Jacobian error |
|---|---:|---:|---:|---:|
| 15 km / 89 / 90 | -0.020% | +22.205% | 0.634% | 150.331% |
| 15 km / 91 / 90 | -0.065% | +17.133% | 1.164% | 94.340% |
| 15 km / 93 / 90 | -0.088% | +0.727% | 1.601% | 27.328% |
| 550 m / 90 / 90 | +0.600% | -0.290% | 12.324% | 6.626% |
| 550 m / 90.5 / 90 | +0.178% | -3.778% | 2.647% | 55.906% |

AC improves radiance, absorption and aerosol-extinction Jacobian agreement in
seven of eight tested geometries. The exactly horizontal 550 m case is an
exception: the layer-averaged source approximation is less accurate than the
simple baseline there. Increasing the number of angular samples does not cure
that vertical-source approximation. These are test results, not a uniform
accuracy guarantee; QS remains preferable for sensitive low-horizon work.

For an eight-wavelength scalar batch at 15 km / VZA 91 / RAA 90, with two
profile-Jacobian families, measured process CPU times are:

| Solver | CPU seconds | AC native FO seconds | AC native MS seconds |
|---|---:|---:|---:|
| QS | 2.734 | — | — |
| AC, 2 points/leg | 1.906 | 0.045 | 1.713 |
| AC, 3 points/leg (default) | 2.313 | 0.047 | 2.093 |
| AC, 5 points/leg | 3.375 | 0.046 | 3.143 |

Without Jacobians, the same eight-wavelength batch takes 0.578 CPU s for AC3
versus 0.719 s for QS (about 20% less); native FO totals 0.019 s. Enabling
Jacobians does not change these radiances.

Thus the default saves about 15% process CPU in the linearized spectral batch. Wall
seconds are 9.37 versus 9.82 for QS because the roughly seven-second execution
delay in this session dominates launches. FO is a direct native calculation;
its compact workspaces use Standard dimensions, and forward-only runs use
the regular FO kernel. Five points requires another MS batch and is slower
than QS here, for little accuracy gain. Standard retains binary128 arithmetic
to protect conservative-scattering Jacobians; QS uses binary64. This is not
a like-for-like precision benchmark.

Polarized runs remain slower: at 15 km / 91 degrees, with delta-M and three
Jacobian families (absorption, extinction, HG asymmetry), QS takes 0.83 CPU s,
AC3 takes 3.41 s, and AC2 takes 2.91 s. AC3 radiance error is -0.066%, extinction
Jacobian error 1.054%, and HG-parameter Jacobian error 0.202%. Reducing to two
points improves runtime but does not beat QS. A focused native profile locates
about 45% of MS time in profile-Jacobian BVP solves and 30% in homogeneous
solutions/derivatives. There are two master calls (one batch per leg), not a
master call per layer: the 544 reported LP BVP calls are internal layer/beam/
Fourier derivative solves. The remaining cost is native MS algebra, not Python
I/O or FO. Sharing homogeneous work between the two legs and a faster stable
linearized core would be needed for a substantial further polarized speedup.
The contribution of precision alone has not been isolated experimentally.

The combined extinction Jacobians pass independent forward finite differences
at two step sizes in four layers, including the split observer parent and
nearly Rayleigh layers (maximum relative discrepancy 2.2e-7). Independent
HG-parameter perturbations also pass with vector RT and delta-M (relative
errors below 2e-5 in the checked cases). Public run_rt() launch-count checks
pass for all eight final-suite geometries: exactly one RT driver invocation,
with the same radiance as direct native execution. The compact FO optimization
changes saved radiances/Jacobians by at most 2e-11 relative; the original two-
executable prototype agrees to 2.4e-12. Forward and linearized radiances agree
for the eight-wavelength batch. See combined_final/verification.json.

The separate public-API regression also passes delta-M aerosol finite
differences, repeated wavelengths, actinic flux and irradiance (including
Jacobians), SS-only runs, and upward/nadir FO comparison with QS. Run
python validation/validate_combined_workflows.py; its outputs are in
combined_workflows/.

Artifacts (all in the ignored Examples/testing/limb_hybrid/ tree):

* combined_final/: final-executable measurements, eight geometries, default
  three points, CSV, plots, arrays, hashes, and one-call API assertions.
* combined_suite/: three-versus-five-point comparison before FO workspace
  reduction; retained as performance-investigation history.
* combined_batch8/: spectral throughput, two/three/five source points.
* combined_fd/: independent extinction finite differences.
* combined_forward8/: forward-only spectral throughput.
* combined_vector/, combined_vector2/: polarization, delta-M and phase-function
  derivatives. Comparisons expose Stokes I, not all Stokes outputs.
* combined_profile/: native polarized timing breakdown; reproducible with
  python validation/profile_combined_ac.py.

## Earlier prototype investigation (historical)

The legacy AC limitation is structural. Its geometry helpers obtain the BOA
viewing angle from `sin(alpha_BOA) = b / (R + h_BOA)`, where `b` is the ray's
impact parameter. When `b > R + h_BOA`, that angle does not exist. The original
AC implementation also interpolates integrated MS source terms between
BOA/aircraft/TOA geometries and traverses a monotonic sequence of layers.
A genuine limb ray instead has two legs and visits some layers twice.
Changing the input angle by a tiny amount cannot fix that case without
changing the physical ray substantially.

The current public `sphericity_correction=0` sparse approximation already
handles upward-looking rays whose backward extension misses BOA. It still
falls back to QS for horizontal and turning rays, because both its old MS
integration and the Standard spherical FO interface assume a monotonic,
surface-referenced ray. Legacy AC (`=2` below TOA) additionally requires an
interior layer-boundary observer. These public rules have not been changed.

## QS defect found and fixed during this investigation

The reference comparison exposed a separate QS geometry bug in
`FOMS_Geometry_Masters.f90`, in the limb branch where camera and tangent are
in the same layer (`nt == nCamera`). The far-side endpoint was assigned the
layer-top altitude but the camera's angle and tangent distance. Consequently
the endpoint did not lie on the ray, and a portion of the path was assigned
to the wrong optical layer. Camera-on-boundary cases hid the error.

The active QS core now computes the far-side angle from the layer-top radius
and updates the tangent distance before traversing the layer above. This is
a geometry correction, not an angle clamp or a convergence adjustment. It
affects both forward and linearized calculations through their shared geometry.

At 550 m, the old QS radiance jumped by 11.1% between VZA 90 and 90.001
degrees; FO alone jumped by 26.9%. Increasing path quadrature from 4 to 16
did not remove it. With the fix and 16-point quadrature, total radiance is
0.06583010205 at 90 degrees and 0.06582804628 at 90.001 degrees, a smooth
decrease of 0.00312%. See the generated `qs_horizon_fix.png`.

`python validation/check_qs_geometry.py` checks constant impact parameter
and independently computed spherical-shell chords through the public native
geometry master. All 27 cases pass. Repeating the same test against the
pre-fix geometry source failed five cases, including interior-layer cameras
at 550 m and 15 km. The corrected QS aerosol Jacobians also pass forward
finite differences at two step sizes, including the camera/tangent parent
layer. Results under `corrected_*` use the repaired QS reference; earlier
result folders retain the investigation history and must not be treated as
current accuracy measurements.

## Implemented research modes

`Fortran/vlidort_standard/sparse_aircraft.f90` now also supports a ray with
impact parameter above BOA. The Standard driver accepts two internal modes:

* `aircraft_limb_raw`: spherical shell segments, direct Cartesian local
  solar/viewing geometry, and separate source sampling for the two legs.
* `aircraft_limb_ms`: the same calculation with internal optical-layer splits
  at observer and tangent altitude. Sublayers retain their parent's optical
  properties; native sublayer derivatives are summed back into the original
  input-layer Jacobians. No atmospheric profile interpolation is performed.

**Both modes return MS only in `radiances.out` and its Jacobian files.** The
driver prints that explicitly. They are research interfaces, not public Python
backend choices. Neither mode calls QS or silently falls back to it.

For each leg, 2–33 local Standard solves sample a layer's effective source
`J = M / (1 - exp(-tau/mu))`. The formal integration uses the actual spherical
shell distance, `I_out = exp(-k ds) I_in + (1-exp(-k ds)) J`, retaining the far
leg's radiance and all derivatives across the tangent. The MS source and its
derivatives are interpolated in `log(mu)`. The zero incoming boundary is at
the far TOA; there is no artificial reflecting surface at the tangent.
The real surface still affects the local diffuse solutions.

The full atmosphere is retained for local MS calculations, including layers
below the tangent. Their influence on diffuse illumination and its Jacobians
must not be discarded. Derivatives of the source normalization, transmission,
and layer mapping are analytic. Finite differences are used only to validate.

The present local-source approximation uses layer-averaged sources rather
than pointwise diffuse emissivities. Source-count convergence therefore does
not establish accuracy. Vertical discretization and indirect illumination
from other layers can affect Jacobians much more than total radiance.

## First order and limitations

The harness makes a separate QS single-scattering call and adds its radiance
and Jacobians to the native Standard MS result. This uses QS **only for FO**,
not QS multiple scattering. It also runs full QS as an independent comparison.
The reported hybrid time includes both FO and MS calls, but excludes the
reference call and Python input generation.

Delta-M is disabled consistently in these comparisons: native QS SS-only
does not apply the same FO truncation used by a full delta-M calculation.
The harness checks that its FO result equals the FO component reported by
the full QS reference. Enabling delta-M in a future public hybrid requires
an explicit interface providing the matching scaled FO and its derivatives.

Research-mode restrictions are explicit: observer at or below TOA, ray misses
BOA, sunlit local geometry throughout the path, and radiance only. Twilight,
refraction, flux corrections, and above-TOA observer conversion are outside
this prototype. Flux requests are rejected rather than mixed with MS-only
radiances. A split can add up to two layers and must fit `MAXLAYERS=101`.
Vector runs currently compare Stokes I; the Python/native file interface does
not expose a full-Stokes comparison here.

## Measurements and recommendation

With corrected QS, the eight-case scalar suite using nine source points per
leg gives the following representative differences. Jacobian errors are
normalized to the largest absolute derivative in the reference profile.

| Observer / VZA / RAA | Radiance error | Absorption Jacobian error | Aerosol extinction Jacobian error |
|---|---:|---:|---:|
| 15 km / 89° / 90° | −0.017% | 0.044% | 0.082% |
| 15 km / 91° / 90° | −0.062% | 0.068% | 1.079% |
| 15 km / 93° / 90° | −0.086% | 0.093% | 1.589% |
| 550 m / 90° / 90° | +0.600% | 0.755% | 12.327% |
| 550 m / 90.5° / 90° | +0.178% | 0.158% | 2.626% |

In particular, the 15 km / 91° aerosol discrepancy previously attributed to
the hybrid was about 17%; it falls to about 1.1% with the corrected QS geometry.
Layer splitting is valuable near the ground: at 550 m / 90.5°, it removes
most of the unsplit prototype's large error. It is not a complete cure for
the horizontal-ray source approximation. On the finer grid, vector Stokes I
errors were about −0.068% at 15 km / 91° and +0.263% at 550 m / 90°, with
aerosol derivative errors about 1.13% and 7.62%, respectively.

The Standard MS derivatives agree with their own forward finite differences
to better than 1e-7 relative in the tested layers. This includes
the parent layer split internally. The remaining model differences should
not be interpreted as a failure of that derivative propagation. Flux,
SS-only, above-TOA, surface-intersecting, and twilight requests are rejected
as intended by the research mode. Existing public sparse workflows passed
their derivative, multi-wavelength, flux, and SS regression checks.

The research hybrid remains roughly twice as slow as QS in these small
four-stream scalar runs. Timings are indicative, not controlled benchmarks:
two executable calls are required, and Standard uses binary128 while QS uses
binary64. A QS FO-only call already takes almost as long as its full call in
this configuration. Increasing source samples from five to nine generally
changes radiance very little and does not cure the remaining Jacobian errors.

### Runtime audit: separate kernel cost from executable startup

The subsequent native profiling audit changes the interpretation of those
whole-executable timings. In this execution environment there is approximately
**7 seconds before the first Fortran statement**, even when both radiance and
flux are disabled. A minimal Fortran program with no model code also exhibits
this delay. Its CPU time is about 0.02–0.03 seconds. This is a pre-entry execution
delay, not FO integration, input export, or atmospheric model initialization.
Its underlying OS/session cause has not been identified; the measurements
should not be assumed to describe launches outside this session.

For the scalar 15 km / 91° case, four streams and two profile-Jacobian families:

| Timed native component | One wavelength |
|---|---:|
| QS geometry | 1–1.5 ms |
| QS FO including profile Jacobians | 11–13 ms |
| QS MS | about 0.26 s |
| Standard source solves, 5 samples/leg | about 0.55 s, 6 full-master calls |
| Standard source solves, 9 samples/leg | about 0.95 s, 10 full-master calls |

The short near leg has one effective source geometry in this case, hence
6 and 10 calls rather than 10 and 18. In an eight-wavelength batch, FO takes
about 0.032 s in total; QS MS takes 2.07 s, and the Standard source solves take
4.09 s (five samples) or 6.91 s (nine samples). Python input export is roughly
0.04–0.11 s for the active calculations; output logging/loading is around
milliseconds. Instrumented radiances match the uninstrumented reference.

There are therefore two independent performance problems in the prototype:

* It launches separate QS-FO and Standard-MS executables, paying the pre-entry
  delay twice. This is validation scaffolding, not an acceptable production
  implementation. **A normal Python SS+MS radiance simulation must launch one
  Fortran driver, with all FO, MS sampling, integration, and combination inside
  that driver.** Public `run_rt()` still dispatches to one executable; it does
  not use the two-executable research harness.
* Each Standard sample invokes a complete master, repeating workspaces,
  homogeneous solutions, and other setup, and calculating unused outgoing FO.
  QS Mk3 instead computes its homogeneous solution once per Fourier order and
  reuses it across segments, batching fine geometries in its particular solves.
  The original AC master also batches its two/three geometries into a single
  Standard call. The prototype does not yet retain that efficiency.

Standard's binary128 versus QS's binary64 is an additional comparison
confounder. The measured difference must not be attributed entirely to sample
count or solely to any one solver component. The proposed implementation
needs native FO plus batched MS sources and reuse of homogeneous/factorization
work, followed by a like-for-like precision benchmark.

The existing Standard spherical FO machinery is applicable to its supported
monotonic rays and is already used by legacy AC. For a surface-missing turning
ray, its current geometry construction starts at BOA and enforces
`b = (R+h_BOA) sin(alpha_BOA)`, so simply enabling its spherical correction
does not represent the two-leg path. Its FO integration can be reused with an
extended geometry adapter; alternatively, the small QS geometry/FO kernels can
be linked directly into the driver. Neither choice requires a second Python
call or a QS MS solve. The existing Standard delta-M FO helper also uses a
full-solve workaround; a direct FO interface should supply consistently scaled
optical inputs without doing discarded MS work.

`profile_limb_hybrid.py` builds instrumented copies under `tmp/`, preserving
production executables, and records parent wall time, time to the first Fortran
statement, process CPU time, native kernel timers, input-export time, and
output-I/O time. Native timers use `system_clock`; short CPU timings on Windows
are quantized and should not be used alone for millisecond-scale kernels.

```sh
python validation/profile_limb_hybrid.py
python validation/profile_limb_hybrid.py --reuse-build --wavelengths 8 --repeats 1
python validation/profile_limb_hybrid.py --startup-control
```

Raw measurements are in `Examples/testing/limb_hybrid/profile/`.

Recommendation: retain QS for production limb calculations, with the geometry
fix, and keep this hybrid explicitly experimental. A useful next iteration
would embed the QS geometry/FO kernels in the Standard executable, then
evaluate/interpolate pointwise MS emissivity rather than a normalized
whole-layer source. Reusing local diffuse solutions would also be necessary
to demonstrate a speed advantage. The current near-horizontal sampling can
still produce a small MS interpolation change when a tiny second leg appears;
the large QS discontinuity is fixed, but this prototype is not a guarantee of
smooth parameter derivatives across every geometry transition.

## Reproducing the checks

Build using the repository makefile, then run, for example:

```sh
python validation/check_qs_geometry.py
python validation/validate_limb_hybrid.py --points 5 9 --output Examples/testing/limb_hybrid/corrected_suite
python validation/validate_limb_hybrid.py --case 15000 91 90 --points 5 --fd --guards --output Examples/testing/limb_hybrid/current_native_checks
python validation/validate_limb_hybrid.py --quick --unsplit --output Examples/testing/limb_hybrid/current_raw
python validation/validate_limb_hybrid.py --quick --fine --vector --points 5 --output Examples/testing/limb_hybrid/corrected_vector_fine
python validation/validate_limb_hybrid.py --case 15000 91 90 --points 5 --qs-fd --output Examples/testing/limb_hybrid/corrected_qs_fd
python validation/validate_limb_hybrid.py --case 550 89.999 90 --case 550 90 90 --case 550 90.001 90 --points 5 --qs-fine 16 --output Examples/testing/limb_hybrid/corrected_horizon
```

Omitting `--quick` exercises eight geometries, including two azimuths in the
solar plane. `--points` sets the samples per leg. `--fine` changes the lower
20 km vertical grid from 2 km to 1 km spacing; this also changes the sampled
atmospheric profiles, so it is not a pure internal quadrature refinement.

All cases use 477 nm, four half-space streams, SZA 60 degrees, surface albedo
0.2, and an optical aerosol with AOD 0.1, exponential 2 km scale height,
SSA 0.95, and HG asymmetry parameter 0.7. Gas absorption is zeroed while
Rayleigh scattering remains. Generic absorption and aerosol extinction
Jacobians are compared. Active upper-atmosphere aerosol values retain the
library's existing extinction floor.

Executables and input/output files are isolated in temporary directories.
`manifest.json` records executable hashes and settings; a rebuild during a
run is detected. Results include component arrays (`*.npz`), native logs,
`comparison.csv`, `comparison.png`, and `jacobian_profiles.png`. The reported
Jacobian error is `max(abs(hybrid-QS)) / max(abs(QS))`, not the largest
pointwise relative error near a zero crossing. FD checks use two step sizes
and include nearly Rayleigh layers and layers below the tangent.

Generated results remain in the already ignored `Examples/testing/` tree.
The harness and these notes are outside that ignored tree so they can be
versioned with the implementation. No changes are committed or pushed by
the validation workflow.

## Standard binary64 versus binary128 (2026-11-25)

The isolated precision benchmark changes only Standard's real-kind promotion:
O2, fno-automatic, dimensions, source code, convergence, and optical inputs are
unchanged. At the time of this benchmark, production Standard was binary128; QS was
binary64/O3. No production executable was replaced. The build helper compiles
all 101 units, including LAPACK and FO, into a separate directory with module
dependencies respected.

Reproduce with:

```sh
python validation/compare_standard_precision.py --build --standard128-exe tmp/standard_binary128/DRIVER.exe
python validation/compare_standard_precision.py --case limb_scalar --wavelengths 128 --repeats 1 --skip-fd --output Examples/testing/standard_precision/batch128
python validation/summarize_standard_precision.py
```

Saved results are under `Examples/testing/standard_precision/`. The original
16-wavelength limb results are in that directory; the satellite continuation
is under `satellite/`. The initial run stopped on an exact-nadir NaN in the
existing binary128 executable; its log is retained. The continuation uses
1 degree off nadir instead. `combined_summary.json` combines the completed
cases, and the individual directories retain raw arrays, logs, timing and
accuracy CSVs, and executable hashes.

The atmosphere has 15 layers, four half-space streams, SZA 60 degrees, RAA
90 degrees, albedo 0.2, optical aerosol AOD 0.1, 2 km exponential scale height,
SSA 0.95 and HG g=0.7. Gas absorption is zero; Rayleigh remains. Wavelengths
span 470-490 nm. Jacobian runs request absorption and aerosol extinction
profiles. Delta-M is enabled for satellite cases and disabled for the limb
baseline, matching its earlier configuration.

Median process CPU seconds, two repetitions for 16 wavelengths:

| Case | Calculation | Standard64 | Standard128 | QS | QS / Standard64 |
|---|---|---:|---:|---:|---:|
| 15 km, VZA 91, scalar | Radiance | 0.1094 | 0.4766 | 1.2344 | 11.3 |
| Same | Radiance + Jacobians | 0.1406 | 1.1719 | 5.0625 | 36.0 |
| TOA, 1 degree off nadir, scalar | Radiance | 0.1250 | 0.5391 | 1.1953 | 9.6 |
| Same | Radiance + Jacobians | 0.2266 | 1.8359 | 6.0703 | 26.8 |
| TOA, 60 degrees off nadir, vector | Radiance | 0.2266 | 6.4688 | 1.5078 | 6.7 |
| Same | Radiance + Jacobians | 0.5313 | 20.7578 | 11.0156 | 20.7 |

The 128-wavelength scalar limb batch (one repetition) gives:

| Calculation | Standard64 | Standard128 | QS | QS / Standard64 |
|---|---:|---:|---:|---:|
| Radiance | 0.3125 | 2.8750 | 8.0000 | 25.6 |
| Radiance + Jacobians | 0.7969 | 9.5156 | 39.2969 | 49.3 |

These are CPU speedups, not wall-clock speedups. The environment still adds
about seven seconds before execution. For the large Jacobian batch, actual
wall times were 7.85 s (Standard64) versus 46.88 s (QS), about 6x. Small CPU
times are quantized on Windows; the repeated small batches and larger batch
are intended to expose fixed costs rather than establish universal ratios.

Satellite cases use the ordinary `standard` driver mode (spherical FO and
pseudo-spherical solar beam, no MS spherical correction). The limb case uses
`standard_native_fo`: native spherical FO plus ordinary local Standard MS,
without AC integration. It avoids the public route's QS fallback for a
surface-missing ray. This is a speed baseline, not an accurate substitute for
QS at limb: its radiance and aerosol Jacobian differences from QS remain
about 18% and 95% in this case.

Precision checks across the 16-wavelength cases:

* Radiances: largest binary64/binary128 difference is 2.12e-11 relative.
* Jacobians: largest difference normalized by the profile peak is 0.0203%.
  The largest local difference is 0.198% for limb absorption near 35 km;
  maximum local aerosol difference is 0.0178%. These are real local errors,
  not merely relative ratios at zero crossings.
* Independent forward aerosol-extinction finite differences at 477 nm,
  1 degree off nadir, use two step sizes and both a nearly Rayleigh layer
  (25-30 km) and an aerosol-rich layer (0-2 km). Standard64 differs by up to
  0.00517% in the nearly Rayleigh layer; QS differs by up to 0.00453%.
  Standard128 differs by less than 1.4e-8 relative across these checks.
  All three agree to better than 8e-8 relative in the aerosol-rich layer.

Thus global binary128 arithmetic is a substantial speed penalty: around
29-39x for the vector Standard tests. Binary64 forward radiances agree
closely, but small conservative-layer Jacobian roundoff persists. The tests
support evaluating a faster precision strategy; they do not certify all
stream counts, optical properties, geometries, fluxes, or Jacobian families.

The larger batch's largest local absorption precision difference is 0.202%,
again near 35 km. A separate exact-nadir reproduction at 470 nm fails with
nonfinite radiance in both binary64 and binary128 (`exact_nadir.json` and
associated logs). This is independent of the precision switch and remains
unfixed by this benchmark.

## Production default switched to binary64

The user accepted the measured local precision tradeoff. Standard/AC now defaults
to binary64; QS and Mie are unchanged. Use `mingw32-make -B PRECISION=quadruple`
in `Fortran/vlidort_standard` for an optional reference build, and rebuild with
`-B PRECISION=double` to restore the default. All units must use the same precision.
The previous production executable was preserved locally at
`tmp/standard_binary128/DRIVER.exe`. Precision comparisons now take an explicit
quadruple reference (`--standard128-exe`) rather than assuming production is
quadruple. Profiling helpers compile their replacement objects in binary64
to match the default production modules.

Production-switch verification passed: all 101 Fortran sources matched the
previously tested isolated binary64 build before its complete object/module
set was installed and relinked with the production makefile. Both precision
flag selections were checked using make dry runs. The installed polarized
Standard radiances and absorption/aerosol Jacobians match the isolated
binary64 results exactly. Combined AC routing, finite-difference aerosol
checks at both step sizes, multi-wavelength actinic flux and irradiance, and
upward/downward single-scattering checks against QS all passed without
relaxing tolerances. Results are under
`Examples/testing/standard_precision/production_double/`.

## Shared-work AC MS solver (2026-11-25)

The combined AC path now calls the supplied QS separated MS kernels directly
from `aircraft_ms_sources.f90`. It does not call the QS master or use QS's
fine LOS quadrature. Both ray legs share one homogeneous eigensystem and BVP
factorization per Fourier order. Each sparse geometry then solves its own
beam and BVP, reusing those solutions for every required source layer. User
homogeneous coefficients are cached across those layer evaluations.

The old AC master returned a plane-layer integrated MS source. The new adapter
reconstructs that same integral analytically from the homogeneous and beam
solutions, including derivatives of eigenvalues, optical depth, beam
attenuation and BVP constants. It handles real and complex eigenpairs and the
removable eigenvalue/viewing-secant coincidence. The existing normalization,
log-mu interpolation, split-layer Jacobian accumulation, surface reflection
and spherical chord integration remain unchanged. This is a performance
change, not a new AC accuracy approximation. More points reduce interpolation
error but need not remove the remaining whole-layer-source approximation.

Three points per leg remain the default; 2 through 33 are supported. The
shared workspace accommodates both legs together. All diffuse Fourier orders
are retained so weak Jacobian terms are not dropped by a geometry-dependent
convergence decision. Initialization touches active layers/beams rather than
maximum array capacity; that optimization also applies to the two supplied
QS setup modules. Native FO and all wavelengths remain inside one driver
invocation. Actinic flux and irradiance retain their existing Standard path.

Reproduce the checks with a saved pre-change binary64 AC executable:

```sh
python validation/check_ac_layer_integrals.py
python validation/validate_shared_ac.py --reference tmp/ac_shared_reference/DRIVER.exe --case 15000 91 90 --points 3 5 9 --wavelengths 64
python validation/validate_shared_ac.py --case 15000 91 90 --vector --delta-m --property-jacobian --property-fd --wavelengths 1 --points 3
python validation/validate_combined_workflows.py --output Examples/testing/limb_hybrid/shared_workflows
```

The comparison saves old/new/QS arrays, executable hashes, kernel timing logs
and peak-normalized radiance/Jacobian differences. Polarized tests exercise
vector transport but compare the exported Stokes I and its derivatives.
Use `--forward` for radiance-only, `--streams` for stream count and
`--profile` for setup/homogeneous/user/beam/linearized-BVP/source timers.
The standalone integral check compares 72 real/complex, thin/thick, upward/
downward cases and analytic derivatives with independent Gaussian quadrature.

### Final shared-work measurements

For 64 wavelengths (470–490 nm), 15 km / VZA 91 / RAA 90, four
half-space streams and the AOD 0.1 optical aerosol described above:

| Calculation | Points/leg | Previous AC MS seconds | Shared AC MS seconds | MS speedup |
|---|---:|---:|---:|---:|
| Scalar radiance + two profile Jacobians | 3 | 0.731 | 0.564 | 1.30x |
| Same | 5 | 1.036 | 0.742 | 1.40x |
| Same | 9 | 1.502 | 1.065 | 1.41x |
| Polarized transport + two profile Jacobians, delta-M | 3 | 2.977 | 2.372 | 1.26x |
| Same | 5 | 4.486 | 3.258 | 1.38x |
| Same | 9 | 7.387 | 5.478 | 1.35x |
| Scalar radiance only | 3 | 0.222 | 0.164 | 1.35x |
| Same | 5 | 0.300 | 0.203 | 1.48x |
| Same | 9 | 0.398 | 0.284 | 1.40x |

These are native MS wall timers, excluding FO and process startup; they are
single long batches, not statistically estimated universal speedups. Full
process CPU at three points was 0.828 s versus 5.766 s for optimized QS
(scalar Jacobians), 2.531 versus 12.547 s (polarized Jacobians), and 0.344
versus 4.078 s (scalar forward). The roughly seven-second launch delay remains.
A separate 33-point, 15 km / VZA 93 case passed: MS fell from 0.1600 to
0.1006 s (1.59x); its geometry and one-wavelength batch differ from the table.

The four-stream comparisons changed radiances by at most 1.23e-12 relative
and profile Jacobians by 1.12e-6 of their peak. Independent extinction checks
in three layers at two step sizes passed; the HG-property derivative agreed
with FD to about 9.7e-6 relative. Ground-looking, horizontal, upward and both
solar-plane limb cases passed. All 72 analytic-integral checks passed, as did
multi-wavelength actinic/irradiance Jacobians, SS-only comparisons and the
public API's single-launch assertion. QS setup optimization was bitwise
identical in the tested scalar and polarized three-wavelength cases.

The eight-stream polarized test exceeded the strict default old/new Jacobian
threshold: 3.39e-5 of the absorption-profile peak (0.00339%) and 2.85e-5 for
aerosol extinction, localized to the 30–40 km layer. Radiances differed by
1.69e-12. An independent extinction FD check with steps 1e-6 and 1e-7 found
about 0.156% local error in the new derivative versus 0.159% in the previous
one. Thus the new result retains a small derivative-accuracy limitation
already present in the reference; it is not exact agreement at eight streams.
The specifically documented eight-stream check uses
`--jacobian-tolerance 5e-5`; the default remains `1e-5`. No solver convergence
or precision settings were relaxed. The raw initial threshold failure and
independent FD results are preserved under `shared_streams8/`.

Results are in `Examples/testing/limb_hybrid/shared_scalar64/`,
`shared_vector64/`, `shared_forward64/`, `shared_maxpoints/`,
`shared_streams8_verified/`, `shared_property/`, `shared_geometries/`,
`shared_workflows/`, `shared_qs_setup/` and `shared_public/`.

The active-workspace initialization change also benefits QS itself. A separate
64-wavelength scalar Jacobian comparison against the saved QS executable gave
19.266 s before versus 5.766 s after in process CPU (3.34x), with bitwise
identical radiances and both Jacobian profiles. Wall time fell from 26.665 to
13.134 s. The AC/QS ratios above therefore use the faster QS reference, not
the older build. See `shared_qs_setup/scalar64_timing.json`.

## Near-horizon accuracy audit before curved-source integration (2026-11-25)

A dedicated check with the current shared-work binary confirms **two distinct
problems** at low observer altitude. Three angular source points remain the
default; no production solver or configuration was changed by this audit.
The scene is the existing scalar 477 nm test, observer 550 m, SZA 60 degrees,
RAA 90 degrees, four streams, delta-M off.

### The exact-horizontal bias is a vertical-source approximation error

At VZA 90 degrees, native AC first order agrees with the FO reported by QS
to 7.2e-12 relative. The multiple-scattering radiances are:

| Solver | MS radiance |
|---|---:|
| QS | 0.04092971289 |
| Standard local MS | 0.04073861819 |
| AC | 0.04132445582 |

Thus Standard is low while AC overshoots. AC converts each plane-layer
integrated source to an effective source by dividing by
`1-exp(-tau/mu)`, then applies it to a curved chord. Within a thick layer,
this does not reproduce the source's variation along the curved ray.
Accurate chord transmission alone does not guarantee a more accurate result.

To isolate this, the exported 0–2 km layer was subdivided while **copying its
pressure, temperature, gas concentrations, aerosol extinction, SSA and phase
function unchanged**. This preserves the piecewise-uniform physical optical
profiles rather than resampling a different exponential atmosphere.
Aerosol derivatives of the child layers were summed back to the original
parent layers. All cases retained three angular samples.

| Lowest-layer spacing | AC radiance error | AC aerosol-Jacobian peak error | Standard radiance error | Standard aerosol-Jacobian peak error |
|---|---:|---:|---:|---:|
| 2000 m | +0.5998% | 12.324% | -0.2904% | 6.626% |
| 500 m | +0.2427% | 4.499% | -0.3071% | 6.789% |
| 125 m | +0.03893% | 0.828% | -0.3074% | 6.792% |
| 31.25 m | +0.00565% | 0.411% | -0.3074% | 6.792% |

QS itself changes by about 0.0274% between the coarse and converged grids,
because its within-layer numerical integration also improves. The reference
is recomputed on each grid. This small reference change does not explain the
large reduction in AC's error. Increasing angular points from 3 to 33 on the
original coarse grid changes AC radiance by only about 0.000117% and leaves
the approximately 0.6% bias.

### A separate sampling discontinuity appears just above 90 degrees

| VZA | AC radiance error | AC aerosol-Jacobian peak error | Standard aerosol-Jacobian peak error |
|---|---:|---:|---:|
| 89.99999 | +0.59979% | 12.324% | 6.627% |
| 90.00000 | +0.59980% | 12.324% | 6.626% |
| 90.00001 | +0.45618% | 39.355% | 6.627% |
| 90.00010 | +0.38709% | 34.639% | 6.626% |
| 90.00100 | +0.30839% | 27.923% | 6.625% |
| 90.01000 | +0.25609% | 17.735% | 6.608% |

QS and the Standard baseline vary smoothly. The AC radiance changes by about
-0.143% between 90 and 90.00001 degrees, whereas QS changes by -0.0000299%.
This is finite output, not the former exact-angle NaN failure.

In `sparse_aircraft.f90`, the global sample range is taken from the minimum
and maximum effective viewing cosines of every traversed layer. Just above
90 degrees, the emerging tangent-to-observer segment is extremely thin, but
its almost-zero effective cosine resets the entire far-leg sampling grid.
The three cosines change from approximately

- at 90 degrees: [0.01067, 0.03839, 0.13819];
- at 90.00001 degrees: [8.74e-8, 1.10e-4, 0.13819].

An almost negligible new path therefore moves two samples away from the
main path. The 1e-9 km split tolerance can also leave this tiny tangent
interval inside its unsplit parent; no guard makes the sampling limit
continuous. Increasing angular samples reduces this artifact but is not a
reliable fix with the required default of three.

The next implementation should keep the three-node sampling continuous when
a vanishing second leg appears, and improve source evaluation along each
curved segment. Reusing the existing diffuse solutions for within-layer
source quadrature is a plausible way to address the second issue without
adding complete MS solves. A fixed VZA offset would merely move or mask the
problem and would not remove the exact-horizontal vertical-source bias.

Reproduction:

```sh
python validation/diagnose_ac_horizon.py
python validation/diagnose_ac_horizon_refinement.py
python validation/plot_ac_horizon_diagnosis.py
```

Data, full component/Jacobian arrays, logs and the figure are under
`Examples/testing/limb_hybrid/horizon_accuracy_audit/`.

## Curved-source AC integration (2026-11-25)

The production AC radiance path now fixes the two issues identified in the
preceding audit. The default remains three angular source points per leg,
double precision, four FO fine points and one Python-to-Fortran invocation.
No extra full MS solve is introduced by the within-layer quadrature.

The new `aircraft_chord_quadrature.f90` evaluates the existing homogeneous
and beam source expansions along the actual spherical chord. It integrates
their exponentials against the outgoing attenuation and supplies analytic
derivatives with respect to eigenvalues and layer optical depth. Angular
source values use local quadratic interpolation (linear with two points).
The geometry samples span actual local viewing cosines, including zero at
a tangent, so an emerging vanishing ray leg does not reset the far-leg nodes.

The original optical-layer grid is retained. Partial observer/tangent layers
are handled by the integration bounds, eliminating the artificial thin
layers and summation of child-layer Jacobians used previously. Homogeneous
solutions and BVP factorizations are still shared; each angular geometry
still requires its beam-dependent solve. A short second leg confined to one
original layer uses one representative angular solution as before.

Quadrature is eight-point Gauss-Legendre on subintervals selected from the
optical variation, with boundaries at angular interpolation stencil changes.
The outgoing optical-depth tail beyond 36 is negligible at binary64 accuracy.
The quadrature also resolves the beam's average solar secant. Optical
Jacobians differentiate the integrand; they do not use finite differences.
Radiance-only runs skip weight-derivative accumulation.

The independent `check_ac_chord_integrals.py` checks 144 real/complex,
thin/thick, vertical/partial/tangent, upward/downward integral cases and
their derivatives against independent high-order quadrature. It also
requires bitwise equality of forward and linearized integral values.
Maximum absolute difference was 8.4e-12. The older
`check_ac_layer_integrals.py` command delegates to the current check.

To compare the changed approximation with a saved previous AC executable,
use `validate_shared_ac.py --allow-model-change --reference PATH`; otherwise
that harness retains its strict old/new equivalence assertions for pure
optimization changes. Case names now retain enough digits to distinguish
90 degrees from 90.00001 degrees.

### Measured accuracy and cost

Against QS at 550 m, 477 nm, SZA 60 degrees and relative azimuth 90 degrees,
with three angular points and six FO/QS fine points:

| VZA | Previous radiance error | Corrected radiance error | Previous aerosol-J peak error | Corrected aerosol-J peak error |
|---|---:|---:|---:|---:|
| 89.99999 | 0.58256% | 0.0003666% | 12.3511% | 0.015164% |
| 90 | 0.58257% | 0.0003665% | 12.3513% | 0.015164% |
| 90.00001 | 0.43898% | 0.0003665% | 39.6576% | 0.015163% |

Jacobian errors are maximum absolute differences divided by the maximum
absolute QS Jacobian, not local relative errors near zero crossings.
The saved-data regression `check_ac_horizon_regression.py` asserts accuracy
and continuity; the largest adjacent-angle changes divided by each field's
peak are 3.13e-7 (radiance), 3.98e-7 (absorption), and 2.15e-6 (aerosol).
The plot is reproduced with `plot_ac_chord_regression.py`.

The six-fine-point result is a reference-resolution check, **not a changed
default**. With the default four fine points, the 64-wavelength horizon
comparison gives maximum radiance error 0.0229% and aerosol-J peak error
0.820%. The much smaller six-point differences indicate sensitivity to
within-layer quadrature in the reference comparison. At 15 km / VZA 91,
the default-four-point scalar errors are 0.0000415% and 0.0909%, respectively;
with polarization and delta-M they are 0.000109% and 0.0713%.
These are selected validation scenes, not general accuracy guarantees.

Paired 64-wavelength measurements (470–490 nm, four streams, three angular
points) used the saved pre-fix binary64 AC executable. Native process CPU
seconds, excluding Python preparation and process-launch delay:

| Case | Previous AC | Corrected AC | QS |
|---|---:|---:|---:|
| 550 m / 90, scalar with Jacobians | 1.531 | 1.672 | 13.047 |
| 15 km / 91, scalar with Jacobians | 1.688 | 1.594 | 9.047 |
| 15 km / 91, polarized + delta-M with Jacobians | 5.594 | 5.125 | 21.641 |
| 15 km / 91, scalar radiance only | 0.641 | 0.578 | 4.875 |

These initial measurements precede the final optimization that skips unused
weight-derivative sums in forward-only runs; the forward integral values are
bitwise unchanged. They show about 9% extra CPU in the horizon case and
6–10% less in the other cases. Timings are machine/load dependent. No broad
claim of equal cost at large angular sample counts is intended: the 33-point
horizon diagnostic adds appreciable quadrature cost. Three remains the default.

Artifacts: `Examples/testing/limb_hybrid/chord_horizon6`,
`chord_horizon64`, `chord_limb64`, `chord_vector64`, `chord_forward64`,
`chord_geometries`, `chord_property`, `chord_streams8`, and `chord_points`.

Additional checks passed for polarized delta-M scattering, eight streams,
forward/backward azimuths, upward and ground-intersecting views, and 2/5/33
angular points. Six independent aerosol-extinction finite differences
(three layers, two step sizes) passed. The HG asymmetry-property derivative
was -0.00389028656 analytically versus -0.00389024887 by finite differences.
Mixed radiance/actinic-flux/irradiance workflows, wavelength changes, and
single-scattering-only comparisons against QS also passed. Their artifacts
are in `chord_workflows`.

The final public-interface regression at 550 m / VZA 90.00001 passed:
exactly one Standard Fortran-driver invocation and radiance equality with
the direct combined-driver run to rtol 1e-12. With default four fine points,
AC radiance error is +0.01744% versus Standard native-FO/local-MS -0.29039%;
aerosol-J peak error is 0.6506% versus 6.6265%. See `chord_public`.
All checks in the curved-source validation sequence completed successfully.

## Default-grid limb altitude sweep (before surface-tangent fix)

Run `python validation/compare_limb_altitudes.py`, followed by
`python validation/summarize_limb_altitudes.py`. Outputs are in
`Examples/testing/limb_altitude_comparison/`: comparison.png/PDF,
differences_above_surface.png, comparison.csv, summary.json, manifest.json,
README.md and per-case logs/NPZ arrays.

This uses 20 existing layer boundaries from 0 to 15000 m on the unmodified
82-layer grid, VZA 90, SZA 60, RAA 90, 477 nm, scalar/four streams/delta-M,
default fine quadrature, and three AC angular points. Aerosol AOD is 0.1 at
477 nm with a 2 km exponential scale height, SSA 0.95 and HG g=0.7; surface
albedo is 0.2. Default gas absorption is retained. O4 SCDs use absorption
Jacobians, verified against the public get_scds calculation. Standard uses
the actual standard backend (FO spherical correction only).

All 60 runs completed, but **the AC surface-horizontal endpoint is invalid**:
the strict impact-parameter test chooses the surface branch, whose path
clipping at VZA 90 and height zero discards every atmospheric layer, so MS
is zero. Radiance is 65.803% low and O4 SCD 18.262% low. This comparison
records the current failure without changing production code.

For the 19 elevated cases, maximum absolute differences from QS are
0.000436% radiance / 0.000470% O4 SCD for AC, versus 27.804% / 55.216%
for Standard. Mean native process CPU time over those same cases is
3.088 s (QS), 0.136 s (Standard), 0.184 s (AC). Mean input-generation plus
native wall time is 10.291 / 7.294 / 7.324 s, respectively; launch overhead
dominates the latter. Measurements are single runs with rotated model order,
not repeated statistical benchmarks. One-time scene initialization and
interpreter startup are excluded from those per-case times.

## Surface-tangent AC fix and continuity validation

The AC driver now routes equality of ray impact radius and BOA radius into
the limb calculation. The limb guard accepts this tangency. No observer-height
or viewing-angle offset is introduced.

Nearby checks also exposed two related boundary cases, now corrected in
sparse_aircraft.f90: upward views use Cartesian local geometry to avoid
near-horizontal cancellation in the older BOA/TOA conversion; a downward
view at BOA still evaluates diffuse surface reflection when its outgoing
atmospheric path has zero length. Both changes reuse existing diffuse
solutions. Three angular source points, binary64 and one driver invocation
are retained.

Run `python validation/check_ac_surface_continuity.py` to reproduce 13 AC/QS
geometry pairs. The `--check-existing` option rechecks their saved results.
Cases include height 0, 0.001, 0.01, 0.1 and 1 m at VZA 90; offsets of
0.00001 degrees at heights 0 and 0.001 m; offsets of 0.001 degrees at 1 m;
and surface nadir / 1 m zenith checks. These deliberately off-grid heights
are continuity diagnostics only: the main comparison retains its 20 exact
default-grid layer boundaries.

At the exact surface horizon, corrected AC gives normalized radiance
0.0549553283241 versus QS 0.0549551696266 (0.000289% difference), and
O4 SCD 8.28021025980e43 versus 8.28018755362e43 (0.000274%).
All 13 pairs pass the 0.01% radiance/SCD agreement bound. Explicit angular
and altitude increment checks also pass (AC minus QS increment discrepancy
below 1e-5 of the reference field).

At exactly ground level, moving below the horizon sees the surface rather
than the sky. This physical discontinuity remains in both solvers; it is
not treated as a numerical continuity failure. The downward surface
radiance and O4 SCD agree with QS to about 3e-11 relative or better.

Updated main plots/data and the surface_continuity subfolder live under
`Examples/testing/limb_altitude_comparison/`. The new main summary includes
all 20 heights in its maximum-error and mean timing statistics; the previous
summary excluded the broken surface endpoint. The earlier audit above is
retained as historical evidence.

## Full-profile limb absorption and aerosol Jacobian comparison

Reproduce with `python validation/compare_limb_jacobians.py` followed by
`python validation/summarize_limb_jacobians.py`. Output directory:
`Examples/testing/limb_altitude_jacobians/`.

The scene, 20 aircraft heights on existing boundaries, 82-layer grid and
three backends match the corrected radiance/O4 altitude comparison.
Both absorption and haze extinction Jacobians are requested in one native
invocation per case. All 60 runs completed with finite full-layer profiles;
radiance and absorption derivatives reproduced the prior absorption-only runs.
No production solver code was changed.

Absorption is d(normalized radiance)/d(layer absorption optical depth).
Haze is d(normalized radiance)/d(layer extinction in km^-1), at 477 nm,
holding aerosol SSA and phase function fixed. These parameterizations have
different units; the extinction derivative includes the layer thickness.
Errors use 100*max_layer|J_model-J_QS|/max_layer|J_QS| independently at each
aircraft height, avoiding divisions by local Jacobian zeros.

Worst profile-peak differences across the full altitude sweep:
- Standard: 202.041% absorption (14 km aircraft, 14–14.5 km layer);
  46.9722% aerosol (15 km aircraft, 15–16 km layer).
- AC: 0.0006602% absorption (3.2 km aircraft, 3.2–3.4 km layer);
  0.0225461% aerosol (15 km aircraft, 15–16 km layer).

These are comparisons with QS, not independent proofs of absolute accuracy.
The large Standard absorption discrepancy occurs at the QS profile peak,
so it is not a near-zero relative-error artifact.

Mean process CPU times with both Jacobian families: QS 4.5555 s,
Standard 0.1492 s, AC 0.2234 s. Mean input-plus-native wall times:
11.7252 / 7.2990 / 7.3722 s. Each geometry was measured once with rotated
solver order; launch overhead remains substantial.

Outputs include:
- absorption/haze heatmaps: model values, common-scale differences and
  magnified AC differences; 0–20 km and full 0–70 km versions;
- selected_profiles.png and a 20-page all_profiles.pdf covering all heights;
- comparison.png/PDF for accuracy and CPU time, plus timings.png;
- all_profiles.csv (all layer values and differences), errors.csv
  (including worst-layer bounds), timings.csv, per-case NPZ/logs,
  summary.json, manifest.json and README.md.

## Off-horizontal radiance, O4 and full-profile Jacobian comparison

Reproduce with `python validation/compare_limb_angles.py` followed by
`python validation/summarize_limb_angles.py`. Outputs are under
`Examples/testing/limb_angle_comparison/`. The runner saves per-case
checkpoints; `--angles` accepts a subset, and the summarizer supports
`--metrics-only` for inspecting an incomplete sweep.

The eight VZAs are 80, 85, 88, 89, 91, 92, 95 and 100 degrees, each at
the same 20 boundary-aligned heights and with the same atmosphere as above.
Each invocation returns radiance and both full Jacobian profiles; O4 SCD
is computed from the absorption Jacobian. No production solver changes
were made for this comparison.

Actual Standard results are omitted where its public API would fall back
to QS: downward rays whose tangent point is above the surface. These are
marked unsupported, never plotted as Standard/QS agreement. For this grid,
Standard supports only 0 and 800 m at VZA 91, and 0 through 3200 m at
VZA 92; all tested heights are supported at the other requested angles.

The overview plots separate accuracy of radiance/O4 alone from accuracy
including absorption and aerosol Jacobians. Profile errors are normalized
by the peak magnitude of the corresponding QS profile, avoiding division
by local zeros. Each angle has raw values, timing/error curves, Jacobian
heatmaps, all 20 full-height profiles in a PDF, and layer-by-layer CSV data.

Native processes for the three solvers ran concurrently in isolated folders
for each geometry, with no simultaneous QS geometries. CPU times are logged
per process; the concurrent wall times are not serial speed benchmarks.
The manifest records executable hashes and the complete settings/grid.

All 447 supported runs succeeded; 33 Standard cases were unsupported.
Across all eight angles and 20 heights, AC maximum differences from QS
are 0.000309% radiance, 0.000305% O4, 0.000464% absorption Jacobians,
and 0.016741% aerosol Jacobians.

For radiance and O4 together at a 5% tolerance, Standard passes all sampled
heights at 80, 95 and 100 degrees; heights through 13.5 km at 85 degrees;
through 3.2 km at 88 and 89 degrees; and all of its supported sampled
heights at 91 and 92 degrees. Requiring both Jacobian profiles as well
reduces passing heights to none at 80–89 degrees, surface only at 91/92,
0 and 800 m at 95, and 0 through 3200 m at 100 degrees. These are sampled
scene-specific comparisons to QS, not universal or absolute accuracy bounds.

AC native CPU was approximately 17.5–24.1 times lower than QS by angle,
with both Jacobian families enabled. See angle_summary.csv for matched
CPU ratios and full height-specific threshold results.

## AC angular-sampling provenance audit

`python validation/audit_ac_sampling.py` builds an isolated logging-only variant
and repeats all 160 AC geometries. All saved arrays were reproduced exactly.
Observed samples: 123 cases used three, 33 turning rays used six (three per
leg), and four downward surface endpoints used one. No full QS MS master
was linked or called. Two-/five-point controls changed counts and results.
Evidence: `Examples/testing/limb_angle_comparison/sampling_audit/README.md`.
Comparison labels now explicitly say three angular points per leg.
