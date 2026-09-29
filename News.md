# GBS Release Notes
v0.5.0 — 29 sep 2026

Behaviour changes
* loft: the skin-direction (v) parameters now default to the averaged chord-length of NURBS Book §10.3 instead of a uniform parametrization (#71) — lofted surfaces differ from 0.4.x for the same input
* `find_span` now follows NURBS Book A2.1 exactly (half-open convention `k[s] <= u < k[s+1]`, #42/#46): correct one-sided derivatives at knots of multiplicity > 1

Performance
* evaluation: allocation-free O(p²) evaluator and one-pass multi-order derivatives (#30); span reduction enabled for derivatives (~12x faster curve derivatives); per-call eval ~4x faster than OCCT
* bulk evaluators are parallel above a size threshold (#88, 12–30x over OCCT's scalar loop); grid sampling (`offset_points`, surface `discretize`) pre-sized and parallel (#89)
* Python: bulk evaluators return numpy arrays without copy (#97, ~10x at 1M points)
* interpolation: banded / sparse / separable solves (#34), then direct band assembly + no-pivot band LU (#96): ~10x on an 800-point build, now faster than OCCT and scipy at every tested size
* approximation: banded one-pass assembly and structured least-squares solve (#65, 2.5–9x vs OCCT)
* loft: v-system factorized once, batched pole solves (#40)
* new `build_batch` (`gbs/execution.h`): builds N independent interpolations/approximations in parallel, exceptions propagated to the caller (#91)

New features
* rational surface derivatives (NURBS Book A4.4, #33)
* spine-guided rational loft, and `loft_approx` (well-posed loft by approximation, #58/#63)
* vectorized derivatives w.r.t. curvilinear abscissa: curves `d_dms` / `d_dm2s`, surfaces `d_dmus` / `d_dmvs` / `d_dmu2s` / `d_dmv2s`
* Python bindings: `d_dm` / `d_dm2` (curves) and `d_dmu` / `d_dmv` / `d_dmu2` / `d_dmv2` (surfaces) accept a list of parameters
* `pygbs.gbs.__version__`; the installed CMake package now carries its version (`find_package(GBS 0.5)`)

Fixes
* `remove_knot` faithful to NURBS Book A5.8 for `num > 1` (#59)
* `refine_approx` stop criteria and pole-count guards (#65)
* rational `loft` with spine compiles and interpolates; rational `flat_v` overload returns a surface (#40/#58)
* hardened Bowyer-Watson cavity construction in the 2D Delaunay mesher (#48)
* gbs-occt: OpenCASCADE 8 support (7.x still supported), 2 runtime divergences fixed (#47)

Build / tooling
* tests migrated from GoogleTest to doctest; test executables are no longer installed
* tuning constants centralized in `gbs/gbsconstants.h` (#38)
* optional precompiled headers (`GBS_USE_PCH`), coverage build (`GBS_COVERAGE`)
* CI on Linux, Windows, macOS arm64 (+ Intel, non-blocking) and OCCT 7/8
* the in-repo conda recipe is removed: gbs is packaged by conda-forge (gbs-feedstock)
* audit reports and benchmarks under `bench/` (NURBS Book fidelity, interpolation, approximation, parallelization, OCCT/scipy/geomdl comparison)
11 oct 2022
* surfaces iso parametric curves
08 mar 2022
* Gordon surface
* constrained surface approximation/interpolation
* add gcc support
* tfi mesh
* increase python wrapping
10 sep 2021
* add name support for iges export
* add function object
* extend Python biding
* add jupyter render

22 may 2021
* add curve and surfaces extension
* add surfaces of revolution
* add igest export support
* many fixes and improvements
![Screencast](img/SurfaceExtension.png)
11 apr 2021
* curves creation from json
* transfinite mesh
* curvature base discretization
* curvilinear abscissa
* curve reparametrization
* cn curves connections
* surface approximation
* curves 2d offset
30 nov 2020
* Add many utilities to connect and join curves.
* Improve render configuration.
* Put everything together to build advanced shape
![Screencast](img/blade1.png)
31 oct 2020
* Loft surface creation

12 oct 2020
* Start python biding

08 oct 2020
* Geometries have now a 3d visualization

01 oct 2020
* Initial public release.