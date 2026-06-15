# Test-coverage audit — core evaluation / interpolation (issue #50)

Baseline of source-based test coverage for the library core, with the untested
branches ranked and the cheap boundary gaps closed in this same PR. Motivation:
the last two bugs (#46 OOB at `u == last knot`, #48 non-robust `in_circle`) were
edge cases no test exercised — this enumerates the rest *before* they bite and
gives a safety net for the loft refactor (#40).

Tracked under epic #44. Analysis + test-adding only: no core algorithm logic is
touched, so it is parallel-safe with the code PRs.

## How it is measured

Clang source-based coverage (`-fprofile-instr-generate -fcoverage-mapping`),
off by default behind the `GBS_COVERAGE` CMake option so a normal build is
untouched. Instrumentation rides on the `gbs_testing` interface target, so only
the test executables are instrumented (not the python bindings / render libs).

Reproduce from scratch:

```sh
scripts/coverage.sh          # configure build-cov/ (clang-21) + build + run + report
scripts/coverage.sh --report # re-run + re-aggregate only (binaries already built)
```

or, in an already-configured `-DGBS_COVERAGE=ON` tree:

```sh
cmake --build build-cov --target coverage
```

Outputs land in `build-cov/coverage/`: `report.txt` (all `gbs/` headers),
`report_core.txt` (the priority perimeter below), `annotated_core.txt`
(line/branch-annotated source), `merged.profdata`.

### Two things the harness does, and why

- **Runs are pinned to one core** (`taskset -c 0`). The evaluators use
  `std::execution::par` (TBB); with instrumentation the per-region profile
  counters become cache-line-contended across all workers, so a millisecond
  parallel sweep balloons to *minutes* of 60-core thrash. Coverage is
  timing-independent, so serializing gives identical coverage in seconds.
- **Micro-benchmarks are excluded** (`-tce=*perf*,cpp_algo.par_vs_seq,tests_bscurve.curve_reparametrized`).
  They only re-run already-covered paths millions of times (10M curve evals,
  100M-element `std::transform`, a 1080-sample arc-length integration). Under
  instrumentation they dominate wall-clock and add no new coverage. Their one
  real cost to this baseline: `curve_reparametrized` was the only caller of some
  `bsctools.h` arc-length branches, which is why `bsctools.h` branch coverage
  (73.5%) lags its line coverage (97.3%).

Scope of the run: the eval/interp/loft test executables
(`tests_basisfunctions`, `tests_bscurve`, `tests_bssurf`,
`tests_curve_derivatives`, `tests_surface_derivatives`, `tests_curves`,
`tests_surfaces`, `tests_bsinterp`, `tests_knotsfunctions`, `tests_bssbuild`,
`tests_bsctools`). Coverage of a header is the union over all of these.

## Baseline — priority perimeter

`before` = master at the start of this PR; `after` = with the edge-case tests
added below. Region / line / branch %.

| File | regions | lines | branches | before→after (region) |
|------|--------:|------:|---------:|----------------------|
| `basisfunctions.ixx` | 81.0% | 88.2% | 76.1% | 80.3 → **81.0** |
| `bscurve.h`          | 92.9% | 95.5% | 73.3% | 89.4 → **92.9** |
| `bssurf.h`           | 88.0% | 91.2% | 70.5% | 85.4 → **88.0** |
| `bscinterp.h`        | 97.7% | 93.1% | 95.5% | 96.6 → **97.7** |
| `bssinterp.h`        | 86.2% | 86.0% | 74.3% | 84.5 → **86.2** |
| `knotsfunctions.ixx` | 90.6% | 93.7% | 83.1% | 90.6 (—) |
| `bssbuild.h`         | 95.8% | 97.2% | 89.3% | 95.8 (—) |
| `bsctools.h`         | 92.1% | 97.3% | 73.5% | 92.1 (—) |
| **TOTAL (perimeter)**| **88.6%** | **93.0%** | **78.3%** | 87.5 → **88.6** |

For context, the whole `gbs/` tree (all headers reachable from the run) sits at
**84.2% region / 81.5% line / 76.0% branch**. The lowest secondary files are the
non-B-spline curve types (`curveonsurface.h` 50%, `curvecomposite.h` 58%,
`curveline.h` 54%, `curvecircle.h` 69%) — out of this audit's eval/interp scope.

## Ranked gaps and disposition

Legend: **[fixed]** test added in this PR · **[intentional]** correct path left
uncovered on purpose, with rationale · **[defer #40]** belongs to the loft work.

### Tier 1 — eval / span / basis

1. **High-degree (p > 24) recursive fallback** — `basisfunctions.ixx`
   L386–391 (curve value), L435–441 (curve all-orders ders), L469–470 (rational
   ders, keyed on derivative *order* n > 24), L548–560 (surface ders), L756–764
   (surface value); mirrored in `bssinterp.h` L254–258 and `knotsfunctions.ixx`
   L872–873 / L899–901 (`fill_basis_row[s]`). **[intentional]** Every fast
   evaluator falls back to the recursive `basis_function` above
   `bspline_stack_max_degree = 24`. That kernel is **O(2^p)**, so a degree-25
   test costs ~26·2²⁵ ≈ 10⁹ ops *per evaluation* — it would re-introduce exactly
   the kind of multi-minute test this audit removed. Real CAD/CAM degrees are
   ≤ ~7; this branch is a "never silently break" guard, not a likely bug source.
   Best closed (if ever) by a single low-`u`-count smoke test gated out of the
   default run, not added here.
2. **Degenerate weight → origin** in `weight_projection` (L808–816, the zero-
   weight branch). **[fixed]** `weight_projection_zero_weight_is_origin`.
3. **`find_span` out-of-domain "no progress" break** (L306–308). **[fixed]**
   `find_span_out_of_domain_below` (u below the first knot returns span `p`).
4. **Degree 0 / degree 1 evaluation** never went through the fast path.
   **[fixed]** `eval_degree_0_and_1` (piecewise-constant + control polyline,
   value and derivative).
5. **Out-of-domain evaluation throws** — `bscurve.h` L603 (`BSCurve::derivatives`),
   L643/L651 (rational value/derivatives). Only `BSCurve::value` was pinned.
   **[fixed]** `tests_bscurve.eval_out_of_bounds_throws`.
6. **`check_curve` validation** — `bscurve.h` L213–214 (unordered knots),
   L219–220 (size equation). The `np ≥ p+1` guard was already tested
   (`ctor_rejects_too_few_poles`); these two siblings were not. **[fixed]**
   `ctor_rejects_unordered_and_wrong_size_knots`.
7. **Pole / knot mutation guards** — `bscurve.h` L472/L487 (`pole(id)` OOB),
   L458 (`movePoles`), L499 (`copyKnots`), L556–557 (`changeBounds(k1,k2)`).
   **[fixed]** `pole_and_buffer_guards`.
8. **Surface OOB / validation** — `bssurf.h` L438/L440/L442 (ctor: pole count,
   U knots, V knots), L1093/L1095 (`value` U/V OOB). **[fixed]**
   `tests_bssurf.ctor_rejects_invalid`, `eval_out_of_bounds_throws`.

   Still open (cheap, low value — secondary): surface **`derivatives()` OOB**
   (L829/L831), **iso-curve OOB** (L1022/L1043) and the **`value(...,du,dv>0)`
   "not implemented"** path (L1103); the base-class `Curve`/`Surface`
   `derivatives` defaults (`bscurve.h` L121, `bssurf.h` L249–251) are only
   reached by non-B-spline curve types outside this scope.

### Tier 2 — interpolation

9. **Interpolation argument guards** — `bscinterp.h` L98 (`build_3pt_tg_dir` < 3
   pts), L301 (`interpolate(Q,p,mode)` p ≥ #pts), L356 (`interpolate(bound,
   bound,cstr,p)` n < p+1). **[fixed]** `interpolation_argument_guards`.
10. **Surface-interpolation guards** — `bssinterp.h` L303 (out-of-bounds
    constraint), L308 (size not multiple of `n_polesv`). **[fixed]**
    `surf_interpolation_argument_guards`. Still open: L315/L319 (mu/mv knot-shape
    throws), L27/L111 (`extract_U` / direct `build_poles` size errors) — narrow
    lower-level guards.
11. Dense vs sparse threshold (`interp_sparse_threshold = 100`), banded vs
    scattered, separable surface path — **already covered**
    (`large_point_interpolation_banded`, `scattered_constr_points_sparse_path`,
    `grid_interp_separable_matches_dense`). No gap.

### Tier 3 — loft (feeds #40)

12. **Loft input validation** — `bssbuild.h` L162 (< 2 curves), L254 (a curve not
    touching the spine), L313 (U/V curves not intersecting). **[defer #40]** The
    `< 2 curves` guard needs a spine `BSCurve`; the other two need geometric
    near-miss inputs. Loft is being refactored in #40 (correctness → perf →
    dedup); these validation tests are best written against the post-refactor
    API rather than pinned to the current one now.

## Tests added in this PR

10 tests / 38 assertions, all on the priority perimeter, all fast (≤ 1 s each):

- `tests_basisfunctions.cpp`: `eval_degree_0_and_1`,
  `weight_projection_zero_weight_is_origin`, `find_span_out_of_domain_below`.
- `tests_bscurve.cpp`: `ctor_rejects_unordered_and_wrong_size_knots`,
  `eval_out_of_bounds_throws`, `pole_and_buffer_guards`.
- `tests_bsinterp.cpp`: `interpolation_argument_guards`.
- `tests_bssurf.cpp`: `ctor_rejects_invalid`, `eval_out_of_bounds_throws`,
  `surf_interpolation_argument_guards`.

No new **bugs** were found — every uncovered path is a correct guard or the
deliberate high-degree fallback, so nothing is filed as a follow-up bug. The
remaining open items above are either intentional (high-degree fallback),
low-value secondary guards, or deferred to #40.
