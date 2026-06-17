# Parallelization audit — runtime data-parallelism opportunities

Scope: where data-parallel execution (`std::execution` via the `gbs/execution.h`
macros, SIMD, CUDA) can speed up hot paths in the library, and where it is
**unsafe** (sequential dependency / determinism) or **not worth it**. This audit is
**read-only on the algorithms** — it adds only this report and the bench harness
`bench/bench_parallel.cpp`. Fixes land via the focused issues listed at the end.
Mirrors the `INTERP_AUDIT` / `APPROX_AUDIT` playbook. Pre-v1.0 perf line (epic #44).

Harness: `bench/bench_parallel.cpp`, header-only, gcc-15 `-O3` (Release), `double`,
`dim=3`, libstdc++ + TBB. Best-of-N wall time. Machine: 64 hardware threads (the
scaling table also reports 1–32 threads, the realistic dev/CI range).

## TL;DR

1. **Multi-point evaluation is the headline win and it is free.** Every bulk
   evaluator (`values`, `d_dms`, `d_dm2s`, surface `d_dmus`/`d_dmvs`/…) is a
   `std::transform` over an independent parameter list writing to a pre-sized
   output — the textbook parallel case. They currently run **serially**; one of
   them even carries the comment *"Parallelization doesn't seems to be worth"*
   (`gbs/bscurve.h:67`). **Measured, that comment is wrong**: adding the policy
   gives **2× at N=1 000, 6–12× at N≥10 000, up to 21×** on the heavier offset
   kernel — and the result is **bit-identical** to the serial output (independent
   writes, no reduction). Drop-in `GBS_PAR_EXEC`, fully portable. → **issue #1**.

2. **Grid sampling leaves the win on the table by growing the output.**
   `offset_points` (`gbs/offsets.h:7`) and `discretize(Surface)`
   (`gbs/bssanalysis.h:8`) build their result with `push_back` / `insert` into a
   shared container, which forces a serial (outer) loop. Pre-sizing to `nu·nv` and
   a single `GBS_PAR_EXEC` `std::transform` over the flattened grid unlocks the
   same eval speedup (offset kernel measured **6–21×**). → **issue #2**.

3. **A lot is already parallel and correct** — `make_points`, the four
   `transform()` pole maps, the `extrema_*` bracketing `min_element`s,
   `extrema_curve_curve_lst`, `closest_curve`, the arc-length `std::transform`s,
   TFI grid + projection, mesh-area `std::reduce`. These are fine; **do not
   double-parallelize** from callers (e.g. `make_points` is already `par`, so a
   caller looping over curves must keep the inner eval as-is or parallelize only
   the outer level).

4. **The meshing core stays serial, by design.** Boyer–Watson / Delaunay
   insertion and edge-recovery flips mutate a shared triangulation in place and are
   determinism-sensitive (#48). **Do not parallelize.** Laplacian smoothing
   (`gbs-mesh/smoothing.h`) is the one meshing exception — a Jacobi double-buffer
   that is genuinely safe per-vertex (→ optional **issue #3**).

5. **SIMD is a weak lever here, measured.** `-march=native -mfma` leaves the scalar
   eval path unchanged within noise — the evaluator is branchy (`find_span`) and
   dependency-bound (de Boor), and the AoS `dim=3` layout doesn't fill a register.
   The only real SIMD lever is build flags letting Eigen vectorize the interp/approx
   *solves* (the dominant cost post #34/#65); a SoA evaluator is a speculative
   refactor. See the dedicated section below.

6. **Portability is the easy part.** The chosen mechanism everywhere is
   `std::execution` behind the `GBS_PAR_EXEC` macro (`gbs/execution.h`): on
   libstdc++ it needs TBB linked (the project already does, `if(UNIX)`), on MSVC it
   uses the native parallel STL, and on Apple libc++ — which ships no parallel
   policies — the macro **degrades gracefully to serial**. So every recommendation
   below is source-portable; the *speedup* is realized on Linux/MSVC and is a
   correctness-preserving no-op on libc++.

## Measured (bench_parallel)

Serial (`std::execution::seq`) vs parallel (`std::execution::par`) on the real
kernel shapes. `determinism` compares the two output arrays element-by-element.

### [1] Curve `value(u)` — the `make_points` / `values` kernel (deg 3, 60 poles)
| N | seq ms | par ms | speedup | determinism |
|---|--------|--------|---------|-------------|
| 1 000     | 0.076   | 0.038  | 2.0×  | bit-identical |
| 10 000    | 0.735   | 0.099  | 7.4×  | bit-identical |
| 100 000   | 3.686   | 0.497  | 7.4×  | bit-identical |
| 1 000 000 | 44.45   | 4.62   | 9.6×  | bit-identical |
| 4 000 000 | 177.2   | 39.63  | 4.5×  | bit-identical |

### [2] Curve `d_dm(u)` — the `d_dms` kernel `C'(u)/|C'(u)|`
| N | seq ms | par ms | speedup | determinism |
|---|--------|--------|---------|-------------|
| 1 000     | 0.088  | 0.034  | 2.6×  | bit-identical |
| 10 000    | 0.881  | 0.125  | 7.1×  | bit-identical |
| 100 000   | 8.801  | 0.757  | 11.6× | bit-identical |
| 1 000 000 | 89.32  | 11.59  | 7.7×  | bit-identical |

### [3] Surface `value(u,v)` (bidegree (3,3), 40×40 poles)
| N | seq ms | par ms | speedup | determinism |
|---|--------|--------|---------|-------------|
| 1 000     | 0.086  | 0.036  | 2.4×  | bit-identical |
| 10 000    | 0.777  | 0.128  | 6.1×  | bit-identical |
| 100 000   | 7.650  | 0.755  | 10.1× | bit-identical |
| 1 000 000 | 77.97  | 10.93  | 7.1×  | bit-identical |

### [4] Surface offset point — `offsets.h` kernel (value + ∂u + ∂v + normal)
| N | seq ms | par ms | speedup | determinism |
|---|--------|--------|---------|-------------|
| 1 000     | 0.374  | 0.056  | 6.7×  | bit-identical |
| 10 000    | 3.632  | 0.312  | 11.6× | bit-identical |
| 100 000   | 36.68  | 2.116  | 17.3× | bit-identical |
| 1 000 000 | 363.8  | 17.42  | 20.9× | bit-identical |

The heavier the per-element kernel, the better it scales (more compute per byte
moved): the offset kernel reaches **20.9×**, the light curve-value kernel saturates
memory bandwidth earlier (the 4M row).

### [5] Thread scaling — surface offset kernel, N = 1 000 000
| threads | ms | speedup | efficiency |
|---------|------|---------|-----------|
| 1  | 360.5 | 1.00×  | 100% |
| 2  | 182.7 | 1.97×  | 99%  |
| 4  | 92.41 | 3.90×  | 98%  |
| 8  | 48.98 | 7.36×  | 92%  |
| 16 | 32.83 | 11.0×  | 69%  |
| 32 | 20.56 | 17.5×  | 55%  |
| 64 | 15.01 | 24.0×  | 38%  |

Near-linear to 8 threads; the realistic **4–8 core dev/CI box gets 4–7×**. Beyond 16
threads the light kernels are memory-bandwidth-bound, so efficiency tapers — there
is no point oversubscribing.

## SIMD / vectorization — a weak lever here (measured)

Thread-level parallelism is the win; **SIMD auto-vectorization buys ~nothing on the
hot evaluation kernels**. Re-running `bench_parallel` built with `-O3 -march=native
-mfma` vs plain `-O3`, the **serial** (scalar-path) columns are unchanged within
noise:

| seq ms | `-O3` | `-march=native -mfma` |
|--------|-------|------------------------|
| curve value, N=100k  | 3.69 | 3.67 |
| curve value, N=1M    | 44.4 | 44.1 |
| surface value, N=1M  | 78.0 | 84.9 |
| offset point, N=1M   | 363.8 | 358.6 |

Structural reasons the auto-vectorizer can't help the evaluator:
- **`find_span` is a branchy binary search** — not vectorizable.
- **de Boor / `ders_basis_funs` is a triangular recurrence** of length `p+1` (= 4 at
  degree 3) with a level→level dependency — too short and too serial for intra-point
  SIMD.
- **AoS layout** (`points_vector = vector<array<T,dim>>`): `dim=3` does not fill a
  vector register (3 doubles ≈ ¾ of an AVX-256 register), so vectorizing the
  `+`/`*`/`cross`/`norm` array ops is counterproductive (remainder handling).

The natural SIMD axis ("many points in lockstep") is exactly the axis the
**thread**-level `GBS_PAR_EXEC` already captures, at far lower code cost. The genuine
SIMD levers are narrow:
- **Eigen solves (interp/approx)** — the dominant cost post #34/#65 — vectorize
  *internally*, gated on build flags. So the real lever is **compiling with
  `-march=native` / `-mavx2 -mfma`** (a flag, not code); header-only ⇒ the consumer
  owns that choice, and distributed binaries need runtime CPU dispatch or an AVX2
  baseline.
- **`par_unseq` instead of `par`** on the eval transforms: marginal — the per-element
  body (one B-spline point) is itself not vectorizable, so it ≈ `par`. Not worth it.
- **A structure-of-arrays evaluator** would let intra-batch SIMD pay *only* when many
  points share a span (dense sampling of one curve), at the cost of divergent
  pole gathers otherwise. Large refactor, speculative ROI — not recommended pre-v1.0.

## Determinism (guard-rail, cf. #48)

- **Eval / transform kernels — bit-identical, proven by the bench.** They write to
  distinct, pre-indexed output slots; there is **no reduction**, so the result does
  not depend on the policy, thread count, or scheduling. Parallelizing these
  introduces **zero** fp-order risk. This is the safe core of the audit.
- **`std::min_element` with `GBS_PAR_EXEC`** (already live in `extrema.h:77,145,293`,
  `bscanalysis.h:661,683`): on a *tie* the parallel version may return a different
  equal-minimum element than the serial scan, which can change a Newton/NLopt
  **warm-start seed**. The refined result still converges, but the seed choice is
  not guaranteed identical across thread counts. This is pre-existing, not
  introduced here — flagged so it is on the record.
- **`std::reduce` with `GBS_PAR_EXEC`** (mesh area, `halfEdgeMeshQuality.h:76`): fp
  addition is non-associative, so the sum differs in the last ULPs across core
  counts. It is a reported **metric**, not an input to a geometric predicate, so it
  is acceptable — but **no parallel reduction must ever feed meshing predicates**.
- **Meshing mutation (Boyer–Watson, edge recovery)** is order-dependent and #48
  pins exact-predicate determinism; parallel insertion would break cross-platform
  reproducibility. Do-not.

## Classified inventory

Verdicts: **PAR-DONE** (already parallel, correct) · **PAR-WIN** (safe, not yet
parallel, high value) · **PAR-MINOR** (safe but low value) · **SIMD/NA** (no
thread-level win; SIMD-only or single-shot) · **DO-NOT** (sequential dependency /
determinism).

| Area | file:line | Verdict | Mechanism | Gain | Portability |
|------|-----------|---------|-----------|------|-------------|
| Curve `values(u_lst,d)` | `gbs/bscurve.h:64` | **PAR-WIN** | add `GBS_PAR_EXEC` to the `transform` | 2–10× (measured) | portable (libc++→serial) |
| Curve `d_dms` / `d_dm2s` | `gbs/bscurve.h:159,177` | **PAR-WIN** | add `GBS_PAR_EXEC` | 2–12× (measured) | portable |
| Surface `d_dmus`/`d_dmvs`/`d_dmu2s`/`d_dmv2s` | `gbs/bssurf.h:316,336,356,376` | **PAR-WIN** | add `GBS_PAR_EXEC` | ~ as [3] (6–10×) | portable |
| `offset_points` (grid) | `gbs/offsets.h:7` | **PAR-WIN** | pre-size `nu·nv` + `GBS_PAR_EXEC` transform | 6–21× (measured) | portable |
| `discretize(Surface)` | `gbs/bssanalysis.h:8` | **PAR-WIN** | flatten to pre-sized grid + one `par` transform | ~ as [3] | portable |
| `make_points` / curve `discretize` | `gbs/bscanalysis.h:505,527` | **PAR-DONE** | `GBS_PAR_EXEC` transform | — | portable |
| `transform()` pole maps ×4 | `gbs/transform.h:16,33,47,64` | **PAR-DONE** | `GBS_PAR_EXEC` | — | portable |
| Arc-length params | `gbs/bscanalysis.h:148,212` | **PAR-DONE** | `GBS_PAR_EXEC` | — | portable |
| `extrema_*` bracketing | `gbs/extrema.h:77,145,293` | **PAR-DONE** | `min_element` `par` | — | tie → seed (see Determinism) |
| `extrema_curve_curve_lst` | `gbs/extrema.h:241` | **PAR-DONE** | `par` transform; **inner solver stays serial** | — | portable; don't double-parallelize |
| `closest_curve` / `closest_p_curve` | `gbs/bscanalysis.h:661,683` | **PAR-DONE** | `min_element` `par` | — | tie → seed |
| TFI grid + projection | `gbs-mesh/tfi.h:616,663` | **PAR-DONE** | `for_each`/`transform` `par` | — | portable |
| Mesh area | `inc/topology/halfEdgeMeshQuality.h:76` | **PAR-DONE** | `reduce` `par` (`…AreaPar`) | — | fp-noise metric only |
| Laplacian smoothing | `gbs-mesh/smoothing.h:35` | **PAR-MINOR** | Jacobi double-buffer per-vertex; `err_max` via `par` reduce | medium (iterative mesh loop) | portable; reduction is a scalar tol |
| Interp/approx row assembly | `knotsfunctions.ixx`, `bssinterp.h` | **PAR-MINOR / not worth** | banded fills write distinct rows | low — #34/#65 made assembly cheap; the **solve** (Eigen) dominates | n/a |
| Interp/approx **solve** (Eigen) | `bscinterp.h`, `bssinterp.h`, `bscapprox.h` | **SIMD-via-flags** | Eigen auto-vectorizes internally | gated on `-march`/AVX (the dominant cost post #34/#65) | build-flag; runtime dispatch for binaries |
| Per-point evaluator (`find_span` + de Boor) | `gbslib.ixx`, `bscurve.h` | **SIMD/NA** | branchy search + short triangular recurrence | ~0% from `-march` (measured) | SoA refactor only, speculative |
| `transformpoints` elementwise | `gbs/transformpoints.h` | **SIMD/NA** | single-point affine; SIMD-only | low | n/a |
| Single-point evaluators (`value`, `d_dm`, iso `isoU`/`isoV`, `CurveOnSurface`) | various | **SIMD/NA** | one point per call; parallelize at the caller | — | n/a |
| Param generation (`uniform_distrib_params`, `deviation_based_params`) | `gbs/bscanalysis.h:339,431` | **DO-NOT** | stateful `generate` / recursive list refine | — | sequential |
| `dev_from_points` | `gbs/bscanalysis.h:43` | **DO-NOT** | `u0` warm-start chain + scalar reductions | — | sequential dependency |
| OCCT `etrema_CS`/`etrema_PC` | `gbs-occt/intersections.cpp:16,43` | **PAR-MINOR** | growing `std::list`; pre-size first; OCCT solver per item | low; OCCT-bound | OCCT/Linux |
| Boyer–Watson / Delaunay insertion | `inc/topology/tessellations.h:265,492` | **DO-NOT** | in-place mesh mutation; order-dependent (#48) | — | — |
| Edge-recovery flips | `inc/topology/edgeRecovery.h:120` | **DO-NOT** | iterative topology mutation | — | — |

## Existing-parallelism inventory (one line each)

- `gbs/execution.h` — the only abstraction: `GBS_PAR_EXEC`/`GBS_SEQ_EXEC` expand to
  `std::execution::par,`/`seq,` when `__cpp_lib_execution` is defined, else to
  nothing (serial). Aliases `gbs::par_exec`/`gbs::seq_exec`.
- ~33 active `GBS_PAR_EXEC` call sites (transforms, reduces, min/max_element);
  ~11 **commented-out** ones (`bsctools.h:305`, `curvescheck.h:19,41`,
  `smoothing.h:38`, `fromjson.h:84,112,150`, `edge.h:183`, `bscbuild.h:127`,
  `mshedge.h:78`) — candidates to re-evaluate, mostly low value.
- **No** direct `tbb::`, **no** `#pragma omp`, **no** `std::thread`/`async`/`mutex`.
  All CPU parallelism is opt-in `std::execution`.
- `gbs-cuda/` — optional (`GBS_USE_CUDA=OFF` by default), a single
  `basisfunctions.cuh` (basis + derivative eval, `float`, Thrust). Out of scope for
  the CPU wins above.

## Recommended issues (high-value, safe wins)

1. **Parallelize the bulk evaluators** (curve `values`/`d_dms`/`d_dm2s`, surface
   `d_dmus`/`d_dmvs`/`d_dmu2s`/`d_dmv2s`): add `GBS_PAR_EXEC` to their
   `std::transform` and delete the stale `gbs/bscurve.h:67` comment. Drop-in,
   bit-identical, measured 2–12×. *Lowest risk, highest leverage.*
2. **Pre-size and parallelize grid sampling**: `offset_points` (`gbs/offsets.h:7`)
   and `discretize(Surface)` (`gbs/bssanalysis.h:8`) — replace the growing
   container + serial outer loop with a pre-sized `nu·nv` buffer and one
   `GBS_PAR_EXEC` transform over the flattened grid. Measured 6–21× on the offset
   kernel. Acceptance: output identical to the current order.
3. *(optional, lower priority)* **Parallelize Laplacian mesh smoothing**
   (`gbs-mesh/smoothing.h`): the Jacobi double-buffer is safe per-vertex; reduce
   `err_max` with a `par` reduction (scalar convergence tol, fp-noise acceptable).

## Reproduce
```
scripts/build_bench.sh bench_parallel
```
(Configures `bench/` Release with Ninja, builds `bench_parallel` — which links
`TBB::tbb` for libstdc++ parallel policies — and runs it.)
