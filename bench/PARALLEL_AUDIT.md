# Parallelization audit — runtime data-parallelism opportunities

Scope: where data-parallel execution (`std::execution` via the `gbs/execution.h`
macros, SIMD, CUDA) can speed up hot paths in the library, and where it is
**unsafe** (sequential dependency / determinism) or **not worth it**. This audit is
**read-only on the algorithms** — it adds only this report and the bench harness
`bench/bench_parallel.cpp`. Fixes land via the focused issues listed at the end.
Mirrors the `INTERP_AUDIT` / `APPROX_AUDIT` playbook. Pre-v1.0 perf line (epic #44).

Harness: `bench/bench_parallel.cpp`, header-only, gcc-15 `-O3` (Release), `double`,
`dim=3`, libstdc++ + TBB. Best-of-N wall time. Machine: 64 hardware threads (the
scaling table sweeps 1–64 threads; a 4–8 core box is the realistic dev/CI range).

## TL;DR

1. **Multi-point evaluation is the headline win and it is free.** Every bulk
   evaluator (`values`, `d_dms`, `d_dm2s`, surface `d_dmus`/`d_dmvs`/…) is a
   `std::transform` over an independent parameter list writing to a pre-sized
   output — the textbook parallel case. They currently run **serially**; one of
   them even carries the comment *"Parallelization doesn't seems to be worth"*
   (`gbs/bscurve.h:67`). **Measured, that comment is wrong**: adding the policy
   gives **2× at N=1 000, 6–12× at N≥10 000, up to 21×** on the heavier offset
   kernel — and the result is **bit-identical** to the serial output (independent
   writes, no reduction). Fully portable — but `par` must be **gated behind a size
   threshold** (~1000; see Threshold), since below it the TBB overhead regresses
   small calls up to 25×. → **#88**.

2. **Grid sampling leaves the win on the table by growing the output.**
   `offset_points` (`gbs/offsets.h:7`) and `discretize(Surface)`
   (`gbs/bssanalysis.h:8`) build their result with `push_back` / `insert` into a
   shared container, which forces a serial (outer) loop. Pre-sizing to `nu·nv` and
   a single threshold-guarded `std::transform` over the flattened grid unlocks the
   same eval speedup (offset kernel measured **6–21×**). → **#89**.

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
   that is genuinely safe per-vertex (→ optional **#90**).

5. **Interp/approx parallelize at the BATCH level, not inside one build.** A single
   `interpolate()`/`approx()` is a sequential banded sparse solve — no `std::execution`
   axis worth taking (assembly is already cheap post #34/#65). But building **N
   independent** curves/surfaces (loft sections, multi-patch, fitting sweeps) is
   embarrassingly parallel: **measured 8–28×** (K=100–1000). Not exploited today (no
   batch API; no internal serial K-build loop). → batch helper / caller guidance, **#91**.

6. **SIMD is a weak lever here, measured.** `-march=native -mfma` leaves the scalar
   eval path unchanged within noise — the evaluator is branchy (`find_span`) and
   dependency-bound (de Boor), and the AoS `dim=3` layout doesn't fill a register.
   The only real SIMD lever is build flags letting Eigen vectorize the interp/approx
   *solves* (the dominant cost post #34/#65); a SoA evaluator is a speculative
   refactor. See the dedicated section below.

7. **Portability is the easy part.** The chosen mechanism everywhere is
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

## Interpolation / approximation builders — intra-build vs batch (measured)

The interp/approx *builders* have **two distinct axes**, and only one is worth
parallelizing:

**Intra-build (a single `interpolate()` / `approx()`) — no useful axis.** After
#34/#65 a build is: banded one-pass assembly → a banded **sparse LU/LDLT
factorization** (`bscinterp.h:138-142`) → a per-`dim` loop of back-substitutions
(`bscinterp.h:147-157`). The factorization is an inherently sequential numeric
kernel (no `std::execution` axis); the assembly is now cheap (the very point of
#34/#65) so parallelizing the row loop is Amdahl-capped and not worth it; the
`dim` back-subs (3 for `dim=3`) could be *batched* into one matrix solve
(APPROX_AUDIT A5) but that is a batching tidy-up, not thread parallelism. Eigen
vectorizes the solve internally (build flags, see SIMD above).

**Batch (N independent builds) — embarrassingly parallel, MEASURED.** Building many
curves/surfaces — loft sections, multi-patch models, a fitting sweep — is fully
data-parallel: each `interpolate()`/`approx()` is independent. `bench_parallel`
[7a]/[7b], a `std::transform` over the batch, seq vs par:

| K builds | interpolate (80 pts, deg 3) | approx (400 pts → 50 poles) |
|----------|------------------------------|------------------------------|
| 10   | 1.4× | 2.7× |
| 100  | 15.8× | 8.7× |
| 1000 | 28.1× | 21.4× |

Each build is far above the per-task threshold, so the batch scales close to core
count. **This is a real interp/approx win that the eval-focused part of the audit
missed.** It is *not* currently exploited: the library has **no batch
interp/approx API**, and there is **no internal serial K-build loop** to fix (loft
does a single tensor build, not N independent ones) — so the opportunity is a batch
helper (or documenting that callers building many entities parallelize the outer
loop, with the same threshold guard). Determinism: each build writes its own object,
so the batch result is order-independent and bit-identical to the serial loop.

## Threshold — gate `par` on input size (measured)

The libstdc++/MSVC parallel algorithms **do not cheaply fall back to serial** for
small ranges: `std::transform(par, …)` still spins up the TBB task machinery, which
costs ~9–15 µs of fixed overhead here. So an *unconditional* `GBS_PAR_EXEC` would
**regress** every small call. The crossover sweep (`bench_parallel` [6a]/[6b],
per-call time, M inner reps to beat timer noise):

| N | curve `value` (lightest) | offset point (heaviest) |
|---|--------------------------|--------------------------|
| 10  | 0.04× — par **25× slower** | 0.33× |
| 50  | 0.16× | **1.00× (break-even)** |
| 100 | 0.25× | 1.64× |
| 200 | 0.36× | 2.63× |
| 500 | 0.72× | 4.39× |
| 1000 | **1.10× (break-even)** | 5.89× |
| 2000 | 1.62× | 7.40× |
| 5000 | 2.83× | 7.82× |

The crossover is **kernel-dependent**: ~50–100 elements for the heavy offset kernel,
~700–1000 for the light curve-value kernel (more compute per element ⇒ overhead
amortizes sooner).

**Recommendation — every `PAR-WIN` fix must guard `par` behind a size threshold**,
not call it unconditionally. Concretely, a small dispatch helper in
`gbs/execution.h` (this is the form the #88/#89/#91 fixes should take, not raw
`GBS_PAR_EXEC`; #88 lands the helper, #89/#91 reuse it):

```cpp
namespace gbs {
    inline size_t parallel_min_size = 1000;   // tunable; conservative for any kernel
    template <class It, class Out, class F>
    Out transform_threshold(It first, It last, Out out, F f) {
        if (size_t(std::distance(first, last)) >= parallel_min_size)
            return std::transform(GBS_PAR_EXEC first, last, out, f);
        return std::transform(first, last, out, f);          // seq below threshold
    }
}
```

- A single **conservative default of ~1000** guarantees no regression on *any* eval
  kernel (the lightest only breaks even at ~1000).
- Heavy kernels (offset, `d_dm2`) cross over near ~64–128, so a lower threshold —
  or a per-kernel value — captures the 100–1000 band (up to ~6× there) at a
  negligible absolute cost on light kernels (<10 µs in the 500–1000 band). Exposing
  `parallel_min_size` as a tunable lets a known-heavy workload opt lower.
- On libc++ `GBS_PAR_EXEC` is empty, so both branches are serial — the guard is a
  correct no-op there.

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
| **Batch** interp/approx (N independent builds) | `bscinterp.h`, `bscapprox.h`, `gbs/execution.h` | **PAR-DONE (#91)** | `build_batch` gated transform over the entity list | **3× @K=8 → 35× @K=1000 measured** | portable; batch overloads shipped, gate `parallel_batch_min_size`=8 |
| Single-build interp/approx (assembly + solve) | `bscinterp.h:122,138`, `bscapprox.h` | **not worth (intra)** | banded assembly cheap (#34/#65); sparse factorization is sequential | low — Amdahl-capped; `dim` back-subs could batch (A5) | n/a |
| Interp/approx **solve** (Eigen) | `bscinterp.h`, `bssinterp.h`, `bscapprox.h` | **SIMD-via-flags** | Eigen auto-vectorizes internally | gated on `-march`/AVX | build-flag; runtime dispatch for binaries |
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

Filed against epic #44. Sequence **#88 first** — it lands the `transform_threshold`
helper that #89 and #91 reuse.

1. **#88 — Parallelize the bulk evaluators** (curve `values`/`d_dms`/`d_dm2s`,
   surface `d_dmus`/`d_dmvs`/`d_dmu2s`/`d_dmv2s`): route their `std::transform`
   through a **size-threshold helper** (see "Threshold" above) — *not* an
   unconditional `GBS_PAR_EXEC`, which regresses small calls up to 25× — and delete
   the stale `gbs/bscurve.h:67` comment. Drop-in, bit-identical, measured 2–12×
   above threshold. *Lowest risk, highest leverage; adds the shared helper.*
2. **#89 — Pre-size and parallelize grid sampling**: `offset_points`
   (`gbs/offsets.h:7`) and `discretize(Surface)` (`gbs/bssanalysis.h:8`) — replace
   the growing container + serial outer loop with a pre-sized `nu·nv` buffer and one
   threshold-guarded transform over the flattened grid. Measured 6–21× on the
   offset kernel. Acceptance: output identical to the current order.
3. **#90 *(optional, lower priority)* — Parallelize Laplacian mesh smoothing**
   (`gbs-mesh/smoothing.h`): the Jacobi double-buffer is safe per-vertex; reduce
   `err_max` with a `par` reduction (scalar convergence tol, fp-noise acceptable).
4. **#91 — Batch interp/approx helper.** Add a threshold-guarded `par` batch builder
   (or document that callers building many curves/surfaces parallelize the outer
   loop) — each build is independent, measured **8–28×** (K=100–1000). *Not*
   intra-build: a single solve has no useful parallel axis.

## #88 — post-implementation validation (gated library path)

#88 (PR #93) landed the shared `gbs::transform_threshold` helper (gate
`gbs::parallel_min_size = 1000`, `gbs/gbsconstants.h`) and routed every bulk
evaluator through it (curve `values`/`d_dms`/`d_dm2s`, surface `values` ×2 /
`d_dmus`/`d_dmvs`/`d_dmu2s`/`d_dmv2s`). This section archives a **fresh
`bench_parallel` run on the post-#88 tree** (gcc-15 `-O3` Release, libstdc++ +
TBB; 64-hw-thread machine; best-of-N wall time) — both to re-confirm the audit's
seq-vs-par numbers still hold and to validate the **shipped, gated** code path.

**Re-confirmation — sections [1]–[4] track the original audit within noise** (same
machine), e.g. curve `value` 7.6× at N=100k / 9.1× at N=1M, surface `value` 13.6× at
N=100k, offset point 16–18× at N≥100k, all `bit-identical`. The crossover sweeps
[6a]/[6b] reproduce too: the light curve-value kernel flips to `par` at **N=1000**
(0.70× at 500 → 1.12× at 1000), the heavy offset kernel near **N=50–100** — which is
exactly why the conservative default gate is **1000** and heavy workloads may lower it.

**New — [8] drives the actual shipped `crv.values()` (gated), not a kernel replica.**
It is timed per call against an explicit-serial reference, across N straddling the
gate. Below the gate the library call takes the **serial** branch (≈ seq cost — *no*
regression, vs the 0.04–0.26× an *un*gated `par` inflicts at the same sizes per [6a]);
at/above 1000 it takes `par` and wins; every row is bit-identical to serial:

| N | seq µs/call | gated µs/call | gated/seq | policy / result |
|---|-------------|---------------|-----------|-----------------|
| 100   | 3.78   | 4.38   | 0.86× | seq / identical |
| 500   | 19.10  | 22.05  | 0.87× | seq / identical |
| 999   | 37.61  | 43.52  | 0.86× | seq / identical |
| 1000  | 37.55  | 35.54  | 1.06× | par / identical |
| 2000  | 74.82  | 48.23  | 1.55× | par / identical |
| 10000 | 369.9  | 113.7  | 3.25× | par / identical |

Reading: the ~0.86× below the gate is the generic-evaluator overhead (virtual
`value()` + the `std::distance` check), **not** a parallel regression — the gate has
correctly kept those calls serial. The win switches on exactly at
`parallel_min_size`, confirming the measured threshold is placed at the break-even.

The `[8]` section was added to `bench/bench_parallel.cpp` for this validation; run it
with the same command below.

## #91 — post-implementation validation (gated batch builders)

#91 ships the **batch axis** the audit identified: building **N independent**
curves/surfaces is embarrassingly parallel, so two threshold-guarded batch
overloads were added and the outer loop routed through a new shared helper:

- `gbs::build_batch(inputs, build, min_size)` (`gbs/execution.h`) — the batch
  sibling of `transform_threshold`. It pre-sizes the output and runs the per-entity
  builder under `par` only when the **build count** is `>= min_size` (default
  `gbs::parallel_batch_min_size`), serial below. Each build owns its solver / matrix
  / poles (the #96 band-LU is a per-call stack local — **no shared mutable state**),
  so the result is **bit-identical** to the serial loop at any thread count.
- `gbs::interpolate(std::vector<points_vector>, p, mode)` and
  `gbs::approx(std::vector<points_vector>, p, n_poles, mode)` — batch overloads that
  wrap the single-entity builders through `build_batch`. (No change to `bssbuild.h`:
  the loft is a single tensor build, not an N-build loop — confirmed, so it is left
  untouched and there is no overlap with #79/#80.)

**Gate placement.** The batch gate is **separate and much lower** than the eval gate
(`parallel_min_size = 1000`): each "element" is a whole banded solve, so the TBB
spin-up (~9-15 µs) is amortized after a handful of builds. The fresh `[9a]` crossover
sweep (gcc-15 `-O3` Release, libstdc++ + TBB, 64-hw-thread machine; the shipped
overloads timed with the gate forced to MAX vs 0):

| K builds | interpolate (80 pts, deg 3) | approx (400 pts → 50 poles) | result |
|----------|------------------------------|------------------------------|--------|
| 2   | 0.91× | 1.65× | bit-identical |
| 4   | 1.78× | 2.58× | bit-identical |
| 8   | 3.08× | 4.26× | bit-identical |
| 16  | 5.30× | 9.20× | bit-identical |
| 32  | 7.18× | 15.2× | bit-identical |
| 64  | 9.78× | 18.8× | bit-identical |
| 128 | 11.4× | 18.7× | bit-identical |

Break-even is at **K ≈ 3–4** for both builders; the only sub-break-even row is K=2
interpolate (0.91×). So `gbs::parallel_batch_min_size` is set to **8** — a safe margin
above the measured crossover that still captures the 3×+ band (the heavier `approx`
crosses even earlier, so 8 is conservative for it). The original `[7a]/[7b]` replicas
confirm the win keeps scaling: **34–35×** at K=1000. A known-heavy workload may lower
the gate toward ~4.

**Production-gate check `[9b]`** drives the shipped overloads at the default gate (8):
K=4 takes the **serial** branch (no regression), K≥8 takes `par`, every row
bit-identical to the per-entity serial loop — the gate switches exactly at
`parallel_batch_min_size`.

Unit coverage: `tests/tests_batch_build.cpp` pins `gbs::parallel_batch_min_size` to 0
(always par) and SIZE_MAX (always serial) and checks the batch `interpolate`/`approx`
output **pole- and knot-exact** against the per-entity serial loop (independent of the
default gate), plus the generic `build_batch` helper and its empty-input case.

The `[9a]/[9b]` sections were added to `bench/bench_parallel.cpp`; run with the command
below.

## Reproduce
```
scripts/build_bench.sh bench_parallel
```
(Configures `bench/` Release with Ninja, builds `bench_parallel` — which links
`TBB::tbb` for libstdc++ parallel policies — and runs it.)
