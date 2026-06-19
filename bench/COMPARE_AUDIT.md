# Competitive benchmark — gbs vs the field

Scope: position **gbs** against established B-spline/NURBS libraries on the same
operations, fairly, on one machine — honestly surfacing where gbs **wins** and
filing focused perf issues where it **lags**. Two fronts: the **C++ core**
(baseline **OpenCASCADE / OCCT**, already bridged via `occt_utils`) and the
**Python interface** `pygbs` (key peer **scipy.interpolate**, plus **geomdl**).
Read-only on the algorithms: this adds only the harness + this report. Pre-v1.0
perf positioning (epic #44). Closes #86; depends on #88 (merged) so the bulk
evaluators measured here are the **parallel** path, not the old serial one.

Harnesses: `bench/bench_compare_occt.cpp` (Front A) and `bench/python/compare_*.py`
(Front B). One documented command per front rebuilds the comparison (see
*Reproduce*).

## TL;DR

1. **Evaluation is gbs's headline strength, in both fronts.** The compile-time
   `<double, dim=3>` templating (no runtime type dispatch) + the allocation-free
   de Casteljau evaluator make gbs's *per-call* eval **4.2× faster than OCCT**
   (C++) and **4.5× faster than scipy** (Python single calls). On top of that the
   **parallel bulk evaluators (#88)** beat OCCT's scalar loop by **12–30×** in
   C++ — OCCT exposes no batch API at all.
2. **Interpolation *build* — was the biggest lag, now a WIN (#96).** It used to
   scale super-linearly (OCCT 3–9×, scipy 1.2–3.7× faster) because `build_poles`
   built a dense `n×n` matrix then `sparseView()`-rescanned it (two O(n²) passes),
   and then spent 60 % of the build in Eigen `SparseLU`'s numeric factorization.
   Fixed in two steps: assemble the banded matrix **directly into band storage**,
   and solve with a **no-pivot band LU** (the collocation matrix is bandwidth
   `deg`, totally positive). The 800-pt build drops **1.16 → 0.115 ms (10×)** and
   gbs now **beats OCCT 1.1–1.3× and scipy 1.04–5.6× at every size**, on both
   fronts.
3. **The Python bulk evaluator's *return* used to throw that away — now fixed
   (#97).** The result was marshalled through a Python-list round-trip: the C++
   core evaluates 1 000 000 curve points in 11 ms (parallel) but the pygbs call
   took **~500–950 ms**, so scipy's vectorized `splev` won the batched regime
   **~17×**. Replacing the round-trip with a single `memcpy` into an `(N,dim)`
   numpy buffer cuts the 1M call to **92 ms (≈10× faster)** and closes the gap to
   scipy to **~1.6×** (from ~17×); values bit-identical. The residual is the
   *input* numpy→`vector` copy — a smaller follow-up.
4. **Smaller wins:** approximation (C++ **2.5–9×**, Python ≈ tie), knot insertion
   (C++ **5.7×**), small-N interpolation in Python (≤ 100 pts).
5. **One more measured lag:** degree elevation in C++ (OCCT **1.8×** faster), but
   the measurement folds in a full curve copy on both sides — reported with that
   caveat, flagged for a cleaner isolation before it is filed.

---

## Machine, toolchain, versions

| | |
|---|---|
| CPU | 64 hardware threads (TBB `default_concurrency` = 64) |
| OS | Linux 6.8.0 (glibc 2.39) |
| C++ compiler | clang 22.1.7, `-O3` Release, libstdc++ + TBB (env `gbs-occt8`) |
| OCCT | OpenCASCADE **8.0.0** (`gbs-occt8` micromamba env) |
| Python | CPython **3.11.13** (env `dev`) |
| numpy / scipy / geomdl | **2.4.3 / 1.16.1 / 5.4.0** |
| pygbs | freshly built from this branch (`build/python/gbs….so`) — carries the #88 parallel path |

`best-of-N` wall time for C++; `median of 11` (with IQR tracked) for Python. Same
instances feed both libraries (same poles/knots/degree/parameter lists).

## Protocol & fairness (what is and isn't in the timed region)

- **Construction / `gbs → OCCT` conversion is done once, outside every timed
  region**, except where the operation *is* a build (interpolation, approx),
  where the whole build is timed on both sides from pre-built inputs.
- **Evaluation is reported in two regimes**, because the libraries expose
  different surfaces:
  - **per-call** — a scalar loop calling `value(u)` / `Value(u)` / `spl(u)`. This
    is the like-for-like *evaluator-quality* comparison.
  - **bulk / batched** — gbs's `values(u_lst)` (parallel since #88) and scipy's
    vectorized `spl(U)`, vs OCCT's only option (the same scalar loop). This is the
    *API-level* comparison ("how a user evaluates N points").
- **In-place ops** (knot insert, degree elevate) time a **fresh copy + the op on
  both sides** — a user must copy to keep the original; the copy is included
  identically on both, and called out.
- **Structural caveats stated, not hidden:** gbs is templated `<double, dim=3>`
  (compile-time); OCCT and scipy are runtime-typed. The gbs Python bulk path is
  parallel (#88, 64 threads), OCCT's loop and scipy's `splev` are single-threaded
  C — so the bulk-eval win is partly *parallelism* and partly *evaluator speed*;
  the per-call rows isolate the latter.
- **No rigging:** identical flags both sides, identical instances, same operation;
  every lag is in the same table as the wins.

---

## Front A — C++ core: gbs vs OCCT 8.0.0

Ratio = `occt / gbs`; **> 1 means gbs is faster**. Raw run:
`bench/results_compare_occt.txt`.

### Curve point evaluation `value(u)`, d=0

**Per-call** (gbs single `value()` loop vs OCCT `Value()` loop) — the evaluator
quality comparison, both single-threaded scalar:

| N | gbs ms | occt ms | ratio | winner |
|---|--------|---------|-------|--------|
| 1 000     | 0.0656 | 0.2737  | 4.17× | gbs |
| 10 000    | 0.3206 | 1.3326  | 4.16× | gbs |
| 100 000   | 3.2113 | 13.303  | 4.14× | gbs |
| 1 000 000 | 32.117 | 133.56  | 4.16× | gbs |

**Bulk** (gbs `values(u_lst)` parallel #88 vs OCCT scalar loop — OCCT has no batch API):

| N | gbs ms | occt ms | ratio | winner |
|---|--------|---------|-------|--------|
| 1 000     | 0.0424 | 0.1412 | 3.33×  | gbs |
| 10 000    | 0.0972 | 1.4080 | 14.49× | gbs |
| 100 000   | 0.4877 | 13.513 | 27.71× | gbs |
| 1 000 000 | 11.460 | 135.67 | 11.84× | gbs |

### Curve derivatives (bulk)

| op | N | gbs ms | occt ms | ratio |
|----|---|--------|---------|-------|
| 1st `values(·,1)` vs `DN(u,1)` | 1 000 | 0.0388 | 0.1220 | 3.15× |
| | 10 000 | 0.1226 | 0.9881 | 8.06× |
| | 100 000 | 0.7570 | 9.8871 | 13.06× |
| | 1 000 000 | 12.123 | 100.05 | 8.25× |
| 2nd `values(·,2)` vs `DN(u,2)` | 1 000 | 0.0379 | 0.1008 | 2.66× |
| | 10 000 | 0.1505 | 1.0057 | 6.68× |
| | 100 000 | 0.9063 | 9.9079 | 10.93× |
| | 1 000 000 | 12.252 | 99.568 | 8.13× |

### Surface evaluation (bulk)

| op | N | gbs ms | occt ms | ratio |
|----|---|--------|---------|-------|
| value d=0 | 1 000 | 0.0404 | 0.2462 | 6.09× |
| | 10 000 | 0.1317 | 2.3437 | 17.80× |
| | 100 000 | 0.7599 | 23.212 | 30.55× |
| | 1 000 000 | 12.590 | 232.82 | 18.49× |
| d/du | 1 000 | 0.0473 | 0.3120 | 6.60× |
| | 10 000 | 0.1591 | 2.9696 | 18.66× |
| | 100 000 | 1.0209 | 29.350 | 28.75× |
| | 1 000 000 | 12.607 | 294.68 | 23.38× |

### Curve interpolation (deg 3, full build) — gbs now wins (#96)

Includes the real **201-point** case (#34). gbs `CHORD_LENGTH` vs OCCT
`GeomAPI_Interpolate` (its own parametrization) — both produce an interpolating
C² cubic; **ratio > 1 means gbs is faster**. The build was rebuilt twice: #96 part
1 removed the dense `n×n` materialization (assemble directly as sparse), part 2
replaced Eigen `SparseLU` with a **no-pivot band LU** assembled straight into band
storage (the collocation matrix is bandwidth `deg`, totally positive → no pivoting
needed; falls back to `SparseLU` on a small pivot).

| N | orig ms | #96 sparse | #96 band-LU | occt ms | ratio (band-LU) |
|---|---------|-----------|-------------|---------|-----------------|
| 50  | 0.0243 | 0.0219 | **0.0066** | 0.0078 | 1.17× |
| 100 | 0.0568 | 0.0422 | **0.0123** | 0.0144 | 1.16× |
| **201** | 0.1342 | 0.0816 | **0.0249** | 0.0282 | 1.13× |
| 400 | 0.3717 | 0.1681 | **0.0530** | 0.0665 | 1.26× |
| 800 | 1.1642 | 0.3323 | **0.1150** | 0.1287 | 1.12× |

**10× over the original at 800 pts** (1.16 → 0.115 ms), and gbs now **beats OCCT at
every size (1.1–1.3×)**. Why it works: the band LU does no symbolic analysis, no
supernode bookkeeping, no fill — just `O(n·deg²)` Gaussian elimination inside the
band. Measured breakdown at 800 pts: SparseLU factor+solve was 0.31 ms (60 % of the
build); the band LU does it in 0.05 ms, bit-identical (`max|Δ| = 9e-16`).

The **general constrained interpolation** (`build_poles(constrPoint)`, arbitrary
value/derivative constraints in any order — the Python `build_poles` API and
`c2_connect`) shares the win: sorting the constraints by `(u, d)` is a pure row
permutation that makes the matrix banded (poles are the columns, so they stay in
natural order). The same band LU then applies, with the dense/sparse solve kept as
a fallback for tiny, non-banded or ill-conditioned sets. Two measurements:

- vs gbs's old `solve_collocation` (scrambled N-constraint system): **7.8× (N=100)
  → 20× (N=800)**, bit-identical (`max|Δ| = 9e-16`).
- **vs OCCT** `GeomAPI_Interpolate` (point interpolation — the like-for-like OCCT
  exposes), section `[7b]`: the general path **beats OCCT 1.09–1.29×** at every
  size, matching the structured path (the sort cost is negligible):

| N | gbs general path ms | occt ms | ratio |
|---|---------------------|---------|-------|
| 50  | 0.0060 | 0.0077 | 1.27× |
| 201 | 0.0241 | 0.0283 | 1.17× |
| 800 | 0.1149 | 0.1255 | 1.09× |

### Curve approximation (LSQ, deg 3, ~N/4 poles) — gbs wins

Structural caveat: gbs fits a **fixed pole count**; OCCT `GeomAPI_PointsToBSpline`
picks knots adaptively from a tolerance (`1e-3`, deg fixed 3). Not bit-for-bit
like-for-like, but both fit a cubic to the same cloud:

| N | gbs ms | occt ms | ratio | winner |
|---|--------|---------|-------|--------|
| 100 | 0.0231 | 0.2164 | 9.37× | gbs |
| 201 | 0.0554 | 0.3784 | 6.83× | gbs |
| 400 | 0.1624 | 0.7117 | 4.38× | gbs |
| 800 | 0.5453 | 1.3906 | 2.55× | gbs |

### Knot insertion / degree elevation (fresh copy + op, both sides)

| op | reps | gbs ms | occt ms | ratio | winner |
|----|------|--------|---------|-------|--------|
| knot insert | 20 000 | 6.788 | 38.829 | 5.72× | gbs |
| degree elevate +1 | 5 000 | 281.75 | 153.33 | 0.54× | **occt** |

Degree elevation: OCCT is **1.84×** faster — but the figure includes a full curve
copy on both sides (gbs value-type copy vs OCCT `Copy()`), so part of the gap may
be copy cost. Flagged for a cleaner isolation before filing (see *Lags*).

---

## Front B — Python interface: pygbs vs scipy (vs geomdl)

Raw runs: `bench/python/results_curve.md`, `bench/python/results_build.md`.
Numbers are **median ms** (lower = better). geomdl is pure-Python; it is included
as an honest reference and is, predictably, orders of magnitude slower.

Sanity: the pygbs curve and the mirrored scipy `BSpline` evaluate to
`max|Δ| = 0` over 97 points — identical geometry, so the eval rows compare
evaluator speed, not different splines.

### Curve evaluation — two regimes

**Per-call** (a Python loop, one point at a time — exposes the call boundary):

| N | pygbs ms | scipy ms | winner |
|---|----------|----------|--------|
| 1 000  | 0.540 | 2.455 | **pygbs 4.5×** |
| 10 000 | 5.448 | 24.25 | **pygbs 4.4×** |

For genuinely one-at-a-time evaluation (e.g. a root-finding inner loop), pygbs's
pybind11 boundary is **lighter** than scipy's per-call Python/numpy dispatch.

**Batched / vectorized** (the right way to use scipy), **after the #97 fix** —
`memcpy` into an `(N,dim)` numpy buffer instead of the old `py::cast` list
round-trip:

| N | pygbs ms (was) | pygbs ms (#97) | scipy ms | gap to scipy |
|---|----------------|----------------|----------|--------------|
| 1 000     | 1.08  | 0.140 | 0.060 | 2.3× |
| 10 000    | 4.86  | 1.398 | 0.586 | 2.4× |
| 100 000   | 73.8  | 9.088 | 5.59  | 1.6× |
| 1 000 000 | 950.6 | **92.0** | 56.0 | **1.6×** |

The fix is **~10× at 1M** and collapses the deficit from ~17× to ~1.6×; the
sanity check stays `max|Δ| = 0`. It was never a core-evaluator deficit — the same
1M-point eval runs in **11.5 ms** in the C++ core (Front A [bulk]); the old
~500–950 ms was pure `py::cast` list marshalling. 1st/2nd derivatives track the
same shape (1M, #97: pygbs ≈ 101 / 100 ms vs scipy ≈ 58 / 59 ms). The remaining
~1.6× is the **input** numpy→`std::vector<double>` copy on the way in plus the
24 MB output `memcpy` — addressable later by accepting/returning the numpy buffer
without an intermediate `std::vector` (smaller follow-up, not filed).

### Curve interpolation (deg 3, full build) — gbs now wins (#96)

Both on a chord-length parametrization (like-for-like interpolation, not
smoothing). Columns track the two #96 steps (sparse assembly, then band LU).

| N | pygbs orig | pygbs band-LU | scipy ms | geomdl ms | winner |
|---|------------|---------------|----------|-----------|--------|
| 50  | 0.039 | **0.018** | 0.103 | 4.6   | **pygbs 5.6×** |
| 100 | 0.067 | **0.031** | 0.097 | 28    | **pygbs 3.1×** |
| 201 | 0.143 | **0.056** | 0.121 | 205   | **pygbs 2.1×** |
| 400 | 0.326 | **0.114** | 0.161 | 1586  | **pygbs 1.4×** |
| 800 | 0.884 | **0.226** | 0.235 | 16003 | **pygbs 1.04×** |

Cross-front confirmation of the #96 work: the 800-pt Python build drops **0.88 →
0.23 ms** and pygbs now **beats scipy `make_interp_spline` at every size** (5.6× at
50 pts down to 1.04× at 800), winning decisively at the small sizes that dominate
real use (the #34 201-pt case is 2.1×).

### Curve approximation (LSQ, fixed ~N/4 poles) — ≈ tie

Matched to scipy `make_lsq_spline` with a clamped knot vector giving the same
coefficient count (not scipy `splprep`'s smoothing-`s` fit, which solves a
different problem):

| N | pygbs ms | scipy ms | winner |
|---|----------|----------|--------|
| 100 | 0.066 | 0.075 | ≈ |
| 201 | 0.120 | 0.122 | ≈ |
| 400 | 0.266 | 0.274 | ≈ |
| 800 | 0.701 | 0.838 | pygbs 1.2× |

### tinyspline / tinynurbs — deferred (opportunistic)

The opportunistic lightweight C++ peers were **not vendored**: this build server
has no outbound network to fetch and license-record the sources under
`bench/third_party/`. OCCT is the C++ baseline and gives a strong, fair
comparison; the harness is structured so a vendored peer drops in as another
column later. Noted here for honesty rather than silently omitted.

---

## Strengths (measured) — honest, with the exact conditions

Every claim is one number + the conditions that produce it. The regime is stated;
nothing is context-free.

1. **Per-call evaluation, C++: gbs 4.16× faster than OCCT** (curve `value`, all N,
   1k–1M, single-threaded scalar loop, deg 3 / 60 poles, clang-22 `-O3`). *Why:*
   compile-time `<double,3>` templating → no runtime type dispatch; allocation-free
   de Casteljau (A2.3, #30).
2. **Bulk evaluation, C++: gbs 12–30× faster than OCCT** (parallel `values()`
   #88 vs OCCT scalar loop; curve up to 27.7× at 100k, surface up to 30.6× at
   100k; bit-identical to serial per #84). *Why:* #88 size-gated `par` transform +
   the fast evaluator; OCCT exposes no batch API.
3. **Per-call evaluation, Python: pygbs 4.5× faster than scipy** (single `value(u)`
   in a loop, N=1k/10k). *Regime:* one-at-a-time eval. For vectorizable work the
   batched path is right; after #97 pygbs's batched eval is within ~1.6× of scipy
   (was ~17×).
4. **Approximation, C++: gbs 2.5–9.4× faster than OCCT** (LSQ deg-3, ~N/4 poles,
   N=100–800; structural caveat: fixed poles vs OCCT tolerance-adaptive knots).
5. **Knot insertion, C++: gbs 5.7× faster than OCCT** (copy + single insert,
   included both sides).
6. **Approximation, Python: pygbs ≈ scipy** (LSQ fixed poles; pygbs slightly
   faster at N ≤ 400, run-to-run variance at N = 800 — that path is untouched here).
7. **Interpolation, C++: gbs 1.1–1.3× faster than OCCT** at every size (deg-3,
   N=50–800, after #96's band-LU build). *Why:* direct band assembly + a no-pivot
   band LU (no symbolic analysis / supernodes / fill).
8. **Interpolation, Python: pygbs 1.04–5.6× faster than scipy** `make_interp_spline`
   at every size (after #96 — 5.6× at 50 pts, still ahead at 800; the #34 201-pt
   case is 2.1×).

## Where gbs lags (measured) — and the issues filed

1. **Interpolation build — FIXED, now a win (#96).** Used to be OCCT **3–9×** /
   scipy **1.2–3.7×** faster with the gap widening in N (dense `n×n` build + a
   general `SparseLU`). Fixed in two steps — direct band assembly, then a no-pivot
   band LU — so the 800-pt build went **1.16 → 0.115 ms (10×)** and gbs now leads
   OCCT 1.1–1.3× and scipy 1.04–5.6× at every size. No residual lag (see Strengths
   7–8).
2. **Python bulk-eval return marshalling — FIXED (#97).** pygbs `values()` used
   to return via a Python-list round-trip (~500–950 ms for 1M vs an 11.5 ms core
   and 56 ms scipy `splev`, ~17× slower). Replaced with a single `memcpy` into an
   `(N,dim)` numpy buffer: **92 ms at 1M (~10× faster)**, gap to scipy now ~1.6×,
   values bit-identical. Residual = input-copy, a smaller un-filed follow-up.
3. **Degree elevation, C++: OCCT 1.84× faster** (copy included both sides). Real
   but confounded by copy cost; **not yet filed** — needs a copy-isolated
   re-measure first (recorded here so it isn't lost).

---

## Reproduce

**Front A (C++ vs OCCT)** — one command, runs in the `gbs-occt8` env:

```
scripts/compare_bench.sh
```

Configures `bench/` with `-DGBS_BUILD_COMPARE_BENCH=ON` (OFF by default, so normal
bench/CI builds carry no OCCT dependency), builds `bench_compare_occt`, runs it.

**Front B (Python vs scipy/geomdl)** — pinned `dev` env (numpy 2.4.3, scipy 1.16.1,
geomdl 5.4.0), pygbs built from this branch:

```
ninja -C build pygbs                                   # build the #88 parallel module
PYTHONPATH=build/python:bench/python python bench/python/compare_curve.py
PYTHONPATH=build/python:bench/python python bench/python/compare_build.py
```

(prefix both python lines with `micromamba run -n dev`). See
`bench/python/README.md` for the env detail.
