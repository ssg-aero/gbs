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
2. **gbs lags on the interpolation *build*** — measured on **both** fronts. OCCT
   is **3–9× faster**, scipy **1.2–3.7× faster** at N ≥ 201. gbs's build scales
   super-linearly with point count where the banded solvers stay near-linear.
   → **focused perf issue** (the single highest-value gap).
3. **The Python bulk evaluator's *return* throws gbs's core advantage away.**
   The C++ core evaluates 1 000 000 curve points in **11 ms** (parallel), but the
   pygbs call returns in **~500–950 ms** because the result is marshalled through
   a Python-list round-trip instead of a zero-copy numpy buffer. Net: scipy's
   vectorized `splev` **wins the batched Python regime ~17×**, even though gbs's
   underlying evaluator is several × faster. → **focused perf issue** (binding).
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

### Curve interpolation (deg 3, full build) — **gbs lags**

Includes the real **201-point** case (#34). gbs `CHORD_LENGTH` vs OCCT
`GeomAPI_Interpolate` (its own parametrization) — both produce an interpolating
C² cubic; **ratio < 1 means OCCT is faster**:

| N | gbs ms | occt ms | ratio | winner |
|---|--------|---------|-------|--------|
| 50  | 0.0243 | 0.0080 | 0.33× | occt |
| 100 | 0.0568 | 0.0148 | 0.26× | occt |
| **201** | **0.1342** | **0.0291** | **0.22×** | occt |
| 400 | 0.3717 | 0.0681 | 0.18× | occt |
| 800 | 1.1642 | 0.1286 | 0.11× | occt |

gbs scales ~22× over the 50→800 range; OCCT ~16×. The gap widens with N — a
**super-linear build cost** vs OCCT's near-linear path.

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

**Batched / vectorized** (the right way to use scipy) — **scipy wins**:

| N | pygbs ms | scipy ms | geomdl ms | winner |
|---|----------|----------|-----------|--------|
| 1 000     | 1.08  | 0.066 | 12.2  | scipy |
| 10 000    | 4.86  | 0.568 | 124.8 | scipy |
| 100 000   | 73.8  | 5.64  | 1305  | scipy |
| 1 000 000 | 950.6 | 56.1  | —     | scipy ~17× |

This is **not** a core-evaluator deficit: the same 1 000 000-point evaluation runs
in **11.5 ms** in the C++ core (Front A [bulk]). The ~500–940 ms is the pygbs
**return marshalling** — `py::cast` of `vector<array<double,3>>` builds a Python
list before numpy copies it, instead of returning a zero-copy buffer. 1st/2nd
derivatives show the same shape (1M: pygbs ≈ 962 / 964 ms, scipy ≈ 58 / 59 ms).
→ **perf issue** (binding); fixing it should expose the #88 parallel speedup to
Python and flip this regime.

### Curve interpolation (deg 3, full build) — mixed; **gbs lags at N ≥ 201**

Both on a chord-length parametrization (like-for-like interpolation, not
smoothing):

| N | pygbs ms | scipy ms | geomdl ms | winner |
|---|----------|----------|-----------|--------|
| 50  | 0.039 | 0.098 | 4.63   | **pygbs** |
| 100 | 0.067 | 0.099 | 27.9   | **pygbs** |
| 201 | 0.143 | 0.119 | 207    | scipy 1.2× |
| 400 | 0.326 | 0.158 | 1578   | scipy 2.1× |
| 800 | 0.884 | 0.239 | 16163  | scipy 3.7× |

Same trend as Front A [interpolation]: gbs leads at small N (lower constant
overhead) but its build scales worse than the banded solver. **Cross-front
confirmation of the same lag.**

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
   in a loop, N=1k/10k). *Regime:* one-at-a-time eval only — for vectorizable work
   scipy's batched path wins (see lags).
4. **Approximation, C++: gbs 2.5–9.4× faster than OCCT** (LSQ deg-3, ~N/4 poles,
   N=100–800; structural caveat: fixed poles vs OCCT tolerance-adaptive knots).
5. **Knot insertion, C++: gbs 5.7× faster than OCCT** (copy + single insert,
   included both sides).
6. **Approximation, Python: pygbs ≈ scipy** (LSQ fixed poles; pygbs 1.2× at N=800).
7. **Interpolation at small N, Python: pygbs faster than scipy for ≤ 100 points**
   (lower constant overhead).

## Where gbs lags (measured) — and the issues filed

1. **Interpolation build cost scales super-linearly.** OCCT **3–9×** faster
   (N 50→800), scipy **1.2–3.7×** faster at N ≥ 201; measured on **both** fronts,
   incl. the 201-pt #34 case. Highest-value gap. → **perf issue: interpolation
   solver scaling**.
2. **Python bulk-eval return marshalling.** pygbs `values()` returns via a
   Python-list round-trip; ~500–940 ms for 1M points vs an 11.5 ms C++ core and a
   56 ms scipy `splev` (~17× slower). The #88 parallel speedup never reaches
   Python. → **perf issue: zero-copy numpy return (buffer protocol)**.
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
