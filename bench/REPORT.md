# Performance study — arc-length derivative evaluators

Scope: the kernels added on `feat/vectorized-curvilinear-derivatives`
(`d_dm`, `d_dm2`, `d_dmu/v`, `d_dmu2/v2` and their vectorized `*s` overloads).

Harness: `bench/bench_arc_length_derivatives.cpp` — header-only, built with
clang-21 `-O3 -march=haswell`, `double`, `dim=3`. Best-of-7 timings, DCE-guarded.
Curve: `BSCurve<double,3>`, degree 3, 60 poles. Surface: `BSSurface<double,3>`,
bidegree (3,3), 30×30 poles.

## Findings

### 1. The vectorized wrappers add no measurable overhead
`d_dms` / `d_dmus` / … are thin `std::transform` loops over the scalar overload.
Measured difference vs a hand-written scalar loop is within run-to-run noise
(±5–11%, no consistent sign) at every N from 1e3 to 1e6. **The PR's design is
sound**: throughput is dictated entirely by the underlying scalar `value()`
kernel, not by the wrapper. There is therefore no SIMD/batched win to chase in
the wrappers themselves — the lever is the per-point kernel below.

### 2. Curve derivative evaluation is O(n_poles), not O(degree) — the dominant cost
`BSCurve::value(u, d)` → `eval_value_decasteljau` (basisfunctions.ixx:209).
Span reduction (restricting the sum to the `p+1` non-zero basis functions) is
gated on `d == 0` (line 215, `// TODO fix for d != 0`). Since `d_dm` uses `d=1`
and `d_dm2` uses `d=1` **and** `d=2`, **the derivative path always sums over ALL
poles** — 60 here — calling the recursive `basis_function` 60× when only 4 are
non-zero. Cost scales with pole count, which is why the *curve* (~2.2 µs) was
slower than the bidegree-(3,3) *surface* (~0.79 µs): the surface kernel
(basisfunctions.ixx:394) span-reduces unconditionally.

The non-zero support of a B-spline basis function and of all its derivatives is
identical, so derivative span reduction is valid in principle. **But removing
the guard alone is NOT safe**: `find_span` was implemented with
`std::lower_bound`, which at a knot of multiplicity > 1 returns the *left* span
— consistent with `basis_function`'s half-open support for values (continuity
hides it) but wrong for one-sided derivatives. The correct fix is to make
`find_span` use the Piegl & Tiller A2.1 convention (`k[s] <= u < k[s+1]`), then
remove the guard. (Verified: removing the guard alone broke
`tests_bsctools.extend_to_point`, which evaluates a derivative exactly at a
`C^0` extension joint; the `find_span` fix makes it pass.) Measured effect of
the complete fix:

| kernel (N=1e6)        | before    | after    | speedup |
|-----------------------|-----------|----------|---------|
| curve `d_dm`          | 2196 ns   | 182 ns   | **12.0×** |
| curve `d_dm2`         | 5063 ns   | 332 ns   | **15.2×** |

Degree sweep (curve `d_dm`, N=2e5, 60 poles) — the win is largest at low degree
(more poles eliminated relative to `p+1`):

| degree | before     | after    | speedup |
|--------|------------|----------|---------|
| 2      | 1183 ns    | 68 ns    | 17.5×   |
| 3      | 2207 ns    | 163 ns   | 13.5×   |
| 5      | 8433 ns    | 871 ns   | 9.7×    |
| 7      | 35161 ns   | 4668 ns  | 7.5×    |

Surface kernels are unchanged by this (already span-reduced) — they act as the
control: `d_dmu` stays ~0.79 µs before and after.

### 3. The basis function is recursive (exponential in degree)
`basis_function` (basisfunctions.ixx:86) is a textbook recursive Cox–de Boor:
each call spawns two recursive calls (lines 106/110) with no memoization, i.e.
~2^p work per evaluation. Even after fix #1, the degree sweep still grows
super-linearly (deg-7 `d_dm` ≈ 29× deg-3). The library already ships an O(p²)
one-pass routine, `basis_funcs` (line 117), but it (a) is not used by the curve
eval and (b) does not compute derivatives. Replacing the per-pole recursive
calls with a `ders_basis_funs`-style routine (all `p+1` functions and all
derivative orders up to `d` in a single O(p·(p+d)) pass) would:
  * collapse the high-degree blow-up to polynomial, and
  * let `d_dm2` get `d=1` and `d=2` from **one** pass instead of two `value()`
    calls (it currently pays ~2.0× `d_dm`, confirmed by the measured ratio).

This is a larger change than #1 and touches the shared curve eval path.

Measured, curve VALUE (d=0), 60 poles, N=2e5 — current per-pole recursive
`eval_value_decasteljau` vs the library's own one-pass `eval_value_deboor_cox`
(`basis_funcs`, O(p^2)):

| degree | recursive | deboor-cox | speedup |
|--------|-----------|------------|---------|
| 2      | 44 ns     | 95 ns      | 0.46x (recursive wins) |
| 3      | 125 ns    | 99 ns      | 1.27x   |
| 5      | 816 ns    | 125 ns     | 6.5x    |
| 7      | 4750 ns   | 154 ns     | 30.8x   |

The recursive path is ~2^p (44 -> 125 -> 816 -> 4750); deboor-cox is ~flat
(O(p^2)). Two takeaways: (a) the better algorithm already exists in the repo but
is not wired into `BSCurve::value`; (b) `eval_value_deboor_cox` allocates a
`std::vector<T> N(p+1)` per call, which is why it LOSES at p=2 — replacing that
heap alloc with a stack buffer (max-degree `std::array`) would make it win at
every degree. The derivative counterpart is `DersBasisFuns` (Piegl A2.3): all
orders 0..d in one O(p^2) pass, which also removes the double `value()` call in
`d_dm2`.

## Recommendations (in priority order)
1. **Fix `find_span` (half-open Piegl convention) + drop the `&& d == 0` guard**
   → 7–17× on every curve derivative (not just the new kernels), and also fixes
   a latent correctness bug for surface derivatives on multiple interior knots.
   Shipped on `perf/curve-derivative-span-reduction` (PR #29).
2. Longer term, replace recursive `basis_function` with a one-pass
   `ders_basis_funs` and have `d_dm2` reuse a single basis-function evaluation.
   This collapses the remaining super-linear degree scaling (deg-7 still ~29×
   deg-3 after fix #1, because `basis_function` is recursive/exponential).
   **Done — see issue #30 below.**
3. The vectorized `*s` overloads need no change.

## Update — issue #30: allocation-free O(p²) A2.3 evaluator (this branch)

Both `eval_value_decasteljau` overloads (curve and surface) were rewritten to
use allocation-free Piegl & Tiller A2.2 (values) / A2.3 (derivatives) on stack
buffers (`ders_basis_funs`/`basis_ders`, `gbs/basisfunctions.ixx`). The per-pole
recursive `basis_function` is gone from the hot path (kept only as a fallback for
`p > bspline_stack_max_degree = 24`). `use_span_reduction` is now a no-op: only
the span's `p+1` basis functions are ever touched, so the sum is always reduced.

### Curve VALUE (d=0), 60 poles, N=2e5, clang-21 -O3 -march=haswell
True before/after (recursive copy kept in `bench_eval_compare.cpp`):

| degree | recursive (old) | A2.3 (new) | speedup | deboor-cox (allocates) |
|--------|-----------------|------------|---------|------------------------|
| 2      | 49 ns           | 9.5 ns     | 5.1×    | 88 ns                  |
| 3      | 143 ns          | 14.8 ns    | 9.6×    | 99 ns                  |
| 5      | 848 ns          | 25.4 ns    | 33.4×   | 121 ns                 |
| 7      | 4971 ns         | 38.1 ns    | 130.6×  | 173 ns                 |

The new path also beats the repo's own `eval_value_deboor_cox` at every degree —
that routine is O(p²) but heap-allocates `std::vector<T> N(p+1)` per call, which
the stack buffer removes. The exponential blow-up is gone: 9.5→38 ns across
p=2..7 (was 49→4971 ns).

### Curve arc-length derivatives (degree sweep, N=2e5, 60 poles)

| degree | d_dm (new) | d_dm (after #29) | d_dm2 (new) | ratio d_dm2/d_dm |
|--------|------------|------------------|-------------|------------------|
| 2      | 44 ns      | 68 ns            | 73 ns       | 1.64             |
| 3      | 50 ns      | 163 ns           | 95 ns       | 1.88             |
| 5      | 65 ns      | 871 ns           | 147 ns      | 2.25             |
| 7      | 90 ns      | 4668 ns          | 211 ns      | 2.33             |

`d_dm` is now ~linear in degree (44→90 ns) instead of exponential; **52× at p=7**
vs the span-reduced-but-recursive #29 path. `d_dm2` still pays ~2× `d_dm` because
it calls `value(u,1)` and `value(u,2)` as two passes — a future one-pass API
(return derivative orders 0..n together) would fold that to ~1×, but it needs a
new method on the curve/surface (the current `value(u,d)` returns a single
order, and the rational path projects through it). Left as a follow-up.

### Surface (bidegree (3,3), 30×30 poles, N=1e6)
`d_dmu` ≈ 85 ns, `d_dmu2` ≈ 159 ns scalar. The surface eval now computes the
`du`-th and `dv`-th derivative basis rows once each (A2.3) and forms the
`(p+1)×(q+1)` outer product, replacing `(p+1)·(q+1)` recursive `basis_function`
calls.

### Correctness
`tests_bscurve.eval_matches_recursive_basis_ground_truth` checks the new A2.3
evaluator against the recursive Cox-de Boor definition (summed over all poles) to
1e-10, for p∈{1,2,3,5,7}, every derivative order d∈[0, p+1] (including d>p ⇒ 0),
sampled on a grid plus every distinct knot — with an interior knot of
multiplicity p (a C⁰ joint). The surface counterpart
`tests_surface_derivatives.deriv_span_reduction_matches_full_span` does the same
against a full-span recursive reference.

## Update — one-pass derivatives (follow-up branch `perf/one-pass-derivatives`)

The `d_dm2` follow-up flagged above is now done. A new one-pass multi-order API
`derivatives(u, n, CK)` (curve) / `derivatives(u, v, nu, nv, SKL)` (surface) is
added as a `virtual` on the `Curve`/`Surface` bases (default: loop `value()`),
overridden on the B-spline classes with single-pass evaluators:
`eval_ders_decasteljau` (curve & surface, one A2.3 basis pass shared across all
orders) and `eval_rational_ders` (rational curve, Piegl & Tiller **A4.2**). The
arc-length helpers `d_dm2`, `d_dmu2`, `d_dmv2` now take **one** basis pass
instead of two `value()` calls.

Curve degree sweep (N=2e5, 60 poles) — `d_dm2/d_dm` ratio collapses from ~2.0:

| degree | d_dm2 before | d_dm2 now | ratio d_dm2/d_dm (before → now) |
|--------|--------------|-----------|----------------------------------|
| 2 | 73 ns  | 44 ns  | 1.64 → **1.04** |
| 3 | 95 ns  | 59 ns  | 1.88 → **1.21** |
| 5 | 147 ns | 89 ns  | 2.25 → **1.33** |
| 7 | 211 ns | 126 ns | 2.33 → **1.54** |

Surface `d_dmu2` (bidegree (3,3), 30×30): ~150 ns → **~110 ns** (≈1.35×). The
residual >1.0 ratio is the unavoidable extra work of the higher derivative order
itself (the A2.3 `a[]`-table recurrence), not a second pass.

Rational **surface** derivatives remain unimplemented (their `value(du>0)` throws
today); the surface override keeps that behaviour via the base default. Rational
**curve** derivatives are now one-pass (A4.2).

Correctness: `tests_curve_derivatives.derivatives_match_per_order_value_{nonrational,rational}`
and `tests_surface_derivatives.derivatives_match_per_order_value` check the
one-pass output against per-order `value()` to 1e-9, over a degree sweep and at
multiplicity-p/q C⁰ joints. All 219 tests pass.

## Reproduce
```
cmake -S bench -B bench/build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=<env>/bin/clang++-21 -DCMAKE_CXX_FLAGS="-march=haswell"
cmake --build bench/build
./bench/build/bench_arc_length_derivatives
```
