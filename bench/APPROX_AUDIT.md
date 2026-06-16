# Performance + correctness audit — B-spline approximation

Scope: the least-squares **approximation** builders in `gbs/bscapprox.h` (curves)
and `gbs/bssapprox.h` (surfaces), plus the shared assembly machinery in
`gbs/knotsfunctions.ixx` (`fill_basis_row` / `fill_basis_band` /
`build_poles_matrix`) and `gbs/bssinterp.h` (`fill_poles_matrix`, `build_poles`).
These turn a point set into control poles by assembling a collocation matrix `N`
and solving the least-squares system `N · P ≈ Q`.

This audit is **read-only on the algorithms** (no fix in this PR); it adds the
bench harness + this report. Fixes are sequenced as PRs in the umbrella issue.

Harness: `bench/bench_approx.cpp`, header-only, clang-21 `-O3`, `double`, `dim=3`.
Best-of-N total `approx()` wall time. Same machine/flags as `bench/INTERP_AUDIT.md`.

## TL;DR

Approximation has the **same three structural inefficiencies as interpolation
(#34)**, but here they are only *partially* fixed, and they sit on a hotter path
(`refine_approx` re-approximates in a loop):

1. **Dense factorization of a structured matrix (all paths).** Every approx solve
   is a dense `colPivHouseholderQr`. For sorted parameters the curve collocation
   is *banded* and `NᵀN` is *banded SPD*; the surface grid collocation is a
   *separable Kronecker* `Nv ⊗ Nu`. None of this structure is exploited — same
   gap #34 closed for interpolation, untouched for approximation.
2. **Surface assembly is dense + recursive.** `fill_poles_matrix`
   (`bssinterp.h:51`) fills **all** `np_u·np_v` columns per row (only
   `(p+1)(q+1)` are non-zero) and each entry uses the **recursive**
   `basis_function` (`~2^p`, `bssinterp.h:65,69`).
3. **`approx_bound_fixed` assembly is dense + recursive** (`bscapprox.h:55,64,65`)
   — the `~2^p` kernel. The plain curve `approx` was already modernized to the
   banded one-pass `build_poles_matrix` (#36), but `approx_bound_fixed` was **not**
   — and it is exactly the path the refine loop and almost every high-level
   overload use (`fix_bound = true`).

Plus a real **logic bug in `refine_approx`** (the average-deviation stop criterion
is computed on the wrong subset) and `refine_approx` **re-assembles + re-factorizes
from scratch every iteration** (up to `n_max = 200` dense QRs).

## Measured

### Curve — banded `approx` vs recursive `approx_bound_fixed`, vs degree
`n_points = 400`, `n_poles = 50` fixed (assembly-dominated; the solve is constant):

| degree | banded ms | recursive ms | rec / band |
|--------|-----------|--------------|------------|
| 1 | 2.78 | 2.84 | 1.02× |
| 2 | 2.78 | 3.02 | 1.08× |
| 3 | 2.80 | 3.46 | 1.23× |
| 5 | 2.80 | 6.07 | 2.17× |
| 7 | 2.89 | 16.03 | 5.54× |
| 9 | 2.90 | **54.43** | **18.74×** |

The banded path (`build_poles_matrix`, `ders`) is **flat** in degree; the recursive
`approx_bound_fixed` explodes `~2^p` (p3→p9 = 16×). This is inefficiency #3 in
isolation — and `approx_bound_fixed` is the *default* high-level path.

### Curve — degree 3, `n_poles = n/4`, vs #points
| n_points | n_poles | banded ms | recursive ms |
|----------|---------|-----------|--------------|
| 50  | 12  | 0.07 | 0.10 |
| 100 | 25  | 0.45 | 0.49 |
| 200 | 50  | 1.82 | 2.16 |
| 400 | 100 | 10.9 | 11.0 |
| 800 | 200 | **76.1** | **81.2** |

At degree 3 the two paths converge because the **dense QR solve dominates** (the
`2^p` term is small): 50→800 points is ~1000× in time for 16× the points — the
dense `colPivHouseholderQr` on the structured `(n × n_poles)` system. With the
banded normal equations (`NᵀN`, bandwidth `p+1`, SPD → band Cholesky) this drops to
`O(n·p²)` — 800 pts should fall from ~76 ms to sub-millisecond (cf. #34 plan #3).

### Surface — grid approx, bidegree (3,3), `n_poles = (N/2)²`, vs grid side N
| N | pts | poles | time ms | × prev |
|---|-----|-------|---------|--------|
| 8  | 64   | 16  | 0.08 | – |
| 12 | 144  | 36  | 0.86 | 10.3× |
| 16 | 256  | 64  | 4.05 | 4.7× |
| 24 | 576  | 144 | 33.4 | 8.3× |
| 32 | 1024 | 256 | **166.8** | 5.0× |

N=16→32 (poles ×4, points ×4) is **41×** in time — the dense `colPivHouseholderQr`
on a `(N² × np²)` matrix, filled densely (all `np²` columns) with the recursive
kernel. Interpolation got `~N^6 → ~N²` from separability (#35); the *same*
separable structure exists for grid **least squares** (`NᵀN = (MuᵀMu) ⊗ (MvᵀMv)`),
so the analogous win is available here and unclaimed.

## Findings

### A. Performance

- **A1 — `approx_bound_fixed` uses the recursive `basis_function`.**
  `bscapprox.h:55` (interior band), `:64`/`:65` (the two boundary basis columns)
  call the recursive `basis_function`, and the loop spans the **full**
  `(n_params-2) × (n_poles-2)` matrix instead of the `p+1`-wide band. This is the
  `~2^p` kernel removed elsewhere in #30/#36. Fix: assemble with `fill_basis_row`
  (which already clamps to the band and uses `basis_ders` A2.3). Measured 18.7× at
  p9 above. **Hot path:** `refine_approx` and most `approx(...)` overloads call it
  via `fix_bound = true`.

- **A2 — Surface assembly is dense + recursive.** `fill_poles_matrix`
  (`bssinterp.h:51`) loops every `(i, j)` pole pair per parameter (`bssinterp.h:63-70`)
  and each entry is `basis_function(...)·basis_function(...)` (recursive,
  `bssinterp.h:65,69`). Only `(p+1)(q+1)` of the `np_u·np_v` columns per row are
  non-zero. Fix: fill the two 1-D bands with `fill_basis_row` per direction and
  take the tensor product on the support only (mirror `build_poles`,
  `bssinterp.h:121-129`). The scattered `bss_constraint` overload
  (`bssapprox.h:42-43`) has the same recursive/dense fill.

- **A3 — All approx solves are dense `colPivHouseholderQr` on a structured
  matrix.** Curve: `bscapprox.h:70` (bound-fixed), `:114` (plain). Surface:
  `bssapprox.h:23`, `:48`. For sorted parameters the curve system is banded and
  `NᵀN` is banded SPD (→ band `LLT`/`LDLT` or `SparseLU NaturalOrdering`, cf.
  #37/#39); the surface grid system is separable Kronecker (→ factor the two small
  1-D normal-equation matrices once, two sequences of 1-D solves, cf. #35). The
  factorization *is* correctly reused across the `dim` components (`build_poles`,
  `bssinterp.h:83`; the per-`d` loops in `bscapprox.h`), so only the structure is
  missed, not the dim reuse.

- **A4 — `refine_approx` re-approximates from scratch each iteration.**
  `bscapprox.h:214`/`:218` call `approx_bound_fixed`/`approx`, which re-assemble
  **and** re-factorize the entire matrix, inside a loop of up to `n_max = 200`
  iterations (`bscapprox.h:173`) that inserts **one** knot per pass. So an N-knot
  refinement is O(N) full dense factorizations. Once A1–A3 land each pass is cheap;
  a further win is to keep the factorization/assembly and update it for the single
  inserted knot rather than rebuild.

- **A5 — `build_poles` solves one RHS per component** (`bssinterp.h:90`,
  `bscapprox.h:83/127`). Could batch the `dim` right-hand sides into a single
  `solve(B)` (matrix RHS). Minor next to A1–A3.

### B. NURBS-Book fidelity (feeds #43)

- **B1 — Method differs from the book, result is the same (faithful).** The book’s
  global least-squares (A9.6 / Eq 9.63–9.67) forms the **normal equations**
  `(NᵀN) P = Nᵀ Q` and Cholesky-solves. gbs instead QR-solves the **rectangular**
  system `N P ≈ Q` directly (`colPivHouseholderQr`). Both yield the same
  least-squares minimizer; QR is *better* conditioned. Not a deviation in result —
  worth a one-line note in #43. (If A3 moves curves to band `NᵀN` Cholesky, gbs
  would then match the book’s method too.)

- **B2 — `approx_bound_fixed` matches the book’s endpoint-interpolation LS.** It
  fixes `P₀ = Q₀`, `Pₙ = Q_m` (`bscapprox.h:90,91`) and builds the residual
  `R_k = Q_k − N_{0,p}(ū_k)Q₀ − N_{n,p}(ū_k)Q_m` (`bscapprox.h:80`), i.e. book
  Eq 9.63, then LS-solves for the interior poles. Structurally faithful. ✓

- **B3 — No weighted least squares.** The book’s LS admits per-point weights
  `w_k`; gbs is unweighted everywhere. Acceptable, but record the limitation
  (relevant if weighted fitting is ever needed).

### C. Correctness / edge cases

- **C1 — `refine_approx` average-deviation criterion is computed on the wrong
  subset (bug).** `d_avg_ += d` (`bscapprox.h:200`) sits **inside**
  `if (d > d_max_)` (`bscapprox.h:184`). So `d_avg_` accumulates only the
  *record-breaking* distances seen so far in the scan, not all points, then is
  divided by `pts.size()` (`bscapprox.h:205`). The reported "average deviation" is
  therefore a monotone-filtered partial sum — the `d_avg_ < d_avg` stop test
  (`bscapprox.h:222`) is meaningless. Likely the `+= d` belongs outside the
  `if`. **Real logic bug.**

- **C2 — `refine_approx` ignores errors near existing knots.** `d_max_`/`u_max`
  are only updated when the worst point sits "well inside" a span
  (`dul > 0.33 && duh > 0.33`, `bscapprox.h:195`). A large error close to a knot
  never sets `u_max`, so if all large errors are near knots the loop sees
  `u_max < 0` and breaks (`bscapprox.h:206`) reporting convergence with the
  tolerance unmet. The `0.33` is also an unexplained magic constant.

- **C3 — Inconsistent input guards.** The `Curve` overload at `bscapprox.h:310`
  throws if `n_poles < p + 1`; the point-set `approx`/`approx_bound_fixed` entry
  points do **not** — they rely on the downstream `BSCurve` ctor / `check_curve`.
  `n_poles > n_pts` (under-determined LS) is checked for **surfaces**
  (`bssapprox.h:17`) but **not for curves**. Unify the guards (reject
  `n_poles < p+1` and `n_poles > n_pts` consistently, with a clear message).

- **C4 — Degenerate `approx_bound_fixed` with only endpoints.** With
  `n_params == 2` the matrix is `0 × (n_poles-2)`; the path is not guarded.

- **C5 — Default `n_poles = p * 2`** (`bscapprox.h:233,408,433,454`) is arbitrary
  and unrelated to point count or curvature; combined with C3 it can silently
  produce a poor or ill-posed fit.

### D. C++ quality / dead code

- **D1 — Dead code: `removeRow` / `removeColumn`** (`bscapprox.h:10-31`) are
  defined but used nowhere in the repo (`grep` confirms). Each does an O(n²)
  `conservativeResize` copy. Remove.
- **D2 — Dead locals / commented parallelism in `refine_approx`:** `j_ =
  make_range(...)` (`bscapprox.h:172`) is unused; the commented
  `std::for_each(std::execution::par, ...)` and `std::cout` debug lines
  (`bscapprox.h:178-179,203,220,226`) should go.
- **D3 — `std::vector<T> knots{crv_refined.knotsFlats()}` copied every iteration**
  (`bscapprox.h:176`) inside the refine loop.
- **D4 — Strided RHS fill** with the noted `// Pas top au niveau de la
  localisation mémoire` (`bscapprox.h:80,124`; `bssinterp.h:87`) — `b(i) =
  pts[i][d]` walks the point array with stride `dim`. Minor; batched `solve(B)`
  (A5) sidesteps it.
- **D5 — `int` / `size_t` index mixing** in the curve assembly/solve loops
  (`bscapprox.h:51-89,120-132`) — signed/unsigned comparisons.

### E. Duplication / magic numbers (feeds #38)

- **E1 — ~10 near-identical curve `approx` overloads** (`bscapprox.h:231-457`):
  each computes parameters (`deviation_based_params` / `uniform_distrib_params` /
  `curve_parametrization`) + points (`make_points` / `discretize`), then calls
  `approx` + `refine_approx` with slightly different `n_poles` / `tol` defaults.
  Collapsible to a few constrained templates.
- **E2 — `approx_bound_fixed` duplicates `build_poles_matrix`’s assembly** (but
  recursive); the surface `approx` duplicates `fill_poles_matrix` from
  `bssinterp.h`. A2/A1 fixes naturally dedup these onto `fill_basis_row`/
  `fill_basis_band`.
- **E3 — Magic numbers** to recentralize in #38: `d_max = 1e-3`, `d_avg = 1e-4`,
  `n_max = 200` (`bscapprox.h:166`), `0.33` (`:195`), `n_poles = p*2`
  (`:233,408,433,454`), `np = 30` (`:289,359,404,421`), `n_poles*10` (`:314`),
  `np/3` (`:363`), `10.*tol` (`:366,387,411,436`), `n_max_pts = 5000` (`:334`),
  `pts.size() - 5` (`:223`).

## Suggested PR sequencing (correctness → perf → dedup)

- **PR 1 — correctness.** Fix C1 (move `d_avg_ += d` out of the `if`) and C2
  (separate "where to insert a knot" from "what is the max/avg error"; surface the
  `0.33` as a named constant). Unify the C3 guards across curve entry points. Add
  regression tests: a refine that must keep iterating until `d_avg` is met; an
  over-/under-specified `n_poles` rejected consistently. Low risk, no perf change.
- **PR 2 — assembly (A1, A2).** Route `approx_bound_fixed` and the surface
  `fill_poles_matrix` / `bss_constraint` fill through `fill_basis_row` /
  `fill_basis_band` (banded one-pass, A2.3). Acceptance: poles match the current
  dense/recursive solver to ~1e-9; curve assembly flat in degree; bench shows the
  p9 18.7× and the surface dense-fill factor gone.
- **PR 3 — structured solve (A3, A4).** Curve: banded `NᵀN` SPD → band `LLT`/`LDLT`
  (or `SparseLU NaturalOrdering`) — `O(n·p²)`. Surface grid: separable least
  squares `(MuᵀMu) ⊗ (MvᵀMv)` — `~N^? → ~N²`. Reuse the factorization across
  `refine_approx` iterations (A4). Acceptance: poles ~1e-9 vs dense; curve ~linear
  in #points; surface ~`N²`; bench before/after.
- **PR 4 — dedup + cleanup (D, E).** Remove dead code (D1, D2), collapse the
  overload set (E1), consolidate surface assembly with `bssinterp` (E2), and feed
  the constants list (E3) to #38.

## References
- #34 — the analogous interpolation perf work (band assembly, separable surface
  solve, banded curve solve). Same three inefficiencies, already fixed there.
- #43 — NURBS-Book transcription audit; B1/B2/B3 feed the A9.x line items.
- #44 — epic / orchestration hub.
