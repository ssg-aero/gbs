# Performance audit — B-spline interpolation

Scope: the interpolation builders in `gbs/bscinterp.h`, `gbs/bssinterp.h` and the
shared `build_poles_matrix` in `gbs/knotsfunctions.ixx`. These turn a set of
constraints (points / derivatives) into control poles by assembling a collocation
matrix `N` and solving `N · poles = Q`.

Harness: `bench/bench_interp.cpp`, header-only, clang-21 `-O3 -march=haswell`,
`double`, `dim=3`. Best-of-N total `interpolate()` wall time.

## TL;DR

Interpolation has **three compounding inefficiencies**, all in matrix assembly +
solve. None are on the curve *evaluation* path fixed in #30/#32 — this is a
separate, larger win:

1. **Dense factorization of a structured matrix.** The collocation matrix is
   *banded* (curve) / *separable Kronecker* (tensor-product surface), but it is
   stored dense and factorized with `partialPivLu` — `O(n³)` for curves,
   `O((nu·nv)³)` for surfaces. The surface code even says `//TODO solve banded
   system`.
2. **All n² entries assembled, though only `p+1` per row are non-zero.** The
   assembly loops every (row, col) pair instead of the band around the span.
3. **Each entry via the recursive `basis_function`** (`~O(p·2^p)`, the same
   exponential kernel replaced in #30) instead of one `ders_basis_funs` pass per
   row.

## Measured

### Curve, degree 3, vs #points
| n_points | time ms | × per doubling |
|----------|---------|----------------|
| 50       | 0.41    | –     |
| 100      | 0.91    | 2.2×  |
| 200      | 4.13    | 4.5×  |
| 400      | 23.2    | 5.6×  |
| 800      | 269.9   | 11.6× |

Super-cubic growth (a doubling costs 5–12×): the dense `O(n³)` LU dominates at
scale, on top of `O(n²)` assembly. 800 points already costs **0.27 s**.

### Curve, n=200, vs degree
| degree | 1 | 2 | 3 | 5 | 7 | 9 |
|--------|---|---|---|---|---|---|
| time ms| 3.17 | 3.49 | 4.15 | 7.67 | 22.7 | 83.7 |

With `n` fixed the solve is ~constant, so this curve is **pure assembly cost** —
and it explodes `~2^p` (p3→p9 is 20×). This is inefficiency #3 in isolation.

### Surface, bidegree (3,3), vs grid side N (N×N)
| N  | poles | time ms | × prev |
|----|-------|---------|--------|
| 8  | 64    | 0.35    | –      |
| 12 | 144   | 1.90    | 5.4×   |
| 16 | 256   | 7.38    | 3.9×   |
| 24 | 576   | 95.8    | 13.0×  |
| 32 | 1024  | 606.1   | 6.3×   |

N=16→32 (poles ×4) is **82×** in time ≈ `N^6.4` — i.e. the dense
`O((nu·nv)³)` solve. A 32×32 grid already costs **0.6 s**; a 64×64 grid would be
~30–40 s. This is the single worst hotspot in the library's batch paths.

## Root causes (by file)

- `knotsfunctions.ixx:848 build_poles_matrix` and `bscinterp.h:139`,
  `bscinterp.h:332-351`: nested `i,j` over the full `n×n`, `N(i,j) =
  basis_function(...)` (recursive). #2 + #3.
- `bscinterp.h:104,146,353`: `N.partialPivLu()` on a dense `n×n`. #1.
- `bssinterp.h:51 fill_poles_matrix`: full `(nu·nv)×(nu·nv)` assembly as the
  Kronecker product `Nv ⊗ Nu`, then `bssinterp.h:116 build_poles(N.partialPivLu(),
  Q)` — one dense `O((nu·nv)³)` solve. #1 + #2 + #3, and ignores separability.
- `bssinterp.h:191`: scattered-constraint surface builder, same dense full
  assembly.

Note: factorization **is** correctly reused across the `dim` spatial components
(one factorization, `dim` cheap back-substitutions) — that part is already good.
The waste is purely structural: dense, full, recursive.

## Recommended fixes (priority order)

1. **Surface: exploit separability (biggest win).** For grid interpolation the
   system is `(Nv ⊗ Nu) · P = Q`. Solve it as two sets of 1-D banded systems:
   factor `Nu` (size nu) and `Nv` (size nv) once, solve `Nu` across the nv rows,
   then `Nv` across the nu columns. Cost drops from `O((nu·nv)³)` to roughly
   `O(nu·nv·(p²+q²))` — from `~N^6` to `~N²`. Expect 32×32: 606 ms → low single
   ms (>100×), and large grids become feasible.
2. **Band assembly + one-pass basis (curve & surface).** Use `find_span` to fill
   only the `p+1` non-zero entries per row, each row from a single
   `ders_basis_funs` pass (the alloc-free A2.3 from #30, now in `master`) instead
   of `p+1` recursive `basis_function` calls. Removes #2 and #3 together; flattens
   the degree explosion (p9 curve 84 ms → a few ms).
3. **Banded solve for curves.** Store `N` banded (Eigen `BandMatrix` / a sparse
   `SparseLU`, or a hand-rolled LU for the totally-positive band) → `O(n·p²)`
   factorization instead of `O(n³)`. 800 pts/p3: 270 ms → low single ms.

Caveats:
- The band/separable structure holds for the **sequence / grid** interpolators
  (simple knots, monotone params). The **scattered-constraint** builders
  (`constrPoint`, `bsc_constraint`, `bss_constraint`) need the constraints sorted
  first to be banded — there is already a `//TODO sort Q ... to get a block
  system`. Banding those is a follow-up; the one-pass basis (#2/#3) applies to
  them regardless.
- Derivative-constrained interpolation (`nc>1`) yields a *block*-banded matrix —
  still bandable, just with a wider band.

## Reproduce
```
cmake -S bench -B bench/build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3 -march=haswell"
cmake --build bench/build --target bench_interp
./bench/build/bench_interp
```
