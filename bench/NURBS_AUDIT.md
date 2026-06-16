# Transcription-fidelity audit — *The NURBS Book* (Piegl & Tiller) algorithms

Scope: verify that gbs faithfully transcribes the canonical Piegl & Tiller
pseudocode for every NURBS-Book algorithm the library is built on. This is a
**read-only audit** of the algorithms; it adds this cumulative report and files a
focused issue per real deviation. Tracking issue: **#43**; epic: **#44**.

Companion in-repo reports: `bench/INTERP_AUDIT.md`, `bench/APPROX_AUDIT.md`,
`bench/COVERAGE_AUDIT.md`.

## Method

For each algorithm the code is compared **line-by-line** to the canonical
pseudocode — indices, loop bounds, `alpha`/coefficients, stop conditions,
room-checks. Every finding is classed:

- **✓ faithful** — transcribes the book (possibly via a documented *equivalent*
  route that yields identical results, e.g. repeated single-knot insertion in
  place of the batch refinement algorithm, or a separable/banded solve in place
  of a dense one);
- **(b) footgun** — a deviation silently compensated by the callers (as A2.1
  `find_span` was before #42);
- **(c) bug** — a real defect.

Algorithms are prioritised by **actual use** in the library: a missing or
divergent algorithm that nothing calls is recorded as a *gap*, not a deviation.
Usage was grepped before classifying any finding.

## TL;DR

The core evaluation / derivative / knot-operation / interpolation / approximation
arc the library is built on is a **faithful** transcription of the NURBS Book.

- **Chapters 2, 3, 4** (basis, point/derivative evaluation, rational projection):
  faithful, including the freshly-implemented **A4.4 RatSurfaceDerivs** (#33 / PR
  #55).
- **Chapter 5** (knot insertion, Bézier decomposition, degree elevation, knot
  removal): faithful. Batch refinement (A5.3/A5.4) and surface decomposition
  (A5.6) are realised as the equivalent *repeated single-curve* operations.
- **Chapter 9** (global interpolation with end conditions, surface interpolation,
  least-squares approximation): faithful; the tensor-product and least-squares
  paths use the separable Kronecker / banded normal-equation forms of the same
  systems (see `INTERP_AUDIT.md` / `APPROX_AUDIT.md`).

Two real deviations were ever found across the whole sweep, both now tracked:
`find_span`/A2.1 (**#42**, fixed) and `remove_knot`/A5.8 (**#59**, fixed). One
quality deviation remains and is filed by this pass: the generic loft default
parametrization is uniform instead of the §10.3-recommended averaged
chord-length (**→ #71**, see below). No new correctness bug was found.

## Verdict table

`file:line` points at the definition. "≡" marks a faithful *equivalent* route.

### Chapter 2 — basis functions

| Algo | What | gbs symbol | Location | Verdict | Ref |
|------|------|-----------|----------|---------|-----|
| A2.1 | FindSpan | `find_span` | `basisfunctions.ixx:278` | ✓ (fixed) — `n` convention was a footgun | #42 |
| A2.2 | BasisFuns | `basis_funcs` | `basisfunctions.ixx:151` | ✓ faithful | — |
| A2.3 | DersBasisFuns | `ders_basis_funs` | `basisfunctions.ixx:181` | ✓ faithful | — |
| A2.4 | OneBasisFun | — | — | gap — **unused**; single-basis eval done by the recursive Cox-de Boor `basis_function` (Eq 2.5/3.8, verified) | — |
| A2.5 | DersOneBasisFun | — | — | gap — **unused** (as A2.4) | — |
| Eq 2.5 | Cox-de Boor value | `basis_function` | `basisfunctions.ixx:41` | ✓ faithful | — |
| Eq 3.8 | Cox-de Boor deriv | `basis_function` | `basisfunctions.ixx:86` | ✓ faithful (incl. 0/0→0) | — |

### Chapter 3 — point & derivative evaluation (non-rational)

| Algo | What | gbs symbol | Location | Verdict | Ref |
|------|------|-----------|----------|---------|-----|
| A3.1 | CurvePoint | `eval_value_decasteljau` (curve) | `basisfunctions.ixx:354` | ✓ faithful (sum of the `p+1` non-zero `N·P`) | — |
| A3.2 | CurveDerivsAlg1 | `eval_ders_decasteljau` (curve) | `basisfunctions.ixx:408` | ✓ faithful (A2.3 pass × poles) | — |
| A3.3 | CurveDerivCpts | — | — | gap — **unused**; derivatives go via A3.2, not via derivative control points | — |
| A3.4 | CurveDerivsAlg2 | — | — | gap — **unused** (A3.2 route used instead; identical result) | — |
| A3.5 | SurfacePoint | `eval_value_decasteljau` (surface) | `basisfunctions.ixx:788` | ✓ faithful (tensor product, `i_min+ri + n_polesU·(j_min+rj)`) | — |
| A3.6 | SurfaceDerivsAlg1 | `eval_ders_decasteljau` (surface) | `basisfunctions.ixx:504` | ✓ faithful | — |

### Chapter 4 — rational (NURBS) evaluation

| Algo | What | gbs symbol | Location | Verdict | Ref |
|------|------|-----------|----------|---------|-----|
| A4.1 | CurvePoint (rational) | `eval_rational_value_simple` (d=0) | `basisfunctions.ixx:1092` | ✓ faithful (homogeneous A3.1 then ÷w) | — |
| A4.2 | RatCurveDerivs | `eval_rational_ders` (curve) | `basisfunctions.ixx:464` | ✓ faithful (quotient recurrence) | — |
| A4.3 | SurfacePoint (rational) | `BSSurfaceRational::value` (d=0) | `bssurf.h:1086` | ✓ faithful (homogeneous A3.5 then ÷w) | — |
| **A4.4** | **RatSurfaceDerivs** | `eval_rational_ders` (surface) | `basisfunctions.ixx:581` | **✓ faithful** (2-D quotient rule on the homogeneous `Aders` table; line-checked term-for-term) | #33 / PR #55 |

A4.4 detail: the homogeneous derivative table `Aders[k][l]` is built by one A2.3
pass per parametric direction on the weighted (`dim+1`) poles (= A3.6 on weighted
poles), then the quotient rule is applied. The recurrence matches the book term
for term — the `j`-direction lower-order terms (`i==0`), the `i`-loop with
`wders[i][0]` plus the nested `v2` sum over `wders[i][j]`, and the final divide
by `w00`. gbs evaluates the full **rectangular** `nu×nv` table where the book
gives the **triangular** `k+l≤d` region; the rectangular generalisation is valid
and is computed in correct dependency order (each `SKL[k][l]` reads only
`k'≤k, l'≤l`). Already pinned by `tests/tests_surface_derivatives.cpp`
(`rational_all_weights_one_matches_non_rational`,
`rational_weighted_matches_finite_differences`,
`d_dmu2_on_nurbs_circle_is_curvature_vector`).

### Chapter 5 — knot & degree operations

| Algo | What | gbs symbol | Location | Verdict | Ref |
|------|------|-----------|----------|---------|-----|
| A5.1 | CurveKnotIns | `insert_knots` | `knotsfunctions.ixx:366` | ✓ faithful (indices, `alpha`, `R[]`, room-check) | — |
| A5.2 | SurfaceKnotIns | `insertKnotU` / `insertKnotV` | `bssurf.h:590` / `:632` | ✓ ≡ A5.1 applied to each row/column (same `u`, same multiplicity ⇒ rows stay synced) | — |
| A5.3 | RefineKnotVectCurve | — (loop of `insert_knots`) | `bsctools.h` `unify_knots` | ✓ ≡ repeated A5.1 (batch form is an efficiency-only variant; identical result) | — |
| A5.4 | RefineKnotVectSurface | — (per-direction insertion) | `bssurf.h` | ✓ ≡ repeated A5.1 | — |
| A5.5 | DecomposeCurve | `bezier_segments` | `knotsfunctions.ixx:832` | ✓ ≡ raise interior knots to mult `p` via A5.1, then slice into `p+1`-pole segments | — |
| A5.6 | DecomposeSurface | — (per-isocurve) | `bssurf.h` (in elevation) | ✓ ≡ A5.5 per isocurve | — |
| Eq 5.36 | Bézier degree elevation | `increase_bezier_degree` | `knotsfunctions.ixx:806` | ✓ faithful | — |
| A5.8 | RemoveCurveKnot (tol) | `remove_knot(u,p,num,U,Pw,tol)` | `knotsfunctions.ixx:554` | ✓ faithful (fixed: single carried `temp[]`, per-pass tol test, commit only on accept) | #59 / PR #61 |
| A5.8′ | RemoveCurveKnot (no-tol) | `remove_knot(u,num,U,P,p)` | `knotsfunctions.ixx:456` | ✓ faithful (A5.8 with the deviation test elided — exact removal for elevation cleanup; always removes exactly `num`, so surface elevation rows stay synced) | — |
| A5.9 | DegreeElevateCurve | `increase_degree` | `knotsfunctions.ixx:946` | ✓ ≡ Bézier-decompose (A5.5) + elevate (Eq 5.36) + recombine + `remove_knot` interior cleanup (`num = p−M[i]`) | — |
| A5.10 | DegreeElevateSurface | `increaseDegreeU` / `increaseDegreeV` | `bssurf.h:867` / `:899` | ✓ ≡ A5.9 per row/column (no-tol removal ⇒ identical knots across rows) | — |
| A5.11 | DegreeReduceCurve | — | — | gap — **not implemented, unused** (no caller anywhere) | — |

### Chapter 9 — interpolation & approximation

| Algo | What | gbs symbol | Location | Verdict | Ref |
|------|------|-----------|----------|---------|-----|
| A9.1 | GlobalCurveInterp | `interpolate(pt_begin,pt_end,cstr_lst,p)` | `bscinterp.h:399` | ✓ faithful (collocation rows = `fill_basis_row` at given derivative order) | — |
| A9.1 | end-derivative / natural ends | `interpolate(Q,u,add_computed)` + `build_3pt_tg_vec` | `bscinterp.h:277`, `:64` | ✓ faithful (Bessel interior tangents §p.386; natural — zero end 2nd-deriv — clamped at the ends) | — |
| A9.4 | GlobalSurfInterp | `build_poles(Q,ku,kv,u,v,p,q)` | `bssinterp.h:131` | ✓ faithful — separable Kronecker `Mu·P·Mv^T = Q` (= the book's tensor system) | #34 |
| A9.4′ | scattered-constraint surf interp | `build_poles(bss_constraint…)` | `bssinterp.h:246` | ✓ faithful (mixed `(du,dv)` rows = outer product of the two banded basis-derivative vectors) | — |
| Eq 9.63–9.67 | LS curve (free) | `approx` | `bscapprox.h:137` | ✓ faithful (normal equations `NᵀN P = NᵀQ`, banded SPD) | #65 (B1) |
| Eq 9.63–9.67 | LS curve (fixed ends) | `approx_bound_fixed` | `bscapprox.h:81` | ✓ faithful (endpoints moved to RHS; interior solves normal equations) | #65 (B2) |
| Eq 9.63–9.67 | LS surface (grid) | `approx` (surface) | `bssapprox.h:36` | ✓ faithful — separable normal equations `Cu·P·Cv = Muᵀ Q Mv` | #65 |
| — | weighted LS | — | — | gap — **not implemented** (unweighted throughout; documented limitation, not a bug) | #65 (B3) |

### Chapter 10 — surface constructions (context, not re-audited here)

| §/Algo | What | gbs symbol | Verdict | Ref |
|--------|------|-----------|---------|-----|
| §10.3 | Skinned (lofted) surface | `loft` family | ✓ faithful recipe (unify degree/knots → V by Eq 9.8 → per-column A9.1 → tensor assembly) — see #43 | — |
| §10.3 | loft default **parametrization** | `loft(curves_info,q_max)` | **(b) deviation** — uniform `make_range(0,1,n)` instead of averaged chord-length/centripetal (Eq 9.5/9.6); `curve_parametrization` exists but is not wired in | **→ #71** |
| §10.5.1 | Gordon surface | `gordon` (`bssbuild.h:562`) | ✓ faithful (line-checked) — Boolean sum `Lu + Lv − Tuv` (= `L1 + L2 − T`); see note below | — |
| §10.2 / §10.4 / §10.5.2 | swung / swept / Coons | — | not implemented (out of "loft" scope; recorded in #43) | — |

Gordon detail (`gordon`, `bssbuild.h:562`): the §10.5.1 Gordon surface is the
**Boolean sum** of a bidirectional curve network, `S = L1 + L2 − T`. The code is
a faithful transcription:

- `build_intersections` returns the network intersection grid `Q[i][j]` and the
  shared `(u, v)` parameters; **all three** terms are built on those same params.
- `Lu = loft_generic(u-curves, v)` ≡ **L1** (skin the u-curves in `v`);
  `Lv = loft_generic(v-curves, u)` then `invertUV()` ≡ **L2** (skin the v-curves
  in `u`). Each is elevated to the common degree.
- `Tuv = build_poles(Q, …, u, v, p_, q_)` ≡ **T**, the tensor-product
  interpolation of the intersection grid (A9.4) at the same `(u, v)`.
- `unify_knots({Lu, Lv, Tuv})` makes the three pole grids compatible **before**
  the pole arithmetic, then `poles_gordon = poles(Lu) + poles(Lv) − poles(Tuv)`
  — exactly `L1 + L2 − T`. The result carries the unified `(ku, kv)` at degree
  `(p, q)`.

The construction matches the book term-for-term; the only inherited caveat is the
loft default parametrization (**#71**) feeding `Lu`/`Lv`, which is a quality, not
correctness, deviation.

## Deviations

### Real, tracked

1. **A2.1 `find_span`** — wrong `n` convention, was patched by call-site clamps
   (footgun). Filed and fixed: **#42** (PR #49).
2. **A5.8 `remove_knot`** — under-removed for `num>1` (stale poles re-read across
   passes; `first`/`last` advanced before the room test). Filed and fixed: **#59**
   (PR #61); now a faithful A5.8 transcription, pinned by
   `tests_knotsfunctions.remove_knot_algo`.

### Real, filed by this pass

3. **§10.3 loft default parametrization** — `loft(curves_info, q_max)`
   (`gbs/loft/loftBase.h:162`) uses a **uniform** skin-direction parameter
   `v = make_range(0,1,n_curves)` (with a `// TODO: Compute an optimal
   parametrization`) instead of the §10.3-recommended **averaged chord-length /
   centripetal** parameters (Eq 9.5/9.6 averaged over the pole rows). The surface
   still interpolates every section, so this is a *quality* deviation (b), not a
   correctness bug — but it distorts the v-mapping for unevenly-spaced sections.
   `curve_parametrization` (`knotsfunctions.ixx:306`) already implements the
   chord/centripetal parameters and the spine loft
   (`bssbuild.h get_v_aligned_poles`) already averages them correctly; they are
   simply not wired into this generic overload. → focused issue (see #43 / #44).

### Gaps (unimplemented and unused — not deviations)

- A2.4 / A2.5 OneBasisFun / DersOneBasisFun (single-basis eval via recursive
  `basis_function` instead).
- A3.3 / A3.4 CurveDerivCpts / CurveDerivsAlg2 (derivatives via A3.2 instead).
- A5.11 DegreeReduceCurve (no caller).
- Weighted least-squares (Chapter 9) — unweighted only (#65, B3).
- Swung (§10.2), swept (§10.4), Coons (§10.5.2), surface of revolution as a
  book-construction — out of scope, recorded for completeness.

## Status

The full audit perimeter (Chapters 2–5 and 9 algorithms the library uses, plus
the Chapter 10 loft/Gordon constructions) is now covered with **no new
correctness bug**. The Gordon surface blend (§10.5.1) — the last pending
line-check — was verified faithful this pass. The only outstanding real deviation,
the loft default parametrization, is filed as its own focused issue so that every
real deviation found across the audit is tracked by a dedicated issue: A2.1 → #42
(fixed), A5.8 → #59 (fixed), loft-param → #71 (open).
