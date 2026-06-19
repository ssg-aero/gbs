#pragma once

#include <cstddef>
#include <limits>

// =============================================================================
// gbs tuning constants — single source of truth
// -----------------------------------------------------------------------------
// Hard-coded numeric constants that tune the core algorithms (tolerances,
// solver thresholds and default algorithm parameters) live here instead of
// being re-typed at every call site. Each one is named and documented so the
// intent is explicit and the value can be reviewed / tuned in one place.
//
// Tolerances and floating defaults are templated on the scalar type so they
// adapt to `float` / `double`; integer thresholds are plain `inline constexpr`.
//
// NOTE: algorithm-intrinsic constants (e.g. the Shewchuk robust-predicate error
// bounds in inc/topology/robustPredicates.h, or the Gauss-Kronrod point count)
// are NOT tuning knobs and deliberately stay with their algorithm.
// =============================================================================

namespace gbs
{
    // ---- Knot / parameter comparison tolerance -----------------------------
    // Two knots (or curve parameters) closer than this are treated as equal.
    // Scaled off the machine epsilon so it tracks the working precision.
    template <typename T>
    inline constexpr T knot_eps = std::numeric_limits<T>::epsilon() * 100;

    // ---- Allocation-free evaluator capacity --------------------------------
    // Largest B-spline degree the allocation-free evaluators handle on the
    // stack. Real CAD/CAM degrees are well under this; beyond it the evaluators
    // fall back to the (correct but slow) recursive basis_function so behaviour
    // never silently breaks. Stack cost is O((max+1)^2) scalars.
    inline constexpr std::size_t bspline_stack_max_degree = 24;

    // ---- Dense vs sparse linear-solve cutoffs ------------------------------
    // Below this number of unknowns a dense LU/LDLT beats the sparse path: the
    // sparse setup (symbolic analysis / fill-reducing ordering) has a fixed
    // overhead that dominates for small systems, while above it the banded
    // collocation sparsity (p+1 non-zeros per row) makes the sparse solve much
    // cheaper than the dense O(n^3) factorization.
    inline constexpr std::ptrdiff_t interp_sparse_threshold = 100; // interpolation (LU)
    inline constexpr std::ptrdiff_t approx_sparse_threshold = 100; // approximation (Cholesky)

    // The interpolation band LU (gbs/bscinterp.h, #96) runs WITHOUT pivoting:
    // valid because the point-interpolation collocation matrix is totally positive
    // (Schoenberg-Whitney), so no-pivot Gaussian elimination is stable (de Boor)
    // with growth factor ~1 and bounded multipliers (measured <= ~3 even on a
    // strongly clustered parametrization). Derivative-constrained (Hermite) and
    // pathological systems are NOT totally positive and can produce a large
    // multiplier (a small pivot); the factorization then bails above this bound
    // and the caller falls back to the pivoted sparse/dense solve. The bound is
    // far above any value-interpolation multiplier, so the fast path is never lost
    // for the common case, yet any genuine instability is caught.
    template <typename T>
    inline constexpr T interp_band_max_growth = T(1e4);

    // ---- Parallel-execution size gate --------------------------------------
    // Bulk evaluators (values / d_dms / d_dmus ...) route std::transform through
    // gbs::transform_threshold (gbs/execution.h), which dispatches to the
    // parallel execution policy only when the input has at least this many
    // elements and runs serially below it. The parallel STL (libstdc++/MSVC)
    // does not cheaply fall back to serial for small ranges: it still spins up
    // the TBB task machinery (~9-15 us fixed cost, measured), so an unconditional
    // `par` regresses small calls by up to 25x. 1000 is the conservative,
    // *measured* break-even for the lightest kernel (curve value); heavier
    // kernels (offset, d_dm2) cross over near ~64-128, so a known-heavy workload
    // may lower this. See bench/PARALLEL_AUDIT.md (#84/#92). Mutable so it can be
    // tuned at runtime per workload.
    inline std::size_t parallel_min_size = 1000;

    // ---- Curve/surface intersection & extrema tolerance --------------------
    // Default tolerance for extrema / intersection / projection searches
    // (extrema_curve_curve, build_intersections, gordon, closest_curve, ...).
    template <typename T>
    inline constexpr T extrema_tol = T(1e-6);

    // ---- Delaunay / Boyer-Watson mesh tolerance ----------------------------
    // Default tolerance for the in-circle / orientation predicates and the
    // mesh-editing geometry tests (Boyer-Watson, edge recovery, point-in-face).
    template <typename T>
    inline constexpr T delaunay_tol = T(1e-10);

    // ---- Least-squares approximation defaults ------------------------------
    // Max / mean deviation targets of the iterative refinement (refine_approx,
    // approx).
    template <typename T>
    inline constexpr T approx_dev_max = T(1e-3);
    template <typename T>
    inline constexpr T approx_dev_avg = T(1e-4);
    // Cap on the refinement iterations (knot-insertion passes).
    inline constexpr std::size_t approx_max_iter = 200;
    // A refinement knot is only inserted when both surrounding spans keep at
    // least this fraction of their length, so knots never pile up at a boundary.
    template <typename T>
    inline constexpr T approx_min_span_fraction = T(0.33);
    // Default number of points for a preliminary curve discretization.
    inline constexpr std::size_t approx_sample_count = 30;
    // Hard cap on the points produced by deviation-based discretization.
    inline constexpr std::size_t approx_max_sample_points = 5000;
    // default_n_poles heuristic: poles ~= this factor * degree.
    inline constexpr std::size_t default_poles_per_degree = 2;

    // ---- Curve -> B-spline conversion / loft defaults ----------------------
    // Defaults shared by to_bs_curve() and loft() when approximating an
    // arbitrary Curve by a BSCurve.
    template <typename T>
    inline constexpr T bs_conversion_dev = T(0.01); // max deviation target
    inline constexpr std::size_t bs_conversion_np = 100;     // discretization points
    inline constexpr std::size_t bs_conversion_degree = 5;   // approximation degree
    inline constexpr std::size_t loft_v_degree_max = 3;      // default max v-degree of a loft
} // namespace gbs
