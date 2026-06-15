#pragma once
#include <vector>
#include <array>
#include <gbs/bscurve.h>

namespace gbs{

    /**
     * Shared v-direction pole solver — the single core both the plain loft and the
     * spine loft route their pole computation through (issue #40 PR3 / #44).
     *
     * `B` holds the already-assembled right-hand sides: one column per
     * (u-pole, spatial component) pair (`B.cols() == n_poles_u * dim`) and one row
     * per v-direction degree of freedom (`B.rows() == nc * v.size()`). `nc` is the
     * number of constraints interpolated at every section: 1 for the plain loft
     * (positions only), 2 for the spine loft (position + spine tangent).
     *
     * The v-collocation matrix is assembled and factorized ONCE and every column is
     * solved together as multiple right-hand sides (the PR2 batched/hoisted solve).
     * The result is written straight into the surface's v-outer / u-inner pole
     * layout, so neither caller needs a transpose round trip.
     *
     * @return A flattened vector of poles in the layout the surface ctor expects.
     */
    template <std::floating_point T, size_t dim, size_t nc>
    auto build_loft_surface_poles_core(const MatrixX<T> &B, const std::vector<T> &flat_v,
                                       const std::vector<T> &v, size_t q)
        -> std::vector<std::array<T, dim>>
    {
        const size_t nv = v.size();
        const auto n_poles_v = Eigen::Index(nc * nv);
        const size_t nu = size_t(B.cols()) / dim;

        MatrixX<T> N(n_poles_v, n_poles_v);
        build_poles_matrix<T, nc>(flat_v, v, q, size_t(n_poles_v), N);
        MatrixX<T> X = N.partialPivLu().solve(B);

        std::vector<std::array<T, dim>> poles(size_t(n_poles_v) * nu);
        for (Eigen::Index j = 0; j < n_poles_v; ++j)
            for (size_t i = 0; i < nu; ++i)
                for (size_t d = 0; d < dim; ++d)
                    poles[size_t(j) * nu + i][d] = X(j, Eigen::Index(i * dim + d));

        return poles;
    }

    /**
     * Builds the poles for a (plain) loft surface from a set of curve poles, using
     * NURBS mathematics. Thin wrapper over build_loft_surface_poles_core for the
     * positions-only (nc == 1) case.
     *
     * @param poles_curves The poles of the curves used to build the loft surface.
     * @param flat_v The flattened knot vector in the v direction.
     * @param v The knot vector in the v direction.
     * @param flat_u The flattened knot vector in the u direction (unused; kept for API).
     * @param p The degree of the curve in the u direction (unused; kept for API).
     * @param q The degree of the curve in the v direction.
     * @return A flattened vector of poles representing the loft surface.
     */
    template <typename T, size_t dim>
    auto build_loft_surface_poles(const std::vector<std::vector<std::array<T, dim>>> &poles_curves,
                            const std::vector<T> &flat_v, const std::vector<T> &v, const std::vector<T> &flat_u, size_t p,
                            size_t q)
    {
        size_t nu = poles_curves[0].size();
        size_t nv = poles_curves.size();

        // Assemble every right-hand side at once (one column per (u-index, dim)),
        // read straight from poles_curves[j][i][d] — no transpose copy — and let
        // the shared core factorize once and batch-solve all columns.
        MatrixX<T> B(Eigen::Index(nv), Eigen::Index(nu * dim));
        for (size_t i = 0; i < nu; ++i)
            for (size_t d = 0; d < dim; ++d)
                for (size_t j = 0; j < nv; ++j)
                    B(Eigen::Index(j), Eigen::Index(i * dim + d)) = poles_curves[j][i][d];

        return build_loft_surface_poles_core<T, dim, 1>(B, flat_v, v, q);
    }

    /**
     * Creates a lofted surface from a series of B-spline curves.
     *
     * The function performs a loft operation, which is a common technique in computer-aided design
     * and computer graphics for creating a surface that smoothly interpolates between a series of curves.
     * This is achieved by first unifying the degrees and knot vectors of the input B-spline curves.
     * The control points for the lofted surface are then computed based on these unified curves.
     *
     * @tparam T The floating point type used for curve coordinates and knots.
     * @tparam dim The dimensionality of the control points.
     * @param curves_info A vector of BSCurveInfo structures, each representing a B-spline curve
     *                    with control points, a knot vector, and a degree.
     * @param v A vector of parameter values in the 'v' direction for the lofting process.
     * @param flat_v A flattened knot vector in the 'v' direction for the lofted surface.
     * @param q The degree of the lofted surface in the 'v' direction.
     * @return A tuple containing the control points of the lofted surface, the knot vector in the 'u' 
     *         direction, and the degree of the surface in the 'u' direction.
     *
     * Example Usage:
     *     auto surface = loft(curves, v_values, flat_v, degree_v);
     */
    template <std::floating_point T, size_t dim>
    auto loft(std::vector<BSCurveInfo<T, dim>> &curves_info, const std::vector<T> &v, const std::vector<T> &flat_v, size_t q)
    {
        // Unify the degrees and knot vectors of the input curves
        unify_degree(curves_info);
        unify_knots(curves_info);

        // Store the unified degree and the size of the knot vector
        auto p = std::get<2>(curves_info.front());
        auto ncr = std::distance(curves_info.begin(), curves_info.end());
        auto nu = std::get<0>(curves_info.front()).size();
        auto flat_u = std::get<1>(curves_info.front());

        // Extract the control points from each curve
        std::vector<std::vector<std::array<T, dim>>> poles_curves(ncr);
        std::transform(
            curves_info.begin(), curves_info.end(),
            poles_curves.begin(),
            [](const auto &curve_info) { return std::get<0>(curve_info); }
        );

        // Build and return the lofted surface
        return std::make_tuple(build_loft_surface_poles(poles_curves, flat_v, v, flat_u, p, q), flat_u, p);
    }

    /**
     * Generates a lofted surface from a series of B-spline curves.
     *
     * This function creates a lofted surface by connecting a set of B-spline curves specified by the `curves_info`.
     * It first computes a simple multiple flat knot vector in the 'v' direction based on the provided parameters. 
     * The function then uses this knot vector, along with the curve information, to create the lofted surface. 
     * The result includes the control points of the surface, the knot vectors in both the 'u' and 'v' directions, 
     * and the degree in the 'u' direction.
     *
     * @tparam T Floating point type used for curve coordinates and knots.
     * @tparam dim Dimensionality of the control points (poles).
     * @param curves_info A vector of BSCurveInfo<T, dim>, each containing essential information of a B-spline curve.
     * @param v A vector of parameter values in the 'v' direction for the lofting process.
     * @param q The degree of the lofted surface in the 'v' direction.
     * @return A tuple containing the following elements:
     *         - The control points of the lofted surface.
     *         - The knot vector in the 'u' direction.
     *         - The knot vector in the 'v' direction.
     *         - The degree of the surface in the 'u' direction.
     *
     * Example Usage:
     *     std::vector<BSCurveInfo<float, 3>> curves_info = ...; // Curve information
     *     std::vector<float> v = ...; // Parameter values in the 'v' direction
     *     auto [poles, flat_u, flat_v, p] = loft(curves_info, v, degree_v); // Lofted surface
     */
    template <std::floating_point T, size_t dim>
    auto loft(std::vector<BSCurveInfo<T, dim>> &curves_info, const std::vector<T> &v, size_t q)
    {
        // Build a simple multiple flat knot vector in the 'v' direction
        auto flat_v = build_simple_mult_flat_knots(v, q);

        // Delegate to the main loft function to create the lofted surface
        auto [poles, flat_u, p] = loft(curves_info, v, flat_v, q);
        return std::make_tuple(poles, flat_u, flat_v, p);
    }

    template <std::floating_point T, size_t dim>
    auto loft(std::vector<BSCurveInfo<T, dim>> &curves_info, size_t q_max)
    {
        // TODO: Compute an optimal parametrization
        auto n_curves = curves_info.size();
        auto v = make_range(T{}, T{1}, n_curves);
        auto q = std::min(q_max, n_curves-1);
        // Build a simple multiple flat knot vector in the 'v' direction
        auto flat_v = build_simple_mult_flat_knots(v, q);

        // Delegate to the main loft function to create the lofted surface
        auto [poles, flat_u, p] = loft(curves_info, v, flat_v, q);
        return std::make_tuple(poles, flat_u, flat_v, p, q);
    }
}