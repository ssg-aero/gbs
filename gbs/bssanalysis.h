#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <list>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <gbs/execution.h>
#include <gbs/bssurf.h>
namespace gbs
{
    namespace bssanalysis_detail
    {
        /// The parameters of a grid of nu x nv points evenly spaced on the bounds of a surface, u
        /// first, v outer
        template <typename T, size_t dim>
        auto grid_params(const Surface<T, dim> &srf, size_t nu, size_t nv) -> std::vector<std::array<T, 2>>
        {
            auto [u1, u2, v1, v2] = srf.bounds();
            auto u = make_range(u1, u2, nu);
            auto v = make_range(v1, v2, nv);
            std::vector<std::array<T, 2>> uv(nu * nv);
            for (size_t j = 0; j < nv; ++j)
                for (size_t i = 0; i < nu; ++i)
                    uv[j * nu + i] = {u[i], v[j]};
            return uv;
        }
    }

    template <typename T, size_t dim>
    auto discretize(const Surface<T, dim> &srf, size_t nu, size_t nv) -> points_vector<T, dim>
    {
        // Evaluate the flattened (u,v) grid in one size-gated parallel loop (#89),
        // keeping the U-first, V-outer order. This drops the serial outer v-loop with
        // insert() (and its double-parallelization: a par transform nested in a serial
        // for_each) for a single flat loop — each node writes a distinct slot, so the
        // result is bit-identical. A surface that throws at a node (an offset where the
        // basis has no normal) has its first exception, in the order of the grid,
        // rethrown after the loop (build_batch): out of a parallel algorithm it would
        // terminate the process.
        const auto uv = bssanalysis_detail::grid_params(srf, nu, nv);
        return build_batch(
            uv, [&srf](const std::array<T, 2> &uv_) { return srf(uv_[0], uv_[1]); }, parallel_min_size);
    }

    /**
     * @brief The unit normal of a surface at (u, v): S_u x S_v, normed.
     *
     * The derivatives are compared as S_u Δu and S_v Δv, Δu and Δv the lengths of the bounds (1
     * where unbounded), so that the rule depends neither on the units nor on a scaling of the
     * parameters: a point is regular where |S_u Δu x S_v Δv| > 1e-6 max(|S_u| Δu, |S_v| Δv)^2,
     * the rule of the meshes of gbs-brep on these derivatives. Elsewhere, at a pole of a sphere or
     * of a surface of revolution closed on its axis, at the apex of a cone, along a collapsed row
     * of poles, the normal is its limit from inside the domain: taken at 1, 2 and 3 steps of 1e-4
     * of the way from (u, v) to the centre of the bounds, and extrapolated to (u, v) (error in the
     * cube of the step). Where that way stays on the degenerate line (from the apex of a double
     * cone to its centre), the way to a corner is taken. At a pole, that is the normal of the
     * surface, which has one there, whatever the parametrization. At an apex, where the surface
     * has none, it is that of the side the way goes to.
     *
     * A surface more than 1e6 times longer than wide (|S_u| Δu < 1e-6 |S_v| Δv around the point)
     * never meets that rule: its normal is then taken, at the point or as the same limit, where
     * S_u and S_v make an angle of more than 1e-6 rad.
     *
     * @throws std::domain_error if the surface has no normal near (u, v) either (degenerate around
     * it), or (u, v) is not finite
     */
    template <typename T>
    auto normal_direction(const Surface<T, 3> &srf, T u, T v) -> point<T, 3>
    {
        const auto no_normal = [&] {
            return std::domain_error("normal_direction: the surface has no normal near (" + std::to_string(u) + ", " + std::to_string(v) + ")");
        };
        if (!std::isfinite(u) || !std::isfinite(v))
            throw no_normal();
        const auto [u1, u2, v1, v2] = srf.bounds();
        const auto length = [](T a, T b) { return std::isfinite(b - a) && b > a ? b - a : T(1); };
        const T lu = length(u1, u2), lv = length(v1, v2);
        // S_u x S_v normed where regular: |S_u Δu x S_v Δv| beyond 1e-6 max(|S_u| Δu, |S_v| Δv)
        // times the same or, on a thin surface, times the shorter (the sine of the angle)
        const auto at = [&](T u_, T v_, bool thin) -> std::optional<point<T, 3>> {
            const auto su = srf.value(u_, v_, 1, 0), sv = srf.value(u_, v_, 0, 1);
            const auto n = cross(su, sv);
            const T l = norm(n), a = norm(su) * lu, b = norm(sv) * lv;
            if (!(l * lu * lv > T(1e-6) * std::max(a, b) * (thin ? std::min(a, b) : std::max(a, b))) || !std::isfinite(l * lu * lv))
                return std::nullopt;
            return n / l;
        };
        // the quadratic through the normals at 1, 2 and 3 steps of the way from (u, v) to (tu, tv), at 0
        constexpr T h = T(1e-4);
        const auto limit = [&](T tu, T tv, bool thin) -> std::optional<point<T, 3>> {
            const T du = tu - u, dv = tv - v;
            if (!std::isfinite(du) || !std::isfinite(dv) || (du == T(0) && dv == T(0)))
                return std::nullopt;
            std::array<point<T, 3>, 3> n;
            for (int k = 0; k < 3; ++k)
                if (const auto nk = at(u + T(k + 1) * h * du, v + T(k + 1) * h * dv, thin))
                    n[k] = *nk;
                else
                    return std::nullopt;
            const auto m = T(3) * n[0] - T(3) * n[1] + n[2];
            const T l = norm(m);
            if (!(l > T(0)) || !std::isfinite(l))
                return std::nullopt;
            return m / l;
        };
        // the centre of the bounds, then their corners (an unbounded side one unit away)
        const auto side = [](T b, T x, T dx) { return std::isfinite(b) ? b : x + dx; };
        const T a1 = side(u1, u, -1), a2 = side(u2, u, 1), b1 = side(v1, v, -1), b2 = side(v2, v, 1);
        const std::array<std::array<T, 2>, 5> towards{{{(a1 + a2) / 2, (b1 + b2) / 2}, {a1, b1}, {a2, b1}, {a1, b2}, {a2, b2}}};
        for (const bool thin : {false, true})
        {
            if (const auto n = at(u, v, thin))
                return *n;
            for (const auto &[tu, tv] : towards)
                if (const auto n = limit(tu, tv, thin))
                    return *n;
        }
        throw no_normal();
    }

    /**
     * @brief The points of a surface on a grid of nu x nv evenly spaced parameters (discretize) and
     * the unit normals there (normal_direction), in the same order: u first, v outer.
     *
     * One loop on the grid gives the point and the normal of each node; the first exception, in
     * the order of the grid, is rethrown after it (see discretize).
     *
     * @throws std::domain_error where normal_direction does
     */
    template <typename T>
    auto discretize_with_normals(const Surface<T, 3> &srf, size_t nu, size_t nv)
    {
        const auto uv = bssanalysis_detail::grid_params(srf, nu, nv);
        const auto nodes = build_batch(
            uv, [&srf](const std::array<T, 2> &uv_) { return std::array{srf(uv_[0], uv_[1]), normal_direction(srf, uv_[0], uv_[1])}; },
            parallel_min_size);
        points_vector<T, 3> points(nodes.size()), normals(nodes.size());
        for (size_t i = 0; i < nodes.size(); ++i)
            points[i] = nodes[i][0], normals[i] = nodes[i][1];
        return std::make_tuple(std::move(points), std::move(normals));
    }

/**
 * @brief Refines the range of isoV curve parameters recursively based on the deviation.
 * 
 * @tparam T Floating-point type
 * @tparam dim Dimension of the curve
 * @param srf The surface to discretize
 * @param v the v position on surface
 * @param dev_max The maximum allowed deviation
 * @param it_u1 Iterator pointing to the start of the range in the parameter list
 * @param it_u3 Iterator pointing to the end of the range in the parameter list
 * @param u_lst The list of curve parameters
 */
    template <typename T, size_t dim>
    void refine_params_u_recursive(const Surface<T, dim> &srf, T v, T dev_max, typename std::list<T>::iterator it_u1, typename std::list<T>::iterator it_u3, std::list<T> &u_lst)
    {
        auto u_ = 0.5 * (*it_u1 + *it_u3);
        auto v1 = srf(u_, v) - srf(*it_u1, v);
        auto v2 = srf(*it_u3, v) - srf(*it_u1, v);
        auto dev_ = norm(cross(v1, v2)) / (norm(v1) * norm(v2));

        if (dev_ > dev_max)
        {
            auto it_u2 = u_lst.insert(it_u3, u_);
            refine_params_u_recursive(srf, v, dev_max, it_u1, it_u2, u_lst);
            refine_params_u_recursive(srf, v, dev_max, it_u2, it_u3, u_lst);
        }
    }

/**
 * @brief Generates a list of isoV curve parameters based on deviation using a recursive approach.
 * 
 * @tparam T Floating-point type
 * @tparam dim Dimension of the curve
 * @param srf The surface to discretize
 * @param v the v position on surface
 * @param u1 Start value of the parameter range
 * @param u2 End value of the parameter range
 * @param n Initial number of points to create
 * @param dev_max The maximum allowed deviation
 * @param n_max_pts Maximum number of points allowed
 * @return A list of curve parameters with the refined deviation
 */
    template <typename T, size_t dim>
    auto deviation_based_u_params(const Surface<T, dim> &srf, T u1, T u2, T v,size_t n, T dev_max, size_t n_max_pts = approx_max_sample_points) -> std::list<T>
    {
        // Create initial list of parameters
        std::list<T> u_lst;
        for (size_t i{}; i < n; i++)
            u_lst.push_back(u1 + (u2 - u1) * i / (n - 1.0));

        // Refine the list recursively if the maximum number of points is not reached
        if (u_lst.size() < n_max_pts)
        {
            auto it_u1 = u_lst.begin();
            auto it_u3 = std::next(it_u1);
            while (it_u3 != u_lst.end())
            {
                refine_params_u_recursive(srf, v, dev_max, it_u1, it_u3, u_lst);

                it_u1 = it_u3;
                it_u3 = std::next(it_u1);
            }
        }
        return u_lst;
    }

/**
 * @brief Refines the range of isoU curve parameters recursively based on the deviation.
 * 
 * @tparam T Floating-point type
 * @tparam dim Dimension of the curve
 * @param srf The surface to discretize
 * @param u the u position on surface
 * @param dev_max The maximum allowed deviation
 * @param it_v1 Iterator pointing to the start of the range in the parameter list
 * @param it_v3 Iterator pointing to the end of the range in the parameter list
 * @param v_lst The list of curve parameters
 */
    template <typename T, size_t dim>
    void refine_params_v_recursive(const Surface<T, dim> &srf, T u, T dev_max, typename std::list<T>::iterator it_v1, typename std::list<T>::iterator it_v3, std::list<T> &v_lst)
    {
        auto u_ = 0.5 * (*it_v1 + *it_v3);
        auto v1 = srf(u, u_) - srf(u, *it_v1);
        auto v2 = srf(u, *it_v3) - srf(u, *it_v1);
        auto dev_ = norm(cross(v1, v2)) / (norm(v1) * norm(v2));

        if (dev_ > dev_max)
        {
            auto it_u2 = v_lst.insert(it_v3, u_);
            refine_params_v_recursive(srf, u, dev_max, it_v1, it_u2, v_lst);
            refine_params_v_recursive(srf, u, dev_max, it_u2, it_v3, v_lst);
        }
    }

/**
 * @brief Generates a list of isoU curve parameters based on deviation using a recursive approach.
 * 
 * @tparam T Floating-point type
 * @tparam dim Dimension of the curve
 * @param srf The curve to discretize
* @param srf The surface to discretize
 * @param u the u position on surface
 * @param u1 Start value of the parameter range
 * @param u2 End value of the parameter range
 * @param n Initial number of points to create
 * @param dev_max The maximum allowed deviation
 * @param n_max_pts Maximum number of points allowed
 * @return A list of curve parameters with the refined deviation
 */
    template <typename T, size_t dim>
    auto deviation_based_v_params(const Surface<T, dim> &srf, T v1, T v2, T u,size_t n, T dev_max, size_t n_max_pts = approx_max_sample_points) -> std::list<T>
    {
        // Create initial list of parameters
        std::list<T> v_lst;
        for (size_t i{}; i < n; i++)
            v_lst.push_back(v1 + (v2 - v1) * i / (n - 1.0));

        // Refine the list recursively if the maximum number of points is not reached
        if (v_lst.size() < n_max_pts)
        {
            auto it_v1 = v_lst.begin();
            auto it_v3 = std::next(it_v1);
            while (it_v3 != v_lst.end())
            {
                refine_params_v_recursive(srf, u, dev_max, it_v1, it_v3, v_lst);

                it_v1 = it_v3;
                it_v3 = std::next(it_v1);
            }
        }
        return v_lst;
    }

}