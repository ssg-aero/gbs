#pragma once
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
        // Pre-size to nu*nv and evaluate the flattened (u,v) grid in one size-gated
        // parallel transform (#89), keeping the U-first, V-outer order. This drops
        // the serial outer v-loop with insert() (and its double-parallelization: a
        // par transform nested in a serial for_each) for a single flat transform —
        // each node writes a distinct slot, so the result is bit-identical.
        const auto uv = bssanalysis_detail::grid_params(srf, nu, nv);

        points_vector<T, dim> points(nu * nv);
        transform_threshold(
            uv.begin(), uv.end(), points.begin(),
            [&srf](const std::array<T, 2> &uv_) { return srf(uv_[0], uv_[1]); });

        return points;
    }

    /**
     * @brief The unit normal of a surface at (u, v): S_u x S_v, normed.
     *
     * Where S_u x S_v vanishes against max(|S_u|, |S_v|)^2 (1e-6, the rule of the meshes of
     * gbs-brep), at a pole of a sphere or of a surface of revolution closed on its axis, at the
     * apex of a cone, along a collapsed row of poles, the normal is its limit from inside the
     * domain: taken at 1, 2 and 3 steps of 1e-4 of the way from (u, v) to the centre of the
     * bounds, and extrapolated to (u, v) (error in the cube of the step). At a pole, that is the
     * normal of the surface, which has one there, whatever the parametrization. At an apex, where
     * the surface has none, it is that of the side the centre of the domain lies on.
     *
     * @throws std::domain_error if the surface has no normal near (u, v) either (degenerate around
     * it), or (u, v) is not finite
     */
    template <typename T>
    auto normal_direction(const Surface<T, 3> &srf, T u, T v) -> point<T, 3>
    {
        const auto at = [&](T u_, T v_) -> std::optional<point<T, 3>> {
            const auto su = srf.value(u_, v_, 1, 0), sv = srf.value(u_, v_, 0, 1);
            const auto n = cross(su, sv);
            const T l = norm(n), scale = std::max(norm(su), norm(sv));
            if (!(l > T(1e-6) * scale * scale) || !std::isfinite(l))
                return std::nullopt;
            return n / l;
        };
        if (const auto n = at(u, v))
            return *n;
        const auto [u1, u2, v1, v2] = srf.bounds();
        T du = (u1 + u2) / 2 - u, dv = (v1 + v2) / 2 - v;
        if (du == T(0) && dv == T(0)) // at the centre: towards a corner
            du = (u2 - u1) / 2, dv = (v2 - v1) / 2;
        constexpr T h = T(1e-4);
        const auto n1 = at(u + h * du, v + h * dv), n2 = at(u + 2 * h * du, v + 2 * h * dv), n3 = at(u + 3 * h * du, v + 3 * h * dv);
        if (!n1 || !n2 || !n3)
            throw std::domain_error("normal_direction: the surface has no normal near (" + std::to_string(u) + ", " + std::to_string(v) + ")");
        const auto n = T(3) * *n1 - T(3) * *n2 + *n3; // the quadratic through the three, at 0
        return n / norm(n);
    }

    /**
     * @brief The points of a surface on a grid of nu x nv evenly spaced parameters (discretize) and
     * the unit normals there (normal_direction), in the same order: u first, v outer.
     *
     * @throws std::domain_error where normal_direction does
     */
    template <typename T>
    auto discretize_with_normals(const Surface<T, 3> &srf, size_t nu, size_t nv)
    {
        auto points = discretize(srf, nu, nv);
        const auto uv = bssanalysis_detail::grid_params(srf, nu, nv);
        points_vector<T, 3> normals(uv.size());
        const auto normal = [&srf](const std::array<T, 2> &uv_) { return normal_direction(srf, uv_[0], uv_[1]); };
        if (uv.size() >= parallel_min_size) // a domain_error may escape: never out of a parallel algorithm
            transform_rethrow(uv.begin(), uv.end(), normals.begin(), normal);
        else
            std::transform(uv.begin(), uv.end(), normals.begin(), normal);
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