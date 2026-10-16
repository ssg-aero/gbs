#pragma once
#include <gbs/solvers.h>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>
#include <gbs/curves>
#include <gbs/curveline.h>
#include <gbs/bssurf.h>
#include <gbs/intersect/projection.h>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

#include <cmath>
#include <optional>
#include <stdexcept>
#include <gbs/execution.h>

namespace gbs
{
    /// Default of the `solver` parameter of the extrema, which is ignored (deprecated, removed with NLopt).
    static const nlopt::algorithm default_nlopt_algo=nlopt::LN_COBYLA;

    /**
     * @brief Throws if the bounds cannot be sampled to bracket an initial guess
     * (unbounded curve such as a Line: u2 - u1 overflows to infinity)
     */
    template <typename T, size_t n>
    void check_bracket_bounds(const std::array<T, n> &bounds)
    {
        for (size_t i{}; i < n; i += 2)
        {
            if (!std::isfinite(bounds[i + 1] - bounds[i]))
            {
                throw std::invalid_argument("extrema: cannot bracket an unbounded parameter range, provide an initial guess or a trimmed curve");
            }
        }
    }

    // The extrema below are computed by Newton's method (gbs/intersect/projection.h) and return
    // the true distance. Their `solver` parameter is ignored: it is kept until NLopt is removed
    // (docs/sources/design/nlopt_removal.md). `tol_x` is the tolerance on the parameters: Newton
    // stops when every step is below tol_x² (1 + |x|).

    /**
     * @brief Project point on curve, by Newton from u0
     *
     * @param crv    : the curve
     * @param pnt    : the point
     * @param u0     : start value
     * @param tol_x  : tolerance on the parameter
     * @param solver : ignored (deprecated)
     * @return {u, distance}
     */
    template <typename T, size_t dim>
    auto extrema_curve_point(const Curve<T, dim> &crv, const std::array<T, dim> &pnt, T u0,T tol_x,[[maybe_unused]] nlopt::algorithm solver=default_nlopt_algo) -> std::array<T,2>
    {
        return newton_point_curve(crv, pnt, u0, crv.bounds(), tol_x);
    }
    /**
     * @brief Project point on curve: n_bracket samples, Newton from the two nearest
     *
     * @param crv    : the curve, bounded
     * @param pnt    : the point
     * @param tol_u  : tolerance on the parameter
     * @param solver : ignored (deprecated)
     * @return {u, distance}
     */
    template <typename T, size_t dim>
    auto extrema_curve_point(const Curve<T, dim> &crv, const std::array<T, dim> &pnt,T tol_u,[[maybe_unused]] nlopt::algorithm solver=default_nlopt_algo,size_t n_bracket=30) -> std::array<T,2>
    {
        check_bracket_bounds(crv.bounds());
        return project_point_curve(crv, pnt, crv.bounds(), n_bracket, tol_u);
    }
    /**
     * @brief Project point on surface, by Newton from (u0, v0)
     *
     * @param srf    : the surface
     * @param pnt    : the point
     * @param u0, v0 : start values
     * @param tol_x  : tolerance on the parameters
     * @param solver : ignored (deprecated)
     * @return {u, v, distance}
     */
    template <typename T, size_t dim>
    auto extrema_surf_pnt(const Surface<T, dim> &srf, const std::array<T, dim> &pnt, T u0, T v0, T tol_x, [[maybe_unused]] nlopt::algorithm solver=default_nlopt_algo)
    {
        return newton_point_surface(srf, pnt, u0, v0, srf.bounds(), tol_x);
    }
    /**
     * @brief Project point on surface: n_bracket_u x n_bracket_v samples, Newton from the two nearest
     *
     * @param srf    : the surface, bounded
     * @param pnt    : the point
     * @param tol_x  : tolerance on the parameters
     * @param solver : ignored (deprecated)
     * @return {u, v, distance}
     */
    template <typename T, size_t dim>
    auto extrema_surf_pnt(const Surface<T, dim> &srf, const std::array<T, dim> &pnt, T tol_x, [[maybe_unused]] nlopt::algorithm solver=default_nlopt_algo, size_t n_bracket_u = 30, size_t n_bracket_v = 30)
    {
        check_bracket_bounds(srf.bounds());
        return project_point_surface(srf, pnt, srf.bounds(), n_bracket_u, n_bracket_v, tol_x);
    }
    /**
     * @brief Closest points of two curves, by Newton from (u10, u20)
     *
     * @return {u1, u2, distance}
     */
    template <typename T, size_t dim>
    auto extrema_curve_curve(const Curve<T, dim> &crv1, const Curve<T, dim> &crv2, T u10, T u20,T tol_x, [[maybe_unused]] nlopt::algorithm solver=default_nlopt_algo)
    {
        return newton_curve_curve(crv1, crv2, u10, u20, crv1.bounds(), crv2.bounds(), tol_x);
    }
    /**
     * @brief Closest points of two curves: n_bracket_u x n_bracket_v samples, Newton from the two nearest pairs
     *
     * @return {u1, u2, distance}
     */
    template <typename T, size_t dim>
    auto extrema_curve_curve(const Curve<T, dim> &crv1, const Curve<T, dim> &crv2,T tol_x, [[maybe_unused]] nlopt::algorithm solver=default_nlopt_algo, size_t n_bracket_u = 30, size_t n_bracket_v = 30)
    {
        check_bracket_bounds(crv1.bounds());
        check_bracket_bounds(crv2.bounds());
        return closest_curve_curve(crv1, crv2, crv1.bounds(), crv2.bounds(), n_bracket_u, n_bracket_v, tol_x);
    }

/**
 * @brief 
 * 
 * @tparam T 
 * @param crv1 
 * @param crv2 
 * @return tuple (u1,u2,0,0) last 2 zeros are for consitency with overloaded function
 */
    template <typename T>
    auto extrema_curve_curve(const Line<T, 2> &crv1, const Line<T, 2> &crv2) // TODO  raise if // and use tuple for other functions
    {
        auto [P1, N1] = crv1.getAx();
        auto [P2, N2] = crv2.getAx();
        auto P21 = P2 - P1;
        auto d  = det(std::array<T,4>{ N1[0], -N2[0],  N1[1], -N2[1]});
        auto l1 = det(std::array<T,4>{P21[0], -N2[0], P21[1], -N2[1]}) / d;
        auto l2 = det(std::array<T,4>{ N1[0], P21[0],  N1[1], P21[1]}) / d;
        return std::array<T,2>{l1,l2};
    }
/**
 * @brief Closest points of a curve and each curve of a list
 * 
 * @param crv 
 * @param first, last : iterators on pointers to curves
 * @param tol_x 
 * @param solver : ignored (deprecated)
 * @param n_bracket_u, n_bracket_v : samples of crv and of each curve of the list
 * @return std::vector of {u, u_i, distance}
 */
    template <typename T, size_t dim, typename ForwardIt>
    auto extrema_curve_curve_lst(const Curve<T,dim> &crv, ForwardIt first,ForwardIt last,T tol_x, [[maybe_unused]] nlopt::algorithm solver=default_nlopt_algo, size_t n_bracket_u = 30, size_t n_bracket_v = 30)
    {
        auto n = std::distance(first, last);
        std::vector<std::array<T,3>> intersection_info(n);
        std::transform(
            GBS_PAR_EXEC
            first, last,
            intersection_info.begin(),
            [&crv,tol_x,n_bracket_u,n_bracket_v](const auto &crv_){
                return extrema_curve_curve(crv, *crv_, tol_x, default_nlopt_algo, n_bracket_u, n_bracket_v);
            }
        );
        return intersection_info;
    }

    /**
     * @brief Closest points of a surface and a curve, by Newton from (u_c0, u_s0, v_s0)
     *
     * @return {u, v, u_c, distance}
     */
    template <typename T, size_t dim>
    auto extrema_surf_curve(const Surface<T, dim> &srf, const Curve<T, dim> &crv, T u_c0, T u_s0, T v_s0, T tol_x, [[maybe_unused]] nlopt::algorithm solver=default_nlopt_algo)
    {
        const auto [t, u, v, d] = newton_curve_surface(crv, srf, u_c0, u_s0, v_s0, crv.bounds(), srf.bounds(), tol_x);
        return std::array<T,4>{u, v, t, d};
    }

    /**
     * @brief Closest points of a surface and a curve: n_bracket_u x n_bracket_v x n_bracket_w
     * samples (u, v of the surface, parameter of the curve), Newton from the two nearest triples
     *
     * @return {u, v, u_c, distance}
     */
    template <typename T, size_t dim>
    auto extrema_surf_curve(const Surface<T, dim> &srf, const Curve<T, dim> &crv, T tol_x, [[maybe_unused]] nlopt::algorithm solver=default_nlopt_algo, size_t n_bracket_u = 30, size_t n_bracket_v = 30, size_t n_bracket_w = 30) 
    {
        check_bracket_bounds(srf.bounds());
        check_bracket_bounds(crv.bounds());
        const auto [t, u, v, d] = closest_curve_surface(crv, srf, crv.bounds(), srf.bounds(), n_bracket_w, n_bracket_u, n_bracket_v, tol_x);
        return std::array<T,4>{u, v, t, d};
    }

    template <typename T, size_t dim>
    auto solve_coordinate(const Curve<T, dim> &crv, size_t iCoord, T Xi, std::optional<T> u0 = std::nullopt, int digits = 20, size_t maxit = 20)
    {
        if(iCoord>=dim)
        {
            throw std::out_of_range("Invalid size.");
        }

        using namespace boost::math::tools;

        auto [u_min, u_max] = crv.bounds();
        auto f_eq = [&crv, iCoord, Xi](T u) {
            return std::make_tuple(
                crv(u)[iCoord] - Xi,
                crv(u, 1)[iCoord],
                crv(u, 2)[iCoord]);
        };

        T u0_ = u0.value_or(0.5*(u_min+u_max));

        return newton_raphson_iterate(f_eq, u0_, u_min, u_max, digits, maxit);
    }

    // template <typename T, size_t dim>
    // auto max_coordinate(const Curve<T, dim> &crv, size_t iCoord, T Xi, std::optional<T> u0 = std::nullopt, int digits = 20, size_t maxit = 20)
    // {
    //     if(iCoord>=dim)
    //     {
    //         throw std::out_of_range("Invalid size.");
    //     }

    //     using namespace boost::math::tools;

    //     auto [u_min, u_max] = crv.bounds();
    //     auto f_eq = [&crv, iCoord, Xi](T u) {
    //         auto x = crv(u)[iCoord];
    //         auto dx = crv(u, 1)[iCoord];
    //         auto ddx = crv(u, 1)[iCoord];
    //         return std::make_tuple(
    //             1/x*x,
    //             -2*dx/x/x/x,
    //             -2*ddx/x/x/x + 6 * dx*dx /x/x/x/x
    //             );
    //     };

    //     T u0_ = u0.value_or(0.5*(u_min+u_max));

    //     return newton_raphson_iterate(f_eq, u0_, u_min, u_max, digits, maxit);
    // }

    template <typename T, size_t dim>
    auto max_coordinate(const Curve<T, dim> &crv, size_t iCoord, T Xi, int digits = 20, size_t maxit = 20)
    {
        if(iCoord>=dim)
        {
            throw std::out_of_range("Invalid size.");
        }

        using namespace boost::math::tools;

        auto [u_min, u_max] = crv.bounds();
        auto f_eq = [&crv, iCoord, Xi](T u) {
            auto x = crv(u)[iCoord];
            return  1/x/x;
        };

        return brent_find_minima(f_eq, u_min, u_max, digits, maxit);
    }

    template <typename T, size_t dim>
    auto min_coordinate(const Curve<T, dim> &crv, size_t iCoord, T Xi, int digits = 20, size_t maxit = 20)
    {
        if(iCoord>=dim)
        {
            throw std::out_of_range("Invalid size.");
        }

        using namespace boost::math::tools;

        auto [u_min, u_max] = crv.bounds();
        auto f_eq = [&crv, iCoord, Xi](T u) {
            auto x = crv(u)[iCoord];
            return  x*x;
        };

        return brent_find_minima(f_eq, u_min, u_max, digits, maxit);
    }

} // namespace gbs