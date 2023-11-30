#pragma once
#include <gbs/solvers.h>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>
#include <gbs/curves>
#include <gbs/curveline.h>
#include <gbs/bssurf.h>
#include <gbs/vecop.h>

#include <optional>

namespace gbs
{
    static const nlopt::algorithm default_nlopt_algo=nlopt::LN_COBYLA;

    /**
     * @brief Project point on curve
     * 
     * @tparam T 
     * @tparam dim 
     * @param pnt     : the point
     * @param crv     : the curve
     * @param u0      : guess value
     * @param tol_x   : tolerance
     * @param solver : solver type please have look to https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#nomenclature to change this value
     * @return auto
     */
    template <typename T, size_t dim>
    auto extrema_curve_point(const Curve<T, dim> &crv, const std::array<T, dim> &pnt, T u0,T tol_x,nlopt::algorithm solver=default_nlopt_algo) -> std::array<T,2>
    {
        auto f = [&pnt,&crv](const std::vector<double> &x)
        {
            return std::vector<double>{sq_norm(crv(x[0]) - pnt)}; //double req for nlop
            // return std::vector<T>{crv(x[0],1)*(crv(x[0]) - pnt)};
        };

        auto g = [&pnt,&crv](const std::vector<double> &x,const std::vector<double> &r)
        {
            auto dfdu = 2.*crv(x[0],1)*(crv(x[0]) - pnt);
            // auto dfdu = crv(x[0],2)*(crv(x[0]) - pnt) + sq_norm(crv(x[0],1)) ;
            return std::vector<double>{2. * r[0] * dfdu}; //double req for nlop
        };

        std::vector<T> x{u0};
        std::vector<T> lb{crv.bounds()[0]};
        std::vector<T> hb{crv.bounds()[1]};
        auto minf = solve_D_nlop(
            f,g,
            x, lb, hb,
            tol_x*tol_x,
            solver
            );
        // return {T(x[0]),T(sqrt(minf))};
        return std::array<T,2>{T(x[0]),T(sqrt(minf))};
    }
    /**
     * @brief Project point on curve
     * 
     * @tparam T 
     * @tparam dim 
     * @param crv 
     * @param pnt 
     * @param tol_u 
     * @return auto 
     */
    template <typename T, size_t dim>
    auto extrema_curve_point(const Curve<T, dim> &crv, const std::array<T, dim> &pnt,T tol_u,nlopt::algorithm solver=default_nlopt_algo,size_t n_bracket=30) -> std::array<T,2>
    {
        auto u = make_range<T>(crv.bounds(),n_bracket);
        auto comp = [&](auto u1, auto u2){
            return sq_norm(pnt-crv(u1)) < sq_norm(pnt-crv(u2));
        };
        auto u0 = *std::min_element(std::execution::par , u.begin(), u.end(), comp);
        return extrema_curve_point(crv, pnt, u0,tol_u,solver);
    }
    /** Project point on surface
     * @brief 
     * 
     * @tparam T 
     * @tparam dim 
     * @tparam rational 
     * @param srf 
     * @param pnt 
     * @param u0 
     * @param v0 
     * @param tol_x 
     * @param solver 
     * @return auto
     */
    template <typename T, size_t dim>
    auto extrema_surf_pnt(const Surface<T, dim> &srf, const std::array<T, dim> &pnt, T u0, T v0, T tol_x, nlopt::algorithm solver=default_nlopt_algo)
    {
        auto f = [&pnt,&srf](const std::vector<double> &x)
        {
            return std::vector<double>{sq_norm(srf(x[0],x[1]) - pnt)};
        };

        auto g = [&pnt,&srf](const std::vector<double> &x,const std::vector<double> &r)
        {
            // return std::vector<T>{2.*crv(x[0],1)*(crv(x[0]) - pnt)};
            auto f_ = srf(x[0],x[1])  - pnt;
            auto dfdu = 2.*srf(x[0],x[1],1,0)* f_;
            auto dfdv = 2.*srf(x[0],x[1],0,1)* f_;
            return std::vector<double>{2. * r[0] * dfdu,2. * r[0] * dfdv};
        };

        auto [u1,u2,v1,v2] = srf.bounds();
        std::vector<T> x{u0,v0};
        std::vector<T> lb{u1,v1};
        std::vector<T> hb{u2,v2};
        auto minf = solve_D_nlop(
            f,g,
            x, lb, hb,
            tol_x*tol_x,
            solver
        );
        return std::array<T,3>{x[0], x[1], std::sqrt(static_cast<T>(minf))};
    }
    /**
     * @brief Project point on surface, initial value for solver is bracketed
     * 
     * @tparam T 
     * @tparam dim 
     * @param srf 
     * @param pnt 
     * @param tol_x 
     * @param solver 
     * @param n_bracket_u 
     * @param n_bracket_v 
     * @return auto 
     */
    template <typename T, size_t dim>
    auto extrema_surf_pnt(const Surface<T, dim> &srf, const std::array<T, dim> &pnt, T tol_x, nlopt::algorithm solver=default_nlopt_algo, size_t n_bracket_u = 30, size_t n_bracket_v = 30)
    {
        auto uv = make_range<T>(srf.bounds(), n_bracket_u , n_bracket_v);
        auto comp = [&](auto uv1, auto uv2){
            auto [u1, v1] = uv1;
            auto [u2, v2] = uv2;
            return sq_norm(pnt-srf(u1,v1)) < sq_norm(pnt-srf(u2,v2));
        };
        auto [u0, v0 ] = *std::min_element(std::execution::par , uv.begin(), uv.end(), comp);
        return extrema_surf_pnt(srf, pnt, u0, v0, tol_x, solver);
    }

    template <typename T, size_t dim>
    auto extrema_curve_curve(const Curve<T, dim> &crv1, const Curve<T, dim> &crv2, T u10, T u20,T tol_x, nlopt::algorithm solver=default_nlopt_algo)
    {
        auto f = [&crv1,&crv2](const std::vector<double> &x)
        {
            return std::vector<double>{sq_norm(crv1(x[0])-crv2(x[1]))};
        };

        auto g = [&crv1,&crv2](const std::vector<double> &x,const std::vector<double> &r)
        {
            auto f_ = crv1(x[0])-crv2(x[1]);
            auto dfdu = 2.*crv1(x[0],1)* f_;
            auto dfdv = 2.*crv2(x[1],1)* f_;
            return std::vector<double>{2. * r[0] * dfdu,2. * r[0] * dfdv};
        };

        auto [u1,u2] = crv1.bounds();
        auto [v1,v2] = crv2.bounds();
        std::vector<T> x{u10,u20};
        std::vector<T> lb{u1,v1};
        std::vector<T> hb{u2,v2};
        auto minf = solve_D_nlop(
            f,g,
            x, lb, hb,
            tol_x*tol_x,
            solver
        );
        return std::array<T,3>{T(x[0]),T(x[1]),T(sqrt(minf))};
    }

    template <typename T, size_t dim>
    auto extrema_curve_curve(const Curve<T, dim> &crv1, const Curve<T, dim> &crv2,T tol_x, nlopt::algorithm solver=default_nlopt_algo, size_t n_bracket_u = 30, size_t n_bracket_v = 30)// -> extrema_CC_result<T>
    {
        auto samples1 = make_range(crv1.bounds(), n_bracket_u); 
        auto samples2 = make_range(crv2.bounds(), n_bracket_v);

        T minDistance = std::numeric_limits<T>::max();
        std::pair<T, T> bestPair;

        for (T u : samples1) {
            for (T v : samples2) {
                T distance = sq_norm(crv1(u) - crv2(v)); 
                if (distance < minDistance) {
                    minDistance = distance;
                    bestPair = {u, v};
                }
            }
        }

        auto [u10, u20] = bestPair;
        return extrema_curve_curve<T,dim>(crv1,crv2,u10,u20,tol_x,solver);
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
 * @brief 
 * 
 * @tparam T 
 * @tparam dim 
 * @tparam ForwardIt 
 * @param crv 
 * @param first 
 * @param last 
 * @param tol_x 
 * @param solver 
 * @param n_bracket_u 
 * @param n_bracket_v 
 * @return auto 
 */
    template <typename T, size_t dim, typename ForwardIt>
    auto extrema_curve_curve_lst(const Curve<T,dim> &crv, ForwardIt first,ForwardIt last,T tol_x, nlopt::algorithm solver=default_nlopt_algo, size_t n_bracket_u = 30, size_t n_bracket_v = 30)
    {
        auto n = std::distance(first, last);
        std::vector<std::array<T,3>> intersection_info(n);
        std::transform(
            std::execution::par,
            first, last,
            intersection_info.begin(),
            [&crv,tol_x](const auto &crv_){
                return extrema_curve_curve(crv, *crv_, tol_x);
            }
        );
        return intersection_info;
    }

    template <typename T, size_t dim>
    auto extrema_surf_curve(const Surface<T, dim> &srf, const Curve<T, dim> &crv, T u_c0, T u_s0, T v_s0, T tol_x, nlopt::algorithm solver=default_nlopt_algo) //-> extrema_CS_result<T>
    {

        auto f = [&crv,&srf](const std::vector<double> &x)
        {
            return std::vector<double>{sq_norm(crv(x[0])-srf(x[1],x[2]))};
        };

        auto g = [&crv,&srf](const std::vector<double> &x,const std::vector<double> &r)
        {
            auto f_ = crv(x[0])-srf(x[1],x[2]);
            auto dfduc = 2.*crv(x[0],1)  * f_;
            auto dfdus =-2.*srf(x[1],1,0)* f_;
            auto dfdvs =-2.*srf(x[2],0,1)* f_;
            return std::vector<double>{2. * r[0] * dfduc,2. * r[0] * dfdus,2. * r[0] * dfdvs};
        };

        auto [uc1,uc2] = crv.bounds();
        auto [us1,us2,vs1,vs2] = srf.bounds();
        std::vector<T> x{u_c0,u_s0,v_s0};
        std::vector<T> lb{uc1,us1, vs1};
        std::vector<T> hb{uc2,us2, vs2};
        auto minf = solve_D_nlop(
            f,g,
            x, lb, hb,
            tol_x*tol_x,
            solver
        );
        return std::array<T,4>{x[1], x[2],x[0],std::sqrt(static_cast<T>(minf))};
    }

    template <typename T, size_t dim>
    auto extrema_surf_curve(const Surface<T, dim> &srf, const Curve<T, dim> &crv, T tol_x, nlopt::algorithm solver=default_nlopt_algo, size_t n_bracket_u = 30, size_t n_bracket_v = 30, size_t n_bracket_w = 30) 
    {
        auto uvw = make_range( srf.bounds(), crv.bounds(), n_bracket_u, n_bracket_v, n_bracket_w );
        auto comp = [&](auto uvw1, auto uvw2){
            auto [u1, v1, w1] = uvw1;
            auto [u2, v2, w2] = uvw2;
            return sq_norm(srf(u1, v1)-crv(w1)) < sq_norm(srf(u2, v2)-crv(w2));
        };
        auto [u0, v0, w0 ] = *std::min_element(std::execution::par , uvw.begin(), uvw.end(), comp);
        return extrema_surf_curve(srf, crv, u0, v0, w0, tol_x, solver);
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