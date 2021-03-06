#pragma once
#include <gbs/solvers.h>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>
#include <gbs/curves>
#include <gbs/curveline.h>
#include <gbs/bssurf.h>

#include <optional>

namespace gbs
{
    static const nlopt::algorithm default_nlopt_algo=nlopt::LN_COBYLA;

    template <typename T, typename F>
    auto approx_min_loc(const std::vector<T> &u, const F &f)
    {
        auto f_ = f(u.front());
        auto f_min = f_;
        auto u_min = u.front();
        for (auto u_ : u)
        {
            f_ = f(u_);
            if (f_ < f_min)
            {
                f_min = f_;
                u_min = u_;
            }
        }
        return std::make_pair(u_min,f_min);
    }

    template <typename T, typename F>
    auto approx_min_loc(const std::vector<T> &u, const std::vector<T> &v, const F &f)
    {
        auto sq_d = f(u.front(),v.front());
        auto sq_d_min = sq_d;
        auto U_min = std::array<T,2>{u.front(), v.front()};
        for (auto u_ : u)
        {
            // auto f_ = [u_](const auto & v_)
            // {
            //     return f(u_,v_);
            // };
            for (auto v_ : v)
            {
                // sq_d = f_(v_);
                sq_d = f(u_, v_);
                if (sq_d < sq_d_min)
                {
                    sq_d_min = sq_d;
                    U_min = {u_, v_};
                }
            }
        }
        return U_min;
    }
//TODO use variadic
    template <typename T, typename F>
    auto approx_min_loc(const std::vector<T> &u, const std::vector<T> &v, const std::vector<T> &w, const F &f)
    {
        auto sq_d = f(u.front(), v.front(), w.front());
        auto sq_d_min = sq_d;
        auto U_min = std::array<T, 3>{u.front(), v.front(), w.front()};
        for (auto u_ : u)
        {
            for (auto v_ : v)
            {
                for (auto w_ : w)
                {
                    sq_d = f(u_, v_, w_);
                    if (sq_d < sq_d_min)
                    {
                        sq_d_min = sq_d;
                        U_min = {u_, v_, w_};
                    }
                }
            }
        }
        return U_min;
    }
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
            return std::vector<double>{gbs::sq_norm(crv(x[0]) - pnt)}; //double req for nlop
            // return std::vector<T>{crv(x[0],1)*(crv(x[0]) - pnt)};
        };

        auto g = [&pnt,&crv](const std::vector<double> &x,const std::vector<double> &r)
        {
            auto dfdu = 2.*crv(x[0],1)*(crv(x[0]) - pnt);
            // auto dfdu = crv(x[0],2)*(crv(x[0]) - pnt) + gbs::sq_norm(crv(x[0],1)) ;
            return std::vector<double>{2. * r[0] * dfdu}; //double req for nlop
        };

        std::vector<T> x{u0};
        std::vector<T> lb{crv.bounds()[0]};
        std::vector<T> hb{crv.bounds()[1]};
        auto minf = gbs::solve_D_nlop(
            f,g,
            x, lb, hb,
            tol_x,
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
    auto extrema_curve_point(const Curve<T, dim> &crv, const std::array<T, dim> &pnt,T tol_u,nlopt::algorithm solver=default_nlopt_algo,size_t n_barcket=30) -> std::array<T,2>
    {
        auto u = make_range<T>(crv.bounds(),n_barcket);
        auto f = [&](const auto &u_)
        {
            return norm(pnt-crv(u_));
        };
        auto [u0,f0] = approx_min_loc(u,f);
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
            return std::vector<double>{gbs::sq_norm(srf(x[0],x[1]) - pnt)};
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
        auto minf = gbs::solve_D_nlop(
            f,g,
            x, lb, hb,
            tol_x,
            solver
        );
        return std::array<T,3>{x[0], x[1], sqrt(minf)};
    }
    template <typename T, size_t dim>
    auto extrema_surf_pnt(const Surface<T, dim> &srf, const std::array<T, dim> &pnt, T tol_x, nlopt::algorithm solver=default_nlopt_algo)
    {
        auto u = make_range<T>(srf.bounds()[0] , srf.bounds()[1],30);
        auto v = make_range<T>(srf.bounds()[2] , srf.bounds()[3],30);
        auto f = [&](const auto &u_,const auto &v_)
        {
            return norm(pnt-srf(u_,v_));
        };
        auto [u0,v0] = approx_min_loc(u,v,f);
        return extrema_surf_pnt(srf, pnt, u0, v0, tol_x, solver);
    }

    template <typename T, size_t dim>
    auto extrema_curve_curve(const Curve<T, dim> &crv1, const Curve<T, dim> &crv2, T u10, T u20,T tol_x, nlopt::algorithm solver=default_nlopt_algo)
    {
        auto f = [&crv1,&crv2](const std::vector<double> &x)
        {
            return std::vector<double>{gbs::sq_norm(crv1(x[0])-crv2(x[1]))};
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
        auto minf = gbs::solve_D_nlop(
            f,g,
            x, lb, hb,
            tol_x,
            solver
        );
        return std::array<T,3>{T(x[0]),T(x[1]),T(sqrt(minf))};

    }

    template <typename T, size_t dim>
    auto extrema_curve_curve(const Curve<T, dim> &crv1, const Curve<T, dim> &crv2,T tol_x, nlopt::algorithm solver=default_nlopt_algo)// -> extrema_CC_result<T>
    {
        auto u1 = make_range<T>(crv1.bounds()[0] , crv1.bounds()[1],30);
        auto u2 = make_range<T>(crv2.bounds()[0] , crv2.bounds()[1],30);
        auto f = [&](const auto &u1_,const auto &u2_)
        {
            return norm(crv1(u1_)-crv2(u2_));
        };
        auto [u10,u20]=approx_min_loc(u1,u2,f);
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

    template <typename T, size_t dim>
    auto extrema_surf_curve(const Surface<T, dim> &srf, const Curve<T, dim> &crv, T u_c0, T u_s0, T v_s0, T tol_x, nlopt::algorithm solver=default_nlopt_algo) //-> extrema_CS_result<T>
    {

        auto f = [&crv,&srf](const std::vector<double> &x)
        {
            return std::vector<double>{gbs::sq_norm(crv(x[0])-srf(x[1],x[2]))};
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
        auto minf = gbs::solve_D_nlop(
            f,g,
            x, lb, hb,
            tol_x,
            solver
        );
        return std::array<T,4>{x[1], x[2],x[0],sqrt(minf)};
    }

    template <typename T, size_t dim>
    auto extrema_surf_curve(const Surface<T, dim> &srf, const Curve<T, dim> &crv, T tol_x, nlopt::algorithm solver=default_nlopt_algo) 
    {
        auto u = make_range<T>(srf.bounds()[0] , srf.bounds()[1],30);
        auto v = make_range<T>(srf.bounds()[2] , srf.bounds()[3],30);
        auto uc= make_range<T>(crv.bounds()[0] , crv.bounds()[1],30);
        auto f = [&](const auto &u_,const auto &v_,const auto &uc_)
        {
            return norm(crv(uc_)-srf(u_,v_));
        };
        auto [u0,v0,uc0]=approx_min_loc(u,v,uc,f);
        return extrema_surf_curve(srf, crv, u0, v0, uc0, tol_x, solver);
    }

    template <typename T, size_t dim>
    auto solve_coordinate(const Curve<T, dim> &crv, size_t iCoord, T Xi, std::optional<T> u0 = std::nullopt, int digits = 20, size_t maxit = 20)
    {
        if(iCoord>=dim)
        {
            throw std::exception("Invalid size.");
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

} // namespace gbs