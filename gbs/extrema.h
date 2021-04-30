#pragma once
#include <gbs/solvers.h>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>
#include <gbs/curves>
#include <gbs/curveline.h>
#include <gbs/bssurf.h>

#include <optional>
// static const char* default_algo = "GN_ORIG_DIRECT";
static const nlopt::algorithm default_algo=nlopt::GN_ORIG_DIRECT;
namespace gbs
{

    template <typename T>
    struct extrema_PC_result
    {
        T u;
        T d;
    };
    template <typename T>
    struct extrema_CC_result
    {
        T u1;
        T u2;
        T d;
    };
    template <typename T>
    struct extrema_PS_result
    {
        T u;
        T v;
        T d;
    };
    template <typename T>
    struct extrema_CS_result
    {
        T u_s;
        T v_s;
        T u_c;
        T d;
    };
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
     * @return extrema_PC_result<T> 
     */
    template <typename T, size_t dim>
    auto extrema_PC(const Curve<T, dim> &crv, const std::array<T, dim> &pnt, T u0,T tol_x,nlopt::algorithm solver=default_algo) -> extrema_PC_result<T>
    {
        auto f = [&pnt,&crv](const std::vector<double> &x)
        {
            return std::vector<double>{gbs::sq_norm(crv(x[0]) - pnt)};
        };

        auto g = [&pnt,&crv](const std::vector<double> &x,const std::vector<double> &r)
        {
            auto dfdu = 2.*crv(x[0],1)*(crv(x[0]) - pnt);
            return std::vector<double>{2. * r[0] * dfdu};
        };

        std::vector<double> x{u0};
        std::vector<double> lb{crv.bounds()[0]};
        std::vector<double> hb{crv.bounds()[1]};
        auto minf = gbs::solve_D_nlop(
            f,g,
            x, lb, hb,
            tol_x,
            solver
            );
        return {T(x[0]),T(sqrt(minf))};
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
    auto extrema_PC(const Curve<T, dim> &crv, const std::array<T, dim> &pnt,T tol_u) -> extrema_PC_result<T>
    {
        auto u0 = T(0.5) * (crv.bounds()[1] - crv.bounds()[0]);
        return extrema_PC(crv, pnt, u0,tol_u);
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
     * @return extrema_PS_result<T> 
     */
    template <typename T, size_t dim>
    auto extrema_PS(const Surface<T, dim> &srf, const std::array<T, dim> &pnt, T u0, T v0, T tol_x, nlopt::algorithm solver=default_algo) -> extrema_PS_result<T>
    {
        auto f = [&pnt,&srf](const std::vector<double> &x)
        {
            return std::vector<double>{gbs::sq_norm(srf(x[0],x[1]) - pnt)};
        };

        auto g = [&pnt,&srf](const std::vector<double> &x,const std::vector<double> &r)
        {
            // return std::vector<double>{2.*crv(x[0],1)*(crv(x[0]) - pnt)};
            auto f_ = srf(x[0],x[1])  - pnt;
            auto dfdu = 2.*srf(x[0],x[1],1,0)* f_;
            auto dfdv = 2.*srf(x[0],x[1],0,1)* f_;
            return std::vector<double>{2. * r[0] * dfdu,2. * r[0] * dfdv};
        };

        auto [u1,u2,v1,v2] = srf.bounds();
        std::vector<double> x{u0,v0};
        std::vector<double> lb{u1,v1};
        std::vector<double> hb{u2,v2};
        auto minf = gbs::solve_D_nlop(
            f,g,
            x, lb, hb,
            tol_x,
            solver
        );
        return {x[0], x[1], sqrt(minf)};
    }
    template <typename T, size_t dim>
    auto extrema_PS(const Surface<T, dim> &srf, const std::array<T, dim> &pnt, T tol_x, nlopt::algorithm solver=default_algo) -> extrema_PS_result<T>
    {
        auto u0  = 0.5 * (srf.bounds()[1] - srf.bounds()[0]);
        auto v0  = 0.5 * (srf.bounds()[3] - srf.bounds()[2]);
        return extrema_PS(srf, pnt, u0, v0, tol_x, solver);
    }

    template <typename T, size_t dim>
    auto extrema_CC(const Curve<T, dim> &crv1, const Curve<T, dim> &crv2, T u10, T u20,T tol_x, nlopt::algorithm solver=default_algo) -> extrema_CC_result<T>
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
        std::vector<double> x{u10,u20};
        std::vector<double> lb{u1,v1};
        std::vector<double> hb{u2,v2};
        auto minf = gbs::solve_D_nlop(
            f,g,
            x, lb, hb,
            tol_x,
            solver
        );
        return {T(x[0]),T(x[1]),T(sqrt(minf))};

    }

    template <typename T, size_t dim>
    auto extrema_CC(const Curve<T, dim> &crv1, const Curve<T, dim> &crv2,T tol_x, nlopt::algorithm solver=default_algo) -> extrema_CC_result<T>
    {
        auto u10 = 0.5 * (crv1.bounds()[0]+crv1.bounds()[1]);
        auto u20 = 0.5 * (crv2.bounds()[0]+crv2.bounds()[1]);
        return extrema_CC<T,dim>(crv1,crv2,u10,u20,tol_x,solver);
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
    auto extrema_CC(const Line<T, 2> &crv1, const Line<T, 2> &crv2) // TODO  raise if // and use tuple for other functions
    {
        auto [P1, N1] = crv1.getAx();
        auto [P2, N2] = crv2.getAx();
        auto P21 = P2 - P1;
        auto d  = det(std::array<double,4>{ N1[0], -N2[0],  N1[1], -N2[1]});
        auto l1 = det(std::array<double,4>{P21[0], -N2[0], P21[1], -N2[1]}) / d;
        auto l2 = det(std::array<double,4>{ N1[0], P21[0],  N1[1], P21[1]}) / d;
        return std::make_tuple(l1,l2,0.,0.);
    }

    template <typename T, size_t dim>
    auto extrema_CS(const Surface<T, dim> &srf, const Curve<T, dim> &crv, T u_c0, T u_s0, T v_s0, T tol_x, nlopt::algorithm solver=default_algo) -> extrema_CS_result<T>
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
        std::vector<double> x{u_c0,u_s0,v_s0};
        std::vector<double> lb{uc1,us1, vs1};
        std::vector<double> hb{uc2,us2, vs2};
        auto minf = gbs::solve_D_nlop(
            f,g,
            x, lb, hb,
            tol_x,
            solver
        );
        return {x[0], x[1], x[2],sqrt(minf)};
    }

    template <typename T, size_t dim>
    auto extrema_CS(const Surface<T, dim> &srf, const Curve<T, dim> &crv, T tol_x, nlopt::algorithm solver=default_algo) -> extrema_CS_result<T>
    {
        auto u0  = 0.5 * (srf.bounds()[1] - srf.bounds()[0]);
        auto v0  = 0.5 * (srf.bounds()[3] - srf.bounds()[2]);
        auto u0c = 0.5 * (crv.bounds()[1] - crv.bounds()[0]);
        return extrema_CS(srf, crv, u0, v0, u0c, tol_x, solver);
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

        // T u0_ {u0 ? u0.value : 0.5*(u_min+u_max)};
        T u0_ = u0.value_or(0.5*(u_min+u_max));

        return newton_raphson_iterate(f_eq, u0_, u_min, u_max, digits, maxit);
    }

} // namespace gbs