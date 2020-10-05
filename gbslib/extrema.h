#pragma once
#include <nlopt.hpp>
#include <gbslib/bscurve.h>
namespace gbs
{
    template <typename T>
    struct extrema_PC_result
    {
        T u;
        T d;
    };
    /**
     * @brief Project point on curve
     * 
     * @tparam T 
     * @tparam dim 
     * @param crv     : the curve
     * @param pnt     : the point
     * @param u0      : guess value
     * @param tol_x   : tolerance
     * @param solver : solver type please have look to https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#nomenclature to change this value
     * @return extrema_PC_result<T> 
     */
    template <typename T, size_t dim,bool rational>
    auto extrema_PC(const BSCurveGeneral<T, dim,rational> &crv, const std::array<T, dim> &pnt, T u0,T tol_x,const char* solver="LN_COBYLA") -> extrema_PC_result<T>
    {

        class UserData
        {
        public:
            UserData(const gbs::BSCurveGeneral<T, dim, rational> &crv, const std::array<T, dim> &pnt) : p{pnt}
            {
                // c = const_cast<gbs::BSCurveGeneral<T, dim, rational> *>(&crv);
                c = &crv;
            }
            const gbs::BSCurveGeneral<T, dim,rational> *c;
            std::array<T, dim> p;
        };
        UserData data(crv, pnt);

        auto f = [](const std::vector<T> &x, std::vector<T> &grad, void *user_data) {
            auto p_d = (UserData *)(user_data);
            auto c_u = p_d->c->value(x[0]);
            if (!grad.empty())
            {
                auto dc_u = p_d->c->value(x[0], 1);
                grad[0] = 0;
                for (int i = 0; i < 3; i++)
                {
                    grad[0] += 2 * dc_u[i] * (c_u[i] - p_d->p[i]);
                }
            }
            return gbs::sq_norm(c_u - p_d->p);
        };

        nlopt::opt opt(solver, 1);
        std::vector<T> lb(1), hb(1);
        lb[0] = crv.knotsFlats().front();
        hb[0] = crv.knotsFlats().back();

        opt.set_lower_bounds(lb);
        opt.set_upper_bounds(hb);
        opt.set_min_objective(f, &data);
        opt.set_xtol_rel(tol_x);
        std::vector<T> x(1);
        x[0] = u0;
        T minf;

        opt.optimize(x, minf); //can raise

        return {x[0],sqrt(minf)};

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
    template <typename T, size_t dim,bool rational>
    auto extrema_PC(const BSCurveGeneral<T, dim,rational> &crv, const std::array<T, dim> &pnt,T tol_u)
    {
        auto u0 = 0.5 * (crv.knotsFlats().back() - crv.knotsFlats().front());
        return extrema_PC(crv, pnt, u0,tol_u);
    }
} // namespace gbs