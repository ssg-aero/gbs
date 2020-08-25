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

    template <typename T, size_t dim>
    auto extrema_PC(const BSCurve<T, dim> &crv, const std::array<T, dim> &pnt, T u0,T tol_x,const char* solveur="LN_COBYLA") -> extrema_PC_result<T>
    {

        class UserData
        {
        public:
            UserData(const gbs::BSCurve<T, dim> &crv, const std::array<T, dim> &pnt) : c{crv}, p{pnt} {}
            gbs::BSCurve<T, dim> c;
            std::array<T, dim> p;
        };
        UserData data(crv, pnt);

        auto f = [](const std::vector<T> &x, std::vector<T> &grad, void *user_data) {
            auto p_d = (UserData *)(user_data);
            auto c_u = p_d->c.value(x[0]);
            if (!grad.empty())
            {
                auto dc_u = p_d->c.value(x[0], 1);
                grad[0] = 0;
                for (int i = 0; i < 3; i++)
                {
                    grad[0] += 2 * dc_u[i] * (c_u[i] - p_d->p[i]);
                }
            }
            return gbs::sq_norm(c_u - p_d->p);
        };

        nlopt::opt opt(solveur, 1);
        // nlopt::opt opt("LD_MMA", 1);
        // nlopt::opt opt("LN_COBYLA", 1);
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

    template <typename T, size_t dim>
    auto extrema_PC(const BSCurve<T, dim> &crv, const std::array<T, dim> &pnt,T tol_u)
    {
        auto u0 = 0.5 * (crv.knotsFlats().back() - crv.knotsFlats().front());
        return extrema_PC(crv, pnt, u0,tol_u);
    }
} // namespace gbs