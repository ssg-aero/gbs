#pragma once
#include <nlopt.hpp>
#include <vector>
namespace gbs
{
    template <typename F>
    struct N_UserData
    {
        N_UserData(const F &f_eq) : f_eq_{f_eq} {}
        const F &f_eq_;
    };

    template <typename F1,typename F2>
    struct G_UserData
    {
        G_UserData(const F1 &f_eq, const F2 &f_grad) : f_eq_{f_eq}, f_grad_{f_grad} {}
        const F1 &f_eq_;
        const F2 &f_grad_;
    };

    template <typename F, typename T>
    auto solve_N_nlop(const F &f, std::vector<T> &x, const std::vector<T> &lb, const std::vector<T> &hb,T tol_x, const char *solver)
    {

        auto eq_ = [](const std::vector<double> &x, std::vector<double> &grad, void *user_data) {
            if (!grad.empty())
            {
                throw std::exception("Not implemented");
            }
            auto p_d = (N_UserData<decltype(f)> *)(user_data);
            return p_d->f_eq_(x);
        };

        N_UserData data(f);

        nlopt::opt opt(solver, x.size());

        opt.set_lower_bounds(lb);
        opt.set_upper_bounds(hb);
        opt.set_min_objective(eq_, &data);
        opt.set_xtol_rel(tol_x);

        T minf;

        opt.optimize(x, minf); 

        return minf;
    }
}
