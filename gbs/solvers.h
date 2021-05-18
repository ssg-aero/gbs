#pragma once
#include <nlopt.hpp>
#include <vector>
#include <numeric>
#include <gbs/vecop.h>
namespace gbs
{
    template <typename F>
    struct N_UserData
    {
        N_UserData(const F &f_eq) : f_eq_{f_eq} {}
        const F &f_eq_;
    };

    template <typename F, typename T>
    auto solve_N_nlop(const F &f, std::vector<T> &x, const std::vector<T> &lb, const std::vector<T> &hb,T tol_x, nlopt::algorithm solver)
    {

        auto eq_ = [](const std::vector<double> &x, std::vector<double> &grad, void *user_data) {
            if (!grad.empty())
            {
                throw std::exception("Not implemented");
            }
            auto p_d = (N_UserData<decltype(f)> *)(user_data);
            return double(p_d->f_eq_(x));
        };

        N_UserData data(f);

        nlopt::opt opt(solver, x.size());

        auto lb_d = vector_of_doubles(lb);
        auto hb_d = vector_of_doubles(hb);
        auto x_d = vector_of_doubles(x);

        opt.set_lower_bounds(lb_d);
        opt.set_upper_bounds(hb_d);
        opt.set_min_objective(eq_, &data);
        opt.set_xtol_rel(tol_x);

        double minf;

        opt.optimize(x_d, minf); 

        std::transform(x_d.begin(),x_d.end(),x.begin(),[](const auto &v){return static_cast<T>(v);});

        return T(minf);
    }

    template <typename F1,typename F2>
    struct D_UserData
    {
        D_UserData(const F1 &f_eq, const F2 &f_grad) : f_eq_{f_eq}, f_grad_{f_grad} {}
        const F1 &f_eq_;
        const F2 &f_grad_;
    };

    template <typename F1, typename F2, typename T>
    auto solve_D_nlop(const F1 &f, const F2 &g, std::vector<T> &x, const std::vector<T> &lb, const std::vector<T> &hb,T tol_x, nlopt::algorithm solver)
    {

        auto eq_ = [](const std::vector<double> &x, std::vector<double> &grad, void *user_data) {
            auto p_d = (D_UserData<decltype(f),decltype(g)> *)(user_data);
            auto r = p_d->f_eq_(x);
            if (!grad.empty())
            {
                grad = p_d->f_grad_(x,r); // passing f_eq_ evaluation can, I some occasions save some computation steps
            }
            std::transform(r.begin(),r.end(),r.begin(),[](const auto &r_){return r_*r_;});
            return std::reduce(r.begin(),r.end());
        };

        D_UserData data(f,g);

        nlopt::opt opt(solver, x.size());

        auto lb_d = vector_of_doubles(lb);
        auto hb_d = vector_of_doubles(hb);
        auto x_d = vector_of_doubles(x);

        opt.set_lower_bounds(lb_d);
        opt.set_upper_bounds(hb_d);
        opt.set_min_objective(eq_, &data);
        opt.set_xtol_rel(tol_x);

        double minf;

        opt.optimize(x_d, minf); 

        std::transform(x_d.begin(),x_d.end(),x.begin(),[](const auto &v){return static_cast<T>(v);});

        return T(minf);

        // opt.set_lower_bounds(lb);
        // opt.set_upper_bounds(hb);
        // opt.set_min_objective(eq_, &data);
        // opt.set_xtol_rel(tol_x);

        // T minf;

        // opt.optimize(x, minf); 

        // return minf;
    }
}
