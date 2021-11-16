#pragma once
#include <Eigen/Dense>
#include <gbs/gbslib.h>
#include <gbs/bssinterp.h>

namespace gbs
{
    template <typename T, size_t dim>
    auto approx(const points_vector<T, dim> &Q, const std::vector<T> &k_flat_u, const std::vector<T> &k_flat_v, const std::vector<T> &u, const std::vector<T> &v, size_t p, size_t q) -> std::vector<std::array<T, dim>>
    {
        auto n_pt = Q.size();
        auto n_params_u = u.size();
        auto n_params_v = v.size();
        auto n_poles_u = k_flat_u.size() - (p+1); 
        auto n_poles_v = k_flat_v.size() - (q+1); 
        auto n_poles = n_poles_u * n_poles_v;
        if (n_poles > n_pt)
        {
            std::length_error("size error");
        }
        MatrixX<T> N(n_params_u*n_params_v, n_poles);
        fill_poles_matrix(N,k_flat_u,k_flat_v,u,v,p,q);
        return build_poles(N.colPivHouseholderQr(),Q);
    }
}