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

    template <typename T, size_t dim>
    auto approx(const std::vector<bss_constrain<T, dim>> &Q, const std::vector<T> &k_flat_u, const std::vector<T> &k_flat_v, size_t p, size_t q) -> gbs::BSSurface<T, dim>
    {
        auto n_constrains = Q.size();
        auto n = k_flat_u.size() - p - 1;
        auto m = k_flat_v.size() - q - 1;
        auto n_poles = n * m;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> N(n_constrains, n_poles);
        for (int k = 0; k < n_constrains; k++) // Eigen::ColMajor is default
        {
            auto [u, v, x, du, dv] = Q[k];
            for (int j{}; j < m; j++)
            {
                for (int i{}; i < n; i++)
                {
                    auto l = i + n * j;
                    N(k, l) = gbs::basis_function(u, i, p, du, k_flat_u) *
                              gbs::basis_function(v, j, q, dv, k_flat_v);
                }
            }
        }

        auto N_inv = N.colPivHouseholderQr();

        VectorX<T> b(n_constrains);
        std::vector<std::array<T, dim>> poles(n_poles);
        for (int d = 0; d < dim; d++)
        {
            for (int i = 0; i < n_constrains; i++)
            {
                b(i) = std::get<2>(Q[i])[d];
            }

            auto x = N_inv.solve(b);

            for (int i = 0; i < n_poles; i++)
            {
                poles[i][d] = x(i);
            }
        }

        return gbs::BSSurface<T, dim>(poles, k_flat_u, k_flat_v, p, q);   
    }
}