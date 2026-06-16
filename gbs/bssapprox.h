#pragma once
#include <string>
#include <stdexcept>
#include <Eigen/Dense>
#include <gbs/gbslib.h>
#include <gbs/bssinterp.h>

namespace gbs
{
    /**
     * @brief Reject ill-posed surface approximation sizes, consistently for every
     * entry point (mirrors gbs::check_approx_sizes for curves).
     *
     * Each parametric direction needs at least @c p+1 / @c q+1 control points, and the
     * total pole count cannot exceed the number of constraints (least-squares problem).
     */
    inline void check_approx_sizes_surf(size_t n_poles_u, size_t p, size_t n_poles_v, size_t q, size_t n_constraints)
    {
        if (n_poles_u < p + 1 || n_poles_v < q + 1)
        {
            throw std::length_error(
                "approx: a (" + std::to_string(p) + "," + std::to_string(q) +
                ") B-spline surface needs at least (" + std::to_string(p + 1) + "," +
                std::to_string(q + 1) + ") poles (got (" + std::to_string(n_poles_u) + "," +
                std::to_string(n_poles_v) + "))");
        }
        if (n_poles_u * n_poles_v > n_constraints)
        {
            throw std::length_error(
                "approx: n_poles (" + std::to_string(n_poles_u * n_poles_v) +
                ") must not exceed the number of points/constraints (" +
                std::to_string(n_constraints) + ")");
        }
    }
    template <typename T, size_t dim>
    auto approx(const points_vector<T, dim> &Q, const std::vector<T> &k_flat_u, const std::vector<T> &k_flat_v, const std::vector<T> &u, const std::vector<T> &v, size_t p, size_t q) -> std::vector<std::array<T, dim>>
    {
        auto n_pt = Q.size();
        auto n_params_u = u.size();
        auto n_params_v = v.size();
        auto n_poles_u = k_flat_u.size() - (p+1);
        auto n_poles_v = k_flat_v.size() - (q+1);
        auto n_poles = n_poles_u * n_poles_v;
        check_approx_sizes_surf(n_poles_u, p, n_poles_v, q, n_pt);
        MatrixX<T> N(n_params_u*n_params_v, n_poles);
        fill_poles_matrix(N,k_flat_u,k_flat_v,u,v,p,q);
        return build_poles(N.colPivHouseholderQr(),Q);
    }

    template <typename T, size_t dim>
    auto approx(const std::vector<bss_constraint<T, dim>> &Q, const std::vector<T> &k_flat_u, const std::vector<T> &k_flat_v, size_t p, size_t q) -> BSSurface<T, dim>
    {
        auto n_constraints = Q.size();
        auto n = k_flat_u.size() - p - 1;
        auto m = k_flat_v.size() - q - 1;
        auto n_poles = n * m;
        check_approx_sizes_surf(n, p, m, q, n_constraints);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> N(n_constraints, n_poles);
        for (int k = 0; k < n_constraints; k++) // Eigen::ColMajor is default
        {
            auto [u, v, x, du, dv] = Q[k];
            for (int j{}; j < m; j++)
            {
                for (int i{}; i < n; i++)
                {
                    auto l = i + n * j;
                    N(k, l) = basis_function(u, i, p, du, k_flat_u) *
                              basis_function(v, j, q, dv, k_flat_v);
                }
            }
        }

        auto N_inv = N.colPivHouseholderQr();

        VectorX<T> b(n_constraints);
        std::vector<std::array<T, dim>> poles(n_poles);
        for (int d = 0; d < dim; d++)
        {
            for (int i = 0; i < n_constraints; i++)
            {
                b(i) = std::get<2>(Q[i])[d];
            }

            auto x = N_inv.solve(b);

            for (int i = 0; i < n_poles; i++)
            {
                poles[i][d] = x(i);
            }
        }

        return BSSurface<T, dim>(poles, k_flat_u, k_flat_v, p, q);   
    }
}