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

        // Separable least squares. The grid collocation matrix is the Kronecker
        // product Mv (x) Mu, so its normal equations factor as
        //   (MvᵀMv) (x) (MuᵀMu),  i.e.  Cu * P * Cv = Muᵀ Q Mv,
        // with Cu = MuᵀMu (n_poles_u²) and Cv = MvᵀMv (n_poles_v²) the two small 1-D
        // normal-equation matrices. We factor each once (SPD -> LDLT) and solve two
        // sequences of 1-D solves instead of one dense O((nu*nv)*(np_u*np_v)^2) QR on
        // the full Kronecker system (surface analog of the curve banded solve / #35).
        // Poles and Q are stored U-first: index i + n_poles_u * j.
        const Eigen::Index nau = Eigen::Index(n_params_u), nav = Eigen::Index(n_params_v);
        const Eigen::Index npu = Eigen::Index(n_poles_u), npv = Eigen::Index(n_poles_v);

        MatrixX<T> Mu(nau, npu);
        Mu.setZero();
        for (size_t iu = 0; iu < n_params_u; ++iu)
            fill_basis_row(Mu, Eigen::Index(iu), u[iu], k_flat_u, p, size_t{0});
        MatrixX<T> Mv(nav, npv);
        Mv.setZero();
        for (size_t iv = 0; iv < n_params_v; ++iv)
            fill_basis_row(Mv, Eigen::Index(iv), v[iv], k_flat_v, q, size_t{0});

        auto ldlt_u = (Mu.transpose() * Mu).ldlt();
        auto ldlt_v = (Mv.transpose() * Mv).ldlt();

        std::vector<std::array<T, dim>> poles(n_poles);
        MatrixX<T> Qd(nau, nav);
        for (int d = 0; d < int(dim); ++d)
        {
            for (size_t iv = 0; iv < n_params_v; ++iv)
                for (size_t iu = 0; iu < n_params_u; ++iu)
                    Qd(Eigen::Index(iu), Eigen::Index(iv)) = Q[iu + n_params_u * iv][d];

            const MatrixX<T> M = Mu.transpose() * Qd * Mv;     // Muᵀ Q Mv (n_poles_u x n_poles_v)
            const MatrixX<T> Y = ldlt_u.solve(M);              // Cu^{-1} M
            const MatrixX<T> P = ldlt_v.solve(Y.transpose()).transpose(); // Y Cv^{-1}

            for (size_t j = 0; j < n_poles_v; ++j)
                for (size_t i = 0; i < n_poles_u; ++i)
                    poles[i + n_poles_u * j][d] = P(Eigen::Index(i), Eigen::Index(j));
        }
        return poles;
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
        // Banded tensor-product assembly: a constraint row is the outer product of
        // the two 1-D basis (derivative) bands at (u, v). Only the (p+1)(q+1) columns
        // in the joint support are non-zero, so we evaluate each 1-D band once
        // (one A2.3 pass) and write only that block instead of a recursive
        // basis_function over all n*m columns of every row.
        N.setZero();
        for (size_t k = 0; k < n_constraints; k++) // Eigen::ColMajor is default
        {
            auto [u, v, x, du, dv] = Q[k];
            if (du > p || dv > q)
                continue; // derivative order above degree => row is identically zero

            const size_t span_u = find_span(n, p, u, k_flat_u) - k_flat_u.begin();
            const size_t imin = (span_u >= p) ? span_u - p : 0;
            const size_t icnt = (imin + p < n) ? p + 1 : n - imin;
            const size_t span_v = find_span(m, q, v, k_flat_v) - k_flat_v.begin();
            const size_t jmin = (span_v >= q) ? span_v - q : 0;
            const size_t jcnt = (jmin + q < m) ? q + 1 : m - jmin;

            const bool on_stack = (p <= bspline_stack_max_degree) && (q <= bspline_stack_max_degree);
            T Nu[bspline_stack_max_degree + 1];
            T Nv[bspline_stack_max_degree + 1];
            if (on_stack)
            {
                basis_ders(span_u, p, du, u, k_flat_u, Nu);
                basis_ders(span_v, q, dv, v, k_flat_v, Nv);
            }
            for (size_t jb = 0; jb < jcnt; ++jb)
            {
                size_t j = jmin + jb;
                T Nvj = on_stack ? Nv[jb] : basis_function(v, j, q, dv, k_flat_v);
                for (size_t ib = 0; ib < icnt; ++ib)
                {
                    size_t i = imin + ib;
                    T Nui = on_stack ? Nu[ib] : basis_function(u, i, p, du, k_flat_u);
                    N(Eigen::Index(k), Eigen::Index(i + n * j)) = Nui * Nvj;
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