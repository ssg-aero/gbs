#pragma once
#include <vector>
#include <Eigen/Dense>
#include <gbslib/bssinterp.h>

namespace gbs
{
    template <typename T, size_t dim>
    auto approx(const std::vector<std::array<T, dim>> &pts, size_t p, size_t n_poles, const std::vector<T> &u) -> gbs::BSCurve<T,dim>
    {
        // auto k_flat = build_simple_mult_flat_knots<T>(u, n_poles, p);
        auto nk = n_poles + p + 1;
        std::vector<double> k_flat(nk);
        std::fill(k_flat.begin(), std::next(k_flat.begin(), p), 0.);
        std::fill(std::next(k_flat.begin(), nk - 1 - p), k_flat.end(), u.back() - u.front());
        for (int j = 1; j < n_poles - p; j++) // TODO use std algo
        {
            k_flat[j + p] = j / double(n_poles - p);
        }

        Eigen::MatrixX<T> N(u.size(), n_poles);
        build_poles_matix<T, 1>(k_flat, u, p, n_poles, N);

        auto N_inv = N.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);

        auto n_pt = pts.size();

        std::vector<std::array<T, dim>> poles(n_poles);
        Eigen::VectorX<T> b(n_pt);
        for (int d = 0; d < dim; d++)
        {
            for (int i = 0; i < n_pt; i++)
            {
                    b( i ) = pts[i][d]; //Pas top au niveau de la localisation mémoire
            }

            auto x = N_inv.solve(b);

            for (int i = 0; i < n_poles; i++)
            {
                poles[i][d] = x(i);
            }
        }

        return BSCurve(poles,k_flat,p);
    }

    template <typename T, size_t dim>
    auto approx(const std::vector<std::array<T, dim>> &pts, size_t p,const std::vector<T> &u) // -> gbs::BSCurve<T,dim>
    {
        // auto k_flat = build_simple_mult_flat_knots<T>(u,p,p); // init avec mini

    }

    template <typename T, size_t dim>
    auto approx(const std::vector<std::array<T, dim>> &pts, size_t p, gbs::KnotsCalcMode mode ) // -> gbs::BSCurve<T,dim>
    {
        auto u = gbs::curve_parametrization(pts, gbs::KnotsCalcMode::CHORD_LENGTH, true);
        return approx(pts,p,u);
    }

}