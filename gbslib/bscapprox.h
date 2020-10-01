#pragma once
#include <vector>
#include <Eigen/Dense>
#include <gbslib/bssinterp.h>

namespace gbs
{
    template <typename T>
    void removeRow(MatrixX<T> &matrix, unsigned int rowToRemove)
    {
        unsigned int numRows = matrix.rows() - 1;
        unsigned int numCols = matrix.cols();

        if (rowToRemove < numRows)
            matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.bottomRows(numRows - rowToRemove);

        matrix.conservativeResize(numRows, numCols);
    }

    template <typename T>
    void removeColumn(MatrixX<T> &matrix, unsigned int colToRemove)
    {
        unsigned int numRows = matrix.rows();
        unsigned int numCols = matrix.cols() - 1;

        if (colToRemove < numCols)
            matrix.block(0, colToRemove, numRows, numCols - colToRemove) = matrix.rightCols(numCols - colToRemove);

        matrix.conservativeResize(numRows, numCols);
    }
    /**
     * @brief Approximate a set of points with exact matchin at bounds
     * 
     * @tparam T 
     * @tparam dim 
     * @param pts     : point set to approximate
     * @param p       : curve's degree
     * @param n_poles : desired poles number
     * @param u       : point parameter on curve
     * @param k_flat  : curve's parametrization
     * @return gbs::BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx_bound_fixed(const std::vector<std::array<T, dim>> &pts, size_t p, size_t n_poles, const std::vector<T> &u, std::vector<double> k_flat) -> gbs::BSCurve<T, dim>
    {
        auto n_params = int(u.size());
        MatrixX<T> N(n_params-2, n_poles-2);


        for (int i = 0; i < n_params-2; i++)
        {
            for (int j = 0; j < n_poles-2; j++)
            {
                    N(i , j) = gbs::basis_function(u[i+1], j+1, p, 0, k_flat);
            }
        }

        auto p_begin = pts.front();
        auto p_end = pts.back();
        std::vector<T> Nbegin(n_params-2),Nend(n_params-2);
        for (int i = 0; i < n_params-2; i++)
        {
            Nbegin[i] = gbs::basis_function(u[i+1], 0, p, 0, k_flat);
            Nend[i] = gbs::basis_function(u[i+1], n_poles-1, p, 0, k_flat);
        }


        // auto N_inv = N.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
        auto N_inv = N.colPivHouseholderQr();

        auto n_pt = pts.size();

        std::vector<std::array<T, dim>> poles(n_poles);
        VectorX<T> b(n_pt-2);
        for (int d = 0; d < dim; d++)
        {
            for (int i = 0; i < n_pt-2; i++)
            {
                b(i) = pts[i+1][d] - p_begin[d] * Nbegin[i] - p_end[d] * Nend[i]; //Pas top au niveau de la localisation mémoire
            }

            auto x = N_inv.solve(b);

            for (int i = 1; i < n_poles - 1; i++)
            {
                poles[i][d] = x(i - 1);
            }
        }
        poles.front() = pts.front();
        poles.back() = pts.back();

        return BSCurve<T,dim>(poles, k_flat, p);
    }
    /**
     * @brief Approximate a point set
     * 
     * @tparam T 
     * @tparam dim 
     * @param pts : point set to approximate
     * @param p       : curve's degree
     * @param n_poles : desired poles number
     * @param u       : point parameter on curve
     * @param k_flat  : curve's parametrization
     * @return gbs::BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx(const std::vector<std::array<T, dim>> &pts, size_t p, size_t n_poles, const std::vector<T> &u, std::vector<double> k_flat) -> gbs::BSCurve<T, dim>
    {
        MatrixX<T> N(u.size(), n_poles);
        build_poles_matix<T, 1>(k_flat, u, p, n_poles, N);

        // auto N_inv = N.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
        auto N_inv = N.colPivHouseholderQr();

        auto n_pt = pts.size();

        std::vector<std::array<T, dim>> poles(n_poles);
        VectorX<T> b(n_pt);
        for (int d = 0; d < dim; d++)
        {
            for (int i = 0; i < n_pt; i++)
            {
                b(i) = pts[i][d]; //Pas top au niveau de la localisation mémoire
            }

            auto x = N_inv.solve(b);

            for (int i = 0; i < n_poles; i++)
            {
                poles[i][d] = x(i);
            }
        }

        return BSCurve<T,dim>(poles, k_flat, p);
    }
    /**
     * @brief Approximate a point set, the curve's parametrization is automaticaly computed
     * 
     * @tparam T 
     * @tparam dim 
     * @param pts     : point set to approximate
     * @param p       : curve's degree
     * @param n_poles : desired poles number
     * @param u       : point parameter on curve
     * @param fix_bound : force match on bounds
     * @return gbs::BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx(const std::vector<std::array<T, dim>> &pts, size_t p, size_t n_poles, const std::vector<T> &u, bool fix_bound) -> gbs::BSCurve<T, dim>
    {
        auto nk = n_poles + p + 1;
        std::vector<double> k_flat(nk);
        std::fill(k_flat.begin(), std::next(k_flat.begin(), p), 0.);
        std::fill(std::next(k_flat.begin(), nk - 1 - p), k_flat.end(), u.back() - u.front());
        for (int j = 1; j < n_poles - p; j++) // TODO use std algo
        {
            k_flat[j + p] = j / double(n_poles - p);
        }

        if (fix_bound)
        {
            return approx_bound_fixed(pts, p, n_poles, u, k_flat);
        }
        else
        {
            return approx(pts, p, n_poles, u, k_flat);
        }
    }
    /**
     * @brief Approximate a point set, the curve's poles' number and parametrization are automaticaly computed  
     * 
     * @tparam T 
     * @tparam dim 
     * @param pts : point set to approximate
     * @param p   : curve's degree
     * @param mode : curve parametrization mode
     * @param fix_bound : force match on bounds
     * @return gbs::BSCurve<T, dim> 
     */
    template <typename T, size_t dim>
    auto approx(const std::vector<std::array<T, dim>> &pts, size_t p, gbs::KnotsCalcMode mode, bool fix_bound) -> gbs::BSCurve<T, dim>
    {
        auto u = gbs::curve_parametrization(pts, mode, true);
        auto n_poles = p * 2;
        // auto n_poles = pts.size() / 5;
        auto crv = approx(pts, p, n_poles, u, fix_bound);

        for (int i = 0; i < 200; i++)
        {
            auto d_avg = 0., d_max = 0., u_max = -1.;
            std::vector<T> knots{crv.knotsFlats()};
            auto u0 = knots.front();
            std::for_each(
                // std::execution::par,
                pts.begin(),
                pts.end(),
                [&](const auto &pnt) {
                    auto res = gbs::extrema_PC(crv, pnt, u0, 1e-6);
                    u0 = res.u;
                    auto it = std::lower_bound(knots.begin(), knots.end(), res.u);
                    auto uh = *it;
                    if (it == knots.begin())
                        it = std::next(it);
                    auto ul = *(std::next(it, -1));
                    auto dul = (res.u - ul) / (uh - ul);
                    auto duh = (uh - res.u) / (uh - ul);
                    if (
                        res.d > d_max 
                        && dul > 0.33 
                        && duh > 0.33
                        )
                    {
                        d_max = res.d;
                        u_max = res.u;
                    }
                    d_avg += res.d;
                });

            auto it = std::lower_bound(knots.begin(), knots.end(), u_max);
            knots.insert(it, u_max);
            n_poles++;
            if (fix_bound)
            {
                crv = approx_bound_fixed(pts, p, n_poles, u, knots);
            }
            else
            {
                crv = approx(pts, p, n_poles, u, knots);
            }
            d_avg /= pts.size();
            std::cout << "d_avg: " << d_avg << ", d_max: " << d_max << ", u_max:" << u_max << std::endl;
            if (
                u_max < 0 
                || d_max < 1e-3
                || d_avg < 1e-4
                )
                break;
            std::cout << "n poles: " << crv.poles().size() << ", n_flat: " << crv.knotsFlats().size() << std::endl;
        }
        return crv;
    }
    /**
     * @brief  Approximate a point set, the curve's parametrization is automaticaly computed and the bounds are matched
     * 
     * @tparam T 
     * @tparam dim 
     * @param pts  : point set to approximate
     * @param p   : curve's degree
     * @param n_poles : desired poles' number
     * @param mode   : curve parametrization mode
     * @return auto 
     */
    template <typename T, size_t dim>
    auto approx(const std::vector<std::array<T, dim>> &pts, size_t p, size_t n_poles, gbs::KnotsCalcMode mode) // -> gbs::BSCurve<T,dim>
    {
        auto u = gbs::curve_parametrization(pts, mode, true);
        return approx(pts, p, n_poles, u, true);
    }

} // namespace gbs