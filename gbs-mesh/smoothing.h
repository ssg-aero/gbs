#pragma once
#include <gbs/gbslib.h>
#include <gbs/execution.h>

#include <algorithm>
#include <numeric>
#include <vector>

namespace gbs
{
    template <typename T>
    auto elliptic_structured_smoothing( points_vector<T,2> &pts, size_t nj, size_t i1, size_t i2, size_t j1, size_t j2, size_t n_it, T tol = 1e-4)
    {
        points_vector<T,2> pts_{pts};
        // Get sizes and steps
        size_t ni = pts.size() / nj;
        T d_ksi = 1 / ( ni - T(1) );
        T d_eth = 1 / ( nj - T(1) );
        // Mesht traversing function
        auto X = [nj, &pts](size_t i, size_t j, size_t d) -> const T&
        {
            return pts[j+nj*i][d];
        };
        // Jacobi sweep of one interior row i: reads only the previous buffer (pts),
        // writes only its own slots of pts_, returns the row's max correction.
        // Rows are therefore independent (#90): the sweep is race-free and
        // bit-identical whatever the execution order, and so is err_max (max is
        // exact in floating point, unlike a sum).
        auto sweep_row = [&](size_t i) -> T
        {
            T row_err{};
            for (size_t j{j1 + 1}; j < j2; j++)
            {
                auto x_ksi = (X(i + 1, j, 0) - X(i - 1, j, 0)) / (2 * d_ksi);
                auto y_ksi = (X(i + 1, j, 1) - X(i - 1, j, 1)) / (2 * d_ksi);
                auto x_eth = (X(i, j + 1, 0) - X(i, j - 1, 0)) / (2 * d_eth);
                auto y_eth = (X(i, j + 1, 1) - X(i, j - 1, 1)) / (2 * d_eth);

                T a = x_eth * x_eth + y_eth * y_eth;
                T b = x_ksi * x_eth + y_ksi * y_eth;
                T c = x_ksi * x_ksi + y_ksi * y_ksi;

                auto f = [&](size_t d)
                {
                    return
                    ( a / d_ksi / d_ksi * (X(i + 1, j, d) + X(i - 1, j, d))
                    + c / d_eth / d_eth * (X(i, j + 1, d) + X(i, j - 1, d))
                    - b / 2 / d_ksi / d_eth * (X(i + 1, j + 1, d) - X(i + 1, j - 1, d) + X(i - 1, j - 1, d) - X(i - 1, j + 1, d))
                    ) / 2 / (a / d_ksi / d_ksi + c / d_eth / d_eth);
                };

                pts_[j+nj*i][0] = f(0);
                pts_[j+nj*i][1] = f(1);
                auto dx = pts_[j+nj*i][0] - X(i,j,0);
                auto dy = pts_[j+nj*i][1] - X(i,j,1);

                row_err = std::max(row_err, std::abs(dx) + std::abs(dy));
            }
            return row_err;
        };
        // Interior rows, and the size gate: parallel only when the number of
        // interior vertices amortizes the task spin-up (cf. transform_threshold).
        std::vector<size_t> rows(i2 > i1 + 1 ? i2 - i1 - 1 : 0);
        std::iota(rows.begin(), rows.end(), i1 + 1);
        std::vector<T> rows_err(rows.size());
        const bool par = rows.size() * (j2 > j1 + 1 ? j2 - j1 - 1 : 0) >= parallel_min_size;
        // Start interation process
        T err_max{};
        size_t it {};
        do
        {
            if (par)
                std::transform(GBS_PAR_EXEC rows.begin(), rows.end(), rows_err.begin(), sweep_row);
            else
                std::transform(rows.begin(), rows.end(), rows_err.begin(), sweep_row);
            err_max = rows_err.empty() ? T{} : *std::max_element(rows_err.begin(), rows_err.end());
            it++;
            std::swap(pts_, pts);
        }while( (it < n_it) && ( err_max > tol ) );
        return std::make_pair(it,err_max);
    }

}
