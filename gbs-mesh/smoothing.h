#include <gbs/gbslib.h>

namespace gbs
{
    template <typename T>
    auto elliptic_structutred_smoothing( points_vector<T,2> &pts, size_t nj, size_t i1, size_t i2, size_t j1, size_t j2, size_t n_it, T tol = 1e-4)
    {
        // auto j_span = make_range<size_t>(j1,j2);
        points_vector<T,2> pts_{pts};
        // Get sizes and steps
        size_t ni = pts.size() / nj;
        T d_ksi = 1 / ( ni - T(1) );
        T d_eth = 1 / ( nj - T(1) );
        // Mesht traversing function
        auto X = [nj, &pts](size_t i, size_t j, size_t d) -> T&
        {
            return pts[j+nj*i][d];
        };
        // smoothing step
        T a,b,c;
        auto f = [&a,&b,&c,&X,d_ksi, d_eth](size_t i, size_t j, size_t d)
        {
            return 
            ( a / d_ksi / d_ksi * (X(i + 1, j, d) + X(i - 1, j, d)) 
            + c / d_eth / d_eth * (X(i, j + 1, d) + X(i, j - 1, d)) 
            - b / 2 / d_ksi / d_eth * (X(i + 1, j + 1, d) - X(i + 1, j - 1, d) + X(i - 1, j - 1, d) - X(i - 1, j + 1, d)) 
            ) / 2 / (a / d_ksi / d_ksi + c / d_eth / d_eth);
        };
        // Start interation process
        T err_max{};
        size_t it {};
        do
        {
            err_max = T{};
            for (size_t i{i1 + 1}; i < i2; i++)
            {
                for (size_t j{j1 + 1}; j < j2; j++)
                // std::for_each(std::execution::par,j_span.begin(), j_span.end(), [&](size_t j)
                {

                    auto x_ksi = (X(i + 1, j, 0) - X(i - 1, j, 0)) / (2 * d_ksi);
                    auto y_ksi = (X(i + 1, j, 1) - X(i - 1, j, 1)) / (2 * d_ksi);
                    auto x_eth = (X(i, j + 1, 0) - X(i, j - 1, 0)) / (2 * d_eth);
                    auto y_eth = (X(i, j + 1, 1) - X(i, j - 1, 1)) / (2 * d_eth);

                    a = x_eth * x_eth + y_eth * y_eth;
                    b = x_ksi * x_eth + y_ksi * y_eth;
                    c = x_ksi * x_ksi + y_ksi * y_ksi;

                    pts_[j+nj*i][0] = f(i,j,0);
                    pts_[j+nj*i][1] = f(i,j,1);
                    auto dx = pts_[j+nj*i][0] - X(i,j,0);
                    auto dy = pts_[j+nj*i][1] - X(i,j,1);

                    auto err = std::abs(dx) + std::abs(dy);

                    err_max = std::max(err_max,err);
                    // printf("\t\terror: %.3e\n",err);        
                }
                // );
            }
            // printf("\titeration : %i, error max: %.3e\n", int(it), err_max);
            it++;
            std::swap(pts_, pts);
        }while( (it < n_it) && ( err_max > tol ) );
        return std::make_pair(it,err_max);
    }

}