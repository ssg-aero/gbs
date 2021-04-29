#pragma once
#include <vector>
#include <gbs/vecop.h>

namespace gbs
{
    const auto pi   =     acos(-1.);
    const auto x2pi = 2 * acos(-1.);
    template <typename L, typename T>
    auto kronecker(T i, T j) -> L
    {
        return i == j ? 1. : 0.;
    }

    template <typename T>
    auto radians(T angle_deg) -> T
    {
        return angle_deg * x2pi / T{360};
    }

    template <typename T>
    auto degrees(T angle_rad) -> T
    {
        return angle_rad / x2pi * T{360};
    }

    template <typename T>
    auto factorial(T n) -> T
    {
        if (n == 1)
        {
            return 1;
        }
        else
        {
            return n * factorial(n - 1);
        }
    }

    template <typename T,typename L>
    auto binomial_law(L n, L k) -> T
    {
        if (k == n)
        {
            return T(1);
        }
        else
        {
            return T(factorial(n)) / T(factorial(k) * factorial(n - k));
        }
    }

    template <typename T>
    auto make_range(T v1, T v2, size_t n) -> std::vector<T>
    {
        if (n < 2)
        {
            throw std::length_error("2 points a required for a range");
        }
            
        std::vector<T> v(n);
        T step = ( v2 -v1 ) / (n - 1.);
        v.front() = v1;
        std::generate(
            std::next(v.begin(), 1),
            std::next(v.end() - 1),
            [&, v_ = v1]() mutable {
                return v_ += step;
            });
        v.back() = v2;

        return v;
    }

    template<typename T>
    auto make_range(T v1, T v2) -> std::vector<T> 
    {
        auto n = v2-v1 +1;
        if ( n < 2)
        {
            throw std::length_error("2 points a required for a range");
        }
            
        std::vector<T> v(n);
        std::generate(
            v.begin(),
            v.end(),
            [&, v_ = v1]() mutable {
                return v_ ++;
            });

        return v;
    }

    template<typename T>
    auto make_range(const std::array<T,2> &bounds, size_t n ) -> std::vector<T>
    {
        return make_range(bounds[0],bounds[1],n);
    }

    template <typename T>
    auto newton_solve = [](const auto &func,auto p, T u0, T tol_f = 1.e-3, T tol_u = 1.e-4, size_t it_max=100) -> T
    {
        auto delta = tol_u * 10., d0 = tol_f * 10.;
        auto u = u0;
        auto count =0;
        while (delta>tol_u && d0 > tol_f && count < it_max)
        {
            auto d0 = func(u)-p;
            auto d1 = func(u,1);
            auto d2 = func(u,2);
            delta = d1*d0 / (d2*d0+d1*d1);
            u -= delta;
            count++;
        }
        return u;
    };

    /**
     * @brief matrix 2x2 determinant
     *   
     *   | M[0] M[1] |
     *   | M[2] M[3] |
     * 
     * @tparam T 
     * @param M 
     * @return auto 
     */
    template <typename T>
    auto det(const std::array<T,4> &M) -> T
    {
        return M[0]*M[3] - M[1]*M[2];
    }


} // namespace gbs