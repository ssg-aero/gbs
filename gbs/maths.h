#pragma once
#include <vector>
#include <gbs/vecop.h>

namespace gbs
{
    const auto pi   =     acos(-1.);
    const auto x2pi = 2 * acos(-1.);

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
    std::vector<T> make_range(T v1, T v2, size_t n)
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

} // namespace gbs