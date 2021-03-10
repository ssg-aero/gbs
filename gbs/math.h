#pragma once
#include <vector>

namespace gbs
{
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

    template <typename T>
    auto binomial_law(size_t n, size_t k) -> T
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
    std::vector<T> make_range(T v1, T v2, size_t n, bool parallel = false)
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