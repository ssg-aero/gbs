#pragma once

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
} // namespace gbs