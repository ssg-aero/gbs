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

    template <typename T>
    std::vector<T> make_range(T v1, T v2, size_t n, bool parallel = false)
    {
        if (n < 2)
        {
            throw std::length_error("2 points a required for a range");
        }
            
        std::vector<T> v(n);
        T step = ( v2 -v1 ) / T( n - 1);

        // T v2_check = v2-step;
        // std::generate(v.begin(),v.end(),[&,v_ = v1-step] () mutable 
        // { 
        //     return v_ >= v2_check ? v2 : v_+=step ; 
        // });
        // T v2_check = v2-step;
        //TODO check wich is faster
        v.front() = v1;
        if(parallel)
        {
            std::generate(
                std::execution::par,
                std::next(v.begin(), 1),
                std::next(v.end() - 1),
                [&, v_ = v1]() mutable {
                    return v_ += step;
                });
        }
        else
        {
            std::generate(
                std::next(v.begin(), 1),
                std::next(v.end() - 1),
                [&, v_ = v1]() mutable {
                    return v_ += step;
                });
        }
        
        v.back() = v2;

        return v;
    }

    template <typename T>
    std::vector<T> make_range(std::array<T,2> bounds, size_t n, bool parallel = false)
    {
        return make_range(bounds[0],bounds[1],n,parallel);
    }

} // namespace gbs