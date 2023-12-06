module;

#include <numbers>

export module math;
export import vecop; // make range can use overloaded -

import <vector>;
import <algorithm>;
import <stdexcept>;


export namespace gbs
{

    template <typename T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

    template <typename L, typename T>
    auto kronecker(T i, T j) -> L
    {
        return i == j ? 1. : 0.;
    }

    template <typename T>
    auto radians(T angle_deg) -> T
    {
        return angle_deg * std::numbers::pi_v<T> / T{180};
    }

    template <typename T>
    auto degrees(T angle_rad) -> T
    {
        return angle_rad / std::numbers::pi_v<T> * T{180};
    }

    template <typename T>
    auto factorial(T n) -> T
    {
        if (n <= 1) // cover 0 and negatives
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
        T step = ( v2 - v1 ) / (n - 1.);
        v.front() = v1;
        if (n > 2)
        {
            std::generate(
                std::next(v.begin(), 1),
                std::next(v.end() - 1),
                [&, v_ = v1]() mutable
                {
                    return v_ += step;
                });
        }
        v.back() = v2;

        return v;
    }

    template <typename T>
    auto make_range(T u1, T u2, T v1, T v2, size_t nu , size_t nv) -> std::vector<std::pair<T,T>> 
    {
        if (nu < 2 || nv < 2)
        {
            throw std::length_error("2 points a required for a range");
        }
        
        T lu{ static_cast<T>( nu - 1 )};
        T lv{ static_cast<T>( nv - 1 )};
        std::vector<std::pair<T,T>> v(nu*nv);
        for(size_t i{}; i < nu; i++)
        {
            for(size_t j{}; j < nv; j++)
            {
                v[j + nv * i] = std::make_pair<T,T>( u1 + (u2-u1) * i / lu, v1 + (v2-v1) * j / lv  );
            }
        }

        return v;
    }

    template <typename T>
    auto make_range(T u1, T u2, T v1, T v2, T w1, T w2, size_t nu , size_t nv, size_t nw) -> std::vector<std::tuple<T,T,T>> 
    {
        if (nu < 2 || nv < 2 || nw < 2)
        {
            throw std::length_error("2 points a required for a range");
        }
        
        T lu{ static_cast<T>( nu - 1 )};
        T lv{ static_cast<T>( nv - 1 )};
        T lw{ static_cast<T>( nw - 1 )};
        std::vector<std::tuple<T,T,T>> v(nu*nv*nw);
        for(size_t i{}; i < nu; i++)
        {
            for(size_t j{}; j < nv; j++)
            {
                for(size_t k{}; k < nw; k++)
                {
                    v[k + nw * ( j + nv * i )] = std::tuple<T,T,T>( u1 + (u2-u1) * i / lu, v1 + (v2-v1) * j / lv, w1 + (w2-w1) * k / lw  );
                }
            }
        }

        return v;
    }

    template <typename T>
    auto make_range(std::array<T,4> bounds1, std::array<T,2> bounds2, size_t nu , size_t nv, size_t nw) ->  std::vector<std::tuple<T,T,T>> 
    {
        return make_range(bounds1[0], bounds1[1], bounds1[2], bounds1[3], bounds2[0], bounds2[1], nu, nv, nw);
    }

    template <typename T>
    auto make_range(std::array<T,4> bounds, size_t nu , size_t nv) -> std::vector<std::pair<T,T>>
    {
        return make_range(bounds[0], bounds[1], bounds[2], bounds[3], nu, nv);
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