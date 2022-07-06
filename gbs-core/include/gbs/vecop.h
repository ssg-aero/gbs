#pragma once
#include <array>
#include <algorithm>
#include <functional>
#include <execution>
#include <cmath>

namespace gbs
{
    const auto vecop_policy = std::execution::seq;

    template <typename T, size_t dim>
    std::array<T, dim> operator+(const std::array<T, dim> &a, const std::array<T, dim> &b)
    {
        std::array<T, dim> c;
        std::transform(
            vecop_policy,
            a.begin(), a.end(), b.begin(), c.begin(), 
            std::plus<T>()
            );
        return c;
    }

    template <typename T, size_t dim>
    std::array<T, dim> & operator+=(std::array<T, dim> &a, const std::array<T, dim> &b)
    {
        std::transform(
            vecop_policy,
            a.begin(), a.end(), b.begin(), a.begin(), 
            std::plus<T>()
            );
        return a;
    }

    template <typename T, size_t dim>
    std::array<T, dim> operator-(const std::array<T, dim> &a, const std::array<T, dim> &b)
    {
        std::array<T, dim> c;
        std::transform(
            vecop_policy,
            a.begin(), a.end(), b.begin(), c.begin(), 
            std::minus<T>()
            );
        return c;
    }

    template <typename T, size_t dim>
    std::array<T, dim> operator-(const std::array<T, dim> &a, T b)
    {
        std::array<T, dim> c;
        std::transform(
            vecop_policy,
            a.begin(), a.end(), c.begin(), 
            [&b](const auto &a_){return a_-b;}
            );
        return c;
    }

    template <typename T, size_t dim>
    std::array<T, dim> operator-(T a, const std::array<T, dim> &b)
    {
        std::array<T, dim> c;
        std::transform(
            vecop_policy,
            a.begin(), a.end(), c.begin(), 
            [&b](const auto &a_){return b-a_;}
            );
        return c;
    }

    template <typename T>
    std::array<T, 3> operator^(const std::array<T, 3> &a, const std::array<T, 3> &b)
    {
       return  {
           a[1] * b[2] - a[2] * b[1],
           a[2] * b[0] - a[0] * b[2],
           a[0] * b[1] - a[1] * b[0],
        };
    }

    template <typename T>
    std::array<T, 1> operator^(const std::array<T, 2> &a, const std::array<T, 2> &b)
    {
       return {a[0] * b[1] - a[1] * b[0]};
    }

    template <typename T>
    std::array<T, 1> operator^(const std::array<T, 1> &a, const std::array<T, 1> &b)
    {
       return {a[0] - b[0]};
    }

    template <typename T,size_t dim>
    T operator*(const std::array<T, dim> &a, const std::array<T, dim> &b)
    {
       return  std::transform_reduce(
           vecop_policy,
            a.begin(),a.end(),
            b.begin(),
            0.0
       ); //c++20
    }

    template <typename T, typename L, size_t dim>
    std::array<T, dim> operator*(const std::array<T, dim> &a, L b)
    {
        std::array<T, dim> c;
        std::transform(
            vecop_policy,
            a.begin(), a.end(),
            c.begin(),
            [&](const auto &ax) { return static_cast<T>(ax * b); });
        return c;
    }

    template <typename T, typename L, size_t dim>
    std::array<T, dim> operator*(L a, const std::array<T, dim> &b)
    {
        return b * a;
    }

    template <typename T, typename L,size_t dim>
    std::array<T, dim>operator/(const std::array<T, dim> &a, L b)
    {
        return a * (L{1}/b);
    }

    template <typename T,size_t dim>
    // std::vector< std::array<T, dim>> coef_div(const std::vector< std::array<T, dim>> &v1, const std::vector<T> &v2)
    std::vector< std::array<T, dim>>  operator/(const std::vector< std::array<T, dim>> &a, const std::vector<T> &b)
    {
        std::vector<std::array<T, dim>> d(a.size());
        std::transform(
            vecop_policy,
            a.begin(), a.end(),
            b.begin(),
            d.begin(),
            [](const auto &a_, const auto &b_) { return a_ / b_; });
        return d;
    }

    template <typename T, size_t dim>
    T sq_norm(const std::array<T, dim> &c)
    {
        return std::transform_reduce(
                       vecop_policy,
                       c.begin(), c.end(), 
                       c.begin(), 
                       0.0); // c++20
    }

    template <typename T, size_t dim>
    T norm(const std::array<T, dim> &c)
    {
        return std::sqrt( sq_norm(c) );
    }

    template <typename T, size_t dim>
    T sq_distance(const std::array<T, dim> &a, const std::array<T, dim> &b)
    {
        auto c = a - b;
        return sq_norm(c);
    }

    template <typename T, size_t dim>
    T distance(const std::array<T, dim> &a, const std::array<T, dim> &b)
    {
        return std::sqrt(sq_distance(a, b));
    }

    template <typename T, size_t dim>
    T length(const std::vector< std::array<T, dim> > &pts)
    {
        return std::transform_reduce(
            vecop_policy,
            ++pts.begin(),pts.end(),pts.begin(),T(0.),
            std::plus<T>(),
            [](const auto &pt1,const auto &pt2){return distance(pt1,pt2);}
            );
    }

    template <typename T>
    std::vector<T> delta(const std::vector<T> &v)
    {
        std::vector<T> d_v(v.size() - 1);
        std::transform(
            vecop_policy,
            v.begin(), std::next(v.end(), -1),
            std::next(v.begin()),
            d_v.begin(),
            [](const auto &u1, const auto &u2) {
                return u2 - u1;
            });
        return d_v;
    }

    template <typename T>
    std::vector<T> delta(const std::vector<T> &v1,const std::vector<T> &v2)
    {
        std::vector<T> d_v(v1.size());
        std::transform(
            vecop_policy,
            v1.begin(), v1.end(),
            v2.begin(),
            d_v.begin(),
            [](const auto &u1, const auto &u2) {
                return u2 - u1;
            });
        return d_v;
    }

    template <typename T>
    void adim(std::vector<T> &v)
    {
        std::transform(
            vecop_policy,
            v.begin(), std::next(v.end(), -1),
            v.begin(),
            [](const auto &v_) {
                return v_ / norm(v_);
            });
    }

    template <typename T,size_t dim>
    void adim(std::array<T,dim> &v)
    {
        auto n = norm(v);
        std::transform(
            vecop_policy,
            v.begin(),v.end(),
            v.begin(),
            [&n](const auto &x_) {
                return x_ / n;
            });
    }

    // r_k = (1-b_k) * v_k + b_k * v_k+1
    template <typename T,typename L,typename It_v,typename It_b,typename It_r>
    auto blend(const std::vector<T> &v,const std::vector<L> &b,std::vector<T> &r,It_v v_start,It_v v_end,It_b b_start,It_r r_start) -> void
    {
        std::vector<T> v1(v_end-v_start+1);
        std::transform(
            vecop_policy,
            v_start,v_end,
            b_start,
            v1.begin(),
            [](const auto &v_,const auto &b_){return (1-b_)*v_;}
        );

        std::vector<T> v2(v_end-v_start+1);
        std::transform(
            vecop_policy,
            std::next(v_start),std::next(v_end),
            b_start,
            v2.begin(),
            [](const auto &v_,const auto &b_){return b_*v_;}
        );

        std::transform(
            vecop_policy,
            v1.begin(),v1.end(),
            v2.begin(),
            r_start,
            [](const auto &v1_,const auto &v2_){return v1_ +v2_;}
        );
    }

    template<typename T>
    auto vector_of_doubles(const T &l){
        std::vector<double> v(l.size());
        // std::copy(std::begin(l), std::end(l), std::back_inserter(v));
        std::transform(l.begin(),l.end(),v.begin(),[](const auto &v_){return static_cast<double>(v_);});
        return v; 
    }
    // auto move_to_vector_of_doubles(const auto &v){return std::vector<double>{ std::make_move_iterator(std::begin(v)), std::make_move_iterator(std::end(v)) }; }


} // namespace gbs