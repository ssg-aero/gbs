#pragma once
#include <gbslib/vecop.h>

namespace gbs
{
    template <typename T, size_t dim>
    auto translate(std::array<T,dim> &x,const std::array<T,dim> &v) -> void
    {
        x+=v;
    }

// operate on non rational poles
    template <typename T, size_t dim>
    auto translate(std::array<T, dim + 1> &x, const std::array<T, dim> &v) -> void
    {
        std::transform(
            x.begin(), std::next(x.end(), -1),
            v.begin(), x.begin(),
            [&x](const auto x_i, const auto v_i) {
                return x_i + v_i * x.back();
            });
    }

    template <typename T, size_t dim>
    auto translated(const std::array<T,dim> &x,const std::array<T,dim> &v) -> std::array<T,dim>
    {
        return x+v;
    }

    template <typename T, size_t dim>
    auto translated(const std::array<T, dim + 1> &x, const std::array<T, dim> &v) -> std::array<T,dim>
    {
        std::array<T,dim+1> res{x};
        translate(res,v);
        return res;
    }

    // template <typename T, size_t dim>
    // auto scale(std::array<T,dim> &x,T s) -> void
    // {
    //     x*=v;
    // }

    // template <typename T, size_t dim>
    // auto scaled(const std::array<T,dim> &x,const std::array<T,dim> &v) -> std::array<T,dim>
    // {
    //     return x*s;
    // }
}