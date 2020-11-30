#pragma once
#include <gbs/vecop.h>

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

    template <typename T>
    auto rotate(std::array<T,2> &x,T a) -> void
    {
        auto c = cos(a);
        auto s = sin(a);
        auto tmp = x[0];
        x[0]=x[0]*c-x[1]*s;
        x[1]=tmp*s+x[1]*c;
    }
    template <typename T>
    auto rotate(std::array<T, 3> &x, T a, const std::array<T, 3> &ax) -> void
    {
        adim(ax);
        auto c = cos(a);
        auto s = sin(a);
        auto tmp = x[0];
        x1 = x[0] * (ax[0] * ax[0] * (1 - c) + c) + x[1] * (ax[0] * ax[1] * (1 - c) + ax[2] * s) + x[2] * (ax[0] * ax[2] * (1 - c) + ax[1] * s);
        x2 = x[0] * (ax[0] * ax[1] * (1 - c) + ax[2] * s) + x[1] * (ax[1] * ax[1] * (1 - c) + c) + x[2] * (ax[1] * ax[2] * (1 - c) - ax[0] * s);
        x3 = x[0] * (ax[0] * ax[2] * (1 - c) - ax[1] * s) + x[1] * (ax[1] * ax[2] * (1 - c) + ax[0] * s) + x[2] * (ax[1] * ax[2] * (1 - c) + c);
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