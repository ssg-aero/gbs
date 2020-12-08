#pragma once
#include <gbs/vecop.h>
#include <gbs/bscurve.h>

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
    auto rotated(const std::array<T,2> &x,T a) -> std::array<T,2>
    {
        std::array<T,2> res{x};
        rotate(res,a);
        return res;
    }

    template <typename T>
    auto rotate(std::array<T, 3> &x, T a, const std::array<T, 3> &ax) -> void
    {
        std::array<T, 3> ax_ = ax;
        adim(ax_);
        auto c = cos(a);
        auto s = sin(a);
        
        auto x1 = x[0] * (ax_[0] * ax_[0] * (1 - c) + c) + x[1] * (ax_[0] * ax_[1] * (1 - c) + ax_[2] * s) + x[2] * (ax_[0] * ax_[2] * (1 - c) + ax_[1] * s);
        auto x2 = x[0] * (ax_[0] * ax_[1] * (1 - c) + ax_[2] * s) + x[1] * (ax_[1] * ax_[1] * (1 - c) + c) + x[2] * (ax_[1] * ax_[2] * (1 - c) - ax_[0] * s);
        auto x3 = x[0] * (ax_[0] * ax_[2] * (1 - c) - ax_[1] * s) + x[1] * (ax_[1] * ax_[2] * (1 - c) + ax_[0] * s) + x[2] * (ax_[1] * ax_[2] * (1 - c) + c);
        x = { x1, x2, x3};
    }

    template <typename T>
    auto rotate(BSCurve<T,2> &crv, T a) -> void
    {
        std::transform(
            crv.poles_begin(),
            crv.poles_end(),
            crv.poles_begin(),
            [&a](const auto &x_) { return rotated(x_,a);}
        );
    }

    template <typename T, size_t dim>
    auto translate(BSCurve<T,dim> &crv, const std::array<T,dim> &v) -> void
    {
        std::transform(
            crv.poles_begin(),
            crv.poles_end(),
            crv.poles_begin(),
            [&v](const auto &x_) { return translated(x_,v);}
        );
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