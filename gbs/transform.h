#pragma once
#include <gbs/vecop.h>
#include <gbs/bscurve.h>
#include <gbs/bssurf.h>

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
        
        auto x1 = x[0] * (ax_[0] * ax_[0] * (1 - c) + c) + x[1] * (ax_[0] * ax_[1] * (1 - c) - ax_[2] * s) + x[2] * (ax_[0] * ax_[2] * (1 - c) + ax_[1] * s);
        auto x2 = x[0] * (ax_[0] * ax_[1] * (1 - c) + ax_[2] * s) + x[1] * (ax_[1] * ax_[1] * (1 - c) + c) + x[2] * (ax_[1] * ax_[2] * (1 - c) - ax_[0] * s);
        auto x3 = x[0] * (ax_[0] * ax_[2] * (1 - c) - ax_[1] * s) + x[1] * (ax_[1] * ax_[2] * (1 - c) + ax_[0] * s) + x[2] * (ax_[2] * ax_[2] * (1 - c) + c);
        x = { x1, x2, x3};
    }

    template <typename T>
    auto rotated(const point<T, 3> &x, T a, const std::array<T, 3> &ax) -> std::array<T, 3>
    {
        std::array<T,3> res{x};
        rotate(res,a,ax);
        return res;
    }

    template <typename T, size_t dim>
    auto scaled(const point<T, dim> &x, T scale_factor, point<T, dim> center = point<T, dim>{}) -> point<T, dim>
    {
        return (x-center) * scale_factor + center;
    }

    template <typename T, size_t dim>
    auto scale(point<T, dim> &x, T scale_factor, point<T, dim> center = point<T, dim>{}) -> void
    {
        x = scaled(x,scale_factor,center);
    }

    template <typename T, size_t dim>
    auto scale(point<T, dim> &x, T scale_factor,size_t coord, point<T, dim> center = point<T, dim>{})
    {
        x[coord] =  (x[coord]-center[coord]) * scale_factor + center[coord];
    }

    template <typename T, size_t dim>
    auto scaled(const point<T, dim> &x, T scale_factor,size_t coord, point<T, dim> center = point<T, dim>{}) -> point<T, dim>
    {
        std::array<T,dim> res{x};
        scale(res,scale_factor,coord,center);
        return res;
    }

    template <typename T, size_t dim>
    auto transform(BSCurve<T,dim> &crv, auto trf) -> void
    {
        auto poles = crv.poles();
        std::transform(
            std::execution::par,
            poles.begin(),
            poles.end(),
            poles.begin(),
            trf
        );
        crv.movePoles(poles);
    }

    template <typename T, size_t dim>
    auto transform(BSCurveRational<T,dim> &crv, auto trf) -> void
    {
        auto poles_w = crv.poles();
        std::vector<std::array<T, dim>> poles;
        std::vector<T> weights;
        separate_weights(poles_w,poles,weights);
        
        std::transform(
            std::execution::par,
            poles.begin(),
            poles.end(),
            poles.begin(),
            trf
        );
        poles_w = add_weights_coord(poles,weights);
        crv.movePoles(poles_w);
    }


    template <typename T, size_t dim>
    auto transform(BSSurface<T,dim> &crv, auto trf) -> void
    {
        auto poles = crv.poles();
        std::transform(
            std::execution::par,
            poles.begin(),
            poles.end(),
            poles.begin(),
            trf
        );
        crv.movePoles(poles);
    }

    template <typename T, size_t dim>
    auto transform(BSSurfaceRational<T,dim> &crv, auto trf) -> void
    {
        auto poles_w = crv.poles();
        std::vector<std::array<T, dim>> poles;
        std::vector<T> weights;
        separate_weights(poles_w,poles,weights);
        
        std::transform(
            std::execution::par,
            poles.begin(),
            poles.end(),
            poles.begin(),
            trf
        );
        poles_w = add_weights_coord(poles,weights);
        crv.movePoles(poles_w);
    }

    template <typename T>
    auto rotate(BSCurve<T,2> &crv, T a) -> void
    {
        auto trf = [&a](const auto &x_) { return rotated(x_,a);};
        transform<T,2>(crv,trf);
    }

    template <typename T>
    auto rotate(BSCurveRational<T,2> &crv, T a) -> void
    {
        auto trf = [&a](const auto &x_) { return rotated(x_,a);};
        transform<T,2>(crv,trf);
    }

    template <typename T>
    auto rotate(BSCurve<T,3> &crv, T a, const std::array<T, 3> &ax) -> void
    {
        auto trf = [&a,&ax](const auto &x_) { return rotated(x_,a,ax);};
        transform<T,3>(crv,trf);
    }

    template <typename T>
    auto rotate(BSCurveRational<T,3> &crv, T a, const std::array<T, 3> &ax) -> void
    {
        auto trf = [&a,&ax](const auto &x_) { return rotated(x_,a,ax);};
        transform<T,3>(crv,trf);
    }

    template <typename T>
    auto rotate(BSSurface<T,3> &crv, T a, const std::array<T, 3> &ax) -> void
    {
        auto trf = [&a,&ax](const auto &x_) { return rotated(x_,a,ax);};
        transform<T,3>(crv,trf);
    }

    template <typename T>
    auto rotate(BSSurfaceRational<T,3> &crv, T a, const std::array<T, 3> &ax) -> void
    {
        auto trf = [&a,&ax](const auto &x_) { return rotated(x_,a,ax);};
        transform<T,3>(crv,trf);
    }

    template <typename T, size_t dim>
    auto translate(BSCurve<T, dim> &crv, const std::array<T, dim> &v) -> void
    {
        auto trf = [&v](const auto &x_) { return translated(x_, v); };
        transform<T,dim>(crv, trf);
    }

    template <typename T, size_t dim>
    auto translate(BSCurveRational<T, dim> &crv, const std::array<T, dim> &v) -> void
    {
        auto trf = [&v](const auto &x_) { return translated(x_, v); };
        transform<T,dim>(crv, trf);
    }

    template <typename T, size_t dim>
    auto scale(BSCurve<T, dim> &crv, T scale_factor, point<T, dim> center = point<T, dim>{}) -> void
    {
        auto trf = [&](const auto &x_) { return scaled(x_, scale_factor,center); };
        transform<T,dim>(crv, trf);
    }

    template <typename T, size_t dim>
    auto scale(BSCurveRational<T, dim> &crv, T scale_factor, point<T, dim> center = point<T, dim>{}) -> void
    {
        auto trf = [&](const auto &x_) { return scaled(x_, scale_factor,center); };
        transform<T,dim>(crv, trf);
    }

    template <typename T, size_t dim>
    auto scale(BSCurve<T, dim> &crv, T scale_factor,size_t coord, point<T, dim> center = point<T, dim>{}) -> void
    {
        auto trf = [&](const auto &x_) { return scaled(x_, scale_factor,coord,center); };
        transform<T,dim>(crv, trf);
    }

    template <typename T, size_t dim>
    auto scale(BSCurveRational<T, dim> &crv, T scale_factor,size_t coord, point<T, dim> center = point<T, dim>{}) -> void
    {
        auto trf = [&](const auto &x_) { return scaled(x_, scale_factor,coord,center); };
        transform<T,dim>(crv, trf);
    }

}