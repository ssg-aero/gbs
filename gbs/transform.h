#pragma once

#include <execution>

#include <gbs/bscurve.h>
#include <gbs/bssurf.h>
#include <gbs/transformpoints.h>

namespace gbs
{
    template <typename T, size_t dim>
    auto transform(BSCurve<T, dim> &crv, auto trf) -> void
    {
        auto poles = crv.poles();
        std::transform(
            std::execution::par,
            poles.begin(),
            poles.end(),
            poles.begin(),
            trf);
        crv.movePoles(poles);
    }

    template <typename T, size_t dim>
    auto transform(BSCurveRational<T, dim> &crv, auto trf) -> void
    {
        auto poles_w = crv.poles();
        std::vector<std::array<T, dim>> poles;
        std::vector<T> weights;
        separate_weights(poles_w, poles, weights);

        std::transform(
            std::execution::par,
            poles.begin(),
            poles.end(),
            poles.begin(),
            trf);
        poles_w = add_weights_coord(poles, weights);
        crv.movePoles(poles_w);
    }

    template <typename T, size_t dim>
    auto transform(BSSurface<T, dim> &crv, auto trf) -> void
    {
        auto poles = crv.poles();
        std::transform(
            std::execution::par,
            poles.begin(),
            poles.end(),
            poles.begin(),
            trf);
        crv.movePoles(poles);
    }

    template <typename T, size_t dim>
    auto transform(BSSurfaceRational<T, dim> &crv, auto trf) -> void
    {
        auto poles_w = crv.poles();
        std::vector<std::array<T, dim>> poles;
        std::vector<T> weights;
        separate_weights(poles_w, poles, weights);

        std::transform(
            std::execution::par,
            poles.begin(),
            poles.end(),
            poles.begin(),
            trf);
        poles_w = add_weights_coord(poles, weights);
        crv.movePoles(poles_w);
    }

    template <typename T,bool rational>
    auto transform(BSCurveGeneral<T, 3,rational> &crv, const Matrix4<T> &M) -> void
    {
        auto trf = [&M](const auto &pt){return transform<T>(pt,M);};
        transform(crv,trf);
    }
    template <typename T,bool rational>
    auto transform(BSSurfaceGeneral<T, 3,rational> &crv, const Matrix4<T> &M) -> void
    {
        auto trf = [&M](const auto &pt){return transform<T>(pt,M);};
        transform(crv,trf);
    }

    // template <typename G>
    // auto transformed(const G &geo, auto trf)
    auto transformed(const auto &geo, auto trf)
    {
        auto res {geo};
        transform(res,trf);
        return res;
    }

    template <typename G,typename T>
    auto rotate(G &crv, T a) -> void
    {
        auto trf = [&a](const auto &x_) { return rotated(x_, a); };
        transform(crv, trf);
    }
    template <typename G,typename T>
    auto rotate(G &crv, T a, const std::array<T, 3> &ax) -> void
    {
        auto trf = [&a,&ax](const auto &x_) { return rotated(x_,a, ax); };
        transform(crv, trf);
    }
    // template <typename T>
    // auto rotate(BSCurve<T, 2> &crv, T a) -> void
    // {
    //     auto trf = [&a](const auto &x_) { return rotated(x_, a); };
    //     transform<T, 2>(crv, trf);
    // }

    // template <typename T>
    // auto rotate(BSCurveRational<T, 2> &crv, T a) -> void
    // {
    //     auto trf = [&a](const auto &x_) { return rotated(x_, a); };
    //     transform<T, 2>(crv, trf);
    // }



    // template <typename T>
    // auto rotate(Geom<T,2> &geo, T a, const std::array<T, 3> &ax) -> void
    // {
    //     auto trf = [&a, &ax](const auto &x_) { return rotated(x_, a, ax); };
    //     transform<T, 3>(geo, trf);
    // }

    // template <typename T, size_t dim>
    // auto translate(Geom<T,dim> &geo, const std::array<T, dim> &v) -> void
    // {
    //     auto trf = [&v](const auto &x_) { return translated(x_, v); };
    //     transform<T, dim>(geo, trf);
    // }

    // template <typename T, size_t dim>
    // auto scale(Geom<T,dim> &geo, T scale_factor, point<T, dim> center = point<T, dim>{}) -> void
    // {
    //     auto trf = [&](const auto &x_) { return scaled(x_, scale_factor, center); };
    //     transform<T, dim>(geo, trf);
    // }

    // template <typename T, size_t dim>
    // auto scale(Geom<T,dim> &geo, T scale_factor, size_t coord, point<T, dim> center = point<T, dim>{}) -> void
    // {
    //     auto trf = [&](const auto &x_) { return scaled(x_, scale_factor, coord, center); };
    //     transform<T, dim>(geo, trf);
    // }
/////////////////////////////////
    // template <typename T>
    // auto rotate(BSCurve<T, 3> &crv, T a, const std::array<T, 3> &ax) -> void
    // {
    //     auto trf = [&a, &ax](const auto &x_) { return rotated(x_, a, ax); };
    //     transform<T, 3>(crv, trf);
    // }

    // template <typename T>
    // auto rotate(BSCurveRational<T, 3> &crv, T a, const std::array<T, 3> &ax) -> void
    // {
    //     auto trf = [&a, &ax](const auto &x_) { return rotated(x_, a, ax); };
    //     transform<T, 3>(crv, trf);
    // }

    // template <typename T>
    // auto rotate(BSSurface<T, 3> &crv, T a, const std::array<T, 3> &ax) -> void
    // {
    //     auto trf = [&a, &ax](const auto &x_) { return rotated(x_, a, ax); };
    //     transform<T, 3>(crv, trf);
    // }

    // template <typename T>
    // auto rotate(BSSurfaceRational<T, 3> &crv, T a, const std::array<T, 3> &ax) -> void
    // {
    //     auto trf = [&a, &ax](const auto &x_) { return rotated(x_, a, ax); };
    //     transform<T, 3>(crv, trf);
    // }
/////////////////////////////////
    template <typename T, size_t dim>
    auto translate(BSCurve<T, dim> &crv, const std::array<T, dim> &v) -> void
    {
        auto trf = [&v](const auto &x_) { return translated(x_, v); };
        transform(crv, trf);
    }

    template <typename T, size_t dim>
    auto translate(BSCurveRational<T, dim> &crv, const std::array<T, dim> &v) -> void
    {
        auto trf = [&v](const auto &x_) { return translated(x_, v); };
        transform(crv, trf);
    }
/////////////////////////////////
    template <typename T, size_t dim>
    auto scale(BSCurve<T, dim> &crv, T scale_factor, point<T, dim> center = point<T, dim>{}) -> void
    {
        auto trf = [&](const auto &x_) { return scaled(x_, scale_factor, center); };
        transform(crv, trf);
    }

    template <typename T, size_t dim>
    auto scale(BSCurveRational<T, dim> &crv, T scale_factor, point<T, dim> center = point<T, dim>{}) -> void
    {
        auto trf = [&](const auto &x_) { return scaled(x_, scale_factor, center); };
        transform(crv, trf);
    }

    template <typename T, size_t dim>
    auto scale(BSCurve<T, dim> &crv, T scale_factor, size_t coord, point<T, dim> center = point<T, dim>{}) -> void
    {
        auto trf = [&](const auto &x_) { return scaled(x_, scale_factor, coord, center); };
        transform(crv, trf);
    }

    template <typename T, size_t dim>
    auto scale(BSCurveRational<T, dim> &crv, T scale_factor, size_t coord, point<T, dim> center = point<T, dim>{}) -> void
    {
        auto trf = [&](const auto &x_) { return scaled(x_, scale_factor, coord, center); };
        transform(crv, trf);
    }

    template <typename T, size_t dim>
    auto scale(BSSurface<T, dim> &srf, T scale_factor, point<T, dim> center = point<T, dim>{}) -> void
    {
        auto trf = [&](const auto &x_) { return scaled(x_, scale_factor, center); };
        transform(srf, trf);
    }

    
    template <typename T, size_t dim>
    auto scale(BSSurfaceRational<T, dim> &srf, T scale_factor, point<T, dim> center = point<T, dim>{}) -> void
    {
        auto trf = [&](const auto &x_) { return scaled(x_, scale_factor, center); };
        transform(srf, trf);
    }
}