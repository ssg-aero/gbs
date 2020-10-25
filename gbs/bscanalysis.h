#pragma once
#include <nlopt.hpp>
#include <gbs/bscurve.h>
#include <gbs/extrema.h>
namespace gbs
{
    template <typename T>
    struct crv_dev_info
    {
        T u_max;
        T d_max;
        T d_avg;
    };
    /**
     * @brief Computes the deviation from a set of points
     * 
     * @tparam T 
     * @tparam dim 
     * @param points 
     * @param crv 
     * @return crv_dev_info<T> 
     */
    template <typename T, size_t dim>
    auto DeviationFromPoints(const std::vector<std::array<T, dim>> &points, const BSCurve<T, dim> &crv) -> crv_dev_info<T>
    {
        auto d_avg = 0., d_max = 0., u_max = 0.;
        auto u0 = crv.knotsFlats().front();
        std::for_each(
            points.begin(),
            points.end(),
            [&](const auto &pnt) {
                auto res = extrema_PC(crv, pnt, u0, 1e-6);
                u0 = res.u;
                if (res.d > d_max)
                {
                    d_max = res.d;
                    u_max = res.u;
                }
                d_avg += res.d;
            });
        d_avg /= points.size();
        return {u_max, d_max, d_avg};
    }
    /**
     * @brief Very basic discretization
     * 
     * @tparam T 
     * @tparam dim 
     * @param crv 
     * @param n 
     * @return PointArray<T,dim> 
     */
    template <typename T, size_t dim,bool rational>
    auto Discretize(const BSCurveGeneral<T,dim,rational> &crv, size_t n) -> gbs::points_vector<T,dim>
    {
        points_vector<T,dim> points(n);
        auto u1 = crv.knotsFlats().front();
        auto u2 = crv.knotsFlats().back();
        auto du = (u2-u1) / (n-1);
        std::generate(points.begin(),points.end(),[&,u=u1-du]() mutable {return crv.value(u+=du);});
        return points;
    }

} // namespace gbs