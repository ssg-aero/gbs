#pragma once
#include <nlopt.hpp>
#include <gbslib/bscurve.h>
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
    auto dev_from_points(const std::vector<std::array<T, dim>> &points, const BSCurve<T, dim> &crv) -> crv_dev_info<T>
    {
        auto d_avg = 0., d_max = 0., u_max = 0.;
        auto u0 = crv.knotsFlats().front();
        std::for_each(
            points.begin(),
            points.end(),
            [&](const auto &pnt) {
                auto res = gbs::extrema_PC(crv, pnt, u0, 1e-6);
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
} // namespace gbs