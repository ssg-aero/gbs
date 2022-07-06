#pragma once
#include <gbs/bssurf.h>
namespace gbs
{
    template <typename T, size_t dim>
    auto discretize(const Surface<T, dim> &srf, size_t nu, size_t nv) -> gbs::points_vector<T, dim>
    {
        points_vector<T, dim> points;
        auto [u1, u2, v1, v2] = srf.bounds();
        auto u = make_range(u1, u2, nu);
        auto v = make_range(v1, v2, nv);

        std::for_each(
            v.begin(),
            v.end(),
            [&srf, &u, &points, nu](const auto &v_) // mutable
            {
                points_vector<T, dim> points_(nu);
                std::transform(
                    std::execution::par,
                    u.begin(),
                    u.end(),
                    points_.begin(),
                    [&srf, v_](const auto &u_)
                    {
                        return srf(u_, v_);
                    });

                points.insert(
                    points.end(),
                    std::make_move_iterator(points_.begin()),
                    std::make_move_iterator(points_.end()));
            });

        return points;
    }
}