#pragma once
// Analytic and NURBS test geometry shared by the gbs-brep face tests.

#include <gbs/curves>
#include <gbs/surfaces>
#include <gbs/bscbuild.h>

#include <cmath>
#include <functional>
#include <memory>
#include <numbers>
#include <stdexcept>

namespace brep_tests
{
    using namespace gbs;
    using T = double;
    constexpr T pi = std::numbers::pi_v<T>;

    // Analytic surface S(u,v) given by a lambda; first derivatives by central differences.
    class Analytic : public Surface<T, 3>
    {
        std::function<point<T, 3>(T, T)> f_;
        std::array<T, 4> b_;

    public:
        Analytic(std::function<point<T, 3>(T, T)> f, std::array<T, 4> b) : f_{std::move(f)}, b_{b} {}
        auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, 3> override
        {
            const T h = 1e-6;
            if (du == 0 && dv == 0)
                return f_(u, v);
            if (du == 1 && dv == 0)
                return (f_(u + h, v) - f_(u - h, v)) * (T(0.5) / h);
            if (du == 0 && dv == 1)
                return (f_(u, v + h) - f_(u, v - h)) * (T(0.5) / h);
            throw std::runtime_error("not implemented");
        }
        auto bounds() const -> std::array<T, 4> override { return b_; }
    };

    std::shared_ptr<Surface<T, 3>> cylinder(T r = 1., T h = 2.)
    {
        return std::make_shared<Analytic>([r](T u, T v) { return point<T, 3>{r * std::cos(u), r * std::sin(u), v}; },
                                          std::array<T, 4>{0., 2. * pi, 0., h});
    }
    std::shared_ptr<Surface<T, 3>> sphere(T r = 1.)
    {
        return std::make_shared<Analytic>(
            [r](T u, T v) { return point<T, 3>{r * std::cos(v) * std::cos(u), r * std::cos(v) * std::sin(u), r * std::sin(v)}; },
            std::array<T, 4>{0., 2. * pi, -pi / 2., pi / 2.});
    }
    std::shared_ptr<Surface<T, 3>> torus(T R = 3., T r = 1.)
    {
        return std::make_shared<Analytic>(
            [R, r](T u, T v) { return point<T, 3>{(R + r * std::cos(v)) * std::cos(u), (R + r * std::cos(v)) * std::sin(u), r * std::sin(v)}; },
            std::array<T, 4>{0., 2. * pi, 0., 2. * pi});
    }
    std::shared_ptr<Surface<T, 3>> cone()
    {
        return std::make_shared<Analytic>([](T u, T v) { return point<T, 3>{v * std::cos(u), v * std::sin(u), v}; },
                                          std::array<T, 4>{0., 2. * pi, 0., 1.});
    }

    // Exact rational NURBS cylinder: the rational circle extruded along z.
    std::shared_ptr<Surface<T, 3>> nurbs_cylinder(T h = 2.)
    {
        auto circle = build_circle<T, 3>(1.);
        points_vector<T, 4> poles;
        for (T z : {T(0), h})
            for (auto p : circle.poles())
            {
                p[2] += p[3] * z; // weighted homogeneous coordinates
                poles.push_back(p);
            }
        return std::make_shared<BSSurfaceRational<T, 3>>(poles, circle.knotsFlats(), std::vector<T>{0., 0., h, h},
                                                         circle.degree(), 1);
    }

    std::shared_ptr<Surface<T, 3>> plane()
    {
        points_vector<T, 3> poles{{0., 0., 0.}, {2., 0., 0.}, {0., 1., 0.}, {2., 1., 0.}};
        std::vector<T> k{0., 0., 1., 1.};
        return std::make_shared<BSSurface<T, 3>>(poles, k, k, 1, 1);
    }

    std::shared_ptr<Surface<T, 3>> revolution(T angle)
    {
        auto profile = std::make_shared<BSCurve<T, 2>>(build_segment<T, 2>({1., 1.}, {2., 3.}));
        return std::make_shared<SurfaceOfRevolution<T>>(profile, ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}}, 0., angle);
    }

    // Analytic 3D curve C(t) given by a lambda; first derivative by central differences.
    class AnalyticCurve : public Curve<T, 3>
    {
        std::function<point<T, 3>(T)> f_;
        std::array<T, 2> b_;

    public:
        AnalyticCurve(std::function<point<T, 3>(T)> f, std::array<T, 2> b) : f_{std::move(f)}, b_{b} {}
        auto value(T t, size_t d = 0) const -> std::array<T, 3> override
        {
            const T h = 1e-6;
            if (d == 0)
                return f_(t);
            if (d == 1)
                return (f_(t + h) - f_(t - h)) * (T(0.5) / h);
            throw std::runtime_error("not implemented");
        }
        auto bounds() const -> std::array<T, 2> override { return b_; }
    };
}
