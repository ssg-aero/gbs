#include <doctest_gtest.hpp>
#include <gbs/bselementary.h>
#include <gbs/bscbuild.h>

#include <cmath>
#include <numbers>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

using namespace gbs;
using gbs::operator-;
using gbs::operator+;
using gbs::operator*;

namespace
{
    using T = double;
    constexpr T pi = std::numbers::pi;
    constexpr T tol = 1e-13;

    // A tilted frame: {origin, axis, reference direction not orthogonal to the axis}.
    const ax2<T, 3> tilted{{{1., -2., 0.5}, {1., 1., 1.}, {1., 0., 0.}}};

    // Orthonormal X, Y, Z of a frame, computed independently of the library.
    std::array<point<T, 3>, 3> xyz(const ax2<T, 3> &ax)
    {
        auto Z = ax[1] / norm(ax[1]);
        auto X = ax[2] - (ax[2] * Z) * Z;
        X = X / norm(X);
        return {X, cross(Z, X), Z};
    }

    T dist(const point<T, 3> &a, const point<T, 3> &b) { return norm(a - b); }
}

TEST(tests_bselementary, circle_arcs)
{
    const auto [X, Y, Z] = xyz(tilted);
    const T r = 2.5;
    // spans below, at and above the quarter turns, including the full circle and negative angles
    for (auto [t1, t2, n] : std::vector<std::tuple<T, T, size_t>>{
             {0., pi / 6, 1}, {0., pi / 2, 1}, {-pi / 4, pi / 2, 2}, {0.3, 0.3 + 1.5 * pi, 3}, {0., 2 * pi, 4}, {-pi, pi, 4}})
    {
        auto c = build_circle_arc<T>(r, t1, t2, tilted);
        ASSERT_EQ(c.degree(), 2);
        ASSERT_EQ(c.poles().size(), 2 * n + 1);
        ASSERT_NEAR(c.bounds()[0], t1, tol);
        ASSERT_NEAR(c.bounds()[1], t2, tol);
        auto exact = [&](T a) { return tilted[0] + (r * std::cos(a)) * X + (r * std::sin(a)) * Y; };
        ASSERT_LT(dist(c.begin(), exact(t1)), tol);
        ASSERT_LT(dist(c.end(), exact(t2)), tol);
        for (int i = 0; i <= 200; ++i)
        {
            const T a = t1 + (t2 - t1) * i / 200.;
            const T u = ellipse_arc_parameter(t1, t2, a);
            ASSERT_NEAR(ellipse_arc_angle(t1, t2, u), a, 1e-12);
            ASSERT_LT(dist(c(u), exact(a)), tol);                       // the point at the angle
            const T v = t1 + (t2 - t1) * i / 200.;                       // any parameter: on the circle
            ASSERT_NEAR(norm(c(v) - tilted[0]), r, tol);
            ASSERT_NEAR((c(v) - tilted[0]) * Z, 0., tol);
        }
        // the parameter is the angle at the knots, and the tangent is continuous there
        const auto &k = c.knotsFlats();
        for (size_t i = 3; i + 3 < k.size(); i += 2)
        {
            ASSERT_LT(dist(c(k[i]), exact(k[i])), tol);
            ASSERT_LT(norm(c(k[i] - 1e-9, 1) - c(k[i] + 1e-9, 1)), 1e-6);
        }
    }
    // 2D arc
    auto c2 = build_circle_arc<T>(1., pi / 3, pi, point<T, 2>{1., 1.});
    ASSERT_NEAR(c2(ellipse_arc_parameter(pi / 3, pi, 2.)) [0], 1. + std::cos(2.), tol);
    ASSERT_NEAR(c2(ellipse_arc_parameter(pi / 3, pi, 2.)) [1], 1. + std::sin(2.), tol);
}

TEST(tests_bselementary, ellipse_arcs)
{
    const auto [X, Y, Z] = xyz(tilted);
    const T a = 3., b = 1.2;
    auto e = build_ellipse_arc<T>(a, b, -2., 3.5, tilted);
    ASSERT_EQ(e.poles().size(), 2 * 4 + 1);
    for (int i = 0; i <= 300; ++i)
    {
        const T t = -2. + 5.5 * i / 300.;
        // eccentric angle as for STEP: C(t) = O + a cos t X + b sin t Y
        ASSERT_LT(dist(e(ellipse_arc_parameter(-2., 3.5, t)), tilted[0] + (a * std::cos(t)) * X + (b * std::sin(t)) * Y), tol);
        const auto d = e(t) - tilted[0];
        const T x = d * X, y = d * Y;
        ASSERT_NEAR(x * x / (a * a) + y * y / (b * b), 1., tol);
        ASSERT_NEAR(d * Z, 0., tol);
    }
}

TEST(tests_bselementary, invalid_arguments)
{
    ASSERT_THROW(build_circle_arc<T>(1., 1., 1., tilted), std::invalid_argument);       // empty span
    ASSERT_THROW(build_circle_arc<T>(1., 1., 0., tilted), std::invalid_argument);       // reversed
    ASSERT_THROW(build_circle_arc<T>(1., 0., 7., tilted), std::invalid_argument);       // more than a turn
    ASSERT_THROW(build_circle_arc<T>(0., 0., 1., tilted), std::invalid_argument);       // null radius
    ASSERT_THROW(build_ellipse_arc<T>(1., -1., 0., 1., tilted), std::invalid_argument); // negative radius
    ASSERT_THROW(build_circle_arc<T>(1., 0., 1., ax2<T, 3>{{{0., 0., 0.}, {0., 0., 0.}, {1., 0., 0.}}}), std::invalid_argument);
    ASSERT_THROW(build_circle_arc<T>(1., 0., 1., ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {0., 0., 2.}}}), std::invalid_argument);
    ASSERT_THROW(build_cone<T>(1., pi / 2, tilted, 0., 1.), std::invalid_argument);
    ASSERT_THROW(build_sphere<T>(1., tilted, 0., pi, -2., 1.), std::invalid_argument);
    ASSERT_THROW(build_extrusion(build_segment<T, 3>({0., 0., 0.}, {1., 0., 0.}), point<T, 3>{0., 0., 1.}, 1., 1.), std::invalid_argument);
}

TEST(tests_bselementary, revolution)
{
    const ax1<T, 3> axis{{{1., 2., 3.}, {0., 1., 1.}}};
    const auto Z = axis[1] / norm(axis[1]);
    // rotation of p by angle a around the axis, right-handed
    auto rotated = [&](const point<T, 3> &p, T a) {
        const auto O = axis[0] + ((p - axis[0]) * Z) * Z;
        const auto Xr = p - O;
        return O + std::cos(a) * Xr + std::sin(a) * cross(Z, Xr);
    };
    // a rational generatrix (an arc of circle off the axis) and a polynomial one touching the axis
    const auto arc = build_circle_arc<T>(0.7, 0.2, 2.8, ax2<T, 3>{{{3., 0., 1.}, {1., 0., 0.}, {0., 1., 0.}}});
    const BSCurve<T, 3> poly{points_vector<T, 3>{axis[0] + 2. * Z, {2., 1., 0.}, {4., 4., 1.}}, std::vector<T>{0., 0., 0., 1., 1., 1.}, 2};
    for (const auto &[t1, t2] : std::vector<std::pair<T, T>>{{0., 2 * pi}, {-0.5, 1.7}})
    {
        auto check = [&](const auto &g) {
            auto s = build_revolution(g, axis, t1, t2);
            ASSERT_EQ(s.degreeU(), 2);
            ASSERT_EQ(s.degreeV(), g.degree());
            ASSERT_NEAR(s.boundsV()[0], g.bounds()[0], tol);
            ASSERT_NEAR(s.boundsV()[1], g.bounds()[1], tol);
            for (int i = 0; i <= 20; ++i)
                for (int j = 0; j <= 20; ++j)
                {
                    const T a = t1 + (t2 - t1) * i / 20.;
                    const T v = g.bounds()[0] + (g.bounds()[1] - g.bounds()[0]) * j / 20.;
                    ASSERT_LT(dist(s(ellipse_arc_parameter(t1, t2, a), v), rotated(g(v), a)), 1e-12);
                }
        };
        check(arc);
        check(poly);
    }
    // the generatrix pole on the axis gives a degenerate row: all its poles coincide
    auto s = build_revolution(poly, axis);
    auto cartesian = [](const std::array<T, 4> &h) { return point<T, 3>{h[0] / h[3], h[1] / h[3], h[2] / h[3]}; };
    for (size_t i = 1; i < s.nPolesU(); ++i)
        ASSERT_LT(dist(cartesian(s.poles()[i]), axis[0] + 2. * Z), tol);
}

TEST(tests_bselementary, extrusion)
{
    const point<T, 3> V{0.5, -1., 2.};
    auto check = [&](const auto &c) {
        auto s = build_extrusion(c, V, -1., 2.);
        ASSERT_EQ(s.degreeV(), 1);
        for (int i = 0; i <= 20; ++i)
            for (T v : {-1., 0., 0.7, 2.})
            {
                const T u = c.bounds()[0] + (c.bounds()[1] - c.bounds()[0]) * i / 20.;
                ASSERT_LT(dist(s(u, v), c(u) + v * V), tol);
            }
        return s;
    };
    auto rational = check(build_circle_arc<T>(1., 0., 2., tilted));
    auto polynomial = check(build_segment<T, 3>({0., 0., 0.}, {1., 2., 0.}));
    static_assert(std::is_same_v<decltype(rational), BSSurfaceRational<T, 3>>);
    static_assert(std::is_same_v<decltype(polynomial), BSSurface<T, 3>>);
}

TEST(tests_bselementary, elementary_surfaces)
{
    const auto [X, Y, Z] = xyz(tilted);
    const auto &O = tilted[0];
    auto radial = [&](T u) { return std::cos(u) * X + std::sin(u) * Y; };
    constexpr int N = 24;

    // cylinder: S(u, v) = O + R (cos u X + sin u Y) + v Z
    {
        const T R = 1.5;
        auto s = build_cylinder<T>(R, tilted, -1., 2.);
        for (int i = 0; i <= N; ++i)
            for (int j = 0; j <= N; ++j)
            {
                const T u = 2 * pi * i / N, v = -1. + 3. * j / N;
                ASSERT_LT(dist(s(ellipse_arc_parameter(0., 2 * pi, u), v), O + R * radial(u) + v * Z), tol);
                const auto d = s(u, v) - O; // any parameter: on the cylinder
                ASSERT_NEAR(norm(d - (d * Z) * Z), R, tol);
            }
    }
    // cone with its apex at v = -R / tan(alpha): S(u, v) = O + (R + v tan a)(cos u X + sin u Y) + v Z
    {
        const T R = 1., alpha = pi / 6;
        const T apex = -R / std::tan(alpha);
        auto s = build_cone<T>(R, alpha, tilted, apex, 1., 0.5, 2.5);
        for (int i = 0; i <= N; ++i)
            for (int j = 0; j <= N; ++j)
            {
                const T u = 0.5 + 2. * i / N, v = apex + (1. - apex) * j / N;
                ASSERT_LT(dist(s(ellipse_arc_parameter(0.5, 2.5, u), v), O + (R + v * std::tan(alpha)) * radial(u) + v * Z), 1e-12);
            }
        for (size_t i = 1; i < s.nPolesU(); ++i) // degenerate side at the apex
            ASSERT_LT(dist(s.value(s.boundsU()[0] + 0.1 * i, apex), O + apex * Z), tol);
    }
    // sphere: S(u, v) = O + R cos v (cos u X + sin u Y) + R sin v Z
    {
        const T R = 2.;
        auto s = build_sphere<T>(R, tilted);
        for (int i = 0; i <= N; ++i)
            for (int j = 0; j <= N; ++j)
            {
                const T u = 2 * pi * i / N, v = -pi / 2 + pi * j / N;
                ASSERT_LT(dist(s(ellipse_arc_parameter(0., 2 * pi, u), ellipse_arc_parameter(-pi / 2, pi / 2, v)),
                               O + (R * std::cos(v)) * radial(u) + (R * std::sin(v)) * Z), tol);
                ASSERT_NEAR(norm(s(u, v) - O), R, tol);
            }
        ASSERT_LT(dist(s(1., -pi / 2), O - R * Z), tol); // poles
        ASSERT_LT(dist(s(4., pi / 2), O + R * Z), tol);
    }
    // torus: S(u, v) = O + (R + r cos v)(cos u X + sin u Y) + r sin v Z
    {
        const T R = 3., r = 0.8;
        auto s = build_torus<T>(R, r, tilted, 0., pi, -1., 2.);
        for (int i = 0; i <= N; ++i)
            for (int j = 0; j <= N; ++j)
            {
                const T u = pi * i / N, v = -1. + 3. * j / N;
                ASSERT_LT(dist(s(ellipse_arc_parameter(0., pi, u), ellipse_arc_parameter(-1., 2., v)),
                               O + (R + r * std::cos(v)) * radial(u) + (r * std::sin(v)) * Z), 1e-12);
                const auto d = s(u, v) - O;
                const T h = d * Z;
                const T rho = norm(d - h * Z);
                ASSERT_NEAR((rho - R) * (rho - R) + h * h, r * r, 1e-12);
            }
    }
}
