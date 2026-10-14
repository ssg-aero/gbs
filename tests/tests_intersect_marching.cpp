#include <doctest_gtest.hpp>
#include <gbs/intersect/surface_surface.h>
#include <gbs/bssinterp.h>
#include <gbs/primitives.h>

#include <chrono>
#include <cmath>
#include <numbers>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

using namespace gbs;

namespace
{
    using T = double;
    constexpr T pi = std::numbers::pi;
    constexpr T tol = 1e-6;

    const ax2<T, 3> tilted{{{1., -2., 0.5}, {1., 1., 1.}, {1., 0., 0.}}};
    const Frame3<T> F{tilted};

    Frame3<T> frame(const point<T, 3> &o, const point<T, 3> &z, const point<T, 3> &x)
    {
        return intersect_detail::frame_of<T>(o, intersect_detail::unit<T>(z), x);
    }

    void expect_valid(const Surface<T, 3> &s1, const std::array<T, 4> &b1, const Surface<T, 3> &s2, const std::array<T, 4> &b2,
                      const SurfaceSurfaceResult<T> &r)
    {
        REQUIRE(r.complete);
        for (const auto &c : r.curves)
        {
            INFO("curve on [" << c.a << ", " << c.b << "]");
            REQUIRE(c.b > c.a);
            REQUIRE(c.deviation <= tol);
            for (int k = 0; k <= 50; ++k)
            {
                const T t = c.a + (c.b - c.a) * k / 50.;
                const auto p = c.curve->value(t);
                const auto q1 = c.pc1->value(t), q2 = c.pc2->value(t);
                REQUIRE(norm(s1.value(q1[0], q1[1]) - p) <= 2 * tol);
                REQUIRE(norm(s2.value(q2[0], q2[1]) - p) <= 2 * tol);
                REQUIRE(q1[0] >= b1[0] - 1e-6);
                REQUIRE(q1[0] <= b1[1] + 1e-6);
                REQUIRE(q1[1] >= b1[2] - 1e-6);
                REQUIRE(q1[1] <= b1[3] + 1e-6);
                REQUIRE(q2[0] >= b2[0] - 1e-6);
                REQUIRE(q2[0] <= b2[1] + 1e-6);
                REQUIRE(q2[1] >= b2[2] - 1e-6);
                REQUIRE(q2[1] <= b2[3] + 1e-6);
            }
        }
    }

    T length(const IntersectionCurve<T> &c)
    {
        T l{};
        constexpr int n = 4000;
        for (int k = 0; k < n; ++k)
            l += norm(c.curve->value(c.a + (c.b - c.a) * (k + 1.) / n) - c.curve->value(c.a + (c.b - c.a) * T(k) / n));
        return l;
    }

    T total_length(const SurfaceSurfaceResult<T> &r)
    {
        T l{};
        for (const auto &c : r.curves)
            l += length(c);
        return l;
    }

    SurfaceSurfaceResult<T> marching(const Surface<T, 3> &s1, const Surface<T, 3> &s2, IntersectOptions<T> o = {})
    {
        o.closed_forms = false;
        auto r = intersect(s1, s1.bounds(), s2, s2.bounds(), o);
        expect_valid(s1, s1.bounds(), s2, s2.bounds(), r);
        return r;
    }
}

TEST(tests_intersect_marching, against_closed_forms)
{
    // the same pairs of primitives through their closed forms and by marching: same lengths
    const CylindricalSurface<T> cyl{2., F, {0., 2 * pi, -5., 5.}};
    const SphericalSurface<T> sphere{2.5, F, {0., 2 * pi, -pi / 2, pi / 2}};
    Frame3<T> across = F;
    across.o = F.o + 1. * F.z;
    const Plane<T> plane{across, {-5., 5., -5., 5.}};
    const Plane<T> oblique{frame(F.o, F.z + 0.6 * F.y, F.x), {-6., 6., -6., 6.}};
    const Plane<T> generatrices{frame(F.o + 0.7 * F.x, F.x, F.z), {-6., 6., -6., 6.}};
    struct Pair
    {
        const Surface<T, 3> *a, *b;
    };
    for (const auto &[a, b] : {Pair{&cyl, &plane}, Pair{&cyl, &oblique}, Pair{&cyl, &generatrices}, Pair{&sphere, &oblique}, Pair{&cyl, &sphere}})
    {
        auto exact = intersect(*a, *b);
        REQUIRE(exact.complete);
        auto m = marching(*a, *b);
        ASSERT_FALSE(m.singular);
        ASSERT_NEAR(total_length(m), total_length(exact), 1e-5);
    }
}

TEST(tests_intersect_marching, crossed_cylinders)
{
    // radii 1 and 0.5, perpendicular axes: two loops, each cut at the seams
    const CylindricalSurface<T> big{1., F, {0., 2 * pi, -3., 3.}};
    const CylindricalSurface<T> small{0.5, frame(F.o, F.x, F.y), {0., 2 * pi, -3., 3.}};
    auto r = intersect(big, small);
    expect_valid(big, big.bounds(), small, small.bounds(), r);
    ASSERT_FALSE(r.singular);
    // reference: on the small cylinder, x = +-sqrt(1 - r² cos² a): both loops
    T ref{};
    constexpr int n = 200000;
    const T R = 1., rr = 0.5;
    auto P = [&](T a) { return point<T, 3>{std::sqrt(R * R - rr * rr * std::cos(a) * std::cos(a)), rr * std::cos(a), rr * std::sin(a)}; };
    for (int k = 0; k < n; ++k)
        ref += norm(P(2 * pi * (k + 1.) / n) - P(2 * pi * T(k) / n));
    ASSERT_NEAR(total_length(r), 2 * ref, 1e-5);

    // equal radii (Steinmetz): the loops cross at branch points, flagged singular
    const CylindricalSurface<T> twin{1., frame(F.o, F.x, F.y), {0., 2 * pi, -3., 3.}};
    auto s = intersect(big, twin);
    REQUIRE(s.complete);
    ASSERT_TRUE(s.singular);
    for (const auto &c : s.curves)
        ASSERT_LE(c.deviation, tol);
}

TEST(tests_intersect_marching, plane_cuts)
{
    // a parabola on the cone: every point satisfies the cone equation in the plane
    const ConicalSurface<T> cone{1., pi / 6, F, {0., 2 * pi, -1., 3.}};
    const auto pz = intersect_detail::unit<T>(std::cos(pi / 6) * F.x - std::sin(pi / 6) * F.z);
    const Plane<T> pl{frame(F.o + 0.5 * F.x, pz, F.y), {-9., 9., -9., 9.}};
    auto r = intersect(cone, pl);
    expect_valid(cone, cone.bounds(), pl, pl.bounds(), r);
    ASSERT_FALSE(r.curves.empty());
    // an oblique plane on a torus: primitive (sampled patches) and exact NURBS give the same length
    const ToroidalSurface<T> torus{3., 1., F, {0., 2 * pi, 0., 2 * pi}};
    const auto tn = torus.to_nurbs();
    const Plane<T> ob{frame(F.o + 0.3 * F.z, F.z + 0.2 * F.x, F.x - 0.2 * F.z), {-9., 9., -9., 9.}};
    auto rt = intersect(torus, ob);
    expect_valid(torus, torus.bounds(), ob, ob.bounds(), rt);
    auto rn = intersect(*tn, ob);
    expect_valid(*tn, tn->bounds(), ob, ob.bounds(), rn);
    ASSERT_NEAR(total_length(rt), total_length(rn), 1e-5);
    // the NURBS of a cylinder against a plane: a closed circle away from the seams? the circle
    // crosses the seam of the NURBS: same length as the closed form
    const CylindricalSurface<T> cyl{2., F, {0., 2 * pi, -5., 5.}};
    const auto cn = cyl.to_nurbs();
    Frame3<T> across = F;
    across.o = F.o + 1. * F.z;
    const Plane<T> plane{across, {-5., 5., -5., 5.}};
    auto rc = intersect(*cn, plane);
    expect_valid(*cn, cn->bounds(), plane, plane.bounds(), rc);
    ASSERT_NEAR(total_length(rc), 4 * pi, 1e-5);
}

TEST(tests_intersect_marching, free_surfaces)
{
    // a bump z = 1 - x² - y² (B-spline) cut by the plane z = 0.5: a closed loop inside both boxes
    points_vector<T, 3> grid;
    const int n = 9;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
        {
            const T x = -1.5 + 3. * i / (n - 1), y = -1.5 + 3. * j / (n - 1);
            grid.push_back({x, y, 1. - x * x - y * y});
        }
    const auto bump = interpolate<T, 3>(grid, n, 3, 3, KnotsCalcMode::CHORD_LENGTH);
    const Plane<T> plane{Frame3<T>{ax2<T, 3>{{{0., 0., 0.5}, {0., 0., 1.}, {1., 0., 0.}}}}, {-3., 3., -3., 3.}};
    auto r = intersect(bump, plane);
    expect_valid(bump, bump.bounds(), plane, plane.bounds(), r);
    ASSERT_EQ(r.curves.size(), 1u);
    ASSERT_TRUE(r.curves[0].closed);
    ASSERT_NEAR(length(r.curves[0]), 2 * pi * std::sqrt(0.5), 0.03); // the B-spline interpolates the paraboloid
    // two free surfaces: the bump against a tilted wave
    points_vector<T, 3> wave;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
        {
            const T x = -2. + 4. * i / (n - 1), y = -2. + 4. * j / (n - 1);
            wave.push_back({x, y, 0.3 * x + 0.2 * std::sin(2 * y)});
        }
    const auto w = interpolate<T, 3>(wave, n, 3, 3, KnotsCalcMode::CHORD_LENGTH);
    auto rw = intersect(bump, w);
    expect_valid(bump, bump.bounds(), w, w.bounds(), rw);
    ASSERT_FALSE(rw.curves.empty());
    ASSERT_FALSE(rw.singular);

    // a sphere touching the plane from above: one contact point, by marching
    for (const point<T, 3> axis : {point<T, 3>{1., 0., 0.}, point<T, 3>{0., 0., 1.}}) // contact off and on a pole
    {
        const SphericalSurface<T> ball{1., Frame3<T>{ax2<T, 3>{{{0.2, 0.1, 1.5}, axis, {0., 1., 0.}}}}, {0., 2 * pi, -pi / 2, pi / 2}};
        auto rp = marching(ball, plane);
        ASSERT_TRUE(rp.curves.empty());
        ASSERT_FALSE(rp.points.empty());
        ASSERT_LT(norm(rp.points[0].p - point<T, 3>{0.2, 0.1, 0.5}), 1e-3);
    }
}

TEST(tests_intersect_marching, singular_points_and_limits)
{
    // Viviani: a sphere of radius 1 and a cylinder of radius 1/2 through its centre; the curve
    // crosses itself at (1, 0, 0), where the surfaces are tangent
    const SphericalSurface<T> sphere{1., frame(F.o, F.y, F.z), {0., 2 * pi, -pi / 2, pi / 2}};
    const CylindricalSurface<T> cyl{0.5, frame(F.o + 0.5 * F.x, F.z, F.y), {0., 2 * pi, -2., 2.}}; // seam away from (1, 0, 0)
    auto r = intersect(sphere, cyl);
    REQUIRE(r.complete);
    ASSERT_TRUE(r.singular);
    ASSERT_FALSE(r.curves.empty());
    for (const auto &c : r.curves)
        ASSERT_LE(c.deviation, tol);
    // the cap on the pairs of patches
    IntersectOptions<T> capped;
    capped.max_pairs = 1;
    const CylindricalSurface<T> big{1., F, {0., 2 * pi, -3., 3.}};
    const CylindricalSurface<T> other{0.5, frame(F.o, F.x, F.y), {0., 2 * pi, -3., 3.}};
    ASSERT_FALSE(intersect(big, other, capped).complete);
}
