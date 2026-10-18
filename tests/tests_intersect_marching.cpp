#include <doctest_gtest.hpp>
#include <gbs/intersect/surface_surface.h>
#include <gbs/bscbuild.h>
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

    // equal radii (Steinmetz): two ellipses of semi-axes sqrt(2) and 1 crossing at two branch
    // points, stepped over; the seam of the twin passes through them, so turn it
    for (T a0 : {0., 1.})
    {
        INFO("seam of the twin at " << a0);
        const CylindricalSurface<T> twin{1., frame(F.o, F.x, std::cos(a0) * F.y + std::sin(a0) * F.z), {0., 2 * pi, -3., 3.}};
        auto s = intersect(big, twin);
        expect_valid(big, big.bounds(), twin, twin.bounds(), s);
        ASSERT_FALSE(s.singular);
        T ellipse{};
        auto E = [](T a) { return point<T, 3>{std::sqrt(2.) * std::cos(a), std::sin(a), 0.}; };
        for (int k = 0; k < n; ++k)
            ellipse += norm(E(2 * pi * (k + 1.) / n) - E(2 * pi * T(k) / n));
        ASSERT_NEAR(total_length(s), 2 * ellipse, 1e-5);
        std::size_t spans{};
        for (const auto &c : s.curves)
            spans += c.branch_spans.size();
        if (a0 != 0.)
            ASSERT_EQ(spans, 4u); // each ellipse steps over both branch points
    }
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
    expect_valid(sphere, sphere.bounds(), cyl, cyl.bounds(), r);
    ASSERT_FALSE(r.singular); // the double point is stepped over
    T viviani{};
    constexpr int n = 200000;
    auto V = [](T a) { return point<T, 3>{0.5 * (1. + std::cos(a)), 0.5 * std::sin(a), std::sin(a / 2.)}; };
    for (int k = 0; k < n; ++k)
        viviani += norm(V(4 * pi * (k + 1.) / n) - V(4 * pi * T(k) / n));
    ASSERT_NEAR(total_length(r), viviani, 1e-5);
    std::size_t spans{};
    for (const auto &c : r.curves)
        spans += c.branch_spans.size();
    ASSERT_EQ(spans, 2u); // the curve goes twice through its double point
    // a tangency along a curve still stops the marching: a sphere in a cylinder of its radius
    IntersectOptions<T> marched;
    marched.closed_forms = false;
    const CylindricalSurface<T> sleeve{1., frame(F.o, F.y, F.z), {0., 2 * pi, -2., 2.}};
    ASSERT_TRUE(intersect(sphere, sleeve, marched).singular);
    // the cap on the pairs of patches
    IntersectOptions<T> capped;
    capped.max_pairs = 1;
    const CylindricalSurface<T> big{1., F, {0., 2 * pi, -3., 3.}};
    const CylindricalSurface<T> other{0.5, frame(F.o, F.x, F.y), {0., 2 * pi, -3., 3.}};
    ASSERT_FALSE(intersect(big, other, capped).complete);
}

TEST(tests_intersect_marching, branch_points_on_a_seam)
{
    // Steinmetz with the seam of the twin through the branch points, or very close to them: the
    // branches land on the seam at the branch point (rounding differs across platforms)
    const CylindricalSurface<T> big{1., F, {0., 2 * pi, -3., 3.}};
    for (T a0 : {0., 1e-15, -1e-15, 1e-13, -1e-13, 1e-11, -1e-11, 1e-9, -1e-9, 1e-7, -1e-7, 1e-5, -1e-5, 1e-3, -1e-3})
        for (T dz : {0., 1e-13, 1e-10})
        {
            INFO("seam at " << a0 << ", shift " << dz);
            const CylindricalSurface<T> twin{1., frame(F.o + dz * F.z, F.x, std::cos(a0) * F.y + std::sin(a0) * F.z), {0., 2 * pi, -3., 3.}};
            auto s = intersect(big, twin);
            REQUIRE(s.complete);
            REQUIRE_FALSE(s.singular);
        }
}

TEST(tests_intersect_marching, through_poles)
{
    // Viviani with the axis of the sphere along the cylinder: the curve goes through both poles,
    // where the parametrization of the sphere degenerates; the branches end at the poles
    const SphericalSurface<T> sphere{1., frame(F.o, F.z, F.y), {0., 2 * pi, -pi / 2, pi / 2}};
    const CylindricalSurface<T> cyl{0.5, frame(F.o + 0.5 * F.x, F.z, F.y), {0., 2 * pi, -2., 2.}};
    auto r = intersect(sphere, cyl);
    expect_valid(sphere, sphere.bounds(), cyl, cyl.bounds(), r);
    ASSERT_FALSE(r.singular);
    T viviani{};
    constexpr int n = 200000;
    auto V = [](T a) { return point<T, 3>{0.5 * (1. + std::cos(a)), 0.5 * std::sin(a), std::sin(a / 2.)}; };
    for (int k = 0; k < n; ++k)
        viviani += norm(V(4 * pi * (k + 1.) / n) - V(4 * pi * T(k) / n));
    ASSERT_NEAR(total_length(r), viviani, 1e-5);
    std::size_t at_poles{};
    for (const auto &c : r.curves)
        for (T t : {c.a, c.b})
            at_poles += std::abs(std::abs(c.pc1->value(t)[1]) - pi / 2) < 1e-9;
    ASSERT_EQ(at_poles, 4u); // two branches arrive at each pole
}

TEST(tests_intersect_marching, along_a_border)
{
    // the side of a cylinder as a B-spline extruded from z = 0 to 2, against planes that run
    // along its borders, where the parameters go beyond the bounds by rounding errors: closed
    // curves on both surfaces, not a score of branches stopped off the intersection
    auto extruded = [](const auto &circle) {
        using Pt = std::remove_cvref_t<decltype(circle.poles().front())>;
        points_vector<T, std::tuple_size_v<Pt>> poles;
        for (T z : {T(0), T(2)})
            for (auto p : circle.poles())
            {
                p[2] += (p.size() == 4 ? p[3] : T(1)) * z;
                poles.push_back(p);
            }
        return poles;
    };
    const auto circle = build_circle<T, 3>(1.);
    const BSSurfaceRational<T, 3> side{extruded(circle), circle.knotsFlats(), std::vector<T>{0., 0., 2., 2.}, circle.degree(), 1};
    points_vector<T, 3> pts;
    for (int k = 0; k <= 16; ++k)
        pts.push_back(circle.value(circle.bounds()[0] + (circle.bounds()[1] - circle.bounds()[0]) * k / 16.));
    const auto loop = interpolate<T, 3>(pts, 3, KnotsCalcMode::CHORD_LENGTH);
    const BSSurface<T, 3> polynomial{extruded(loop), loop.knotsFlats(), std::vector<T>{0., 0., 2., 2.}, loop.degree(), 1};
    T loop_length{};
    for (int k = 0; k < 4000; ++k)
        loop_length += norm(loop.value(loop.bounds()[0] + (loop.bounds()[1] - loop.bounds()[0]) * (k + 1.) / 4000.) -
                            loop.value(loop.bounds()[0] + (loop.bounds()[1] - loop.bounds()[0]) * k / 4000.));
    IntersectOptions<T> o;
    o.tol = tol;
    const std::array<T, 4> touching{-1., 1., -1., 1.}, wide{-2., 2., -2., 2.};
    auto cut = [&](const Surface<T, 3> &s, T z, const std::array<T, 4> &box) {
        const Plane<T> plane{ax2<T, 3>{{{0., 0., z}, {0., 0., 1.}, {1., 0., 0.}}}, wide};
        auto r = intersect(plane, box, s, s.bounds(), o);
        expect_valid(plane, box, s, s.bounds(), r);
        CHECK_FALSE(r.singular);
        return r;
    };
    // the lower and the upper border, on the box [-1, 1]² that the circle touches (at x = 1,
    // where it crosses the seam of the side: the curve leaves the side there, not the box) and
    // on a wider box
    for (T z : {T(0), T(2)})
        for (const auto &box : {touching, wide})
        {
            INFO("rational side, z = " << z << ", box " << box[0] << ", " << box[1]);
            auto r = cut(side, z, box);
            REQUIRE(r.curves.size() == 1u);
            const auto &c = r.curves[0];
            CHECK(c.closed);
            CHECK(norm(c.curve->value(c.a) - c.curve->value(c.b)) <= 1e-12);
            CHECK(length(c) == doctest::Approx(2 * pi).epsilon(1e-6));
        }
    for (T z : {T(0), T(2)})
    {
        INFO("polynomial side, z = " << z);
        auto r = cut(polynomial, z, wide);
        REQUIRE(r.curves.size() == 1u);
        CHECK(r.curves[0].closed);
        CHECK(length(r.curves[0]) == doctest::Approx(loop_length).epsilon(1e-6));
    }
    // beyond the bound by more than the rounding margin: the plane misses the side, no curve
    // snapped on the border
    for (T z : {T(-1e-6), T(-1e-9), T(2 + 1e-9)})
    {
        INFO("rational side, z = " << z);
        CHECK(cut(side, z, wide).curves.empty());
    }
    // a border in u: the plane y = 0 through the axis cuts the side along its seam (u = 0) and
    // along the opposite line
    const Plane<T> meridian{ax2<T, 3>{{{0., 0., 0.}, {0., 1., 0.}, {1., 0., 0.}}}, {-2., 2., -3., 3.}};
    auto rm = intersect(meridian, meridian.bounds(), side, side.bounds(), o);
    expect_valid(meridian, meridian.bounds(), side, side.bounds(), rm);
    CHECK_FALSE(rm.singular);
    CHECK(total_length(rm) == doctest::Approx(4.).epsilon(1e-6));
}
