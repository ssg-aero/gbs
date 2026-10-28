#include <doctest_gtest.hpp>
#include <gbs/intersect/surface_surface.h>
#include <gbs-brep/offset.h>
#include <gbs/bscbuild.h>
#include <gbs/bssinterp.h>
#include <gbs/primitives.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
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

    /// z = c0 + c1 x + c2 y + c3 x² + c4 x y + c5 y² on [x0, x1] x [y0, y1], exact: a Bézier patch
    /// of degree 2 (poles by blossoming), parameters [0, 1]², moved by `off`
    BSSurface<T, 3> quadratic(T x0, T x1, T y0, T y1, const std::array<T, 6> &c, const point<T, 3> &off = {0., 0., 0.})
    {
        points_vector<T, 3> poles;
        const T xs[3][2] = {{x0, x0}, {x0, x1}, {x1, x1}}, ys[3][2] = {{y0, y0}, {y0, y1}, {y1, y1}};
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 3; ++i)
            {
                const T xa = xs[i][0], xb = xs[i][1], ya = ys[j][0], yb = ys[j][1];
                const T z = c[0] + c[1] * (xa + xb) / 2 + c[2] * (ya + yb) / 2 + c[3] * xa * xb + c[4] * (xa + xb) * (ya + yb) / 4 + c[5] * ya * yb;
                poles.push_back(point<T, 3>{(xa + xb) / 2, (ya + yb) / 2, z} + off);
            }
        const std::vector<T> k{0., 0., 0., 1., 1., 1.};
        return {poles, k, k, 2, 2};
    }

    /// A cubic sheet z = 0 on [-L, L]² of n x n poles, the pole nearest to (x, y) lifted by `lift`;
    /// uniform knots on [a, b]
    BSSurface<T, 3> sheet(T L, int n, T x, T y, T lift, T a = 0., T b = 1.)
    {
        points_vector<T, 3> poles;
        auto at = [&](int i) { return -L + 2 * L * i / (n - 1); };
        const int ib = int(std::lround((x + L) / (2 * L) * (n - 1))), jb = int(std::lround((y + L) / (2 * L) * (n - 1)));
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
                poles.push_back({at(i), at(j), i == ib && j == jb ? lift : 0.});
        std::vector<T> k{a, a, a, a};
        for (int s = 1; s < n - 3; ++s)
            k.push_back(a + (b - a) * s / (n - 3));
        k.insert(k.end(), {b, b, b, b});
        return {poles, k, k, 3, 3};
    }

    /// Length of a curve with fine chords (sharp ends of thin loops).
    T fine_length(const IntersectionCurve<T> &c)
    {
        T l{};
        constexpr int n = 200000;
        for (int k = 0; k < n; ++k)
            l += norm(c.curve->value(c.a + (c.b - c.a) * (k + 1.) / n) - c.curve->value(c.a + (c.b - c.a) * T(k) / n));
        return l;
    }

    /// Perimeter of the ellipse of semi-axes a and b.
    T ellipse(T a, T b)
    {
        T l{};
        constexpr int n = 200000;
        auto E = [&](T t) { return point<T, 3>{a * std::cos(t), b * std::sin(t), 0.}; };
        for (int k = 0; k < n; ++k)
            l += norm(E(2 * pi * (k + 1.) / n) - E(2 * pi * T(k) / n));
        return l;
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
    // a tangency along a curve is traced: a sphere in a cylinder of its radius, along the equator
    IntersectOptions<T> marched;
    marched.closed_forms = false;
    const CylindricalSurface<T> sleeve{1., frame(F.o, F.y, F.z), {0., 2 * pi, -2., 2.}};
    const auto rs = intersect(sphere, sleeve, marched);
    expect_valid(sphere, sphere.bounds(), sleeve, sleeve.bounds(), rs);
    ASSERT_FALSE(rs.singular);
    ASSERT_FALSE(rs.curves.empty());
    ASSERT_TRUE(std::ranges::all_of(rs.curves, [](const auto &c) { return c.tangent; }));
    ASSERT_NEAR(total_length(rs), 2 * pi, 1e-6);
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
    auto cut = [&](const Surface<T, 3> &s, T z, const std::array<T, 4> &box, const std::array<T, 4> &bounds) {
        const Plane<T> plane{ax2<T, 3>{{{0., 0., z}, {0., 0., 1.}, {1., 0., 0.}}}, bounds};
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
            auto r = cut(side, z, box, wide);
            REQUIRE(r.curves.size() == 1u);
            const auto &c = r.curves[0];
            CHECK(c.closed);
            CHECK(norm(c.curve->value(c.a) - c.curve->value(c.b)) <= 1e-12);
            CHECK(length(c) == doctest::Approx(2 * pi).epsilon(1e-6));
        }
    // the plane itself bounded by [-1, 1]²: the circle touches the boundary of the plane at
    // four points, where the branches end (a section has a vertex where a curve touches the
    // boundary of a face, as at the double point of Viviani's curve): four quarters, whose ends
    // are the contacts exactly. At (1, 0, z) the step also leaves the side across its seam, and
    // the landing on the seam is exact although the guess of the step is beyond x = 1, off the
    // plane; elsewhere the landing, tangent, is put on the contact
    for (T z : {T(0), T(2)})
    {
        INFO("rational side, z = " << z << ", plane bounded by the box it touches");
        auto r = cut(side, z, touching, touching);
        REQUIRE(r.curves.size() == 4u);
        CHECK(total_length(r) == doctest::Approx(2 * pi).epsilon(1e-6));
        const std::array<point<T, 3>, 4> contacts{point<T, 3>{1., 0., z}, {0., 1., z}, {-1., 0., z}, {0., -1., z}};
        for (const auto &c : r.curves)
            for (T t : {c.a, c.b})
            {
                const auto p = c.curve->value(t);
                CHECK(std::ranges::any_of(contacts, [&](const auto &x) { return norm(p - x) <= 1e-12; }));
            }
    }
    // the crossing criterion at the contact (1, 0, 0) of the touching box, on the seam: the
    // circle runs along the isos x = 1 of the plane and v = 0 of the side, and crosses y = 0
    // and the seam u = 0 at a right angle
    {
        const Plane<T> plane{ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}}, wide};
        const intersect_detail::Marcher<T> m{plane, touching, side, side.bounds(), o, T(4)};
        const auto l = m.local({1., 0., side.bounds()[0], 0.});
        REQUIRE(norm(l.p1 - l.p2) <= 1e-12);
        CHECK_FALSE(m.crosses(l, 0));
        CHECK(m.crosses(l, 1));
        CHECK(m.crosses(l, 2));
        CHECK_FALSE(m.crosses(l, 3));
        // at a pole of a sphere the isos of the sphere are not curves through the point: crossed
        const SphericalSurface<T> sphere{1., frame({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}), {0., 2 * pi, -pi / 2, pi / 2}};
        const Plane<T> meridian{ax2<T, 3>{{{0., 0., 0.}, {0., 1., 0.}, {1., 0., 0.}}}, wide};
        const intersect_detail::Marcher<T> mp{sphere, sphere.bounds(), meridian, wide, o, T(4)};
        const auto lp = mp.local({0., pi / 2, 0., -1.});
        REQUIRE(norm(lp.p1 - lp.p2) <= 1e-12);
        CHECK(mp.crosses(lp, 0));
        CHECK(mp.crosses(lp, 1));
    }
    for (T z : {T(0), T(2)})
    {
        INFO("polynomial side, z = " << z);
        auto r = cut(polynomial, z, wide, wide);
        REQUIRE(r.curves.size() == 1u);
        CHECK(r.curves[0].closed);
        CHECK(length(r.curves[0]) == doctest::Approx(loop_length).epsilon(1e-6));
        // the interpolated circle bulges out of the box [-1, 1]² between its points: cut by it
        // (the curves are checked on both surfaces and in the boxes)
        r = cut(polynomial, z, touching, wide);
        CHECK(total_length(r) < loop_length);
    }
    // beyond the bound by more than the rounding margin: the plane misses the side, no curve
    // snapped on the border
    for (T z : {T(-1e-6), T(-1e-9), T(2 + 1e-9)})
    {
        INFO("rational side, z = " << z);
        CHECK(cut(side, z, wide, wide).curves.empty());
    }
    // a border in u: the plane y = 0 through the axis cuts the side along its seam (u = 0) and
    // along the opposite line
    const Plane<T> meridian{ax2<T, 3>{{{0., 0., 0.}, {0., 1., 0.}, {1., 0., 0.}}}, {-2., 2., -3., 3.}};
    auto rm = intersect(meridian, meridian.bounds(), side, side.bounds(), o);
    expect_valid(meridian, meridian.bounds(), side, side.bounds(), rm);
    CHECK_FALSE(rm.singular);
    CHECK(total_length(rm) == doctest::Approx(4.).epsilon(1e-6));
}

TEST(tests_intersect_marching, tangent_along_a_curve)
{
    // surfaces tangent along a curve, where n1 x n2 vanishes: traced as a tangency, not singular
    IntersectOptions<T> o;
    o.tol = tol;
    o.closed_forms = false;
    // the side of a cylinder as a rational B-spline, lying on the plane y = 1 along x = 0
    {
        auto circle = build_circle<T, 3>(1.);
        points_vector<T, 4> poles;
        for (T z : {T(0), T(2)})
            for (auto p : circle.poles())
            {
                p[2] += p[3] * z;
                poles.push_back(p);
            }
        const BSSurfaceRational<T, 3> side{poles, circle.knotsFlats(), std::vector<T>{0., 0., 2., 2.}, circle.degree(), 1};
        const Plane<T> plane{ax2<T, 3>{{{0., 1., 0.}, {0., 1., 0.}, {1., 0., 0.}}}, {-2., 2., -3., 3.}};
        auto r = intersect(plane, plane.bounds(), side, side.bounds(), o);
        expect_valid(plane, plane.bounds(), side, side.bounds(), r);
        CHECK_FALSE(r.singular);
        REQUIRE(r.curves.size() == 1u);
        CHECK(r.curves[0].tangent);
        CHECK(length(r.curves[0]) == doctest::Approx(2.).epsilon(1e-6));
        for (T t : {r.curves[0].a, r.curves[0].b, (r.curves[0].a + r.curves[0].b) / 2})
        {
            const auto p = r.curves[0].curve->value(t);
            CHECK(std::abs(p[0]) <= 1e-6);
            CHECK(std::abs(p[1] - 1.) <= 1e-6);
        }
    }
    // a torus in a coaxial cylinder, tangent along the circle of radius 3
    {
        const ToroidalSurface<T> torus{2., 1., F, {0., 2 * pi, -pi / 2, 3 * pi / 2}};
        const CylindricalSurface<T> cyl{3., F, {0., 2 * pi, -2., 2.}};
        auto r = intersect(torus, cyl, o);
        expect_valid(torus, torus.bounds(), cyl, cyl.bounds(), r);
        CHECK_FALSE(r.singular);
        REQUIRE_FALSE(r.curves.empty());
        for (const auto &c : r.curves)
            CHECK(c.tangent);
        CHECK(total_length(r) == doctest::Approx(6 * pi).epsilon(1e-6));
    }
    // along a seam: the outer equator of the torus as its seam v = 0, the generator x = 1 of
    // the B-spline side as its seam u = 0
    {
        const ToroidalSurface<T> torus{2., 1., F, {0., 2 * pi, 0., 2 * pi}};
        const CylindricalSurface<T> cyl{3., F, {0., 2 * pi, -2., 2.}};
        auto r = intersect(torus, cyl, o);
        expect_valid(torus, torus.bounds(), cyl, cyl.bounds(), r);
        CHECK_FALSE(r.singular);
        for (const auto &c : r.curves)
            CHECK(c.tangent);
        CHECK(total_length(r) == doctest::Approx(6 * pi).epsilon(1e-6));
    }
    {
        auto circle = build_circle<T, 3>(1.);
        points_vector<T, 4> poles;
        for (T z : {T(0), T(2)})
            for (auto p : circle.poles())
            {
                p[2] += p[3] * z;
                poles.push_back(p);
            }
        const BSSurfaceRational<T, 3> side{poles, circle.knotsFlats(), std::vector<T>{0., 0., 2., 2.}, circle.degree(), 1};
        const Plane<T> plane{ax2<T, 3>{{{1., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}}}, {-2., 2., -3., 3.}};
        auto r = intersect(plane, plane.bounds(), side, side.bounds(), o);
        expect_valid(plane, plane.bounds(), side, side.bounds(), r);
        CHECK_FALSE(r.singular);
        for (const auto &c : r.curves)
            CHECK(c.tangent);
        MESSAGE("seam side: " << r.curves.size() << " curves");
        CHECK(total_length(r) >= 2. - 1e-6);
    }
}

TEST(tests_intersect_marching, tangency_limits)
{
    IntersectOptions<T> o;
    o.tol = tol;
    o.closed_forms = false;
    // a tangency that ends inside the boxes: the cylinder of radius 1 along x lying on the face
    // z = -bend max(x - a, 0)^3 (a cubic B-spline in x with a knot at a, exact). Tangent for
    // x <= a; beyond, the surfaces stop bending alike, and the tangency ends there, x = a +
    // angular_tol / 6 whatever the tolerance. With bend = 1 the face bends away: the surfaces
    // part, the branch ends there and is not singular. With bend = -1 the face rises into the
    // cylinder: the surfaces cross, two transversal curves start there, the branch is singular.
    // Both orders of the surfaces (the gap and the bend are taken on the normal of the second)
    auto face_under = [](T bend, T a) {
        const std::vector<T> kx{-1., -1., -1., -1., a, 1., 1., 1., 1.}, ky{-1., -1., 1., 1.};
        const std::array<T, 5> gx{-1., (a - 2.) / 3., a / 3., (a + 2.) / 3., 1.};
        const std::array<T, 5> zx{0., 0., 0., 0., -bend * std::pow(1. - a, 3)};
        points_vector<T, 3> poles;
        for (T y : {T(-1), T(1)})
            for (int i = 0; i < 5; ++i)
                poles.push_back({gx[i], y, zx[i]});
        return BSSurface<T, 3>{poles, kx, ky, 3, 1};
    };
    const CylindricalSurface<T> cyl{1., frame({0., 0., 1.}, {1., 0., 0.}, {0., 0., 1.}), {0., 2 * pi, -1., 1.}};
    for (T bend : {T(1), T(-1)})
        for (bool face_first : {true, false})
            for (T tolerance : {T(1e-4), T(1e-6), T(1e-8)})
            {
                INFO("bend " << bend << ", face first " << face_first << ", tol " << tolerance);
                const auto face = face_under(bend, 0.);
                REQUIRE(std::abs(face.value(0.5, 0.)[2] + bend * 0.125) <= 1e-14);
                REQUIRE(std::abs(face.value(-0.5, 0.3)[2]) <= 1e-14);
                IntersectOptions<T> ot = o;
                ot.tol = tolerance;
                const Surface<T, 3> &s1 = face_first ? static_cast<const Surface<T, 3> &>(face) : cyl;
                const Surface<T, 3> &s2 = face_first ? static_cast<const Surface<T, 3> &>(cyl) : face;
                auto r = intersect(s1, s1.bounds(), s2, s2.bounds(), ot);
                CHECK(r.singular == (bend < 0));
                const auto tangent = std::ranges::find_if(r.curves, [](const auto &c) { return c.tangent; });
                REQUIRE(tangent != r.curves.end());
                if (bend > 0)
                    CHECK(r.curves.size() == 1u);
                const auto a = tangent->curve->value(tangent->a), b = tangent->curve->value(tangent->b);
                CHECK(std::min(a[0], b[0]) == doctest::Approx(-1.).epsilon(1e-9));
                // the end: x = 1e-3 / 6, within two minimal steps (10 tol each: the last step
                // rejected is shorter than two)
                const T x_end = std::max(a[0], b[0]);
                CHECK(x_end >= 1e-3 / 6. - 20. * tolerance - 1e-6);
                CHECK(x_end <= 1e-3 / 6. + 20. * tolerance + 1e-6);
                // what follows the end, looking further along x
                const T te = a[0] > b[0] ? tangent->a : tangent->b;
                const auto q1 = tangent->pc1->value(te), q2 = tangent->pc2->value(te);
                const intersect_detail::Marcher<T> m{s1, s1.bounds(), s2, s2.bounds(), ot, T(1)};
                const intersect_detail::TangentMarcher<T> tm{m, s1, s2};
                const intersect_detail::SSPoint<T> end{tangent->curve->value(te), {q1[0], q1[1], q2[0], q2[1]}};
                using Beyond = intersect_detail::TangentMarcher<T>::Beyond;
                CHECK(tm.parting(end, {1., 0., 0.}).first == (bend > 0 ? Beyond::parts : Beyond::crosses));
            }
    // the face bending away 0.005 before its border x = 1: the gap (0.005^3) stays within the
    // tolerance up to the border, where the branch lands, not singular
    {
        const auto face = face_under(1., 0.995);
        auto r = intersect(face, face.bounds(), cyl, cyl.bounds(), o);
        expect_valid(face, face.bounds(), cyl, cyl.bounds(), r);
        CHECK_FALSE(r.singular);
        REQUIRE(r.curves.size() == 1u);
        const auto a = r.curves[0].curve->value(r.curves[0].a), b = r.curves[0].curve->value(r.curves[0].b);
        CHECK(std::max(a[0], b[0]) >= 1. - tol); // (a landing within the tolerance of the last point is not added)
        // (a step lands over the border there; from just past the end of the tangency, the probe
        // leaves the boxes before the gap reaches the tolerance: the boundary, to land on)
        const intersect_detail::Marcher<T> m{face, face.bounds(), cyl, cyl.bounds(), o, T(1)};
        const intersect_detail::TangentMarcher<T> tm{m, face, cyl};
        const T x0 = 0.995 + 2e-4;
        const intersect_detail::SSPoint<T> end{{x0, 0., 0.}, {x0, 0., pi, x0}};
        REQUIRE(norm(cyl.value(pi, x0) - end.p) <= 1e-12);
        REQUIRE(norm(face.value(x0, 0.) - end.p) <= 1e-10);
        const auto [beyond, off] = tm.parting(end, {1., 0., 0.});
        CHECK(beyond == intersect_detail::TangentMarcher<T>::Beyond::boundary);
        REQUIRE(off);
        CHECK_FALSE(m.inside(off->q));
    }
    // starts past a seam that a closed tangency crosses: coming back, the branch meets its start
    // across the seam; it lands on the seam, its parameters monotone, and is not closed back
    // across it (the parameters would jump there by 2 pi)
    {
        const ToroidalSurface<T> torus{2., 1., F, {0., 2 * pi, -pi / 2, 3 * pi / 2}};
        const CylindricalSurface<T> cyl{3., F, {0., 2 * pi, -2., 2.}};
        const intersect_detail::Marcher<T> m{torus, torus.bounds(), cyl, cyl.bounds(), o, T(1)};
        const intersect_detail::TangentMarcher<T> tm{m, torus, cyl};
        for (int i = 1; i <= 20; ++i)
        {
            const T u0 = 0.005 * i;
            INFO("start at u = " << u0);
            const intersect_detail::SSPoint<T> start{cyl.value(u0, 0.), {u0, 0., u0, 0.}};
            REQUIRE(norm(torus.value(u0, 0.) - start.p) <= 1e-12);
            std::vector<intersect_detail::SSPoint<T>> longest;
            for (T sign : {T(1), T(-1)})
            {
                intersect_detail::Marcher<T>::End e{};
                auto pts = tm.march(start, sign, e);
                CHECK(e != intersect_detail::Marcher<T>::End::closed);
                if (pts.size() > longest.size())
                    longest = std::move(pts);
            }
            REQUIRE(longest.size() > 10u);
            for (std::size_t k = 1; k < longest.size(); ++k)
                for (int j = 0; j < 4; ++j)
                    CHECK(std::abs(longest[k].q[j] - longest[k - 1].q[j]) < 1.);
            CHECK(std::abs(longest.back().q[2] - 2 * pi) <= 1e-9);
        }
    }
    // a surface without second derivatives (an offset): no tangency tracing, the intersection
    // stays singular as before, without throwing
    {
        const auto ball = std::make_shared<SphericalSurface<T>>(0.5, F, std::array<T, 4>{0., 2 * pi, -pi / 2, pi / 2});
        const brep::OffsetSurface<T> sphere{ball, 0.5};
        REQUIRE_FALSE(sphere.has_second_derivatives());
        const CylindricalSurface<T> sleeve{1., F, {0., 2 * pi, -2., 2.}};
        SurfaceSurfaceResult<T> r;
        REQUIRE_NOTHROW(r = intersect(sphere, sphere.bounds(), sleeve, sleeve.bounds(), o));
        CHECK(r.singular);
    }
    // an isolated contact, a sphere on a plane: certified as such by the second forms, not traced;
    // the marching reports it as before, a contact point (and singular, its start being tangent)
    {
        const SphericalSurface<T> sphere{1., frame({0., 0., 1.}, {0., 0., 1.}, {1., 0., 0.}), {0., 2 * pi, -pi / 2, pi / 2}};
        const Plane<T> plane{ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}}, {-2., 2., -2., 2.}};
        auto r = intersect(plane, plane.bounds(), sphere, sphere.bounds(), o);
        CHECK(r.curves.empty());
        REQUIRE(r.points.size() == 1u);
        CHECK(norm(r.points[0].p) <= 1e-6);
    }
}

TEST(tests_intersect_marching, start_on_a_border)
{
    // the top z = -(max(x, 0)^3 - 16 max(x - 0.5, 0)^3) over x in 0.5..1 against the vertical
    // face y = -0.8, x in -0.8..0.9, z in 0..1: they cross from (x0, -0.8, 0), x0 = 0.5 c / (c - 1),
    // c = 16^(1/3), on the bottom border of the face, to (0.9, -0.8, 0.295) on its side border. The
    // start from the centres of the patches meets the face beyond its borders: the correction
    // slides along them instead of stalling there
    const std::vector<T> kx{-1., -1., -1., -1., 0., 0.5, 1., 1., 1., 1.}, k1{0., 0., 1., 1.};
    const std::array<T, 6> gx{-1., -2. / 3., -1. / 6., 0.5, 5. / 6., 1.}, zx{0., 0., 0., 0., -0.5, 1.};
    points_vector<T, 3> poles;
    for (T y : {T(-1), T(1)})
        for (int i = 0; i < 6; ++i)
            poles.push_back({gx[i], y, zx[i]});
    const BSSurface<T, 3> top{poles, kx, k1, 3, 1};
    const std::vector<T> k{0., 0., 1., 1.};
    const BSSurface<T, 3> side{points_vector<T, 3>{{-0.8, -0.8, 0.}, {0.9, -0.8, 0.}, {-0.8, -0.8, 1.}, {0.9, -0.8, 1.}}, k, k, 1, 1};
    const std::array<T, 4> box{0.5, 1., 0., 1.};
    IntersectOptions<T> o;
    o.tol = tol;
    auto r = intersect(top, box, side, side.bounds(), o);
    expect_valid(top, box, side, side.bounds(), r);
    CHECK_FALSE(r.singular);
    REQUIRE(r.curves.size() == 1u);
    const T c = std::cbrt(16.), x0 = 0.5 * c / (c - 1.);
    auto a = r.curves[0].curve->value(r.curves[0].a), b = r.curves[0].curve->value(r.curves[0].b);
    if (a[0] > b[0])
        std::swap(a, b);
    CHECK(norm(a - point<T, 3>{x0, -0.8, 0.}) <= 1e-6);
    CHECK(norm(b - point<T, 3>{0.9, -0.8, -(std::pow(0.9, 3) - 16. * std::pow(0.4, 3))}) <= 1e-6);
}

// ---- review of the marching (lot 2 of the review of the BREP kernel)

TEST(tests_intersect_marching, thin_loop)
{
    // the plane z = 0.5 under the ridge z = 1 - x'^2 / a^2 - y'^2 / b^2 (exact), x' along the
    // diagonal of the boxes: an ellipse of semi-axes a / sqrt 2 and b / sqrt 2, closed once and
    // not across its width, where the start is within a step but on the other side
    const T a = 1.5;
    const Plane<T> plane{ax2<T, 3>{{{0., 0., 0.5}, {0., 0., 1.}, {1., 0., 0.}}}, {-3., 3., -3., 3.}};
    for (T b : {0.4, 0.2, 0.1, 0.05})
    {
        INFO("b = " << b);
        const T A = 1 / (a * a), B = 1 / (b * b);
        const auto ridge = quadratic(-2., 2., -2., 2., {1., 0., 0., -(A + B) / 2, B - A, -(A + B) / 2});
        const auto r = marching(ridge, plane);
        REQUIRE(r.curves.size() == 1u);
        CHECK(r.curves[0].closed);
        CHECK(fine_length(r.curves[0]) == doctest::Approx(ellipse(a / std::sqrt(2.), b / std::sqrt(2.))).epsilon(1e-6));
    }
}

TEST(tests_intersect_marching, closing_near_an_extremum)
{
    // the paraboloid z = 1 - x^2 - y^2 (exact) and the plane z = 0.5: from starts near the
    // extremum of x (and u1) at theta = 0, the step past the start hardly changes u1; the loop is
    // closed once, not marched around again
    const auto bump = quadratic(-1.5, 1.5, -1.5, 1.5, {1., 0., 0., -1., 0., -1.});
    const Plane<T> plane{ax2<T, 3>{{{0., 0., 0.5}, {0., 0., 1.}, {1., 0., 0.}}}, {-3., 3., -3., 3.}};
    const intersect_detail::Marcher<T> m{bump, bump.bounds(), plane, plane.bounds(), {}, T(6)};
    const T R = std::sqrt(0.5);
    for (int k = -100; k <= 100; ++k)
    {
        const T th = k * 2e-4;
        const point<T, 3> P{R * std::cos(th), R * std::sin(th), 0.5};
        const auto start = m.correct({(P[0] + 1.5) / 3, (P[1] + 1.5) / 3, P[0], P[1]}, nullptr, nullptr);
        REQUIRE(start);
        for (T sign : {T(1), T(-1)})
        {
            INFO("start at theta " << th << ", direction " << sign);
            intersect_detail::Marcher<T>::End e{};
            std::vector<std::pair<intersect_detail::SSPoint<T>, intersect_detail::SSPoint<T>>> jumps;
            const auto pts = m.march(*start, sign, e, jumps);
            T polyline{};
            for (std::size_t j = 1; j < pts.size(); ++j)
                polyline += norm(pts[j].p - pts[j - 1].p);
            CHECK(e == intersect_detail::Marcher<T>::End::closed);
            CHECK(polyline < 1.01 * 2 * pi * R);
        }
    }
    // through the intersection: paraboloids z = h0 - (x - dx)^2 - (y - dy)^2 marched two or
    // three times around
    for (const auto &[dx, dy, h0] : {std::array{-0.3 + 0.6 * 3 / 39, -0.3, 0.8 + 0.4 * 21 / 39.}, std::array{-0.3 + 0.6 * 6 / 39, -0.3 + 0.6 * 12 / 39, 0.8 + 0.4 * 18 / 39.},
                                     std::array{-0.3 + 0.6 * 6 / 39, -0.3 + 0.6 * 26 / 39, 0.8 + 0.4 * 0 / 39.}})
    {
        INFO("dx " << dx << ", dy " << dy << ", h0 " << h0);
        const auto shifted = quadratic(-1.5, 1.5, -1.5, 1.5, {h0 - dx * dx - dy * dy, 2 * dx, 2 * dy, -1., 0., -1.});
        const auto r = marching(shifted, plane);
        REQUIRE(r.curves.size() == 1u);
        CHECK(length(r.curves[0]) == doctest::Approx(2 * pi * std::sqrt(h0 - 0.5)).epsilon(1e-6));
    }
}

TEST(tests_intersect_marching, close_parallel_curves)
{
    // the plane z = 0 across the valley z = (x^2 - a^2) / a (exact): the lines x = +-a, 2a apart.
    // On straight curves the step grows to a tenth of the size; the second line is no longer
    // taken as traced within 1 % of a step
    const Plane<T> plane{ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}}, {-3., 3., -3., 3.}};
    for (T a : {2e-3, 1e-3, 5e-4})
    {
        INFO("a = " << a);
        const auto valley = quadratic(-2., 2., -2., 2., {-a, 0., 0., 1 / a, 0., 0.});
        const auto r = marching(valley, plane);
        CHECK(r.curves.size() == 2u);
        CHECK(total_length(r) == doctest::Approx(8.).epsilon(1e-9));
    }
}

TEST(tests_intersect_marching, small_loop_any_parametrization)
{
    // a cubic sheet z = 0 on [-20, 20]^2 (81 x 81 poles) with one pole lifted by 1, against the
    // plane z = 0.2 on the same square: a small loop around the bump, found whatever the
    // parameters of the sheet (the start from the centres of a patch of the bump and of the whole
    // plane moved along the surface of the larger parameters, off the bump), and closed once
    const T L = 20.;
    const Plane<T> plane{ax2<T, 3>{{{0., 0., 0.2}, {0., 0., 1.}, {1., 0., 0.}}}, {-L, L, -L, L}};
    T reference{};
    for (const auto &[a, b] : {std::pair{-20., 20.}, std::pair{0., 1.}, std::pair{0., 1000.}})
        for (T c : {2., 5., 10., 15.})
        {
            INFO("parameters [" << a << ", " << b << "], bump at (" << c << ", " << c << ")");
            const auto s = sheet(L, 81, c, c, 1., a, b);
            const auto r = marching(s, plane);
            REQUIRE(r.curves.size() == 1u);
            CHECK(r.curves[0].closed);
            if (reference == 0.)
                reference = length(r.curves[0]);
            CHECK(length(r.curves[0]) == doctest::Approx(reference).epsilon(1e-6));
        }
}

TEST(tests_intersect_marching, narrow_detour)
{
    // the slope z = 0.3 x (cubic, poles 0.05 apart) with the pole at (0, y0) lifted by H, cut by
    // z = 0: the line x = 0 goes round the bump, narrower than the steps. The steps went across
    // (a curve 8e-2 off at H = 0.3), and a start on the detour traced the line again
    const Plane<T> plane{ax2<T, 3>{{{0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}}}, {-3., 3., -3., 3.}};
    for (const auto &[H, y0] : {std::pair{0.3, -0.81}, std::pair{0.01, 0.}, std::pair{0.01, 0.37}, std::pair{0.1, 0.37}})
    {
        INFO("H = " << H << ", y0 = " << y0);
        const int n = 81;
        points_vector<T, 3> poles;
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
            {
                const T x = -2. + 4. * i / (n - 1), y = -2. + 4. * j / (n - 1);
                poles.push_back({x, y, 0.3 * x + (i == n / 2 && std::abs(y - y0) < 0.025 ? H : 0.)});
            }
        std::vector<T> k{0., 0., 0., 0.};
        for (int s = 1; s < n - 3; ++s)
            k.push_back(T(s) / (n - 3));
        k.insert(k.end(), {1., 1., 1., 1.});
        const BSSurface<T, 3> slope{poles, k, k, 3, 3};
        const auto r = marching(slope, plane);
        CHECK(r.curves.size() == 1u);
    }
}

TEST(tests_intersect_marching, close_parallel_faces)
{
    // faces parallel and close, without intersection: the flat patches are apart (the poles of one
    // beyond the slab of the other), not halved down to the cap on the pairs
    const T c = std::sqrt(0.5);
    const point<T, 3> n{0., -c, c}, e1{1., 0., 0.}, e2{0., c, c};
    auto square = [&](T side, T off) {
        points_vector<T, 3> poles;
        for (int j = 0; j < 2; ++j)
            for (int i = 0; i < 2; ++i)
                poles.push_back(side * (T(i) * e1 + T(j) * e2) + off * n);
        return BSSurface<T, 3>{poles, std::vector<T>{0., 0., 1., 1.}, std::vector<T>{0., 0., 1., 1.}, 1, 1};
    };
    for (const auto &[side, d] : {std::pair{1., 1e-3}, std::pair{1., 1e-4}, std::pair{100., 0.1}, std::pair{100., 0.01}})
    {
        INFO("squares of side " << side << ", " << d << " apart");
        const auto a = square(side, 0.), b = square(side, d);
        const auto r = intersect(a, b);
        CHECK(r.complete);
        CHECK(r.curves.empty());
    }
    // concentric NURBS cylinders (refused as NonGeneric by the booleans)
    for (T d : {1e-2, 1e-3})
    {
        INFO("radii 1 and " << 1 + d);
        const CylindricalSurface<T> c1{1., Frame3<T>{}, {0., 2 * pi, 0., 2.}}, c2{1. + d, Frame3<T>{}, {0., 2 * pi, 0., 2.}};
        const auto n1 = c1.to_nurbs(), n2 = c2.to_nurbs();
        const auto r = intersect(*n1, *n2);
        CHECK(r.complete);
        CHECK_FALSE(r.singular);
        CHECK(r.curves.empty());
    }
}

TEST(tests_intersect_marching, far_from_the_origin)
{
    // the paraboloid z = 1 - x^2 - y^2 (exact) and the plane z = 0.5 moved by (D, D, D): the
    // corrector stopped on a move below 1e-4 tol, under the rounding of the coordinates far from
    // the origin (300 m in mm at 1e-7)
    const T reference = [] {
        const auto bump = quadratic(-1.5, 1.5, -1.5, 1.5, {1., 0., 0., -1., 0., -1.});
        const Plane<T> plane{ax2<T, 3>{{{0., 0., 0.5}, {0., 0., 1.}, {1., 0., 0.}}}, {-3., 3., -3., 3.}};
        return total_length(marching(bump, plane));
    }();
    for (const auto &[D, tolerance] : {std::pair{3e5, 1e-7}, std::pair{1e6, 1e-7}, std::pair{1e7, 1e-6}})
    {
        INFO("D = " << D << ", tol " << tolerance);
        const point<T, 3> off{D, D, D};
        const auto bump = quadratic(-1.5, 1.5, -1.5, 1.5, {1., 0., 0., -1., 0., -1.}, off);
        const Plane<T> plane{ax2<T, 3>{{off + point<T, 3>{0., 0., 0.5}, {0., 0., 1.}, {1., 0., 0.}}}, {-3., 3., -3., 3.}};
        IntersectOptions<T> o;
        o.tol = tolerance;
        const auto r = marching(bump, plane, o);
        CHECK_FALSE(r.singular);
        REQUIRE(r.curves.size() == 1u);
        CHECK(r.curves[0].closed);
        CHECK(r.curves[0].branch_spans.empty());
        CHECK(length(r.curves[0]) == doctest::Approx(reference).epsilon(1e-6));
    }
}

TEST(tests_intersect_marching, small_circles)
{
    // a sphere of radius R cut through its centre by a plane, by marching: a circle of a few dozen
    // to a few hundred tol, which the steps of 10 tol at least turning by 0.05 at most could not
    // follow
    for (T R : {3e-4, 1e-4, 5e-5})
    {
        INFO("R = " << R);
        const SphericalSurface<T> ball{R, frame({0.3, 0.2, 0.1}, {1., 0., 0.}, {0., 1., 0.}), {0., 2 * pi, -pi / 2, pi / 2}};
        const Plane<T> plane{frame({0.3, 0.2, 0.1}, {0., 0., 1.}, {1., 0., 0.}), {-1., 1., -1., 1.}};
        const auto r = marching(ball, plane);
        CHECK_FALSE(r.singular);
        CHECK_FALSE(r.curves.empty());
        CHECK(total_length(r) == doctest::Approx(2 * pi * R).epsilon(1e-5));
    }
}

TEST(tests_intersect_marching, piece_beyond_a_seam)
{
    // planes across cylinders whose ellipse is cut by a bound of the cylinder just past its seam:
    // the short piece between the seam and the bound starts at the end of the long piece in 3D,
    // on the other side of the seam; it was taken as traced (two random cases of the review)
    struct Case
    {
        T R;
        point<T, 3> o, z, x, po, pz, px;
    };
    const std::array cases{
        Case{2.9202658626484581,
             {6.0281340806674022, -40.352943810561051, -72.534496655845928},
             {-0.5202466751629129, 0.83617645384099604, 0.17364427725620987},
             {-0.36777810562166313, -0.40286850313296363, 0.83811469036686492},
             {6.749309094660461, -41.008264002296748, -73.137785691569462},
             {0.55861176321432515, -0.17315985835352177, 0.81115261292408536},
             {0.60675420805403724, -0.5814708914803538, -0.54197872040308093}},
        Case{1.8011327128878081,
             {-94.071415020922771, 32.852986746093848, -76.289098657690047},
             {-0.55803040188844388, -0.82959357395229028, -0.019405479259423278},
             {0.73530742228552692, -0.50517654432356551, 0.45179603118787387},
             {-95.131333101781166, 33.749943030911055, -76.985350024572469},
             {-0.18708015601340414, 0.72135838139259967, -0.66681564230351487},
             {0.69041207565775176, 0.5794276541080895, 0.43312210684834718}}};
    for (const auto &c : cases)
    {
        INFO("cylinder of radius " << c.R);
        const CylindricalSurface<T> cyl{c.R, frame(c.o, c.z, c.x), {0., 2 * pi, -6., 6.}};
        const Plane<T> plane{frame(c.po, c.pz, c.px), {-40., 40., -40., 40.}};
        const auto exact = intersect(cyl, plane);
        REQUIRE(exact.complete);
        const auto r = marching(cyl, plane);
        CHECK(r.curves.size() == exact.curves.size());
        CHECK(total_length(r) == doctest::Approx(total_length(exact)).epsilon(1e-6));
    }
}

TEST(tests_intersect_marching, landing_across_an_apex)
{
    // a plane 0.055 from the apex of a cone (v in -2..4, apex at v = -1.5): a hyperbola, a branch on
    // each nappe. The last step of a branch, landing on the bound v = -2, went across the apex
    // onto the other nappe (a curve 8e-2 off, and pieces traced again: 40.2 instead of 13.5); it
    // obeys the same checks as the other steps (a random case of the review)
    const ConicalSurface<T> cone{0.86594183059086172, 0.52309424731992438,
                                 frame({57.893161731294462, 35.104250420039307, -58.462922596045765}, {-0.55562531041052332, 0.58862681313021192, -0.58719586961709824},
                                       {0.14633782693219236, -0.62598570712439461, -0.76598115830922076}),
                                 {0., 2 * pi, -2., 4.}};
    const Plane<T> plane{frame({58.952531118033164, 33.990303218965884, -57.327049568804213}, {0.74977162241212125, 0.35275360344833578, -0.55982801776948432},
                               {0.59041141209958714, 0.025348955782622435, 0.80670427971054071}),
                         {-40., 40., -40., 40.}};
    // reference: on the generatrix u the plane is met at v(u), linear in v; the arcs where v is
    // within the bounds, their ends found by bisection
    const auto &Fc = cone.frame();
    const auto &Fp = plane.frame();
    const T ta = std::tan(cone.semi_angle());
    auto v_of = [&](T u) {
        const auto rad = std::cos(u) * Fc.x + std::sin(u) * Fc.y;
        return -((Fc.o - Fp.o) * Fp.z + cone.radius() * (rad * Fp.z)) / (ta * (rad * Fp.z) + Fc.z * Fp.z);
    };
    auto at = [&](T u) { return cone.value(u, v_of(u)); };
    auto in = [&](T u) { return v_of(u) >= -2. && v_of(u) <= 4.; };
    T reference{};
    constexpr int n = 200000;
    for (int k = 0; k < n; ++k)
    {
        T a = 2 * pi * k / n, b = 2 * pi * (k + 1) / n;
        if (!in(a) && !in(b))
            continue;
        if (in(a) != in(b))
        {
            T x = a, y = b;
            for (int i = 0; i < 100; ++i)
                (in((x + y) / 2) == in(x) ? x : y) = (x + y) / 2;
            (in(a) ? b : a) = in(a) ? x : y;
        }
        for (int j = 0; j < 20; ++j)
            reference += norm(at(a + (b - a) * (j + 1) / 20) - at(a + (b - a) * j / 20));
    }
    const auto r = marching(cone, plane);
    CHECK_FALSE(r.singular);
    T total{};
    for (const auto &c : r.curves)
        total += fine_length(c);
    CHECK(total == doctest::Approx(reference).epsilon(1e-6));
}

TEST(tests_intersect_marching, tangent_sleeve_perturbed)
{
    // the sphere in a cylinder of its radius of `singular_points_and_limits`, its frames, origins
    // and radii moved by 1e-16 to 1e-12 and a few ulps (rounding differs across platforms): the
    // start search on the tangency wandered on rounding errors without converging (pairs halved up
    // to the cap), the landing of the tangency on the seam had a margin below its sag (the tangency
    // traced twice), and a start at the corner of both seams, off the boxes once put on the
    // tangency, was an isolated contact (singular, as on Windows)
    IntersectOptions<T> marched;
    marched.closed_forms = false;
    auto check = [&](const SphericalSurface<T> &sphere, const CylindricalSurface<T> &sleeve) {
        const auto r = intersect(sphere, sleeve, marched);
        CHECK(r.complete);
        CHECK_FALSE(r.singular);
        CHECK(r.points.empty());
        CHECK(std::ranges::all_of(r.curves, [](const auto &c) { return c.tangent; }));
        CHECK(total_length(r) == doctest::Approx(2 * pi * sphere.radius()).epsilon(1e-6));
    };
    // perturbations found failing, exact (frames as given, not normalized again)
    struct Placed
    {
        T r;
        std::array<point<T, 3>, 4> f; // o, x, y, z
    };
    auto exact = [](const std::array<point<T, 3>, 4> &f) {
        Frame3<T> fr;
        fr.o = f[0];
        fr.x = f[1];
        fr.y = f[2];
        fr.z = f[3];
        return fr;
    };
    const std::array<std::pair<Placed, Placed>, 5> cases{{
        {{0x1.0000000000005p+0,
          {{{0x1.ffffffffffffep-1, -0x1.fffffffffffffp+0, 0x1.0000000000002p-1},
            {0x1.279a74590331dp-1, 0x1.279a74590331dp-1, 0x1.279a74590331cp-1},
            {0x1.a20bd700c2c4p-1, -0x1.a20bd700c2c41p-2, -0x1.a20bd700c2c3ep-2},
            {0x1.36d4664105862p-54, 0x1.6a09e667f3bcdp-1, -0x1.6a09e667f3bcep-1}}}},
         {0x1.ffffffffffff6p-1,
          {{{0x1.0000000000001p+0, -0x1p+1, 0x1.ffffffffffffcp-2},
            {0x1.279a74590331cp-1, 0x1.279a74590331cp-1, 0x1.279a74590331ep-1},
            {0x1.a20bd700c2c4p-1, -0x1.a20bd700c2c3bp-2, -0x1.a20bd700c2c41p-2},
            {-0x1.65d015229bec3p-54, 0x1.6a09e667f3bcfp-1, -0x1.6a09e667f3bcbp-1}}}}},
        {{0x1.ffffffffffffep-1,
          {{{0x1.ffffffffffffp-1, -0x1.ffffffffffff6p+0, 0x1.000000000000ap-1},
            {0x1.279a74590337dp-1, 0x1.279a7459032c7p-1, 0x1.279a745903315p-1},
            {0x1.a20bd700c2bfcp-1, -0x1.a20bd700c2d58p-2, -0x1.a20bd700c2c35p-2},
            {0x1.d0d1ad0a54d79p-47, 0x1.6a09e667f3bc2p-1, -0x1.6a09e667f3bd6p-1}}}},
         {0x1.ffffffffffffep-1,
          {{{0x1.fffffffffffe2p-1, -0x1.0000000000004p+1, 0x1.000000000000ep-1},
            {0x1.279a745903346p-1, 0x1.279a7459032ep-1, 0x1.279a74590333p-1},
            {0x1.a20bd700c2c22p-1, -0x1.a20bd700c2d08p-2, -0x1.a20bd700c2bebp-2},
            {0x1.ccbdb26ebf41ep-47, 0x1.6a09e667f3bc5p-1, -0x1.6a09e667f3bd5p-1}}}}},
        {{0x1.ffffffffffff8p-1,
          {{{0x1.000000000044cp+0, -0x1.ffffffffffc0ep+0, 0x1.0000000000bc5p-1},
            {0x1.279a745903f6cp-1, 0x1.279a745902f85p-1, 0x1.279a745902a65p-1},
            {0x1.a20bd700c238bp-1, -0x1.a20bd700c4115p-2, -0x1.a20bd700c3a3cp-2},
            {-0x1.d78d6e71d3422p-49, 0x1.6a09e667f38b9p-1, -0x1.6a09e667f3ee2p-1}}}},
         {0x1.0000000000007p+0,
          {{{0x1.ffffffffffba4p-1, -0x1.000000000035dp+1, 0x1.00000000007eap-1},
            {0x1.279a7459038d4p-1, 0x1.279a745903011p-1, 0x1.279a745903072p-1},
            {0x1.a20bd700c2835p-1, -0x1.a20bd700c3e68p-2, -0x1.a20bd700c2a41p-2},
            {0x1.7e247481a59f9p-43, 0x1.6a09e667f390cp-1, -0x1.6a09e667f3e8ep-1}}}}},
        {{0x1.0000000000002p+0,
          {{{0x1p+0, -0x1.ffffffffffffep+0, 0x1.0000000000003p-1},
            {0x1.279a74590331dp-1, 0x1.279a74590331fp-1, 0x1.279a745903319p-1},
            {0x1.a20bd700c2c3ep-1, -0x1.a20bd700c2c42p-2, -0x1.a20bd700c2c3dp-2},
            {-0x1.be132a9db19f2p-54, 0x1.6a09e667f3bcap-1, -0x1.6a09e667f3bdp-1}}}},
         {0x1.0000000000003p+0,
          {{{0x1.ffffffffffffdp-1, -0x1.0000000000001p+1, 0x1.ffffffffffffdp-2},
            {0x1.279a74590331bp-1, 0x1.279a74590331cp-1, 0x1.279a74590331ep-1},
            {0x1.a20bd700c2c4p-1, -0x1.a20bd700c2c3dp-2, -0x1.a20bd700c2c3dp-2},
            {0x1.0830ec466b5acp-54, 0x1.6a09e667f3bcep-1, -0x1.6a09e667f3bcdp-1}}}}},
        {{0x1.ffffffffffffap-1,
          {{{0x1.0000000000002p+0, -0x1.0000000000002p+1, 0x1.0000000000001p-1},
            {0x1.279a74590331ep-1, 0x1.279a745903322p-1, 0x1.279a745903318p-1},
            {0x1.a20bd700c2c4p-1, -0x1.a20bd700c2c4cp-2, -0x1.a20bd700c2c36p-2},
            {0x1.47c791159a518p-52, 0x1.6a09e667f3bc6p-1, -0x1.6a09e667f3bd4p-1}}}},
         {0x1.0000000000005p+0,
          {{{0x1.0000000000001p+0, -0x1.ffffffffffffbp+0, 0x1.0000000000001p-1},
            {0x1.279a745903319p-1, 0x1.279a745903324p-1, 0x1.279a745903319p-1},
            {0x1.a20bd700c2c42p-1, -0x1.a20bd700c2c3cp-2, -0x1.a20bd700c2c39p-2},
            {-0x1.87febf90afc1ap-52, 0x1.6a09e667f3bc9p-1, -0x1.6a09e667f3bd2p-1}}}}},
    }};
    for (std::size_t k = 0; k < cases.size(); ++k)
    {
        INFO("case " << k);
        const auto &[a, b] = cases[k];
        check(SphericalSurface<T>{a.r, exact(a.f), {0., 2 * pi, -pi / 2, pi / 2}}, CylindricalSurface<T>{b.r, exact(b.f), {0., 2 * pi, -2., 2.}});
    }
    // and a sample of others
    std::uint64_t state = 0x9E3779B97F4A7C15ull;
    auto next = [&] {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        return T(state >> 11) * 0x1p-52 - 1; // in [-1, 1)
    };
    constexpr T e = std::numeric_limits<T>::epsilon();
    for (int k = 0; k < 40; ++k)
    {
        const T s = std::pow(10., -14. + 2. * next()); // 1e-16 to 1e-12
        auto tiny = [&] { return s * point<T, 3>{next(), next(), next()}; };
        const T r1 = 1. + std::round(8. * next()) * e, r2 = 1. + std::round(8. * next()) * e;
        INFO("perturbation " << k << ", scale " << s);
        check(SphericalSurface<T>{r1, frame(F.o + tiny(), F.y + tiny(), F.z + tiny()), {0., 2 * pi, -pi / 2, pi / 2}},
              CylindricalSurface<T>{r2, frame(F.o + tiny(), F.y + tiny(), F.z + tiny()), {0., 2 * pi, -2., 2.}});
    }
}

TEST(tests_intersect_marching, closing_on_a_step_end)
{
    // the circle of the cylinder x^2 + y^2 = 1 and the plane z = 0.5, stepped by 0.04 rad: a start
    // on the end of a step, or a rounding error beyond or before it, is passed by one of the steps
    // around it, and by none further on
    const CylindricalSurface<T> cyl{1., Frame3<T>{}, {0., 2 * pi, -1., 1.}};
    const Plane<T> plane{ax2<T, 3>{{{0., 0., 0.5}, {0., 0., 1.}, {1., 0., 0.}}}, {-2., 2., -2., 2.}};
    const intersect_detail::Marcher<T> m{cyl, cyl.bounds(), plane, plane.bounds(), {}, T(1)};
    auto at = [](T a) { return intersect_detail::SSPoint<T>{{std::cos(a), std::sin(a), 0.5}, {a, 0.5, std::cos(a), std::sin(a)}}; };
    auto tangent = [](T a) { return point<T, 3>{-std::sin(a), std::cos(a), 0.}; };
    const T th = 0.04;
    for (T a0 : {0.5, 1., 2., 3.})
        for (T d : {0., 1e-15, -1e-15, 1e-12, -1e-12, 1e-9, -1e-9})
        {
            INFO("start at " << a0 << " + " << d);
            const auto start = at(a0 + d);
            REQUIRE(norm(cyl.value(start.q[0], start.q[1]) - start.p) <= 1e-15);
            const auto g0 = intersect_detail::Marcher<T>::rates(m.local(start.q));
            auto passes = [&](T a) { return m.passes_start(start, tangent(a0 + d), g0, at(a), tangent(a), at(a + th), th); };
            CHECK((passes(a0 - 2 * th) || passes(a0 - th) || passes(a0)));
            CHECK_FALSE(passes(a0 + th));
            CHECK_FALSE(passes(a0 - 3 * th));
        }
}
