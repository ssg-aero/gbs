#include <doctest_gtest.hpp>
#include <gbs/intersect/surface_surface.h>
#include <gbs/primitives.h>

#include <cmath>
#include <sstream>
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

    // tilted frame, reference direction not orthogonal to the axis
    const ax2<T, 3> tilted{{{1., -2., 0.5}, {1., 1., 1.}, {1., 0., 0.}}};
    const Frame3<T> F{tilted};

    // a frame with the same axis as F, origin moved along it by h and turned by a around it
    Frame3<T> along(T h, T a = 0.)
    {
        Frame3<T> f = F;
        f.o = F.o + h * F.z;
        f.x = std::cos(a) * F.x + std::sin(a) * F.y;
        f.y = cross(f.z, f.x);
        return f;
    }

    // every curve lies on both surfaces through its pcurves, inside the boxes
    void expect_valid(const Surface<T, 3> &s1, const std::array<T, 4> &b1, const Surface<T, 3> &s2, const std::array<T, 4> &b2,
                      const SurfaceSurfaceResult<T> &r)
    {
        REQUIRE(r.complete);
        for (const auto &c : r.curves)
        {
            INFO("curve on [" << c.a << ", " << c.b << "], exact " << c.exact << ", box1 " << b1[0] << " " << b1[1] << " " << b1[2] << " " << b1[3]);
            REQUIRE(c.b > c.a);
            REQUIRE(c.deviation <= tol);
            for (int k = 0; k <= 20; ++k)
            {
                const T t = c.a + (c.b - c.a) * k / 20.;
                const auto p = c.curve->value(t);
                const auto q1 = c.pc1->value(t), q2 = c.pc2->value(t);
                REQUIRE(norm(s1.value(q1[0], q1[1]) - p) <= c.deviation + 1e-9);
                REQUIRE(norm(s2.value(q2[0], q2[1]) - p) <= c.deviation + 1e-9);
                REQUIRE(q1[0] >= b1[0] - 1e-7);
                REQUIRE(q1[0] <= b1[1] + 1e-7);
                REQUIRE(q1[1] >= b1[2] - 1e-7);
                REQUIRE(q1[1] <= b1[3] + 1e-7);
                REQUIRE(q2[0] >= b2[0] - 1e-7);
                REQUIRE(q2[0] <= b2[1] + 1e-7);
                REQUIRE(q2[1] >= b2[2] - 1e-7);
                REQUIRE(q2[1] <= b2[3] + 1e-7);
            }
        }
        for (const auto &h : r.points)
        {
            REQUIRE(norm(s1.value(h.u1, h.v1) - h.p) <= 1e-9);
            REQUIRE(norm(s2.value(h.u2, h.v2) - h.p) <= 1e-9);
        }
    }

    SurfaceSurfaceResult<T> check(const Surface<T, 3> &s1, const Surface<T, 3> &s2)
    {
        auto r = intersect(s1, s2);
        expect_valid(s1, s1.bounds(), s2, s2.bounds(), r);
        // symmetric
        auto q = intersect(s2, s1);
        expect_valid(s2, s2.bounds(), s1, s1.bounds(), q);
        REQUIRE(q.curves.size() == r.curves.size());
        REQUIRE(q.points.size() == r.points.size());
        REQUIRE(q.coincident == r.coincident);
        return r;
    }

    T length(const IntersectionCurve<T> &c)
    {
        T l{};
        constexpr int n = 2000;
        for (int k = 0; k < n; ++k)
            l += norm(c.curve->value(c.a + (c.b - c.a) * (k + 1.) / n) - c.curve->value(c.a + (c.b - c.a) * T(k) / n));
        return l;
    }
}

TEST(tests_intersect_surfaces, planes)
{
    const Plane<T> p1{F, {-5., 5., -5., 5.}};
    Frame3<T> g = F;
    g.z = F.x;
    g.x = F.y;
    g.y = cross(g.z, g.x);
    const Plane<T> p2{g, {-5., 5., -5., 5.}};
    auto r = check(p1, p2);
    ASSERT_EQ(r.curves.size(), 1u);
    const auto &c = r.curves[0];
    ASSERT_TRUE(c.exact);
    ASSERT_FALSE(c.closed);
    ASSERT_NEAR(length(c), 10., 1e-6); // the line y of both planes, clipped to [-5, 5]
    // parallel planes, coincident planes
    Frame3<T> up = F;
    up.o = F.o + 2. * F.z;
    ASSERT_TRUE(intersect(p1, Plane<T>{up, {-5., 5., -5., 5.}}).curves.empty());
    auto co = intersect(p1, Plane<T>{F, {0., 1., 0., 1.}});
    ASSERT_TRUE(co.coincident);
    ASSERT_TRUE(co.curves.empty());
}

TEST(tests_intersect_surfaces, plane_and_cylinder)
{
    const CylindricalSurface<T> cyl{2., F, {0., 2 * pi, -5., 5.}};
    // perpendicular: a closed circle, exact on the cylinder
    const Plane<T> across{along(1.), {-5., 5., -5., 5.}};
    auto r = check(cyl, across);
    ASSERT_EQ(r.curves.size(), 1u);
    ASSERT_TRUE(r.curves[0].closed);
    ASSERT_NEAR(length(r.curves[0]), 4 * pi, 1e-5);
    const auto q = r.curves[0].pc1->value(1.);
    ASSERT_NEAR(q[0], 1., 1e-12); // u = t
    ASSERT_NEAR(q[1], 1., 1e-12);
    // a box missing a quarter around the origin of the circle parameter: one arc across it
    const CylindricalSurface<T> part{2., F, {pi / 2, 2 * pi, -5., 5.}};
    const Plane<T> rotated{along(1., -pi / 4), {-5., 5., -5., 5.}};
    auto pr = check(part, rotated);
    ASSERT_EQ(pr.curves.size(), 1u);
    ASSERT_NEAR(length(pr.curves[0]), 3 * pi, 1e-5);
    // half a cylinder: half a circle
    const CylindricalSurface<T> half{2., F, {0., pi, -5., 5.}};
    auto h = check(half, across);
    ASSERT_EQ(h.curves.size(), 1u);
    ASSERT_FALSE(h.curves[0].closed);
    ASSERT_NEAR(length(h.curves[0]), 2 * pi, 1e-5);

    // parallel to the axis: two generatrices clipped to v in [-5, 5]; tangent: one
    Frame3<T> side = F;
    side.z = F.x;
    side.x = F.z;
    side.y = cross(side.z, side.x);
    side.o = F.o + 1. * F.x;
    auto l = check(cyl, Plane<T>{side, {-10., 10., -10., 10.}});
    ASSERT_EQ(l.curves.size(), 2u);
    for (const auto &c : l.curves)
    {
        ASSERT_TRUE(c.exact);
        ASSERT_NEAR(length(c), 10., 1e-6);
    }
    side.o = F.o + 2. * F.x;
    auto t = check(cyl, Plane<T>{side, {-10., 10., -10., 10.}});
    ASSERT_EQ(t.curves.size(), 1u);
    ASSERT_TRUE(t.curves[0].tangent);

    // oblique: an ellipse of axes 2 / cos(a) and 2, cut at the seam of the cylinder
    Frame3<T> ob = F;
    const T a = pi / 5;
    ob.z = std::cos(a) * F.z + std::sin(a) * F.y;
    ob.x = F.x;
    ob.y = cross(ob.z, ob.x);
    auto e = check(cyl, Plane<T>{ob, {-10., 10., -10., 10.}});
    ASSERT_FALSE(e.curves.empty());
    T len{};
    for (const auto &c : e.curves)
    {
        ASSERT_FALSE(c.exact); // the pcurve on the plane is interpolated
        len += length(c);
    }
    // perimeter of the ellipse (Ramanujan)
    const T A = 2. / std::cos(a), B = 2., hh = (A - B) * (A - B) / ((A + B) * (A + B));
    ASSERT_NEAR(len, pi * (A + B) * (1 + 3 * hh / (10 + std::sqrt(4 - 3 * hh))), 1e-5);
}

TEST(tests_intersect_surfaces, plane_sphere_cone)
{
    const SphericalSurface<T> sphere{2., F, {0., 2 * pi, -pi / 2, pi / 2}};
    Frame3<T> ob = F;
    ob.z = intersect_detail::unit<T>(F.z + F.x);
    ob.x = intersect_detail::unit<T>(F.x - F.z);
    ob.y = cross(ob.z, ob.x);
    ob.o = F.o + 1. * ob.z;
    auto r = check(sphere, Plane<T>{ob, {-5., 5., -5., 5.}});
    T len{};
    for (const auto &c : r.curves)
        len += length(c);
    ASSERT_NEAR(len, 2 * pi * std::sqrt(3.), 1e-5);
    // tangent plane: a contact point; far plane: nothing
    ob.o = F.o + 2. * ob.z;
    auto p = check(sphere, Plane<T>{ob, {-5., 5., -5., 5.}});
    ASSERT_TRUE(p.curves.empty());
    ASSERT_EQ(p.points.size(), 1u);
    ob.o = F.o + 3. * ob.z;
    ASSERT_TRUE(intersect(sphere, Plane<T>{ob, {-5., 5., -5., 5.}}).points.empty());

    // cone: circle (perpendicular), ellipse (steep), generatrices (through the apex), parabola (not yet)
    const T k = std::tan(pi / 6);
    const ConicalSurface<T> cone{1., pi / 6, F, {0., 2 * pi, -1., 3.}};
    auto c = check(cone, Plane<T>{along(2.), {-9., 9., -9., 9.}});
    ASSERT_EQ(c.curves.size(), 1u);
    ASSERT_EQ((dynamic_cast<const BSCurve<T, 2> &>(*c.curves[0].pc1).degree()), 1u); // exact on the cone, interpolated on the plane
    ASSERT_NEAR(length(c.curves[0]), 2 * pi * (1. + 2. * k), 1e-5);
    Frame3<T> st = F;
    st.z = intersect_detail::unit<T>(F.z + 0.3 * F.x);
    st.x = intersect_detail::unit<T>(F.x - 0.3 * F.z);
    st.y = cross(st.z, st.x);
    st.o = F.o + 1. * F.z;
    auto el = check(cone, Plane<T>{st, {-9., 9., -9., 9.}});
    ASSERT_FALSE(el.curves.empty());
    for (const auto &x : el.curves)
        ASSERT_FALSE(x.tangent);
    // a plane containing the axis: through the apex, two generatrices
    Frame3<T> ax = F;
    ax.z = F.y;
    ax.x = F.z;
    ax.y = cross(ax.z, ax.x);
    auto g = check(cone, Plane<T>{ax, {-9., 9., -9., 9.}});
    ASSERT_EQ(g.curves.size(), 2u);
    for (const auto &x : g.curves)
    {
        ASSERT_TRUE(x.exact);
        ASSERT_NEAR(length(x), 4. / std::cos(pi / 6), 1e-6); // v in [-1, 3]
    }
    // parallel to a generatrix: a parabola, by marching
    Frame3<T> pa = F;
    pa.z = intersect_detail::unit<T>(std::cos(pi / 6) * F.x - std::sin(pi / 6) * F.z);
    pa.x = F.y;
    pa.y = cross(pa.z, pa.x);
    pa.o = F.o + 0.5 * F.x;
    auto para = check(cone, Plane<T>{pa, {-9., 9., -9., 9.}});
    ASSERT_FALSE(para.curves.empty());
    for (const auto &x : para.curves)
        ASSERT_FALSE(x.exact);
}

TEST(tests_intersect_surfaces, coaxial)
{
    const CylindricalSurface<T> cyl{3., F, {0., 2 * pi, -5., 5.}};
    const SphericalSurface<T> sphere{5., along(0.5, 0.7), {0., 2 * pi, -pi / 2, pi / 2}};
    const ToroidalSurface<T> torus{3., 1., along(-2., 0.3), {0., 2 * pi, 0., 2 * pi}};
    const ConicalSurface<T> cone{1., pi / 4, along(0.), {0., 2 * pi, -0.5, 4.}};

    // cylinder and sphere: two circles at h = 0.5 +- 4, cut at the seam of the sphere (its frame is
    // turned by 0.7): four arcs
    auto cs = check(cyl, sphere);
    ASSERT_EQ(cs.curves.size(), 4u);
    T len{};
    for (const auto &c : cs.curves)
    {
        ASSERT_FALSE(c.closed);
        len += length(c);
    }
    ASSERT_NEAR(len, 12 * pi, 1e-4); // lengths measured by chords
    // cylinder and torus: R = 3 meets the tube at h = -2 +- 1; cut at the seam of the torus,
    // whose frame is turned by 0.3, and exact on both (u2 = u1 - 0.3)
    auto ct = check(cyl, torus);
    ASSERT_EQ(ct.curves.size(), 4u);
    for (const auto &c : ct.curves)
        ASSERT_TRUE(c.exact);
    // cylinder and cone: rho = 1 + h = 3 at h = 2
    auto cc = check(cyl, cone);
    ASSERT_EQ(cc.curves.size(), 1u);
    ASSERT_TRUE(cc.curves[0].exact);
    // plane normal to the axis and torus: two circles, then one tangent circle on top
    auto pt = check(torus, Plane<T>{along(-2.), {-9., 9., -9., 9.}});
    ASSERT_EQ(pt.curves.size(), 2u);
    auto top = check(torus, Plane<T>{along(-1.), {-9., 9., -9., 9.}});
    ASSERT_EQ(top.curves.size(), 1u);
    ASSERT_TRUE(top.curves[0].tangent);
    // sphere and cone: one circle at h = 3.2 (the other is out of the box of the cone), cut at the
    // seam of the sphere
    auto sc = check(sphere, cone);
    ASSERT_EQ(sc.curves.size(), 2u);
    // the same surface twice
    ASSERT_TRUE(intersect(cyl, CylindricalSurface<T>{3., along(4., 1.), {0., pi, 0., 1.}}).coincident);
}

TEST(tests_intersect_surfaces, spheres_tori_cylinders)
{
    // two spheres: a circle; tangent: a point; concentric and equal: coincident
    const SphericalSurface<T> s1{2., F, {0., 2 * pi, -pi / 2, pi / 2}};
    Frame3<T> f2 = F;
    f2.o = F.o + 2. * F.x;
    auto r = check(s1, SphericalSurface<T>{2., f2, {0., 2 * pi, -pi / 2, pi / 2}});
    T len{};
    for (const auto &c : r.curves)
        len += length(c);
    ASSERT_NEAR(len, 2 * pi * std::sqrt(3.), 1e-5);
    f2.o = F.o + 4. * F.x;
    auto t = check(s1, SphericalSurface<T>{2., f2, {0., 2 * pi, -pi / 2, pi / 2}});
    ASSERT_TRUE(t.curves.empty());
    ASSERT_EQ(t.points.size(), 1u);
    ASSERT_TRUE(intersect(s1, SphericalSurface<T>{2., along(0., 1.), {0., pi, 0., 1.}}).coincident);

    // a plane through the axis of a torus: two meridian circles, exact on the torus
    const ToroidalSurface<T> torus{3., 1., F, {0., 2 * pi, 0., 2 * pi}};
    Frame3<T> ax = F;
    ax.z = F.y;
    ax.x = F.z;
    ax.y = cross(ax.z, ax.x);
    auto m = check(torus, Plane<T>{ax, {-9., 9., -9., 9.}});
    ASSERT_EQ(m.curves.size(), 2u);
    for (const auto &c : m.curves)
    {
        ASSERT_TRUE(c.closed);
        ASSERT_NEAR(length(c), 2 * pi, 1e-5);
    }
    // a plane cutting the torus obliquely: by marching
    Frame3<T> ob = F;
    ob.z = intersect_detail::unit<T>(F.z + 0.2 * F.x);
    ob.x = intersect_detail::unit<T>(F.x - 0.2 * F.z);
    ob.y = cross(ob.z, ob.x);
    ob.o = F.o + 0.3 * F.z;
    ASSERT_FALSE(check(torus, Plane<T>{ob, {-9., 9., -9., 9.}}).curves.empty());

    // cylinders of parallel axes: two lines; tangent: one
    const CylindricalSurface<T> c1{2., F, {0., 2 * pi, -5., 5.}};
    Frame3<T> g = F;
    g.o = F.o + 3. * F.x + 7. * F.z;
    auto l = check(c1, CylindricalSurface<T>{2., g, {0., 2 * pi, -10., 0.}});
    ASSERT_EQ(l.curves.size(), 2u);
    for (const auto &c : l.curves)
    {
        ASSERT_TRUE(c.exact);
        ASSERT_NEAR(length(c), 8., 1e-6); // v1 in [-3, 5]: the second cylinder covers [-3, 7]
    }
    g.o = F.o + 4. * F.x;
    auto tg = check(c1, CylindricalSurface<T>{2., g, {0., 2 * pi, -5., 5.}});
    ASSERT_EQ(tg.curves.size(), 1u);
    ASSERT_TRUE(tg.curves[0].tangent);
    // skew cylinders: by marching
    Frame3<T> sk = F;
    sk.z = F.x;
    sk.x = F.y;
    sk.y = cross(sk.z, sk.x);
    ASSERT_FALSE(check(c1, CylindricalSurface<T>{1., sk, {0., 2 * pi, -5., 5.}}).curves.empty());
}

TEST(tests_intersect_surfaces, short_parts)
{
    // parts of the closed-form curves shorter than a sampling step of the former restriction
    // (1/512 of the range): found by crossing the pcurves with the sides of the boxes
    const CylindricalSurface<T> cyl{1., F, {0., 2 * pi, -2., 2.}};
    // the circle of a plane across the axis, on a cylinder face 0.005 rad wide
    const Plane<T> across{along(0.3), {-5., 5., -5., 5.}};
    const CylindricalSurface<T> sliver{1., F, {1., 1.005, -2., 2.}};
    auto r = intersect(sliver, sliver.bounds(), across, across.bounds());
    expect_valid(sliver, sliver.bounds(), across, across.bounds(), r);
    REQUIRE(r.curves.size() == 1u);
    ASSERT_NEAR(length(r.curves[0]), 0.005, 1e-9);
    // the same circle on plane faces whose strip y >= 1 - 1e-8 it only grazes (an arc of 2.8e-4
    // rad, a fortieth of a former sampling step), the strip turned around the axis; on the seam of
    // the cylinder the arc is cut in two
    const T y0 = 1. - 1e-8;
    for (int k = 0; k < 16; ++k)
    {
        INFO("strip turned by " << k * pi / 8.);
        const Plane<T> strip{along(0.3, k * pi / 8.), {-1., 1., y0, 1.1}};
        auto s = intersect(cyl, cyl.bounds(), strip, strip.bounds());
        expect_valid(cyl, cyl.bounds(), strip, strip.bounds(), s);
        REQUIRE(s.curves.size() == (k == 12 ? 2u : 1u));
        T l{};
        std::ostringstream pieces;
        for (const auto &c : s.curves)
        {
            pieces << " [" << c.a << ", " << c.b << "] " << length(c);
            l += length(c);
        }
        INFO("pieces" << pieces.str());
        // (on a crossing at 1.4e-4 rad the ends slide along the curve by the tolerance of the
        // pcurves over that angle)
        ASSERT_NEAR(l, 2. * std::acos(y0), 1e-5);
    }
    // a seam inside the range of the carriers: a sphere whose seam is turned to the angle 1 and a
    // coaxial band of cylinder around that angle (1 -+ 1e-3); each circle is cut in two parts of
    // 1e-3 at the seam of the sphere, which stay apart with their pcurves in the boxes
    {
        const SphericalSurface<T> turned{std::sqrt(2.), along(0., 1.), {0., 2 * pi, -pi / 2, pi / 2}};
        const CylindricalSurface<T> band{1., F, {1. - 1e-3, 1. + 1e-3, -2., 2.}};
        for (int order = 0; order < 2; ++order)
        {
            INFO("order " << order);
            auto s = order == 0 ? intersect(turned, turned.bounds(), band, band.bounds()) : intersect(band, band.bounds(), turned, turned.bounds());
            if (order == 0)
                expect_valid(turned, turned.bounds(), band, band.bounds(), s);
            else
                expect_valid(band, band.bounds(), turned, turned.bounds(), s);
            REQUIRE(s.curves.size() == 4u);
            for (const auto &c : s.curves)
                ASSERT_NEAR(length(c), 1e-3, 1e-9);
        }
    }
    // a full turn on a full cylinder: one closed curve, cut at the seam only
    auto f = intersect(cyl, cyl.bounds(), across, across.bounds());
    expect_valid(cyl, cyl.bounds(), across, across.bounds(), f);
    REQUIRE(f.curves.size() == 1u);
    ASSERT_TRUE(f.curves[0].closed);
    ASSERT_NEAR(length(f.curves[0]), 2. * pi, 1e-5); // (2000 chords)
}

TEST(tests_intersect_surfaces, restriction_torus_cone)
{
    // a plane through the axis of a torus (R 3, r 1) restricted to half of its tube (v in 0..pi):
    // two half meridians of length pi
    const ToroidalSurface<T> torus{3., 1., F, {0., 2 * pi, 0., pi}};
    const Plane<T> meridian{Frame3<T>{ax2<T, 3>{{F.o, F.x, F.z}}}, {-9., 9., -9., 9.}};
    auto r = intersect(torus, torus.bounds(), meridian, meridian.bounds());
    expect_valid(torus, torus.bounds(), meridian, meridian.bounds(), r);
    T l{};
    for (const auto &c : r.curves)
        l += length(c);
    ASSERT_NEAR(l, 2. * pi, 1e-5);
    // an ellipse of a plane on a cone near its apex (v = -1.732), on a part of the cone
    // (u in 0.5..4, v in -1.7..1): the parts stay in the box although their pcurves bulge
    const ConicalSurface<T> cone{1., pi / 6, F, {0.5, 4., -1.7, 1.}};
    const auto n = intersect_detail::unit<T>(F.z + 0.4 * F.x);
    const Plane<T> cut{Frame3<T>{ax2<T, 3>{{F.o - 1.5 * F.z, n, F.y}}}, {-9., 9., -9., 9.}};
    auto e = intersect(cone, cone.bounds(), cut, cut.bounds());
    expect_valid(cone, cone.bounds(), cut, cut.bounds(), e);
    REQUIRE_FALSE(e.curves.empty());
}

// ---- review of the closed forms (lot 2 of the review of the BREP kernel, brep2_pr03_analytic.md § 7)

namespace
{
    Frame3<T> frame(const point<T, 3> &o, const point<T, 3> &z, const point<T, 3> &x) { return Frame3<T>{ax2<T, 3>{{o, z, x}}}; }

    T total_length(const SurfaceSurfaceResult<T> &r)
    {
        T l{};
        for (const auto &c : r.curves)
            l += length(c);
        return l;
    }
}

TEST(tests_intersect_surfaces, planes_at_a_small_angle)
{
    // planes through P at a small angle: the line was computed with 1 - (n1 . n2)², which rounds to
    // 0 below 1e-8 rad (an exception on a line at infinity), and from absolute offsets (far from
    // the origin, a line off the planes or none at all); within tol across the boxes, the same plane
    for (T far : {0., 1e3, 1e5})
        for (T th : {1e-2, 1e-4, 1e-6, 1e-8, 1e-9})
        {
            INFO("far " << far << ", angle " << th);
            const point<T, 3> P{far, 2 * far, 0.5 * far + 0.3};
            const Plane<T> p1{frame(P, {0., 0., 1.}, {1., 0., 0.}), {-5., 5., -5., 5.}};
            const Plane<T> p2{frame(P, {0., std::sin(th), std::cos(th)}, {1., 0., 0.}), {-5., 5., -5., 5.}};
            const auto r = check(p1, p2);
            if (th >= 1e-6)
            {
                REQUIRE(r.curves.size() == 1u);
                ASSERT_NEAR(length(r.curves[0]), 10., 1e-6);
            }
            else
                REQUIRE(r.coincident);
        }
}

TEST(tests_intersect_surfaces, spindle_torus_and_double_cone)
{
    // a spindle torus (R 1 < r 2, a DEGENERATE_TOROIDAL_SURFACE of STEP) also passes across its
    // axis: its meridian is the tube circle and its mirror at rho = -R
    const Frame3<T> f = frame({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.});
    const ToroidalSurface<T> spindle{1., 2., f, {0., 2 * pi, -pi, pi}};
    Frame3<T> up = f;
    up.o = {0., 0., 1.};
    // the plane at h = 1: circles of radii 1 + sqrt(3) and sqrt(3) - 1
    const auto pl = check(spindle, Plane<T>{up, {-9., 9., -9., 9.}});
    ASSERT_EQ(pl.curves.size(), 2u);
    ASSERT_NEAR(total_length(pl), 2 * pi * 2 * std::sqrt(3.), 1e-5);
    // a coaxial cylinder of radius 0.5 meets both circles twice: four circles, the inner two cut
    // at the seam of the torus (at u + pi), which the seam of the cylinder keeps apart
    const auto cy = check(spindle, CylindricalSurface<T>{0.5, f, {0., 2 * pi, -5., 5.}});
    ASSERT_EQ(cy.curves.size(), 6u);
    ASSERT_NEAR(total_length(cy), 4 * pi, 1e-5);
    // the pcurves on the torus of the inner part (R + r cos v < 0) are at u + pi: exact
    for (const auto &c : cy.curves)
        ASSERT_TRUE(c.exact);

    // a cone on both nappes (v in -3..3, apex at v = 0) and planes through its axis: beyond the apex
    // the pcurves were off by pi (deviation 6)
    for (T semi : {pi / 4, pi / 3})
        for (T rot : {0., 1e-6, 0.3})
        {
            INFO("semi angle " << semi << ", turned by " << rot);
            const ConicalSurface<T> cone{0., semi, f, {0., 2 * pi, -3., 3.}};
            const Plane<T> ax{frame({0., 0., 0.}, {0., 1., 0.}, {std::sin(semi + rot), 0., std::cos(semi + rot)}), {-9., 9., -9., 9.}};
            const auto r = check(cone, ax);
            ASSERT_EQ(r.curves.size(), 4u); // two generatrices cut at the apex
            ASSERT_NEAR(total_length(r), 4 * 3 / std::cos(semi), 1e-9);
        }
}

TEST(tests_intersect_surfaces, planes_through_cone_apex)
{
    // a cone face starting at its apex (v in 0..10) and planes through the apex: two generatrices,
    // complete; the span beyond the apex was counted as holding a part it touched at 1e-32
    const Frame3<T> f = frame({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.});
    for (T semi : {pi / 6, pi / 4, pi / 3})
    {
        const ConicalSurface<T> cone{0., semi, f, {0., 2 * pi, 0., 10.}};
        for (int k = 0; k < 40; ++k)
        {
            INFO("semi angle " << semi << ", plane " << k);
            const T a = 2.39996 * k; // golden angle
            const point<T, 3> n{std::cos(a), std::sin(a), 0.3 * std::sin(3. * a + 1.)};
            const point<T, 3> x{std::sin(1.7 * a + 0.5), std::cos(0.9 * a), std::sin(2.3 * a)};
            const auto r = check(cone, Plane<T>{frame({0., 0., 0.}, n, x), {-30., 30., -30., 30.}});
            ASSERT_EQ(r.curves.size(), 2u);
            for (const auto &c : r.curves)
                ASSERT_NEAR(length(c), 10. / std::cos(semi), 1e-6);
        }
    }
}

TEST(tests_intersect_surfaces, tangent_within_tolerance)
{
    // tangencies and an apex accepted within tol by the closed forms: the carrier is on one surface,
    // up to tol off the other; the restriction required tol / 10 of every pcurve
    const Frame3<T> f = frame({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.});
    const CylindricalSurface<T> cyl{10., f, {0., 2 * pi, 0., 10.}};
    const ConicalSurface<T> cone{0., pi / 6, f, {0., 2 * pi, 0., 10.}};
    const ToroidalSurface<T> torus{5., 1., f, {0., 2 * pi, -pi, pi}};
    for (T g : {0.5, 0.9})
    {
        const T e = g * tol;
        INFO("gap " << e);
        const auto a = check(cyl, Plane<T>{frame({10. + e, 0., 0.}, {1., 0., 0.}, {0., 1., 0.}), {-20., 20., -20., 20.}});
        ASSERT_EQ(a.curves.size(), 1u);
        ASSERT_TRUE(a.curves[0].tangent);
        const auto b = check(cone, Plane<T>{frame({0., e, 0.}, {0., 1., 0.}, {1., 0., 0.}), {-20., 20., -20., 20.}});
        ASSERT_EQ(b.curves.size(), 2u);
        const auto c = check(torus, Plane<T>{frame({0., 0., 1. + e}, {0., 0., 1.}, {1., 0., 0.}), {-20., 20., -20., 20.}});
        ASSERT_EQ(c.curves.size(), 1u);
        ASSERT_TRUE(c.curves[0].tangent);
        const auto d = check(cyl, CylindricalSurface<T>{3., frame({13. + e, 0., 0.}, {0., 0., 1.}, {1., 0., 0.}), {0., 2 * pi, 0., 10.}});
        ASSERT_EQ(d.curves.size(), 1u);
        ASSERT_TRUE(d.curves[0].tangent);
        const auto i = check(torus, CylindricalSurface<T>{4. - e, frame({0., 0., -5.}, {0., 0., 1.}, {1., 0., 0.}), {0., 2 * pi, 0., 10.}});
        ASSERT_EQ(i.curves.size(), 1u);
        ASSERT_TRUE(i.curves[0].tangent);
    }
}

TEST(tests_intersect_surfaces, circle_near_pole)
{
    // a circle passing at d from the poles of a sphere: its pcurve turns by pi over a length of the
    // order of d, which uniform samples (2048 at most) could not follow
    const Frame3<T> f = frame({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.});
    const SphericalSurface<T> sphere{10., f, {0., 2 * pi, -pi / 2, pi / 2}};
    for (T d : {5e-7, 1e-5, 1e-3, 1e-1})
    {
        INFO("distance to the poles " << d);
        const auto r = check(sphere, Plane<T>{frame({d, 0., 0.}, {1., 0., 0.}, {0., 1., 0.}), {-20., 20., -20., 20.}});
        ASSERT_NEAR(total_length(r), 2 * pi * std::sqrt(100. - d * d), 1e-4);
        // two spheres (concentric within tol for the first)
        if (d < 1e-5)
            continue;
        const auto s = check(sphere, SphericalSurface<T>{10., frame({2 * d, 0., 0.}, {0., 0., 1.}, {1., 0., 0.}), {0., 2 * pi, -pi / 2, pi / 2}});
        ASSERT_NEAR(total_length(s), 2 * pi * std::sqrt(100. - d * d), 1e-4);
    }
}

TEST(tests_intersect_surfaces, near_parallel)
{
    // parallel or coaxial within tol across the boxes, not below a fixed 1e-10 rad: the same cylinder
    // twice, its axis written with 9 digits, gave contact points and an incomplete marching
    const Frame3<T> f = frame({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.});
    const CylindricalSurface<T> cyl{10., f, {0., 2 * pi, 0., 10.}};
    for (T e : {1e-9, 1e-8})
    {
        INFO("tilt " << e);
        const auto same = intersect(cyl, CylindricalSurface<T>{10., frame({0., 0., 0.}, {std::sin(e), 0., std::cos(e)}, {1., 0., 0.}), {0., 2 * pi, 0., 10.}});
        ASSERT_TRUE(same.coincident);
        ASSERT_TRUE(same.complete);
        // a coaxial cone: one circle (rho = 5 + h tan 30° = 10)
        const auto c = check(cyl, ConicalSurface<T>{5., pi / 6, frame({0., 0., 0.}, {std::sin(e), 0., std::cos(e)}, {1., 0., 0.}), {0., 2 * pi, 0., 10.}});
        ASSERT_EQ(c.curves.size(), 1u);
        ASSERT_TRUE(c.curves[0].closed);
        // a plane parallel to the axis: two generatrices
        const auto p = check(cyl, Plane<T>{frame({5., 0., 0.}, {std::cos(e), 0., std::sin(e)}, {0., 1., 0.}), {-20., 20., -20., 20.}});
        ASSERT_EQ(p.curves.size(), 2u);
        ASSERT_NEAR(total_length(p), 20., 1e-6);
    }
}

TEST(tests_intersect_surfaces, long_ellipse)
{
    // a plane nearly parallel to the axis of a cylinder: an ellipse of 1e6 or 1e7 of which 20 lie in
    // the boxes; its pcurves were fitted over the whole turn and missed tol / 10
    const CylindricalSurface<T> cyl{10., frame({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}), {0., 2 * pi, 0., 10.}};
    for (T e : {1e-6, 1e-5, 1e-4})
    {
        INFO("tilt " << e);
        const auto r = check(cyl, Plane<T>{frame({5., 0., 0.}, {std::cos(e), 0., std::sin(e)}, {0., 1., 0.}), {-20., 20., -20., 20.}});
        ASSERT_EQ(r.curves.size(), 2u);
        ASSERT_NEAR(total_length(r), 20., 1e-4);
    }
}

TEST(tests_intersect_surfaces, far_plane_origin)
{
    // the frame of a plane far from the part (an origin of the assembly in STEP): the constant term
    // of the conic cancelled (5.7e-6 off the cone at 1e6)
    const point<T, 3> P{0., 0., 0.};
    const auto fc = frame(P, {1., 1., 1.}, {1., 0., 0.});
    const CylindricalSurface<T> cyl{10., fc, {0., 2 * pi, -20., 20.}};
    const ConicalSurface<T> cone{5., 0.3, fc, {0., 2 * pi, -5., 20.}};
    const auto fp = frame(P + point<T, 3>{1., 2., 3.}, {1., -2., 4.}, {0., 1., 0.});
    for (T D : {1e5, 1e6})
    {
        INFO("plane origin " << D << " away");
        Frame3<T> far = fp;
        far.o = fp.o + D * fp.x;
        const Plane<T> pl{far, {-D - 40., -D + 40., -40., 40.}};
        ASSERT_NEAR(total_length(check(cyl, pl)), total_length(check(cyl, Plane<T>{fp, {-40., 40., -40., 40.}})), 1e-6);
        ASSERT_NEAR(total_length(check(cone, pl)), total_length(check(cone, Plane<T>{fp, {-40., 40., -40., 40.}})), 1e-6);
    }
}

TEST(tests_intersect_surfaces, coincident_within_a_distance)
{
    // a coincidence decided at a distance other than the tolerance of the intersection (the
    // confusion radius of the BREP operations, 10 tol)
    const std::array<T, 4> b{-5., 5., -5., 5.};
    const Plane<T> p1{frame({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}), b};
    // parallel planes 5 tol apart
    const Plane<T> p2{frame({0., 0., 5. * tol}, {0., 0., 1.}, {1., 0., 0.}), b};
    ASSERT_FALSE(intersect(p1, p2).coincident);
    ASSERT_TRUE(coincident(p1, b, p2, b, 10. * tol));
    ASSERT_FALSE(coincident(p1, b, p2, b, 4. * tol));
    // turned by a about y: the gap reaches 5 a at x = +-5, by the largest gap across the region
    for (auto [a, same] : {std::pair{1.8 * tol, true}, std::pair{2.2 * tol, false}})
    {
        INFO("a = " << a);
        const Plane<T> turned{frame({0., 0., 0.}, {-std::sin(a), 0., std::cos(a)}, {std::cos(a), 0., std::sin(a)}), b};
        ASSERT_EQ(coincident(p1, b, turned, b, 10. * tol), same);
    }
    // coaxial cylinders whose radii differ by 5 tol, cylinders whose axes are 5 tol apart
    const std::array<T, 4> bc{0., 2 * pi, -5., 5.};
    const CylindricalSurface<T> c1{3., F, bc};
    const CylindricalSurface<T> c2{3. + 5. * tol, F, bc};
    ASSERT_FALSE(intersect(c1, c2).coincident);
    ASSERT_TRUE(coincident(c1, bc, c2, bc, 10. * tol));
    ASSERT_FALSE(coincident(c1, bc, CylindricalSurface<T>{3. + 15. * tol, F, bc}, bc, 10. * tol));
    Frame3<T> off = F;
    off.o = F.o + 5. * tol * F.x;
    ASSERT_TRUE(coincident(c1, bc, CylindricalSurface<T>{3., off, bc}, bc, 10. * tol));
}

TEST(tests_intersect_surfaces, pcurve_deviation_at_the_ends)
{
    // the deviation of a pcurve covers the gap of its curve to the surface at its ends, where it
    // peaks on a slightly turned or offset surface (an edge imprinted on a face a few tolerances
    // from its own): measured between samples, it was under it, and the pcurve was off the
    // tolerance of the edge raised to it
    auto gap = [](const Curve<T, 3> &c, T t, const Surface<T, 3> &s, const PCurveFit<T> &f) {
        const auto q = f.pcurve->value(t);
        return norm(s.value(q[0], q[1]) - c.value(t));
    };
    // a segment on a plane turned by 5 tol: an exact pcurve
    const T a = 5. * tol;
    const std::array<T, 4> b{-2., 2., -2., 2.};
    const Plane<T> turned{frame({0., 0., 0.}, {-std::sin(a), 0., std::cos(a)}, {std::cos(a), 0., std::sin(a)}), b};
    const LineCurve<T> segment{{0., 0., 0.}, {1., 0., 0.}, 0., 1.};
    const auto fs = intersect_detail::primitive_pcurve(segment, T(0), T(1), turned, b, tol);
    ASSERT_TRUE(fs.exact);
    ASSERT_GE(fs.deviation, gap(segment, 1., turned, fs));
    // a circle R 5 about an axis 5 tol from that of the cylinder: an interpolated pcurve
    const std::array<T, 4> bc{0., 2 * pi, -1., 1.};
    const CylindricalSurface<T> cyl{5., frame({0., 0., 0.}, {0., 0., 1.}, {1., 0., 0.}), bc};
    const EllipseCurve<T> circle{5., 5., frame({a, 0., 0.}, {0., 0., 1.}, {1., 0., 0.}), 0., 2 * pi};
    const auto fc = intersect_detail::primitive_pcurve(circle, T(0), 2 * pi, cyl, bc, tol);
    ASSERT_FALSE(fc.exact);
    ASSERT_GE(fc.deviation, gap(circle, 0., cyl, fc));
}
