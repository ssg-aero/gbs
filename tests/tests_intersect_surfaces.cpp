#include <doctest_gtest.hpp>
#include <gbs/intersect/surface_surface.h>
#include <gbs/primitives.h>

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
