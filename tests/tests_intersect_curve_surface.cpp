#include <doctest_gtest.hpp>
#include <gbs/intersect/curve_surface.h>
#include <gbs/bscbuild.h>
#include <gbs/bscinterp.h>
#include <gbs/bssinterp.h>
#include <gbs/bselementary.h>
#include <gbs/primitives.h>

#include <cmath>
#include <memory>
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

    void expect_on_both(const Curve<T, 3> &c, const Surface<T, 3> &s, const std::array<T, 4> &box, const CurveSurfaceResult<T> &r)
    {
        REQUIRE(r.complete);
        for (const auto &h : r.hits)
        {
            REQUIRE(norm(c.value(h.t) - h.p) < tol);
            REQUIRE(norm(s.value(h.u, h.v) - h.p) < tol);
            REQUIRE(h.u >= box[0] - 1e-9);
            REQUIRE(h.u <= box[1] + 1e-9);
            REQUIRE(h.v >= box[2] - 1e-9);
            REQUIRE(h.v <= box[3] + 1e-9);
        }
    }

    // the line through p with direction d, on [-L, L]: as a primitive and as a B-spline with an
    // inner knot, which goes through subdivision instead of the closed form
    struct Lines
    {
        LineCurve<T> analytic;
        BSCurve<T, 3> spline;
        Lines(const point<T, 3> &p, const point<T, 3> &d, T L)
            : analytic{p, d, -L, L}, spline{points_vector<T, 3>{p - L * d, p, p + L * d}, std::vector<T>{-L, -L, 0., L, L}, 1} {}
    };

    std::vector<point<T, 3>> sorted_points(const CurveSurfaceResult<T> &r)
    {
        std::vector<point<T, 3>> p;
        for (const auto &h : r.hits)
            p.push_back(h.p);
        return p;
    }
}

TEST(tests_intersect_curve_surface, polynomial_roots)
{
    // (t - 1)(t - 2)(t - 3)
    auto r = polynomial_roots<T>({-6., 11., -6., 1.}, 0., 4.);
    ASSERT_EQ(r.roots.size(), 3u);
    for (int k = 0; k < 3; ++k)
        ASSERT_NEAR(r.roots[k], k + 1., 1e-12);
    ASSERT_EQ(polynomial_roots<T>({-6., 11., -6., 1.}, 1.5, 2.5).roots.size(), 1u);
    // near double root: an extreme, not a sign change; an exact double root is a root
    auto d = polynomial_roots<T>({1. + 1e-9, -2., 1.}, 0., 3.);
    ASSERT_TRUE(d.roots.empty());
    ASSERT_EQ(d.extremes.size(), 1u);
    ASSERT_NEAR(d.extremes[0], 1., 1e-12);
    ASSERT_EQ(polynomial_roots<T>({1., -2., 1.}, 0., 3.).roots.size(), 1u);
    // quartic (t² - 1)(t² - 4), leading zeros dropped
    auto q = polynomial_roots<T>({4., 0., -5., 0., 1., 0.}, -3., 3.);
    ASSERT_EQ(q.roots.size(), 4u);
    ASSERT_NEAR(q.roots[0], -2., 1e-12);
    ASSERT_NEAR(q.roots[3], 2., 1e-12);
    ASSERT_NEAR(polynomial_value<T>({1., 2., 3.}, 2.), 17., 1e-15);
}

TEST(tests_intersect_curve_surface, segment_and_box)
{
    Box<T, 3> b;
    b.add({0., 0., 0.});
    b.add({1., 1., 1.});
    ASSERT_TRUE((segment_overlaps<T, 3>({-1., 0.5, 0.5}, {2., 0.5, 0.5}, b)));
    ASSERT_TRUE((segment_overlaps<T, 3>({0.5, 0.5, 0.5}, {0.6, 0.6, 0.6}, b))); // inside
    ASSERT_FALSE((segment_overlaps<T, 3>({-1., 2., 0.5}, {2., 2., 0.5}, b)));
    ASSERT_TRUE((segment_overlaps<T, 3>({-1., 2., 0.5}, {2., 2., 0.5}, b, 1.)));
    // the diagonal of a big box misses a corner box
    ASSERT_FALSE((segment_overlaps<T, 3>({-1., 3., -1.}, {3., -1., -1.}, b)));
    ASSERT_FALSE((segment_overlaps<T, 3>({0., 0., 0.}, {1., 1., 1.}, Box<T, 3>{})));
}

TEST(tests_intersect_curve_surface, lines_and_primitives)
{
    const Frame3<T> f{tilted};
    const CylindricalSurface<T> cylinder{2., f, {0., 2 * pi, -5., 5.}};
    const ConicalSurface<T> cone{1., pi / 6, f, {0., 2 * pi, -1., 4.}};
    const SphericalSurface<T> sphere{2., f, {0., 2 * pi, -pi / 2, pi / 2}};
    const ToroidalSurface<T> torus{3., 1., f, {0., 2 * pi, 0., 2 * pi}};
    const Plane<T> plane{f, {-10., 10., -10., 10.}};

    // a line across the axis: through the frame origin along x
    const Lines across{f.o + T(0.5) * f.z, f.x, 8.};
    struct Case
    {
        const ElementarySurface<T> *s;
        std::size_t n;
    };
    for (const auto &[s, n] : {Case{&cylinder, 2}, Case{&cone, 2}, Case{&sphere, 2}, Case{&torus, 4}, Case{&plane, 0}})
    {
        const auto box = s->bounds();
        auto ra = intersect(across.analytic, -8., 8., *s, box);
        auto rs = intersect(across.spline, -8., 8., *s, box);
        expect_on_both(across.analytic, *s, box, ra);
        expect_on_both(across.spline, *s, box, rs);
        if (s == &plane) // the line lies in the plane z = 0.5? no: it is parallel to it, above
        {
            ASSERT_TRUE(ra.hits.empty());
            ASSERT_TRUE(ra.portions.empty());
            continue;
        }
        ASSERT_EQ(ra.hits.size(), n);
        ASSERT_EQ(rs.hits.size(), n);
        for (std::size_t k = 0; k < n; ++k)
        {
            ASSERT_LT(norm(ra.hits[k].p - rs.hits[k].p), 1e-8);
            ASSERT_FALSE(ra.hits[k].tangent);
        }
    }
    // known abscissas: cylinder +-2, sphere +-sqrt(4 - 0.25), torus +-(3 +- sqrt(1 - 0.25)), cone radius 1 + 0.5 tan 30°
    auto xs = [&](const CurveSurfaceResult<T> &r) {
        std::vector<T> x;
        for (const auto &h : r.hits)
            x.push_back(h.t);
        return x;
    };
    auto rc = xs(intersect(across.analytic, -8., 8., cylinder, cylinder.bounds()));
    ASSERT_NEAR(rc[0], -2., 1e-12);
    ASSERT_NEAR(rc[1], 2., 1e-12);
    auto rt = xs(intersect(across.analytic, -8., 8., torus, torus.bounds()));
    const T w = std::sqrt(0.75);
    ASSERT_NEAR(rt[0], -3. - w, 1e-10);
    ASSERT_NEAR(rt[1], -3. + w, 1e-10);
    ASSERT_NEAR(rt[2], 3. - w, 1e-10);
    ASSERT_NEAR(rt[3], 3. + w, 1e-10);
    auto rk = xs(intersect(across.analytic, -8., 8., cone, cone.bounds()));
    ASSERT_NEAR(rk[1], 1. + 0.5 * std::tan(pi / 6), 1e-12);

    // the plane crossed obliquely; a hit outside the box is dropped
    const Lines oblique{f.o + T(3) * f.x, f.x + f.z, 8.};
    auto rp = intersect(oblique.analytic, -8., 8., plane, plane.bounds());
    ASSERT_EQ(rp.hits.size(), 1u);
    ASSERT_LT(norm(rp.hits[0].p - (f.o + T(3) * f.x)), 1e-12);
    ASSERT_NEAR(rp.hits[0].u, 3., 1e-12);
    ASSERT_TRUE(intersect(oblique.analytic, -8., 8., plane, {-2., 2., -2., 2.}).hits.empty());
    ASSERT_EQ(intersect(oblique.spline, -8., 8., plane, plane.bounds()).hits.size(), 1u);

    // tangent lines: a tangent crossing
    const Lines tangent_cyl{f.o + T(2) * f.y, f.x, 8.};
    for (const Curve<T, 3> *c : std::vector<const Curve<T, 3> *>{&tangent_cyl.analytic, &tangent_cyl.spline})
    {
        auto r = intersect(*c, -8., 8., cylinder, cylinder.bounds());
        expect_on_both(*c, cylinder, cylinder.bounds(), r);
        ASSERT_EQ(r.hits.size(), 1u);
        ASSERT_TRUE(r.hits[0].tangent);
        ASSERT_TRUE(r.portions.empty());
        ASSERT_LT(norm(r.hits[0].p - (f.o + T(2) * f.y)), 1e-5);
    }
    const Lines tangent_sphere{f.o + T(2) * f.z, f.x + f.y, 8.};
    auto ts = intersect(tangent_sphere.analytic, -8., 8., sphere, sphere.bounds());
    ASSERT_EQ(ts.hits.size(), 1u);
    ASSERT_TRUE(ts.hits[0].tangent);
}

TEST(tests_intersect_curve_surface, curves_on_surfaces)
{
    const Frame3<T> f{tilted};
    const CylindricalSurface<T> cylinder{2., f, {0., 2 * pi, -5., 5.}};
    const Plane<T> plane{f, {-1., 1., -2., 3.}};
    // a generatrix of the cylinder: the portion inside the box of the surface
    const Lines generatrix{f.o + T(2) * f.x, f.z, 8.};
    for (const Curve<T, 3> *c : std::vector<const Curve<T, 3> *>{&generatrix.analytic, &generatrix.spline})
    {
        auto r = intersect(*c, -8., 8., cylinder, cylinder.bounds());
        ASSERT_TRUE(r.hits.empty());
        ASSERT_EQ(r.portions.size(), 1u);
        ASSERT_NEAR(r.portions[0].ta, -5., 1e-4);
        ASSERT_NEAR(r.portions[0].tb, 5., 1e-4);
    }
    // a line in the plane, clipped by its box [-1, 1] along x
    const Lines in_plane{f.o, f.x, 4.};
    auto rp = intersect(in_plane.analytic, -4., 4., plane, plane.bounds());
    ASSERT_TRUE(rp.hits.empty());
    ASSERT_EQ(rp.portions.size(), 1u);
    ASSERT_NEAR(rp.portions[0].ta, -1., 1e-4);
    ASSERT_NEAR(rp.portions[0].tb, 1., 1e-4);
    ASSERT_NEAR(rp.portions[0].ua, -1., 1e-4);
    ASSERT_NEAR(rp.portions[0].ub, 1., 1e-4);
    // a parallel of the cylinder (circle) lies on it
    const auto parallel = cylinder.iso_v(1., 0., 2 * pi);
    auto rc = intersect(*parallel, 0., 2 * pi, cylinder, cylinder.bounds());
    ASSERT_TRUE(rc.hits.empty());
    ASSERT_EQ(rc.portions.size(), 1u);
    ASSERT_NEAR(rc.portions[0].ta, 0., 1e-4);
    ASSERT_NEAR(rc.portions[0].tb, 2 * pi, 1e-4);
}

TEST(tests_intersect_curve_surface, circles_and_planes)
{
    const Frame3<T> f{tilted};
    const Plane<T> plane{f, {-10., 10., -10., 10.}};
    // a circle in a plane containing z: crosses the plane z = 0 twice, at +-x
    Frame3<T> g = f;
    g.y = f.z;
    g.z = cross(g.x, g.y);
    const EllipseCurve<T> circle{2., 2., g, 0., 2 * pi};
    auto r = intersect(circle, 0., 2 * pi, plane, plane.bounds());
    expect_on_both(circle, plane, plane.bounds(), r);
    ASSERT_EQ(r.hits.size(), 2u); // t = 0 (start point, once) and t = pi
    ASSERT_NEAR(r.hits[0].t, 0., 1e-12);
    ASSERT_NEAR(r.hits[1].t, pi, 1e-12);
    // through subdivision (the exact NURBS of the circle): the same points
    const auto nurbs = circle.to_nurbs();
    auto rn = intersect(*nurbs, plane);
    expect_on_both(*nurbs, plane, plane.bounds(), rn);
    ASSERT_EQ(rn.hits.size(), 2u);
    ASSERT_LT(norm(rn.hits[0].p - r.hits[0].p), 1e-9);
    ASSERT_LT(norm(rn.hits[1].p - r.hits[1].p), 1e-9);
    // an ellipse on a range that misses the plane; a circle touching it
    const EllipseCurve<T> arc{2., 1., g, 0.2, 3.};
    ASSERT_TRUE(intersect(arc, 0.2, 3., plane, plane.bounds()).hits.empty());
    Frame3<T> h = g;
    h.o = g.o + T(2) * f.z;
    const EllipseCurve<T> touching{2., 2., h, 0., 2 * pi};
    auto rt = intersect(touching, 0., 2 * pi, plane, plane.bounds());
    ASSERT_EQ(rt.hits.size(), 1u);
    ASSERT_TRUE(rt.hits[0].tangent);
    // a circle lying in the plane: one portion
    const EllipseCurve<T> flat{1., 1., f, 0., 2 * pi};
    auto rf = intersect(flat, 0., 2 * pi, plane, plane.bounds());
    ASSERT_TRUE(rf.hits.empty());
    ASSERT_EQ(rf.portions.size(), 1u);
}

TEST(tests_intersect_curve_surface, bsplines)
{
    // a sine wave in the plane (x, z) against the NURBS of a cylinder of axis x... the curve
    // crosses the half-cylinder z = sqrt(1 - y²) where its height equals the radius
    const Frame3<T> f;
    const CylindricalSurface<T> cylinder{1., Frame3<T>{ax2<T, 3>{{{0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}}}}, {0., pi, 0., 4 * pi}};
    const auto nurbs = cylinder.to_nurbs();
    points_vector<T, 3> q;
    for (int i = 0; i <= 80; ++i)
    {
        const T x = 4 * pi * i / 80.;
        q.push_back({x, 0., 1. + 0.5 * std::sin(x)});
    }
    const auto wave = interpolate<T, 3>(q, 3, KnotsCalcMode::CHORD_LENGTH);
    for (const Surface<T, 3> *s : std::vector<const Surface<T, 3> *>{&cylinder, nurbs.get()})
    {
        auto r = intersect(wave, *s);
        expect_on_both(wave, *s, s->bounds(), r);
        ASSERT_EQ(r.hits.size(), 5u); // x = 0, pi, 2 pi, 3 pi, 4 pi
        for (std::size_t k = 0; k < 5; ++k)
            ASSERT_NEAR(r.hits[k].p[0], k * pi, 1e-3);
    }
    // a cubic through a tilted circle of radius 3 crosses the torus tube twice per turn of the
    // tube it meets: against the primitive (sampled) and its exact NURBS (Bézier patches)
    const Frame3<T> ft{tilted};
    const ToroidalSurface<T> torus{3., 1., ft, {0., 2 * pi, 0., 2 * pi}};
    const auto tnurbs = torus.to_nurbs();
    points_vector<T, 3> ring;
    for (int i = 0; i <= 40; ++i)
    {
        const T a = 2 * pi * i / 40.;
        ring.push_back(ft.o + T(3) * ft.radial(a) + T(1.5) * std::sin(a) * ft.z);
    }
    const auto loop = interpolate<T, 3>(ring, 3, KnotsCalcMode::CHORD_LENGTH);
    std::size_t n_hits{};
    for (const Surface<T, 3> *s : std::vector<const Surface<T, 3> *>{&torus, tnurbs.get()})
    {
        auto r = intersect(loop, *s);
        expect_on_both(loop, *s, s->bounds(), r);
        if (n_hits == 0)
            n_hits = r.hits.size();
        ASSERT_EQ(r.hits.size(), n_hits);
    }
    ASSERT_EQ(n_hits, 4u); // leaves the tube at a = asin(2/3) and pi - asin(2/3), comes back below

    // a B-spline surface against a line through its patches
    points_vector<T, 3> grid;
    for (int j = 0; j <= 6; ++j)
        for (int i = 0; i <= 6; ++i)
            grid.push_back({T(i), T(j), 0.3 * std::sin(T(i)) * std::cos(T(j))});
    const auto patch = interpolate<T, 3>(grid, 7, 3, 3, KnotsCalcMode::CHORD_LENGTH);
    const Lines vertical{{2.5, 3.5, -1.}, {0., 0., 1.}, 3.};
    auto rv = intersect(vertical.analytic, -3., 3., patch, patch.bounds());
    expect_on_both(vertical.analytic, patch, patch.bounds(), rv);
    ASSERT_EQ(rv.hits.size(), 1u);
    ASSERT_NEAR(rv.hits[0].p[0], 2.5, 1e-9);
    ASSERT_NEAR(rv.hits[0].p[1], 3.5, 1e-9);
}
