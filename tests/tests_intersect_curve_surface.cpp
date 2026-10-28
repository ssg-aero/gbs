#include <doctest_gtest.hpp>
#include <gbs/intersect/curve_surface.h>
#include <gbs/bscbuild.h>
#include <gbs/bscinterp.h>
#include <gbs/bssinterp.h>
#include <gbs/bselementary.h>
#include <gbs/primitives.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numbers>
#include <random>
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

// ---- review of the curve-surface intersection (lot 2 of the review of the BREP kernel)

TEST(tests_intersect_curve_surface, lying_curves_on_tilted_surfaces)
{
    // a generatrix of a cylinder: the subdivision culled its pairs by axis-aligned boxes only and
    // seeded a piece only if a single patch box held it: max_pairs reached, a truncated portion
    // or thousands of hits on a tilted frame
    const Frame3<T> f{tilted};
    for (T R : {1., 2.})
    {
        const CylindricalSurface<T> cylinder{R, f, {0., 2 * pi, -5., 5.}};
        const auto nurbs = cylinder.to_nurbs();
        for (T u : {0., pi / 2, 0.3, 1.})
        {
            const auto o = f.o + R * f.radial(u);
            const BSCurve<T, 3> generatrix{points_vector<T, 3>{o - T(8) * f.z, o, o + T(8) * f.z}, std::vector<T>{-8., -8., 0., 8., 8.}, 1};
            for (const Surface<T, 3> *s : std::vector<const Surface<T, 3> *>{&cylinder, nurbs.get()})
            {
                auto r = intersect(generatrix, -8., 8., *s, s->bounds());
                ASSERT_TRUE(r.complete);
                ASSERT_TRUE(r.hits.empty());
                ASSERT_EQ(r.portions.size(), 1u);
                ASSERT_NEAR(r.portions[0].ta, -5., 1e-4);
                ASSERT_NEAR(r.portions[0].tb, 5., 1e-4);
            }
        }
    }
    const CylindricalSurface<T> tall{1., f, {0., 2 * pi, -20., 20.}};
    const auto tall_nurbs = tall.to_nurbs();
    const auto o = f.o + f.radial(1.);
    const BSCurve<T, 3> long_generatrix{points_vector<T, 3>{o - T(23) * f.z, o, o + T(23) * f.z}, std::vector<T>{-23., -23., 0., 23., 23.}, 1};
    auto rl = intersect(long_generatrix, -23., 23., *tall_nurbs, tall_nurbs->bounds());
    ASSERT_TRUE(rl.complete);
    ASSERT_TRUE(rl.hits.empty());
    ASSERT_EQ(rl.portions.size(), 1u);
    ASSERT_NEAR(rl.portions[0].ta, -20., 1e-4);
    ASSERT_NEAR(rl.portions[0].tb, 20., 1e-4);
    // a cubic running a few tens of tolerances above a tilted B-spline plane: no crossing, complete
    const Plane<T> plane{f, {-10., 10., -10., 10.}};
    const auto bplane = plane.to_nurbs();
    for (T delta : {1e-4, 3e-5})
    {
        const auto d = f.x + T(0.3) * f.y;
        const auto p = f.o + delta * f.z;
        const BSCurve<T, 3> line{points_vector<T, 3>{p - T(4) * d, p - T(2) * d, p + T(2) * d, p + T(4) * d}, std::vector<T>{0, 0, 0, 0, 1, 1, 1, 1}, 3};
        auto r = intersect(line, *bplane);
        ASSERT_TRUE(r.complete);
        ASSERT_TRUE(r.hits.empty());
        ASSERT_TRUE(r.portions.empty());
    }
}

TEST(tests_intersect_curve_surface, tangent_and_grazing_contacts)
{
    // a tangency through subdivision: its fragments were portions, its near misses dozens of hits
    const Frame3<T> f{tilted};
    const Plane<T> plane{f, {-10., 10., -10., 10.}};
    const auto bplane = plane.to_nurbs();
    const CylindricalSurface<T> cylinder{1., f, {0., 2 * pi, -5., 5.}};
    const auto ncylinder = cylinder.to_nurbs();
    for (T gap : {0., 5e-7, -5e-7, 1e-8}) // touching, missing, dipping under the tolerance
    {
        INFO("gap ", gap);
        Frame3<T> h = f;
        h.y = f.z;
        h.z = cross(h.x, h.y);
        h.o = f.o + T(1. + gap) * f.z; // vertical circle of radius 1 touching z = 0 at t = 3 pi / 2
        const auto circle = EllipseCurve<T>{1., 1., h, 0., 2 * pi}.to_nurbs();
        for (const Surface<T, 3> *s : std::vector<const Surface<T, 3> *>{&plane, bplane.get()})
        {
            auto r = intersect(*circle, *s);
            ASSERT_TRUE(r.portions.empty());
            ASSERT_EQ(r.hits.size(), 1u);
            ASSERT_TRUE(r.hits[0].tangent);
            ASSERT_LT(norm(r.hits[0].p - f.o), 2e-3);
        }
        const auto p = f.o + T(1. + gap) * f.y;
        const BSCurve<T, 3> line{points_vector<T, 3>{p - T(3) * f.x, p - f.x, p + f.x, p + T(3) * f.x}, std::vector<T>{-3, -3, -3, -3, 3, 3, 3, 3}, 3};
        for (const Surface<T, 3> *s : std::vector<const Surface<T, 3> *>{&cylinder, ncylinder.get()})
        {
            auto r = intersect(line, *s);
            ASSERT_TRUE(r.portions.empty());
            ASSERT_EQ(r.hits.size(), 1u);
            ASSERT_TRUE(r.hits[0].tangent);
        }
    }
    // a curve within 0.4 tol of the plane over its length, touching it at x = k pi: a portion
    points_vector<T, 3> q;
    for (int i = 0; i <= 200; ++i)
    {
        const T x = 8. * i / 200.;
        q.push_back(f.o + x * f.x + T(0.4e-6) * std::pow(std::sin(x), 2) * f.z);
    }
    const auto wavy = interpolate<T, 3>(q, 3, KnotsCalcMode::CHORD_LENGTH);
    for (const Surface<T, 3> *s : std::vector<const Surface<T, 3> *>{&plane, bplane.get()})
    {
        auto r = intersect(wavy, *s);
        ASSERT_TRUE(r.hits.empty());
        ASSERT_EQ(r.portions.size(), 1u);
        ASSERT_NEAR(r.portions[0].ta, wavy.bounds()[0], 1e-9);
        ASSERT_NEAR(r.portions[0].tb, wavy.bounds()[1], 1e-9);
    }
    // a segment ending on a sphere, crossing it at a small angle (sine 0.01): one crossing at its
    // end, not a portion
    const SphericalSurface<T> sphere{1.7, f, {0., 2 * pi, -pi / 2, pi / 2}};
    const T u = 2.1, v = 0.4;
    const auto p = sphere.value(u, v), n = cross(sphere.value(u, v, 1, 0), sphere.value(u, v, 0, 1));
    const auto w = sphere.value(u, v, 1, 0);
    const auto dir = T(std::sqrt(1. - 1e-4)) * (w / norm(w)) - T(0.01) * (n / norm(n));
    const BSCurve<T, 3> end{points_vector<T, 3>{p - T(3) * dir, p - T(1.5) * dir, p}, std::vector<T>{0., 0., 0.5, 1., 1.}, 1};
    auto re = intersect(end, 0., 1., sphere, sphere.bounds());
    ASSERT_TRUE(re.portions.empty());
    ASSERT_EQ(re.hits.size(), 1u);
    ASSERT_LT(norm(re.hits[0].p - p), tol);
}

TEST(tests_intersect_curve_surface, lines_far_from_the_origin)
{
    // the closed forms depended on the parameter origin of the line: coefficients of a ray starting
    // far away cancelled, the leading ones were dropped against the others
    const Frame3<T> f{tilted};
    const ToroidalSurface<T> torus{3., 1., f, {0., 2 * pi, 0., 2 * pi}};
    const T w = std::sqrt(0.75);
    for (T D : {0., 1e3, 1e4, 3e4})
    {
        for (const auto &dir : {f.x, (f.x + T(0.2) * f.y) / norm(f.x + T(0.2) * f.y)})
        {
            INFO("D ", D);
            const auto c = f.o + T(0.5) * f.z; // crosses z = 0.5 of the tube four times
            const LineCurve<T> ray{c - (D + 10.) * dir, dir, 0., D + 20.};
            auto r = intersect(ray, 0., D + 20., torus, torus.bounds());
            ASSERT_EQ(r.hits.size(), 4u);
            const std::array<T, 4> x{-3. - w, -3. + w, 3. - w, 3. + w};
            if (dir == f.x)
                for (std::size_t k = 0; k < 4; ++k)
                    ASSERT_NEAR(r.hits[k].t - D - 10., x[k], 1e-6);
        }
    }
    // 1e-15 t² - 1 on [0, 1e8]: the leading coefficient is small, not its term (root sqrt(1e15))
    const auto pr = polynomial_roots<T>({-1., 0., 1e-15}, 0., 1e8);
    ASSERT_EQ(pr.roots.size(), 1u);
    ASSERT_NEAR(pr.roots[0], std::sqrt(1e15), 1e-3);
    // a degree-1 B-spline on knots far from 0 across a unit cylinder
    const CylindricalSurface<T> cylinder{1., Frame3<T>{}, {0., 2 * pi, -5., 5.}};
    for (T k0 : {0., 1e7, 3e7})
    {
        const BSCurve<T, 3> line{points_vector<T, 3>{{-1.5, 0.2, 0.}, {1.5, 0.2, 0.}}, std::vector<T>{k0, k0, k0 + 2, k0 + 2}, 1};
        auto r = intersect(line, k0, k0 + 2, cylinder, cylinder.bounds());
        ASSERT_EQ(r.hits.size(), 2u);
        ASSERT_NEAR(r.hits[0].p[0], -std::sqrt(0.96), 1e-8); // (t is known to 2e-9 at 1e7)
        ASSERT_NEAR(r.hits[1].p[0], std::sqrt(0.96), 1e-8);
    }
}

TEST(tests_intersect_curve_surface, closed_forms_within_tolerance)
{
    // the closed forms ignored the tolerance: a curve on the surface up to the rounding (or within
    // the tolerance) gave nothing instead of a portion, a near miss nothing instead of a contact
    const Frame3<T> f{tilted};
    const Plane<T> plane{f, {-10., 10., -10., 10.}};
    for (T delta : {0., 1e-9, 1e-7})
    {
        const LineCurve<T> line{f.o + delta * f.z, f.x, -4., 4.};
        auto r = intersect(line, -4., 4., plane, plane.bounds());
        ASSERT_TRUE(r.hits.empty());
        ASSERT_EQ(r.portions.size(), 1u);
        ASSERT_NEAR(r.portions[0].ta, -4., 1e-9);
        ASSERT_NEAR(r.portions[0].tb, 4., 1e-9);
    }
    for (T R : {1., 50.})
    {
        const CylindricalSurface<T> cylinder{R, f, {0., 2 * pi, -5., 5.}};
        const LineCurve<T> line{f.o + (R + 1e-9) * f.radial(0.3), f.z, -8., 8.};
        auto r = intersect(line, -8., 8., cylinder, cylinder.bounds());
        ASSERT_EQ(r.portions.size(), 1u);
        ASSERT_NEAR(r.portions[0].ta, -5., 1e-4);
        ASSERT_NEAR(r.portions[0].tb, 5., 1e-4);
    }
    // the iso u lines of a cylinder of radius 100, built by the surface itself
    std::mt19937 rng(5);
    std::uniform_real_distribution<T> U(-1., 1.);
    for (int k = 0; k < 20; ++k)
    {
        const point<T, 3> z{U(rng), U(rng), U(rng)}, x{U(rng), U(rng), U(rng)}, o{100 * U(rng), 100 * U(rng), 100 * U(rng)};
        const CylindricalSurface<T> cylinder{100., Frame3<T>{ax2<T, 3>{o, z, x}}, {0., 2 * pi, -100., 100.}};
        const auto iso = cylinder.iso_u(pi * (U(rng) + 1), -100., 100.);
        auto r = intersect(*iso, -100., 100., cylinder, cylinder.bounds());
        ASSERT_TRUE(r.hits.empty());
        ASSERT_EQ(r.portions.size(), 1u);
    }
    // a circle parallel to the plane at 1e-8: a portion; a circle missing it by 1e-8: a contact
    Frame3<T> g = f;
    g.o = f.o + T(1e-8) * f.z;
    auto rf = intersect(EllipseCurve<T>{1., 1., g, 0., 2 * pi}, 0., 2 * pi, plane, plane.bounds());
    ASSERT_TRUE(rf.hits.empty());
    ASSERT_EQ(rf.portions.size(), 1u);
    Frame3<T> h = f;
    h.y = f.z;
    h.z = cross(h.x, h.y);
    h.o = f.o + T(1. + 1e-8) * f.z;
    auto rn = intersect(EllipseCurve<T>{1., 1., h, 0., 2 * pi}, 0., 2 * pi, plane, plane.bounds());
    ASSERT_TRUE(rn.portions.empty());
    ASSERT_EQ(rn.hits.size(), 1u);
    ASSERT_TRUE(rn.hits[0].tangent);
    ASSERT_NEAR(rn.hits[0].t, 1.5 * pi, 1e-6);
}

TEST(tests_intersect_curve_surface, crossings_at_range_ends)
{
    // segments ending on a surface: the closed forms dropped the crossing at t = b about 40 % of the
    // time (strict bounds, rounding of the polynomial)
    std::mt19937 rng(42);
    std::uniform_real_distribution<T> U(-1., 1.);
    for (int k = 0; k < 100; ++k)
    {
        const point<T, 3> z{U(rng), U(rng), U(rng)}, x{U(rng), U(rng), U(rng)}, o{10 * U(rng), 10 * U(rng), 10 * U(rng)};
        const Frame3<T> f{ax2<T, 3>{o, z, x}};
        const Plane<T> plane{f, {-10., 10., -10., 10.}};
        const CylindricalSurface<T> cylinder{1.7, f, {0., 2 * pi, -10., 10.}};
        const SphericalSurface<T> sphere{1.7, f, {0., 2 * pi, -pi / 2, pi / 2}};
        const T u = 0.3 + 5 * (U(rng) + 1) / 2, v = 0.9 * U(rng);
        const point<T, 3> dir{U(rng), U(rng), U(rng)};
        for (const Surface<T, 3> *s : std::vector<const Surface<T, 3> *>{&plane, &cylinder, &sphere})
        {
            const auto p = s->value(u, v);
            const LineCurve<T> line{p - T(3) * dir, T(3) * dir, 0., 1.};
            const BSCurve<T, 3> spline{points_vector<T, 3>{p - T(3) * dir, p - T(1.5) * dir, p}, std::vector<T>{0., 0., 0.5, 1., 1.}, 1};
            for (const Curve<T, 3> *c : std::vector<const Curve<T, 3> *>{&line, &spline})
            {
                auto r = intersect(*c, 0., 1., *s, s->bounds());
                ASSERT_TRUE(std::ranges::any_of(r.hits, [&](const auto &h) { return norm(h.p - p) < tol; }));
            }
        }
    }
}

TEST(tests_intersect_curve_surface, crossings_near_poles)
{
    // cubics crossing a NURBS sphere a few 1e-6 from its pole: Gauss-Newton stopped on the pole
    // (Su = 0), where its steps cannot turn u, and the crossing was lost
    const SphericalSurface<T> sphere{1., Frame3<T>{}, {0., 2 * pi, -pi / 2, pi / 2}};
    const auto nurbs = sphere.to_nurbs();
    struct Case
    {
        T eps, a;
        point<T, 3> dir, bend;
    };
    for (const auto &[eps, a, dir, bend] :
         {Case{2.9765690130303959e-06, 1.7312369483311913, {-0.045226251629314373, -0.10251708638720436, 0.11652534965088179},
               {-0.045882782934088415, 0.024367746539819671, 0.017007774713931955}},
          Case{1.1885551968140967e-06, 1.9619117849795269, {0.27770906543218787, -0.91645902554302072, 0.49052144164841982},
               {-0.0077686373408950796, -0.019811681486300133, -0.049999348071474423}},
          Case{1.4514602903034925e-06, 1.7394167046019093, {0.63849059465224811, -0.34856970646551411, 0.066378724351603635},
               {0.046742392117526524, 0.023581561516546959, -0.0081560655001304671}},
          Case{2.8839046729659392e-08, 2.1467368986311386, {-0.62871644103164237, 0.74712505899454196, -0.061048698205485974},
               {0.046821861397898693, 0.009807563618170867, -0.049991634720572262}}})
    {
        INFO("eps ", eps);
        const auto target = sphere.value(a, pi / 2 - eps);
        const BSCurve<T, 3> c{points_vector<T, 3>{target - T(0.6) * dir, target - T(0.2) * dir + bend, target + T(0.2) * dir - bend, target + T(0.6) * dir},
                              std::vector<T>{0, 0, 0, 0, 1, 1, 1, 1}, 3};
        for (const Surface<T, 3> *s : std::vector<const Surface<T, 3> *>{&sphere, nurbs.get()})
        {
            auto r = intersect(c, *s);
            expect_on_both(c, *s, s->bounds(), r);
            ASSERT_TRUE(std::ranges::any_of(r.hits, [&](const auto &h) { return norm(h.p - target) < 1e-5; }));
        }
    }
}

TEST(tests_intersect_curve_surface, closed_curve_crossing_at_start)
{
    // a closed curve crossing the surface at its start point was reported twice unless both
    // refinements landed exactly on t = a and t = b
    std::mt19937 rng(3);
    std::uniform_real_distribution<T> U(-1., 1.);
    for (int k = 0; k < 50; ++k)
    {
        const point<T, 3> z{U(rng), U(rng), U(rng)}, x{U(rng), U(rng), U(rng)}, o{5 * U(rng), 5 * U(rng), 5 * U(rng)};
        const Frame3<T> g{ax2<T, 3>{o, z, x}};
        const T r = 1. + U(rng) * 0.5;
        const EllipseCurve<T> circle{r, r, g, 0., 2 * pi};
        const auto p0 = circle.value(0.);
        const point<T, 3> n = g.z + T(0.7) * U(rng) * g.y + T(0.5) * g.x;
        const point<T, 3> ref{U(rng), U(rng), U(rng)};
        const Plane<T> plane{Frame3<T>{ax2<T, 3>{p0, n, ref}}, {-5., 5., -5., 5.}};
        const auto pn = plane.to_nurbs();
        const auto cn = circle.to_nurbs();
        auto at_start = [&](const CurveSurfaceResult<T> &res) { return std::ranges::count_if(res.hits, [&](const auto &h) { return norm(h.p - p0) < tol; }); };
        ASSERT_EQ(at_start(intersect(circle, 0., 2 * pi, plane, plane.bounds())), 1);
        ASSERT_EQ(at_start(intersect(*cn, plane)), 1);
        ASSERT_EQ(at_start(intersect(*cn, *pn)), 1);
    }
}

TEST(tests_intersect_curve_surface, poles_apexes_and_nappes)
{
    // the closed forms filtered the parameters of the candidates by the box before refining: at a
    // pole u is arbitrary, on the other nappe of a cone it is off by pi
    const Frame3<T> w{};
    const SphericalSurface<T> lune{1., w, {1., 2., -pi / 2, pi / 2}};
    const LineCurve<T> axis{{0., 0., -2.}, {0., 0., 1.}, 0., 4.};
    auto ra = intersect(axis, 0., 4., lune, lune.bounds());
    expect_on_both(axis, lune, lune.bounds(), ra);
    ASSERT_EQ(ra.hits.size(), 2u);
    ASSERT_NEAR(ra.hits[0].t, 1., 1e-9);
    ASSERT_NEAR(ra.hits[1].t, 3., 1e-9);
    // a line across the nappe v < 0 of a cone of apex O (the radial direction reversed): crossings
    // at t = 3.05942 and 3.12241
    const ConicalSurface<T> cone{0., pi / 6, w, {0., 2 * pi, -2., 2.}};
    const LineCurve<T> across{{-3., 0.3, 1.}, {1., 0., -0.5}, 0., 6.};
    auto rc = intersect(across, 0., 6., cone, cone.bounds());
    expect_on_both(across, cone, cone.bounds(), rc);
    ASSERT_EQ(rc.hits.size(), 2u);
    ASSERT_NEAR(rc.hits[0].t, 3.05942, 1e-5);
    ASSERT_NEAR(rc.hits[1].t, 3.12241, 1e-5);
}
