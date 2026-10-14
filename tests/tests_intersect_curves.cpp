#include <doctest_gtest.hpp>
#include <gbs/intersect/curves.h>
#include <gbs/bscbuild.h>
#include <gbs/bscinterp.h>
#include <gbs/bselementary.h>
#include <gbs/curves>

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

    template <std::size_t dim>
    void expect_on_both(const Curve<T, dim> &c1, const Curve<T, dim> &c2, const CurveCurveResult<T, dim> &r)
    {
        REQUIRE(r.complete);
        for (const auto &h : r.hits)
        {
            REQUIRE(norm(c1.value(h.t1) - h.p) < tol);
            REQUIRE(norm(c2.value(h.t2) - h.p) < tol);
        }
    }
}

TEST(tests_intersect_curves, boxes)
{
    Box<T, 2> b;
    ASSERT_TRUE(b.empty());
    ASSERT_EQ(b.diagonal(), 0.);
    b.add({0., 0.});
    b.add({3., 4.});
    ASSERT_NEAR(b.diagonal(), 5., 1e-15);
    Box<T, 2> c;
    c.add({3.5, 0.});
    ASSERT_FALSE(b.overlaps(c));
    ASSERT_TRUE(b.overlaps(c, 0.5));
    ASSERT_TRUE(b.contains(point<T, 2>{1., 1.}));
    c.add({1., 1.});
    ASSERT_FALSE(b.contains(c));
    ASSERT_TRUE(b.contains(c, 0.5));
    ASSERT_FALSE(b.overlaps(Box<T, 2>{}));
}

TEST(tests_intersect_curves, segments)
{
    const auto s1 = build_segment<T, 2>({0., 0.}, {2., 2.});
    const auto s2 = build_segment<T, 2>({0., 2.}, {2., 0.});
    auto r = intersect(s1, s2);
    expect_on_both(s1, s2, r);
    ASSERT_EQ(r.hits.size(), 1u);
    ASSERT_LT(norm(r.hits[0].p - point<T, 2>{1., 1.}), 1e-12);
    ASSERT_FALSE(r.hits[0].tangent);
    ASSERT_TRUE(r.overlaps.empty());

    // parallel and apart; meeting at an end (T junction)
    ASSERT_TRUE(intersect(s1, build_segment<T, 2>({0., 1.}, {2., 3.})).hits.empty());
    auto t = intersect(s1, build_segment<T, 2>({2., 2.}, {3., 0.}));
    ASSERT_EQ(t.hits.size(), 1u);
    ASSERT_NEAR(t.hits[0].t1, s1.bounds()[1], 1e-12);

    // collinear: common portion, in the same and in the opposite direction
    const auto a = build_segment<T, 2>({0., 0.}, {2., 0.}, true);
    for (bool reversed : {false, true})
    {
        const auto b = reversed ? build_segment<T, 2>({3., 0.}, {1., 0.}, true) : build_segment<T, 2>({1., 0.}, {3., 0.}, true);
        auto o = intersect(a, b);
        ASSERT_TRUE(o.hits.empty());
        ASSERT_EQ(o.overlaps.size(), 1u);
        const auto &ov = o.overlaps[0];
        ASSERT_NEAR(ov.t1a, 0.5, 1e-9);
        ASSERT_NEAR(ov.t1b, 1., 1e-9);
        ASSERT_NEAR(ov.t2a, reversed ? 1. : 0., 1e-9);
        ASSERT_NEAR(ov.t2b, 0.5, 1e-9);
    }
}

TEST(tests_intersect_curves, circles_and_lines)
{
    // rational circle (B-spline path) and a generic circle (sampled path)
    const auto bs_circle = build_circle<T, 2>(1.);
    const Circle2d<T> circle{1.};
    for (const Curve<T, 2> *c : std::vector<const Curve<T, 2> *>{&bs_circle, &circle})
    {
        const auto line = build_segment<T, 2>({-2., 0.5}, {2., 0.5});
        auto r = intersect(*c, line);
        expect_on_both(*c, line, r);
        ASSERT_EQ(r.hits.size(), 2u);
        const T x = std::sqrt(0.75);
        std::vector<T> xs{r.hits[0].p[0], r.hits[1].p[0]};
        std::ranges::sort(xs);
        ASSERT_NEAR(xs[0], -x, 1e-9);
        ASSERT_NEAR(xs[1], x, 1e-9);

        // tangent line: one tangent crossing
        const auto tangent = build_segment<T, 2>({-2., 1.}, {2., 1.});
        auto rt = intersect(*c, tangent);
        expect_on_both(*c, tangent, rt);
        ASSERT_EQ(rt.hits.size(), 1u);
        ASSERT_TRUE(rt.hits[0].tangent);
        ASSERT_LT(norm(rt.hits[0].p - point<T, 2>{0., 1.}), 1e-5);

        // a line through the start point of the closed circle: (1, 0) once
        const auto diameter = build_segment<T, 2>({-2., 0.}, {2., 0.});
        auto rd = intersect(*c, diameter);
        expect_on_both(*c, diameter, rd);
        ASSERT_EQ(rd.hits.size(), 2u);
    }

    // two circles: two crossings at x = 1/2
    const auto c2 = build_circle<T, 2>(1., {1., 0.});
    auto rc = intersect(bs_circle, c2);
    expect_on_both(bs_circle, c2, rc);
    ASSERT_EQ(rc.hits.size(), 2u);
    for (const auto &h : rc.hits)
        ASSERT_NEAR(h.p[0], 0.5, 1e-9);

    // an arc lying on the circle: a common portion
    const auto arc = build_circle_arc<T>(1., pi / 4, 3 * pi / 4);
    auto ra = intersect(bs_circle, arc);
    ASSERT_TRUE(ra.hits.empty());
    ASSERT_EQ(ra.overlaps.size(), 1u);
    const auto &ov = ra.overlaps[0];
    ASSERT_LT(norm(bs_circle.value(ov.t1a) - arc.value(ov.t2a)), 1e-6);
    ASSERT_LT(norm(bs_circle.value(ov.t1a) - point<T, 2>{std::sqrt(0.5), std::sqrt(0.5)}), 1e-5);
    ASSERT_LT(norm(bs_circle.value(ov.t1b) - point<T, 2>{-std::sqrt(0.5), std::sqrt(0.5)}), 1e-5);

    // a curve with itself: one overlap on the whole range
    auto self = intersect(c2, c2);
    ASSERT_TRUE(self.hits.empty());
    ASSERT_EQ(self.overlaps.size(), 1u);
    ASSERT_NEAR(self.overlaps[0].t1a, c2.bounds()[0], 1e-9);
    ASSERT_NEAR(self.overlaps[0].t1b, c2.bounds()[1], 1e-9);
}

TEST(tests_intersect_curves, bsplines)
{
    // a cubic through a sine wave crosses the axis at each multiple of pi
    points_vector<T, 2> q;
    for (int i = 0; i <= 60; ++i)
    {
        const T x = 4 * pi * i / 60.;
        q.push_back({x, std::sin(x)});
    }
    const auto wave = interpolate<T, 2>(q, 3, KnotsCalcMode::CHORD_LENGTH);
    const auto axis = build_segment<T, 2>({-1., 0.}, {14., 0.});
    auto r = intersect(wave, axis);
    expect_on_both(wave, axis, r);
    ASSERT_EQ(r.hits.size(), 5u); // 0, pi, 2 pi, 3 pi, 4 pi
    for (std::size_t k = 0; k < 5; ++k)
        ASSERT_NEAR(r.hits[k].p[0], k * pi, 1e-3);
    for (std::size_t k = 1; k < 5; ++k)
        ASSERT_LT(r.hits[k - 1].t1, r.hits[k].t1);

    // the same wave against a shifted wave: crossings on a sub-range only
    points_vector<T, 2> q2;
    for (const auto &p : q)
        q2.push_back({p[0], -p[1] + 0.5});
    const auto wave2 = interpolate<T, 2>(q2, 3, KnotsCalcMode::CHORD_LENGTH);
    auto r2 = intersect(wave, wave2);
    expect_on_both(wave, wave2, r2);
    ASSERT_EQ(r2.hits.size(), 4u); // sin x = 1/4, twice per period
    for (const auto &h : r2.hits)
        ASSERT_NEAR(h.p[1], 0.25, 1e-3);
    // horizontal lines across the wave: two crossings per period
    for (T y : {-0.95, -0.5, 0.1, 0.7, 0.99})
    {
        const auto l = build_segment<T, 2>({-1., y}, {14., y});
        auto ry = intersect(wave, l);
        expect_on_both(wave, l, ry);
        ASSERT_EQ(ry.hits.size(), 4u);
    }
    const auto b = wave.bounds();
    auto half = intersect(wave, b[0], (b[0] + b[1]) / 2, wave2, wave2.bounds()[0], wave2.bounds()[1]);
    ASSERT_EQ(half.hits.size(), 2u);
}

TEST(tests_intersect_curves, curves_3d)
{
    const auto l1 = build_segment<T, 3>({0., 0., 0.}, {1., 1., 1.});
    const auto l2 = build_segment<T, 3>({1., 0., 0.}, {0., 1., 1.});
    auto r = intersect(l1, l2);
    expect_on_both(l1, l2, r);
    ASSERT_EQ(r.hits.size(), 1u);
    ASSERT_LT(norm(r.hits[0].p - point<T, 3>{0.5, 0.5, 0.5}), 1e-12);
    // skew lines do not meet
    ASSERT_TRUE(intersect(l1, build_segment<T, 3>({1., 0., 0.}, {0., 1., 0.})).hits.empty());
    // circle and its diameter, in a tilted plane
    const ax2<T, 3> ax{{{1., 2., 3.}, {1., 1., 1.}, {1., -1., 0.}}};
    const auto circle = build_circle_arc<T>(2., 0., 2 * pi, ax);
    const auto d = build_segment<T, 3>(circle.value(0.3), circle.value(2.1));
    auto rc = intersect(circle, d);
    expect_on_both(circle, d, rc);
    ASSERT_EQ(rc.hits.size(), 2u);
}

TEST(tests_intersect_curves, robustness)
{
    // the generic circle (sampled path) on the NURBS circle: one overlap, whole range
    const auto bs_circle = build_circle<T, 2>(1.);
    const Circle2d<T> circle{1.};
    auto o = intersect(circle, bs_circle);
    ASSERT_TRUE(o.hits.empty());
    ASSERT_EQ(o.overlaps.size(), 1u);
    ASSERT_NEAR(o.overlaps[0].t1a, circle.bounds()[0], 1e-6);
    ASSERT_NEAR(o.overlaps[0].t1b, circle.bounds()[1], 1e-6);

    // grazing crossings: a line at a small angle through the circle
    for (T y : {0.9, 0.99, 0.999, 0.9999})
    {
        const auto line = build_segment<T, 2>({-2., y}, {2., y});
        auto r = intersect(bs_circle, line);
        expect_on_both(bs_circle, line, r);
        ASSERT_EQ(r.hits.size(), 2u);
        ASSERT_TRUE(r.overlaps.empty());
        const T x = std::sqrt(1 - y * y);
        ASSERT_NEAR(std::abs(r.hits[0].p[0]), x, 1e-9);
        ASSERT_NEAR(std::abs(r.hits[1].p[0]), x, 1e-9);
    }
    // a square crossed by a circle, as the boundary of a face in its parameter space
    const std::vector<point<T, 2>> sq{{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}};
    const auto hole = build_circle<T, 2>(0.6, {0.5, 0.5});
    std::size_t n{};
    for (std::size_t k = 0; k < 4; ++k)
    {
        const auto side = build_segment<T, 2>(sq[k], sq[(k + 1) % 4]);
        auto r = intersect(side, hole);
        expect_on_both(side, hole, r);
        n += r.hits.size();
    }
    ASSERT_EQ(n, 8u);
}

TEST(tests_intersect_curves, limits)
{
    const auto s1 = build_segment<T, 2>({0., 0.}, {2., 2.});
    const auto s2 = build_segment<T, 2>({0., 2.}, {2., 0.});
    // ranges: the crossing is outside the first half of s1
    const auto b = s1.bounds();
    ASSERT_TRUE(intersect(s1, b[0], b[0] + 0.4 * (b[1] - b[0]), s2, s2.bounds()[0], s2.bounds()[1]).hits.empty());
    // a looser tolerance finds near misses
    const auto s3 = build_segment<T, 2>({0., -1e-4}, {2., -1e-4});
    const auto s4 = build_segment<T, 2>({1., 0.}, {1., 1.});
    ASSERT_TRUE(intersect(s3, s4).hits.empty());
    ASSERT_EQ(intersect(s3, s4, IntersectOptions<T>{.tol = 1e-3}).hits.size(), 1u);
    // the cap on the pairs examined
    ASSERT_FALSE(intersect(s1, s2, IntersectOptions<T>{.max_pairs = 0}).complete);
}
