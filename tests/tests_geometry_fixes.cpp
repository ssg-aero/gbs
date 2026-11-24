// Defects of the geometry found while writing its documentation, each checked against its fix.
#include <doctest_gtest.hpp>
#include <gbs/bscurve.h>
#include <gbs/bssurf.h>
#include <gbs/bsstools.h>
#include <gbs/bscanalysis.h>
#include <gbs/bscapprox.h>
#include <gbs/bscinterp.h>
#include <gbs/bscshaping.h>
#include <gbs/bsctools.h>
#include <gbs/curveline.h>
#include <gbs/curveonsurface.h>
#include <gbs/extrema.h>
#include <gbs/loft/loftCurves.h>
#include <gbs/math.ixx>
#include <cmath>
#include <list>
#include <stdexcept>

using namespace gbs; // the operators on points

namespace
{
    // a planar cubic, not symmetric, of y max 2.1-ish inside the span
    gbs::BSCurve2d_d wavy()
    {
        return {{{0., 0.}, {1., 3.}, {2., -1.}, {3., 2.5}, {4., 0.5}}, {0., 0.5, 1.}, {4, 1, 4}, 3};
    }
} // namespace

// binomial_law computed through factorials in int overflowed beyond 12!: degree elevation to 13
// and above gave a wrong curve
TEST(tests_geometry_fixes, binomial_law_large)
{
    EXPECT_EQ((gbs::binomial_law<int, int>(4, 2)), 6);
    EXPECT_DOUBLE_EQ((gbs::binomial_law<double, int>(30, 15)), 155117520.);
    EXPECT_DOUBLE_EQ((gbs::binomial_law<double, int>(13, 0)), 1.);
    EXPECT_DOUBLE_EQ((gbs::binomial_law<double, int>(13, 13)), 1.);

    gbs::BSCurve2d_d c{{{0., 0.}, {1., 2.}, {3., -1.}, {4., 1.}}, {0., 1.}, {4, 4}, 3};
    for (size_t p : {12, 13, 14, 20})
    {
        auto e = c;
        e.increaseDegree(p - 3);
        EXPECT_EQ(e.degree(), p);
        for (double u : {0., 0.17, 0.5, 0.83, 1.})
            EXPECT_LT(gbs::distance(e(u), c(u)), 1e-12) << "degree " << p << " u " << u;
    }
}

// the second derivative of a curve on a surface used d2u where d2v was due
TEST(tests_geometry_fixes, curve_on_surface_second_derivative)
{
    gbs::BSSurface<double, 3> srf{
        {{0., 0., 0.}, {1., 0., 1.}, {2., 0., 0.}, {0., 1., 1.}, {1., 1., 3.}, {2., 1., 1.}, {0., 2., 0.}, {1., 2., 1.}, {2., 2., 0.}},
        {0., 0., 0., 1., 1., 1.},
        {0., 0., 0., 1., 1., 1.},
        2,
        2};
    gbs::BSCurve2d_d uv{{{0.1, 0.2}, {0.5, 0.9}, {0.9, 0.4}}, {0., 1.}, {3, 3}, 2};
    gbs::CurveOnSurface<double, 3> cs{uv, srf};
    const double u = 0.5, h = 1e-5;
    auto fd = (cs(u + h, 1) - cs(u - h, 1)) / (2 * h);
    EXPECT_LT(gbs::norm(cs(u, 2) - fd), 1e-6);
}

// moved_to_constraints moved the curve by the target instead of onto it
TEST(tests_geometry_fixes, moved_to_constraints_reaches_target)
{
    auto c = wavy();
    std::vector<gbs::bsc_constraint<double, 2>> cs{{0.5, {3., 2.5}, 0}};
    auto m = gbs::moved_to_constraints(c, cs);
    EXPECT_LT(gbs::distance(m(0.5), gbs::point<double, 2>{3., 2.5}), 1e-12);
    gbs::move_to_constraints(c, cs);
    EXPECT_LT(gbs::distance(c(0.5), gbs::point<double, 2>{3., 2.5}), 1e-12);
}

// pole(i, j) read the pole j + i nu: i is the rank along u, as the poles are stored
TEST(tests_geometry_fixes, surface_pole_indices)
{
    std::vector<std::array<double, 3>> poles;
    for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 3; ++i)
            poles.push_back({double(i), double(j), double(10 * j + i)});
    gbs::BSSurface<double, 3> s{poles, {0., 0., 0., 1., 1., 1.}, {0., 0., 1., 1.}, 2, 1};
    for (size_t j = 0; j < 2; ++j)
        for (size_t i = 0; i < 3; ++i)
        {
            EXPECT_EQ(s.pole(i, j), (gbs::point<double, 3>{double(i), double(j), double(10 * j + i)}));
            EXPECT_EQ(std::as_const(s).pole(i, j), s.poles()[i + 3 * j]);
        }
    EXPECT_THROW(s.pole(0, 2), std::out_of_range);
}

// BSCurve2d_f followed an orphan template line and named no type
TEST(tests_geometry_fixes, float_alias)
{
    gbs::BSCurve2d_f c{{{0.f, 0.f}, {1.f, 1.f}}, {0.f, 1.f}, {2, 2}, 1};
    EXPECT_FLOAT_EQ(c(0.5f)[0], 0.5f);
}

// max_coordinate minimised 1/x^2 and min_coordinate x^2, locally: neither the extremum nor its value
TEST(tests_geometry_fixes, extreme_coordinates)
{
    auto c = wavy();
    double y_max = -1e9, y_min = 1e9;
    for (int i = 0; i <= 100000; ++i)
    {
        auto y = c(i / 100000.)[1];
        y_max = std::max(y_max, y);
        y_min = std::min(y_min, y);
    }
    auto [u1, y1] = gbs::max_coordinate(c, 1, 0.);
    EXPECT_NEAR(y1, y_max, 1e-9);
    EXPECT_NEAR(c(u1)[1], y1, 1e-12);
    auto [u2, y2] = gbs::min_coordinate(c, 1, 0.);
    EXPECT_NEAR(y2, y_min, 1e-9);
    EXPECT_NEAR(c(u2)[1], y2, 1e-12);
    // an extremum at an end of the curve
    auto [u3, x3] = gbs::max_coordinate(c, 0, 0.);
    EXPECT_DOUBLE_EQ(u3, 1.);
    EXPECT_DOUBLE_EQ(x3, 4.);
    EXPECT_THROW(gbs::max_coordinate(c, 2, 0.), std::out_of_range);
}

// the length of an unbounded curve was a number (Gauss) or inf (adaptive)
TEST(tests_geometry_fixes, length_of_unbounded_curve)
{
    gbs::Line<double, 2> l{gbs::ax1<double, 2>{{{0., 0.}, {1., 0.}}}};
    EXPECT_THROW(gbs::length(l), std::domain_error);
    EXPECT_THROW(gbs::length(l, 0, false), std::domain_error);
    EXPECT_NEAR(gbs::length(l, 0., 2.), 2., 1e-12);
}

// the points alone interpolated with the computed tangents: the end tangents were zero
TEST(tests_geometry_fixes, interpolate_c1_end_tangents)
{
    std::vector<gbs::constrType<double, 2, 1>> Q{{{{0., 0.}}}, {{{1., 1.}}}, {{{2., 0.}}}, {{{3., 1.}}}};
    auto c = gbs::interpolate<double, 2, 1>(Q, gbs::KnotsCalcMode::CHORD_LENGTH, true);
    auto [u1, u2] = c.bounds();
    EXPECT_GT(gbs::norm(c(u1, 1)), 0.1);
    EXPECT_GT(gbs::norm(c(u2, 1)), 0.1);
    for (const auto &q : Q)
        EXPECT_LT(gbs::extrema_curve_point(c, q[0], 1e-10)[1], 1e-9);
}

// sections of large parameters: the bounds of the second changed to those of the first rounded
// its knots beyond the absolute knot_eps, inserted then a second time
TEST(tests_geometry_fixes, loft_of_large_parameters)
{
    const std::vector<double> k{0., 0.13, 0.29, 0.51, 0.77, 1.};
    const std::vector<size_t> m{4, 1, 1, 1, 1, 4};
    std::vector<std::array<double, 3>> poles;
    for (int i = 0; i < 8; ++i)
        poles.push_back({1000. * std::cos(i * 0.4), 1000. * std::sin(i * 0.4), 0.});
    gbs::BSCurve3d_d c1{poles, k, m, 3}, c2{poles, k, m, 3};
    c1.changeBounds(0., 3141.46);
    c2.changeBounds(0., 2718.28);
    gbs::translate(c2, {0., 0., 500.});
    auto s = gbs::loft(std::list<gbs::BSCurve3d_d>{c1, c2}, 1);
    EXPECT_EQ(s.nPolesU(), 8);
    for (double t : {0., 0.2, 0.5, 0.9, 1.})
        EXPECT_LT(gbs::distance(s(3141.46 * t, s.bounds()[2]), c1(3141.46 * t)), 1e-9);
}

// approx(crv, u1, u2, ...) sampled the whole curve whatever u1 and u2
TEST(tests_geometry_fixes, approx_between_parameters)
{
    auto c = wavy();
    auto a = gbs::approx(c, 0.2, 0.7, 1e-3, 1e-6, 5);
    EXPECT_NEAR(a.bounds()[0], 0.2, 1e-12);
    EXPECT_NEAR(a.bounds()[1], 0.7, 1e-12);
    for (double u : {0.2, 0.35, 0.5, 0.7})
        EXPECT_LT(gbs::distance(a(u), c(u)), 1e-4) << u;
    EXPECT_THROW(gbs::approx(c, 0.7, 0.2, 1e-3, 1e-6, 5), std::invalid_argument);
}

// join averaged the ends of curves that did not meet, silently: an optional tolerance checks it
TEST(tests_geometry_fixes, join_with_tolerance)
{
    gbs::BSCurve2d_d c1{{{0., 0.}, {1., 0.}}, {0., 1.}, {2, 2}, 1};
    gbs::BSCurve2d_d c2{{{1., 0.}, {2., 1.}}, {0., 1.}, {2, 2}, 1};
    gbs::BSCurve2d_d c3{{{1.1, 0.}, {2., 1.}}, {0., 1.}, {2, 2}, 1};
    EXPECT_LT(gbs::distance(gbs::join(c1, c2, 1e-9)(1.), gbs::point<double, 2>{1., 0.}), 1e-15);
    EXPECT_THROW(gbs::join(c1, c3, 1e-3), std::invalid_argument);
    EXPECT_NO_THROW(gbs::join(c1, c3, 0.2));
    // without a tolerance, as before: the mean of the ends
    EXPECT_LT(gbs::distance(gbs::join(c1, c3)(1.), gbs::point<double, 2>{1.05, 0.}), 1e-15);
    EXPECT_THROW(gbs::join(&c1, &c3, 1e-3), std::invalid_argument);
}

// extended at the start of u or v put the surface after the band, shifted to start at u1 (v1):
// it keeps its parameters now, the band before them, as at the end
TEST(tests_geometry_fixes, extended_at_start_keeps_parameters)
{
    gbs::BSSurface<double, 3> srf{
        {{0., 0., 0.}, {1., 0., 1.}, {2., 0., 0.}, {0., 1., 1.}, {1., 1., 3.}, {2., 1., 1.}, {0., 2., 0.}, {1., 2., 1.}, {2., 2., 0.}},
        {0., 0., 0., 1., 1., 1.},
        {0., 0., 0., 1., 1., 1.},
        2,
        2};
    for (auto b : {gbs::SurfaceBound::U_START, gbs::SurfaceBound::V_START, gbs::SurfaceBound::U_END, gbs::SurfaceBound::V_END})
    {
        auto e = gbs::extended(srf, 0.3, b, true);
        auto [e1, e2, f1, f2] = e.bounds();
        if (b == gbs::SurfaceBound::U_START)
            EXPECT_LT(e1, 0.);
        if (b == gbs::SurfaceBound::V_START)
            EXPECT_LT(f1, 0.);
        for (double u : {0., 0.3, 0.8, 1.})
            for (double v : {0., 0.4, 1.})
                EXPECT_LT(gbs::distance(e(u, v), srf(u, v)), 1e-9) << int(b) << " " << u << " " << v;
    }
}
