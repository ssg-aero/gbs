#include <doctest_gtest.hpp>
#include <fstream>
#include <numbers>
#include <gbs/bscanalysis.h>
#include <gbs/bscbuild.h>
#include <gbs/curvecircle.h>
#include <gbs-render/vtkGbsRender.h>
using gbs::operator-;

#ifdef TEST_PLOT_ON
    const bool PLOT_ON = true;
#else
    const bool PLOT_ON = false;
#endif
namespace {
    const double tol = 1e-10;
    const double PI = acos(-1.);
}

template<typename T, size_t dim>
auto print(const gbs::point<T,dim> &p)
{
    std::cout << "[ ";
    std::for_each(p.begin(),std::next(p.end(),-1),[](const auto &v){std::cout << " " << v << ",";});
    std::cout << p.back() << "]" << std::endl;
}

TEST(tests_bscanalysis, discretize_basic)
{
    auto c = gbs::build_circle<double,3>(1.,{0.,0.,0.});
    auto points = gbs::discretize(c,36);

    ASSERT_LT(gbs::norm(points.front()-points.back()),tol);
    std::for_each(
        points.begin(),
        points.end(),
        [&](const auto &pt_)
        {
            auto r = gbs::norm(pt_);
            ASSERT_NEAR(r,1.,tol);
        }
        );

}

TEST(tests_bscanalysis, abs_curv)
{
    auto c = gbs::build_circle<double,3>(1.,{0.,0.,0.});

    auto f_u = gbs::abs_curv(c);

    ASSERT_NEAR(f_u(  PI),0.5,1e-6);
    ASSERT_NEAR(f_u(2*PI),1.0,1e-6);

    std::vector<double> k1 = {0., 0., 0., 0., 1., 1., 1., 1.};
    std::vector<std::array<double,3> > poles1 =
    {
        {0.,1.,0.},
        {1.,2.,0.},
        {2.,2.,0.},
        {3.,0.,0.},
    };

    size_t p1 = 3;

    gbs::BSCurve3d_d c1(poles1,k1,p1);

    f_u = gbs::abs_curv<double,3>(c1, 0.3, 0.7);

    auto points = gbs::discretize(c1,5);

    std::for_each(points.begin(),points.end(),[](const auto &pt){print(pt);});
}

TEST(tests_bscanalysis, abs_curv_adaptive)
{
    auto c = gbs::build_circle<double,3>(1.,{0.,0.,0.});

    auto f_u = gbs::abs_curv_adaptive(c);

    ASSERT_NEAR(f_u(  PI),0.5,1e-6);
    ASSERT_NEAR(f_u(2*PI),1.0,1e-6);

    std::vector<double> k1 = {0., 0., 0., 0., 1., 1., 1., 1.};
    std::vector<std::array<double,3> > poles1 =
    {
        {0.,1.,0.},
        {1.,2.,0.},
        {2.,2.,0.},
        {3.,0.,0.},
    };

    size_t p1 = 3;

    gbs::BSCurve3d_d c1(poles1,k1,p1);

    f_u = gbs::abs_curv_adaptive<double,3>(c1, 0.3, 0.7);

    auto points = gbs::discretize(c1,5);

    std::for_each(points.begin(),points.end(),[](const auto &pt){print(pt);});
}

TEST(tests_bscanalysis, abs_curv_d)
{
    std::vector<double> k1 = {0., 0., 0., 0., 1., 1., 1., 1.};
    std::vector<std::array<double,3> > poles1 =
    {
        {0.,1.,0.},
        {1.,2.,0.},
        {2.,2.,0.},
        {3.,0.,0.},
    };

    size_t p1 = 3;

    gbs::BSCurve3d_d c1(poles1,k1,p1);

    auto f_u = gbs::abs_curv<double,3>(c1,30);

    // auto points = gbs::discretize(c1,5);
}

TEST(tests_bscanalysis, integrate_kronrod_any_width)
{
    // Boost's adaptive rule bisects an interval narrower than ~3e-8 down to its depth limit (983 025 evaluations of a
    // constant); taken on [-1, 1], one pass of 15 points whatever the width
    for (double h : {1e-12, 1e-8, 1., 1e6})
    {
        size_t n{};
        const double a = 0.5, b = a + h;
        const double r = gbs::integrate_kronrod([&](double) { ++n; return 2.; }, a, b);
        ASSERT_NEAR(r / (2. * (b - a)), 1., 1e-15);
        ASSERT_EQ(n, 15u);
    }
    // the bounds either way, the estimated error and the integral of |f|
    double error{}, l1{};
    const double s = gbs::integrate_kronrod([](double u) { return std::sin(u); }, std::numbers::pi, 0., 15, 1e-12, &error, &l1);
    ASSERT_NEAR(s, -2., 1e-14);
    ASSERT_NEAR(l1, 2., 1e-14);
    ASSERT_LT(error, 1e-10);

    // the length of a curve: few evaluations over a narrow range (it took 2 031 585 over a range of 1e-8), the same
    // length and about the same evaluations whatever its bounds (a range of 1e3 stopped at 2e-8 of it)
    struct Counted : gbs::BSCurve<double, 2>
    {
        using gbs::BSCurve<double, 2>::BSCurve;
        mutable size_t evaluations{};
        auto value(double u, size_t d = 0) const -> std::array<double, 2> override
        {
            ++evaluations;
            return gbs::BSCurve<double, 2>::value(u, d);
        }
        auto derivatives(double u, size_t n, std::array<double, 2> *CK) const -> void override
        {
            ++evaluations;
            gbs::BSCurve<double, 2>::derivatives(u, n, CK);
        }
    };
    Counted c{gbs::points_vector<double, 2>{{1., 0.}, {0.7, 0.07}, {0.3, 0.12}, {0.05, 0.06}, {0., 0.}},
              {0., 0., 0., 0., 0.5, 1., 1., 1., 1.}, 3};
    const double narrow = gbs::length(c, 0.2, 0.2 + 1e-8);
    ASSERT_LT(c.evaluations, 100u);
    ASSERT_NEAR(narrow / (gbs::norm(c.value(0.2, 1)) * 1e-8), 1., 1e-7);
    c.evaluations = 0;
    const double whole = gbs::length(c);
    const size_t whole_evaluations = c.evaluations;
    for (double w : {1e-9, 1e3, 1e6})
    {
        Counted scaled = c;
        scaled.changeBounds(0., w);
        scaled.evaluations = 0;
        ASSERT_NEAR(gbs::length(scaled) / whole, 1., 1e-14);
        ASSERT_LE(scaled.evaluations, 2 * whole_evaluations);
    }
}

TEST(tests_bscanalysis, discretize)
{
    auto c = gbs::build_circle<double,3>(1.,{0.,0.,0.});
    size_t n =10;
    auto u = gbs::uniform_distrib_params(c,0.,0.5,n,100);
    for(size_t i {}; i < n ; i++)
    {
        auto s = i / (n-1.) * std::numbers::pi;
        gbs::point<double,3> pt {std::cos(s), std::sin(s),0.};
        auto u_ = *std::next(u.begin(),i);
        ASSERT_LT(gbs::norm(c(u_) - pt ), 1e-6);
    }
}

TEST(tests_bscanalysis, arc_length_distrib_params)
{
    // Half unit circle: arc length = angle, so the point at normalized arc
    // length s is at angle pi * s.
    auto c = gbs::build_circle<double,3>(1.,{0.,0.,0.});
    auto check = [&](const std::list<double> &u, const std::vector<double> &s)
    {
        ASSERT_EQ(u.size(), s.size());
        auto it = u.begin();
        for(size_t i {}; i < s.size(); i++, ++it)
        {
            gbs::point<double,3> pt {std::cos(PI * s[i]), std::sin(PI * s[i]),0.};
            ASSERT_LT(gbs::norm(c(*it) - pt ), 1e-6);
        }
    };

    // Prescribed normalized arc lengths, ends exact
    std::vector<double> s {0., 0.05, 0.2, 0.5, 0.9, 1.};
    auto u = gbs::arc_length_distrib_params(c, 0., 0.5, s, 200);
    check(u, s);
    ASSERT_EQ(u.front(), 0.);
    ASSERT_EQ(u.back(), 0.5);

    // Law s(xi) = xi^2: nodes clustered at the start
    size_t n = 11;
    auto u_law = gbs::arc_length_distrib_params(c, 0., 0.5, n, [](double xi){ return xi * xi; }, 200);
    std::vector<double> s_law(n);
    for(size_t i {}; i < n; i++) s_law[i] = std::pow(i / (n - 1.), 2);
    check(u_law, s_law);

    // Identity law: same as the uniform distribution
    auto u_id = gbs::arc_length_distrib_params(c, 0., 0.5, n, [](double xi){ return xi; }, 200);
    auto u_unif = gbs::uniform_distrib_params(c, 0., 0.5, n, 200);
    for(auto it_id = u_id.begin(), it_unif = u_unif.begin(); it_id != u_id.end(); ++it_id, ++it_unif)
        ASSERT_NEAR(*it_id, *it_unif, 1e-12);

    // Whole curve (bounds overloads); over a full turn the arc length law
    // (abs_curv) converges more slowly with n_law
    auto u_full = gbs::arc_length_distrib_params(c, std::vector<double>{0., 0.25, 1.}, 200);
    auto [u1, u2] = c.bounds();
    ASSERT_EQ(u_full.front(), u1);
    ASSERT_EQ(u_full.back(), u2);
    ASSERT_LT(gbs::norm(c(*std::next(u_full.begin())) - gbs::point<double,3>{0., 1., 0.}), 1e-4);
    ASSERT_EQ(gbs::arc_length_distrib_params(c, 5, [](double xi){ return xi; }, 200).size(), 5);

    // Invalid inputs
    ASSERT_THROW(gbs::arc_length_distrib_params(c, 0., 0.5, std::vector<double>{0., 0.6, 0.4, 1.}), std::invalid_argument);
    ASSERT_THROW(gbs::arc_length_distrib_params(c, 0., 0.5, std::vector<double>{0., 1.2}), std::invalid_argument);
    ASSERT_THROW(gbs::arc_length_distrib_params(c, 0., 0.5, 1, [](double xi){ return xi; }), std::invalid_argument);
}

TEST(tests_bscanalysis, discretize_refined)
{
    auto c = gbs::build_circle<double,3>(1.,{0.,0.,0.});
    auto points = gbs::discretize(c,5,0.01);

    ASSERT_LT(gbs::norm(points.front()-points.back()),tol);
    std::for_each(
        points.begin(),
        points.end(),
        [&](const auto &pt_)
        {
            auto r = gbs::norm(pt_);
            ASSERT_NEAR(r,1.,tol);
        }
        );

}

TEST(tests_bscanalysis, max_curvature_pos)
{
    // ellipse of radii 1 and 0.2: curvature a / b^2 = 25 at the vertices (+-1, 0), u = 0, 0.5, 1
    auto e = gbs::build_ellipse<double,2>(1.,0.2);
    auto [u, cur] = gbs::max_curvature_pos(e, 0.1, 0.9, 1e-63);
    ASSERT_NEAR(u, 0.5, 1e-7);
    ASSERT_NEAR(cur, 25., 1e-9); // the curvature, not |C''|^2
    {
        const auto [u_, cur_] = gbs::max_curvature_pos(e, 1e-6);
        ASSERT_NEAR(cur_, 25., 1e-9);
        ASSERT_NEAR(std::abs(e(u_)[0]), 1., 1e-12);
    }
    // a circle has the curvature 1 / R everywhere
    {
        const auto c = gbs::build_circle<double, 3>(2.);
        ASSERT_NEAR(gbs::max_curvature_pos(c, 1e-6)[1], 0.5, 1e-12);
    }
    // parabola y = x^2 on [-1, 2] (exact quadratic B-spline): curvature 2 at its vertex x = 0
    {
        const gbs::BSCurve<double, 2> parabola{{{-1., 1.}, {0.5, -2.}, {2., 4.}}, {-1., -1., -1., 2., 2., 2.}, 2};
        const auto [up, kp] = gbs::max_curvature_pos(parabola, 1e-8);
        ASSERT_NEAR(parabola(up)[0], 0., 1e-7);
        ASSERT_NEAR(kp, 2., 1e-10);
    }
    if(PLOT_ON)
        gbs::plot(
            gbs::crv_dsp<double,2,true>{
                .c =&e,
                .col_crv = {1.,0.,0.},
                .poles_on = true,
                .col_poles = {0.,1.,0.},
                .col_ctrl = {0.,0.,0.},
                .show_curvature=true,
                },
            gbs::points_vector<double,2>{e(u)});
}
namespace
{
    // Degree-p open curve with interior knots, one of them of multiplicity 2, and
    // a non-planar control net.
    gbs::BSCurve<double, 3> make_hodograph_test_curve(size_t p)
    {
        std::vector<double> k(p + 1, 0.);
        for (double u : {0.2, 0.45, 0.45, 0.7, 0.9})
            k.push_back(u);
        k.insert(k.end(), p + 1, 1.);
        const size_t n = k.size() - p - 1;
        gbs::points_vector<double, 3> poles(n);
        for (size_t i = 0; i < n; i++)
            poles[i] = {double(i), std::sin(1.3 * i), std::cos(0.7 * i) + 0.1 * i * i};
        return gbs::BSCurve<double, 3>(poles, k, p);
    }
}

// #83: derivative_curve(crv, k) is the hodograph: evaluating it reproduces the
// point derivative crv.value(u, k), for every order up to the degree.
TEST(tests_bscanalysis, derivative_curve_matches_point_derivatives)
{
    for (size_t p : {1, 2, 3, 5})
    {
        auto crv = make_hodograph_test_curve(p);
        const auto &U = crv.knotsFlats();
        std::vector<double> params = gbs::make_range<double>(0., 1., 101);
        params.insert(params.end(), U.begin(), U.end()); // knots too (one-sided)
        for (size_t k = 1; k <= p; k++)
        {
            auto dcrv = gbs::derivative_curve(crv, k);
            CAPTURE(p);
            CAPTURE(k);
            ASSERT_EQ(dcrv.degree(), p - k);
            ASSERT_EQ(dcrv.poles().size(), crv.poles().size() - k);
            ASSERT_TRUE(dcrv.knotsFlats() == std::vector<double>(U.begin() + k, U.end() - k));
            for (double u : params)
            {
                auto d_ref = crv.value(u, k);
                auto d = dcrv.value(u);
                ASSERT_LT(gbs::norm(d - d_ref), 1e-9 * (1. + gbs::norm(d_ref)));
            }
        }
    }
}

// #83: derivatives compose (d/du of the hodograph is the next hodograph), and
// orders beyond the degree give the zero curve.
TEST(tests_bscanalysis, derivative_curve_composition_and_high_order)
{
    auto crv = make_hodograph_test_curve(3);
    auto d2 = gbs::derivative_curve(crv, 2);
    auto d1d1 = gbs::derivative_curve(gbs::derivative_curve(crv), 1);
    ASSERT_TRUE(d2.knotsFlats() == d1d1.knotsFlats());
    for (size_t i = 0; i < d2.poles().size(); i++)
        ASSERT_LT(gbs::norm(d2.poles()[i] - d1d1.poles()[i]), 1e-12);

    auto d0 = gbs::derivative_curve(crv, 0);
    ASSERT_TRUE(d0.poles() == crv.poles());
    ASSERT_TRUE(d0.knotsFlats() == crv.knotsFlats());

    auto d4 = gbs::derivative_curve(crv, 4);
    ASSERT_EQ(d4.degree(), 0u);
    for (double u : gbs::make_range<double>(0., 1., 11))
        ASSERT_EQ(gbs::norm(d4.value(u)), 0.);
}

// #83: the hodograph of a straight segment is its constant velocity.
TEST(tests_bscanalysis, derivative_curve_line)
{
    gbs::points_vector<double, 2> poles{{0., 0.}, {1., 2.}, {2., 4.}, {3., 6.}};
    gbs::BSCurve<double, 2> line(poles, {0., 0., 0., 0., 2., 2., 2., 2.}, 3);
    auto v = gbs::derivative_curve(line);
    for (const auto &q : v.poles())
    {
        ASSERT_NEAR(q[0], 1.5, 1e-14);
        ASSERT_NEAR(q[1], 3., 1e-14);
    }
}

TEST(tests_bscanalysis, compute_rolling_ball)
{
    using P = gbs::point<double, 2>;
    // lines y = 0 and y = x: the ball tangent to y = 0 at (a, 0) has its centre at (a, r),
    // r = a / (1 + sqrt(2)), and touches y = x at ((a + r) / 2, (a + r) / 2)
    const auto c1 = gbs::build_segment<double, 2>({0., 0.}, {4., 0.});
    const auto c2 = gbs::build_segment<double, 2>({0., 0.}, {4., 4.});
    for (const double a : {0.5, 1., 2.5})
    {
        const auto u1 = a; // parameter = abscissa on c1
        const auto [u2, centre, parallel] = gbs::compute_rolling_ball<double>(c1, c2, u1);
        const double r = a / (1. + std::sqrt(2.));
        ASSERT_FALSE(parallel);
        ASSERT_LT(gbs::norm(centre - P{a, r}), 1e-9);
        ASSERT_LT(gbs::norm(c2(u2) - P{(a + r) / 2., (a + r) / 2.}), 1e-9);
    }
    // concentric circles: the normals are parallel, the centre is in the middle
    {
        const auto inner = gbs::build_circle<double, 2>(1.), outer = gbs::build_circle<double, 2>(3.);
        const auto [u2, centre, parallel] = gbs::compute_rolling_ball<double>(inner, outer, 0.3);
        ASSERT_TRUE(parallel);
        ASSERT_NEAR(gbs::norm(centre), 2., 1e-9);
    }
    // a ball between a circle of radius 1 and the line y = 3 above it: radius 1, centre at
    // distance 2 from the centre of the circle and 1 from the line, so at height 2
    {
        const auto circle = gbs::build_circle<double, 2>(1.);
        const auto line = gbs::build_segment<double, 2>({-5., 3.}, {5., 3.});
        const double u1 = 0.2; // a point of the upper half of the circle
        const auto p = circle(u1);
        const auto [u2, centre, parallel] = gbs::compute_rolling_ball<double>(circle, line, u1);
        ASSERT_FALSE(parallel);
        const double rho = gbs::norm(centre), r = rho - 1.;
        ASSERT_NEAR(3. - centre[1], r, 1e-9);                     // tangent to the line
        ASSERT_NEAR(gbs::norm(centre - p), r, 1e-9);              // tangent to the circle at p
        ASSERT_NEAR(centre[0] * p[1] - centre[1] * p[0], 0., 1e-9); // on the normal of the circle at p
    }
}

TEST(tests_bscanalysis, compute_rolling_ball_review)
{
    using P = gbs::point<double, 2>;
    // two circles facing each other: scanning crv2 from the projection, the normal of crv2 turns
    // parallel to the normal of crv1 (a pole of the radii, where their difference jumps) before
    // the ball is reached; that jump was taken for a root (centre at about 1e16). On the back of
    // c1 (x < 0) no ball touches c2 without holding points of c1: it throws
    {
        const auto c1 = gbs::build_circle<double, 2>(1.), c2 = gbs::build_circle<double, 2>(0.5, P{4., 0.});
        for (int i = 1; i < 40; ++i)
        {
            const double u1 = i / 40.;
            const auto p1 = c1(u1);
            if (p1[0] < 1e-12)
            {
                ASSERT_THROW(gbs::compute_rolling_ball<double>(c1, c2, u1), std::runtime_error);
                continue;
            }
            const auto [u2, centre, parallel] = gbs::compute_rolling_ball<double>(c1, c2, u1);
            const double r1 = gbs::norm(centre - p1), r2 = gbs::norm(centre - c2(u2));
            ASSERT_LT(r1, 1e3);
            ASSERT_NEAR(r1, r2, 1e-9 * (1. + r1));                            // tangent to both
            ASSERT_NEAR(centre[0] * p1[1] - centre[1] * p1[0], 0., 1e-9 * (1. + r1)); // on the normal of c1
            ASSERT_GE(gbs::norm(centre), 1. + r1 - 1e-9 * (1. + r1));          // outside c1
            ASSERT_GE(gbs::norm(centre - P{4., 0.}), 0.5 + r1 - 1e-9 * (1. + r1)); // outside c2
        }
    }
    // an unbounded crv2 cannot be searched: it throws (it used to hang)
    {
        const auto c1 = gbs::build_circle<double, 2>(1.);
        const gbs::Line<double, 2> line{P{-5., 3.}, P{5., 3.}};
        ASSERT_THROW(gbs::compute_rolling_ball<double>(c1, line, 0.2), std::invalid_argument);
    }
    // no circle tangent to both: it throws instead of returning the closest approach
    {
        const auto c1 = gbs::build_circle<double, 2>(1.);
        const auto c2 = gbs::build_segment<double, 2>({3., -0.1}, {3., 0.1});
        ASSERT_THROW(gbs::compute_rolling_ball<double>(c1, c2, 0.25), std::runtime_error);
    }
}

TEST(tests_bscanalysis, compute_rolling_ball_corner)
{
    using P = gbs::point<double, 2>;
    // y = 0 and the half line y = 1, x >= 2: left of x = 2 the ball tangent to y = 0 at (a, 0)
    // rests on the end (2, 1), r = ((2 - a)² + 1) / 2; beyond, it touches the line, r = 1/2.
    // It used to take the midpoint of (a, 0) and of the projection (2, 1) (parallel normals),
    // a circle tangent to neither
    const auto c1 = gbs::build_segment<double, 2>({0., 0.}, {4., 0.});
    const auto c2 = gbs::build_segment<double, 2>({2., 1.}, {5., 1.});
    for (const double a : {0., 1., 1.5, 1.99})
    {
        const auto [u2, centre, parallel] = gbs::compute_rolling_ball<double>(c1, c2, a);
        ASSERT_LT(gbs::norm(centre - P{a, ((2. - a) * (2. - a) + 1.) / 2.}), 1e-12);
        ASSERT_LT(gbs::norm(c2(u2) - P{2., 1.}), 1e-12);
    }
    for (const double a : {2.5, 3.})
    {
        const auto [u2, centre, parallel] = gbs::compute_rolling_ball<double>(c1, c2, a);
        ASSERT_TRUE(parallel);
        ASSERT_LT(gbs::norm(centre - P{a, 0.5}), 1e-12);
        ASSERT_LT(gbs::norm(c2(u2) - P{a, 1.}), 1e-9);
    }
}

TEST(tests_bscanalysis, compute_rolling_ball_cases)
{
    using P = gbs::point<double, 2>;
    const auto c1 = gbs::build_segment<double, 2>({-2., 0.}, {2., 0.}); // y = 0, p = (u1 - 2, 0)
    // the half line y = 1 - x crosses the tangent of y = 0 next to the contact: rho is +inf
    // there, among the points Brent keeps. The ball at (a, 0): r = (1 - a) / (1 + sqrt(2))
    {
        const auto c2 = gbs::build_segment<double, 2>({-1., 2.}, {2., -1.});
        for (const double a : {0., 0.9, 0.99, 0.999})
        {
            const auto [u2, centre, parallel] = gbs::compute_rolling_ball<double>(c1, c2, a + 2.);
            ASSERT_LT(gbs::norm(centre - P{a, (1. - a) / (1. + std::sqrt(2.))}), 1e-12);
        }
    }
    // a circle, not a B-spline (64 samples): radius 1 around (0, 3); the ball at (a, 0) has
    // |centre - (0, 3)| = 1 + r, r = (a² + 8) / 8
    {
        const gbs::Circle2d<double> c2{1., P{0., 3.}};
        for (const double a : {0., 0.5, 1.5})
        {
            const auto [u2, centre, parallel] = gbs::compute_rolling_ball<double>(c1, c2, a + 2.);
            const double r = (a * a + 8.) / 8.;
            ASSERT_LT(gbs::norm(centre - P{a, r}), 1e-12);
            ASSERT_NEAR(gbs::norm(c2(u2) - P{0., 3.}), 1., 1e-12);
            ASSERT_NEAR(gbs::norm(c2(u2) - centre), r, 1e-12);
        }
    }
    // crv2 passing at h = 1e-4 over p, between two samples, then back at height 0.6 (rho = 0.3
    // there, the smallest sample): the first ball, of radius 0.3, holds the close pass; the
    // projection of its centre on crv2 finds it, and the ball is r = h / 2
    {
        const double h = 1e-4;
        const gbs::BSCurve<double, 2> c2({{-1., h}, {1., h}, {1., 0.6}, {-1., 0.6}}, {0., 0., 1., 2., 3., 3.}, 1);
        for (const double x0 : {0.013, 0.031})
        {
            const auto [u2, centre, parallel] = gbs::compute_rolling_ball<double>(c1, c2, x0 + 2.);
            ASSERT_LT(gbs::norm(centre - P{x0, h / 2.}), 1e-12);
            ASSERT_LT(gbs::norm(c2(u2) - P{x0, h}), 1e-12);
        }
    }
    // the ball on the second side only: p = (1, 0) on the unit circle, crv2 nearest at
    // (-2.5, 0) inside, where the ball (r = 1.75) would hold the unit circle; outside, it rests
    // on the end (3, 3.2) of crv2: r = (2² + 3.2²) / 4, centre (1 + r, 0)
    {
        const gbs::Circle2d<double> unit{1.};
        const gbs::BSCurve<double, 2> c2({{-2.5, 0.}, {-2.5, 4.}, {3., 4.}, {3., 3.2}}, {0., 0., 1., 2., 3., 3.}, 1);
        ASSERT_LT(gbs::norm(unit(0.) - P{1., 0.}), 1e-15);
        const auto [u2, centre, parallel] = gbs::compute_rolling_ball<double>(unit, c2, 0.);
        const double r = (4. + 3.2 * 3.2) / 4.;
        ASSERT_LT(gbs::norm(centre - P{1. + r, 0.}), 1e-12);
        ASSERT_LT(gbs::norm(c2(u2) - P{3., 3.2}), 1e-12);
    }
}

TEST(tests_bscanalysis, compute_rolling_ball_splitter)
{
    // the two sides of a tbs blade passage: near the leading edge the ball rests on the end of
    // side2. There it used to throw (no root of the difference of the radii), or to return a
    // centre at 1e15 or the midpoint of two points (parallel normals at the end of side2)
    std::ifstream in(std::string(GBS_TESTS_IN_DIR) + "/rolling_ball_splitter.txt");
    ASSERT_TRUE(in.good());
    for (std::string line; in.peek() == '#';)
        std::getline(in, line);
    auto read = [&in] {
        std::size_t p, nk, np;
        in >> p >> nk >> np;
        std::vector<double> k(nk);
        for (auto &x : k)
            in >> x;
        gbs::points_vector<double, 2> poles(np);
        for (auto &q : poles)
            in >> q[0] >> q[1];
        return gbs::BSCurve<double, 2>(poles, k, p);
    };
    const auto side1 = read(), side2 = read();
    const auto [a1, b1] = side1.bounds();
    const auto [a2, b2] = side2.bounds();
    int on_corner{};
    for (int i = 0; i <= 100; ++i)
    {
        const double u1 = a1 + (b1 - a1) * i / 100.;
        const auto [u2, centre, parallel] = gbs::compute_rolling_ball<double>(side1, side2, u1);
        const double r = gbs::norm(centre - side1(u1));
        ASSERT_LT(r, 1.);
        ASSERT_NEAR(gbs::norm(centre - side2(u2)), r, 1e-12 * (1. + r)); // touches side2 at u2
        // no point of either side inside the ball
        ASSERT_GE(gbs::extrema_curve_point(side1, centre, 1e-10)[1], r * (1. - 1e-9));
        ASSERT_GE(gbs::extrema_curve_point(side2, centre, 1e-10)[1], r * (1. - 1e-9));
        on_corner += u2 == a2;
    }
    ASSERT_GE(on_corner, 10); // the leading part rests on the end of side2
}

TEST(tests_bscanalysis, max_curvature_pos_narrow_peak)
{
    // a quadratic B-spline of 300 spans along y = 0 with a corner rounded over 2 spans: the peak
    // of curvature spans 1/150 of the range, between samples taken uniformly 64 times
    gbs::points_vector<double, 2> poles;
    for (int i = 0; i <= 150; ++i)
        poles.push_back({double(i), 0.});
    for (int i = 1; i <= 151; ++i)
        poles.push_back({150., double(i)});
    const auto knots = gbs::build_simple_mult_flat_knots<double>(0., 1., poles.size(), 2);
    const gbs::BSCurve<double, 2> crv{poles, knots, 2};
    double c_max{}, u_max{};
    for (int k = 0; k <= 400000; ++k)
        if (const double u = k / 400000., c = gbs::curvature(crv, u); c > c_max)
            c_max = c, u_max = u;
    const auto [u, c] = gbs::max_curvature_pos(crv, 1e-8);
    ASSERT_NEAR(c, c_max, 1e-6 * c_max);
    ASSERT_NEAR(u, u_max, 1e-4);
    ASSERT_THROW(gbs::max_curvature_pos(crv, 0.5, 0.5, 1e-8), std::invalid_argument);
    ASSERT_THROW(gbs::max_curvature_pos(crv, 0.6, 0.5, 1e-8), std::invalid_argument);
}

TEST(tests_bscanalysis, curvature_exact)
{
    // nearly straight quadratic Bézier P0 P1 P2 with C' and C'' nearly parallel (accelerating
    // along x): its curvature is 4 |A x B| / |C'|^3 with A = P1 - P0, B = P2 - P1. Lagrange's
    // identity cancels to about 1e-8 there; the 2D determinant is exact
    const double delta = 1e-9;
    const gbs::BSCurve<double, 2> crv{{{0., 0.}, {0.1, 0.}, {1., delta}}, {0., 0., 0., 1., 1., 1.}, 2};
    for (const double t : {0., 0.3, 0.5, 1.})
    {
        const auto d1 = crv(t, 1);
        const double exact = 4. * 0.1 * delta / std::pow(gbs::norm(d1), 3);
        ASSERT_NEAR(gbs::curvature(crv, t), exact, 1e-9 * exact);
    }
}

TEST(tests_bscanalysis, closest_curve)
{
    std::vector<gbs::BSCurve<double, 3>> curves{gbs::build_segment<double, 3>({0., 0., 0.}, {1., 0., 0.}),
                                                gbs::build_segment<double, 3>({0., 1., 0.}, {1., 1., 0.}),
                                                gbs::build_segment<double, 3>({0., 2., 0.}, {1., 2., 0.})};
    const gbs::point<double, 3> p{0.5, 1.3, 0.};
    ASSERT_EQ(std::distance(curves.begin(), gbs::closest_curve(p, curves.begin(), curves.end())), 1);
    std::vector<const gbs::Curve<double, 3> *> ptrs{&curves[2], &curves[0], &curves[1]};
    ASSERT_EQ(std::distance(ptrs.begin(), gbs::closest_p_curve(p, ptrs.begin(), ptrs.end())), 2);
    // a curve that cannot be projected (unbounded line): the exception reaches the caller
    std::vector<const gbs::Curve<double, 3> *> with_line{&curves[0], nullptr};
    gbs::Line<double, 3> line{gbs::point<double, 3>{0., 0., 0.}, gbs::point<double, 3>{1., 0., 0.}};
    with_line[1] = &line;
    ASSERT_THROW(gbs::closest_p_curve(p, with_line.begin(), with_line.end()), std::invalid_argument);
    // an empty range: the end
    std::vector<gbs::BSCurve<double, 3>> none;
    ASSERT_TRUE(gbs::closest_curve(p, none.begin(), none.end()) == none.end());
}

TEST(tests_bscanalysis, signed_curvature)
{
    // a circle of radius 2 run counterclockwise: +1/2; clockwise: -1/2; the unsigned curvature its absolute value
    auto circle = gbs::build_circle<double, 2>(2., {1., -1.});
    auto clockwise = circle;
    clockwise.reverse();
    for (double t : {0.1, 0.4, 0.8})
    {
        const auto [a, b] = circle.bounds();
        const double u = a + t * (b - a);
        ASSERT_NEAR(gbs::signed_curvature(circle, u), 0.5, 1e-9);
        ASSERT_NEAR(gbs::signed_curvature(clockwise, u), -0.5, 1e-9);
        ASSERT_EQ(gbs::curvature(clockwise, u), std::abs(gbs::signed_curvature(clockwise, u)));
        // the left normal points to the centre on the counterclockwise circle
        const auto n = gbs::normal_direction(circle, u);
        const auto p = circle(u);
        ASSERT_NEAR(p[0] + 2. * n[0], 1., 1e-9);
        ASSERT_NEAR(p[1] + 2. * n[1], -1., 1e-9);
    }
    ASSERT_EQ(gbs::signed_curvature<double>({0., 0.}, {1., 0.}), 0.);
}

namespace
{
    // a B-spline that says it has no second derivatives (as a curve on such a surface)
    struct FirstOrderCurve : gbs::BSCurve<double, 2>
    {
        using gbs::BSCurve<double, 2>::BSCurve;
        auto has_second_derivatives() const -> bool override { return false; }
    };

    auto tooth(const gbs::CurvatureCombPiece<double> &piece, size_t i) -> double
    {
        const auto &p = piece.teeth[3 * i], &tip = piece.envelope[i];
        return std::hypot(tip[0] - p[0], tip[1] - p[1]);
    }
}

TEST(tests_bscanalysis, curvature_comb)
{
    // a circle run counterclockwise: every tooth outwards, as long as the longest
    auto circle = gbs::build_circle<double, 2>(2., {1., -1.});
    const auto single = gbs::curvature_comb(circle, 0.5, 40);
    ASSERT_GE(single.envelope.size(), 41u);
    ASSERT_EQ(single.teeth.size(), 3 * single.envelope.size());
    for (const auto &tip : single.envelope)
        ASSERT_NEAR(std::hypot(tip[0] - 1., tip[1] + 1.), 2.5, 1e-9);

    // y = x^2 run left to right: teeth below, the longest at the vertex; refined where the curve turns, at the
    // vertex, sparser on the flanks
    gbs::BSCurve<double, 2> parabola{gbs::points_vector<double, 2>{{-1., 1.}, {0., -1.}, {1., 1.}}, {0., 0., 0., 1., 1., 1.}, 2};
    const auto comb = gbs::curvature_comb(parabola, 0.5, 20);
    double at_vertex{};
    for (size_t i = 0; i < comb.envelope.size(); ++i)
    {
        const auto &p = comb.teeth[3 * i];
        ASSERT_LT(comb.envelope[i][1], p[1]);
        if (std::abs(p[0]) < 1e-12)
            at_vertex = tooth(comb, i);
    }
    ASSERT_NEAR(at_vertex, 0.5, 1e-12);
    const auto coarse = gbs::curvature_comb(parabola, 0.5, 4);
    ASSERT_GT(coarse.envelope.size(), 5u);
    auto gap_near = [&](double x) {
        double best = 1e300;
        for (size_t i = 0; i + 1 < coarse.envelope.size(); ++i)
        {
            const auto &p = coarse.teeth[3 * i], &q = coarse.teeth[3 * (i + 1)];
            if (std::abs(0.5 * (p[0] + q[0]) - x) < 0.15)
                best = std::min(best, std::hypot(q[0] - p[0], q[1] - p[1]));
        }
        return best;
    };
    ASSERT_LT(gap_near(0.), 0.5 * gap_near(0.9));
}

TEST(tests_bscanalysis, curvature_comb_chain)
{
    // an S then a circle of radius 0.1, tighter than any point of the S: one scale for the chain, the circle's
    // teeth `length` long, the S's shorter in the ratio of the square roots of the curvatures
    gbs::BSCurve<double, 2> s{gbs::points_vector<double, 2>{{0., 0.}, {1., 1.}, {2., -1.}, {3., 0.}},
                              {0., 0., 0., 0., 1., 1., 1., 1.}, 3};
    auto circle = gbs::build_circle<double, 2>(0.1, {3., 0.1});
    const auto comb = gbs::curvature_comb({&s, &circle}, 1.); // T deduced from a braced list
    ASSERT_EQ(comb.size(), 2u);
    for (size_t i = 0; i < comb[1].envelope.size(); ++i)
        ASSERT_NEAR(tooth(comb[1], i), 1., 1e-9);
    const auto [a, b] = s.bounds();
    double longest{};
    for (size_t i = 0; i < comb[0].envelope.size(); ++i)
        longest = std::max(longest, tooth(comb[0], i));
    double k_max{};
    for (int i = 0; i <= 2000; ++i)
        k_max = std::max(k_max, gbs::curvature(s, a + (b - a) * i / 2000.));
    ASSERT_NEAR(longest, std::sqrt(k_max / 10.), 1e-3);
    ASSERT_LT(longest, 0.9);

    // the S turns right then left (both ends along (1, 1)): its teeth inwards (to the left of the run, k < 0) at
    // the start, outwards (to the right, k > 0) at the end
    const auto &first = comb[0].envelope.front(), &last = comb[0].envelope.back();
    ASSERT_GT(first[1] - first[0], 0.);
    ASSERT_LT(last[1] - (last[0] - 3.), 0.);

    // the same chain from a vector of pointers, an array of curves, shared pointers: the same comb
    const std::vector<const gbs::Curve<double, 2> *> pointers{&s, &circle};
    ASSERT_EQ(gbs::curvature_comb(pointers, 1.)[0].envelope, comb[0].envelope);
    const std::array<gbs::BSCurve<double, 2>, 2> by_value{s, s};
    ASSERT_EQ(gbs::curvature_comb(by_value, 1.).size(), 2u);
    const std::vector<std::shared_ptr<gbs::BSCurve<double, 2>>> shared{std::make_shared<gbs::BSCurve<double, 2>>(s)};
    ASSERT_EQ(gbs::curvature_comb(shared, 1.)[0].envelope.size(), gbs::curvature_comb(s, 1.).envelope.size());
    ASSERT_THROW(gbs::curvature_comb(std::vector<const gbs::Curve<double, 2> *>{&s, nullptr}, 1.), std::invalid_argument);

    // a point gets no tooth, nor a piece without second derivatives; segments, teeth of zero length
    gbs::BSCurve<double, 2> point{gbs::points_vector<double, 2>{{3., 0.}, {3., 0.}}, {0., 0., 1., 1.}, 1};
    const FirstOrderCurve first_order{gbs::points_vector<double, 2>{{0., 0.}, {1., 1.}, {2., 0.}}, {0., 0., 0., 1., 1., 1.}, 2};
    const auto skipped = gbs::curvature_comb({&s, &point, &first_order}, 1.);
    ASSERT_TRUE(skipped[1].envelope.empty());
    ASSERT_TRUE(skipped[2].envelope.empty());
    ASSERT_FALSE(skipped[0].envelope.empty());
    gbs::BSCurve<double, 2> seg1{gbs::points_vector<double, 2>{{0., 0.}, {1., 0.}}, {0., 0., 1., 1.}, 1};
    gbs::BSCurve<double, 2> seg2{gbs::points_vector<double, 2>{{1., 0.}, {1., 2.}}, {0., 0., 2., 2.}, 1};
    const auto flat = gbs::curvature_comb({&seg1, &seg2}, 1.);
    for (const auto &piece : flat)
    {
        ASSERT_FALSE(piece.envelope.empty());
        for (size_t i = 0; i < piece.envelope.size(); ++i)
            ASSERT_EQ(tooth(piece, i), 0.);
    }

    // the threshold on the speed is relative: the S scaled down by 1e-9, its teeth unchanged once scaled back
    auto tiny_poles = s.poles();
    for (auto &p : tiny_poles)
        p = {p[0] * 1e-9, p[1] * 1e-9};
    gbs::BSCurve<double, 2> tiny{tiny_poles, s.knotsFlats(), 3};
    const auto small = gbs::curvature_comb(tiny, 1e-9);
    const auto big = gbs::curvature_comb(s, 1.);
    ASSERT_EQ(small.envelope.size(), big.envelope.size());
    for (size_t i = 0; i < big.envelope.size(); ++i)
        ASSERT_NEAR(small.envelope[i][0] * 1e9, big.envelope[i][0], 1e-6);
}

TEST(tests_bscanalysis, area_properties)
{
    // a circle far from the origin, run counterclockwise: its area and its centre
    const double pi = std::numbers::pi;
    auto circle = gbs::build_circle<double, 2>(2., {1e3, -5e2});
    const auto disk = gbs::area_properties(circle);
    ASSERT_NEAR(disk.area, 4. * pi, 1e-9);
    ASSERT_NEAR(disk.centroid[0], 1e3, 1e-9);
    ASSERT_NEAR(disk.centroid[1], -5e2, 1e-9);

    // the unit square as four segments, run clockwise: the area negative, the same centroid
    auto seg = [](gbs::point<double, 2> a, gbs::point<double, 2> b) {
        return gbs::BSCurve<double, 2>{gbs::points_vector<double, 2>{a, b}, {0., 0., 1., 1.}, 1};
    };
    const std::vector<gbs::BSCurve<double, 2>> square{seg({0., 0.}, {0., 1.}), seg({0., 1.}, {1., 1.}),
                                                      seg({1., 1.}, {1., 0.}), seg({1., 0.}, {0., 0.})};
    const auto sq = gbs::area_properties(square);
    ASSERT_NEAR(sq.area, -1., 1e-14);
    ASSERT_NEAR(sq.centroid[0], 0.5, 1e-14);
    ASSERT_NEAR(sq.centroid[1], 0.5, 1e-14);

    // half a circle of radius 1 (an exact rational arc through (1, 0), (0, 1), (-1, 0)), closed by its diameter:
    // the centroid at 4 / (3 pi) above it
    auto half = gbs::build_circle<double, 2>(1.);
    const auto [u0, u1] = half.bounds();
    gbs::BSCurveRational<double, 2> arc = half;
    arc.trim(u0, 0.5 * (u0 + u1));
    const auto hd = gbs::area_properties(arc);
    ASSERT_NEAR(hd.area, 0.5 * pi, 1e-9);
    ASSERT_NEAR(hd.centroid[0], 0., 1e-9);
    ASSERT_NEAR(hd.centroid[1], 4. / (3. * pi), 1e-9);
    // the same half, trimmed without dropping its knots: integrated within its bounds only
    gbs::BSCurveRational<double, 2> soft = half;
    soft.trim(u0, 0.5 * (u0 + u1), false);
    const auto hs = gbs::area_properties(soft);
    ASSERT_NEAR(hs.area, hd.area, 1e-14);
    ASSERT_NEAR(hs.centroid[0], hd.centroid[0], 1e-14);
    ASSERT_NEAR(hs.centroid[1], hd.centroid[1], 1e-14);
    // the two halves both run from (1, 0) to (-1, 0), as the sides of a section from its leading edge: the second
    // one runs backwards, from where the first one ends
    gbs::BSCurveRational<double, 2> lower = half;
    lower.trim(0.5 * (u0 + u1), u1);
    const auto lower_back = lower.reversed();
    const auto whole = gbs::area_properties({&arc, &lower_back});
    ASSERT_NEAR(whole.area, pi, 1e-14);
    ASSERT_NEAR(whole.centroid[0], 0., 1e-14);
    ASSERT_NEAR(whole.centroid[1], 0., 1e-14);

    // a chain with gaps (two sides of a triangle and a point): closed by segments; a braced list, pointers
    const auto a = seg({0., 0.}, {3., 0.}), b = seg({3., 0.}, {0., 3.});
    const auto point = seg({0., 3.}, {0., 3.});
    const auto tri = gbs::area_properties({&b, &point, &a});
    ASSERT_NEAR(tri.area, 4.5, 1e-14);
    ASSERT_NEAR(tri.centroid[0], 1., 1e-14);
    ASSERT_NEAR(tri.centroid[1], 1., 1e-14);
    const std::vector<std::shared_ptr<gbs::BSCurve<double, 2>>> shared{std::make_shared<gbs::BSCurve<double, 2>>(a),
                                                                       std::make_shared<gbs::BSCurve<double, 2>>(b)};
    ASSERT_NEAR(gbs::area_properties(shared).area, 4.5, 1e-14);
    // the pieces either way: the same chain; the order of the pieces gives its direction
    const auto a_back = a.reversed(), b_back = b.reversed();
    ASSERT_NEAR(gbs::area_properties({&b_back, &point, &a_back}).area, 4.5, 1e-14);
    const auto clockwise = gbs::area_properties({&a, &point, &b});
    ASSERT_NEAR(clockwise.area, -4.5, 1e-14);
    ASSERT_NEAR(clockwise.centroid[0], 1., 1e-14);
    ASSERT_NEAR(clockwise.centroid[1], 1., 1e-14);
    // max_gap: the segment closing three sides of a square is its fourth side
    const auto bottom = seg({0., 0.}, {1., 0.}), right = seg({1., 0.}, {1., 1.}), top = seg({1., 1.}, {0., 1.});
    ASSERT_NEAR(gbs::area_properties({&bottom, &right, &top}).area, 1., 1e-14);
    ASSERT_NEAR(gbs::area_properties({&bottom, &right, &top}, 1.5).area, 1., 1e-14);
    ASSERT_THROW(gbs::area_properties({&bottom, &right, &top}, 0.5), std::invalid_argument);
    // a view that is not a range when const: the pieces kept by a filter
    auto kept = square | std::views::filter([](const auto &) { return true; });
    ASSERT_NEAR(gbs::area_properties(kept).area, -1., 1e-14);
    // a null curve
    ASSERT_THROW(gbs::area_properties(std::vector<const gbs::Curve<double, 2> *>{&a, nullptr}), std::invalid_argument);
    ASSERT_THROW(gbs::area_properties(std::vector<std::shared_ptr<gbs::BSCurve<double, 2>>>{nullptr}),
                 std::invalid_argument);

    // the triangle as a CurveComposite, a join off the dyadic points: integrated curve by curve
    auto b_fast = std::make_shared<gbs::BSCurve<double, 2>>(b);
    b_fast->changeBounds(0., 0.37);
    const gbs::CurveComposite<double, 2> contour{{std::make_shared<gbs::BSCurve<double, 2>>(a), b_fast,
                                                  std::make_shared<gbs::BSCurve<double, 2>>(seg({0., 3.}, {0., 0.}))}};
    const auto ct = gbs::area_properties(contour);
    ASSERT_NEAR(ct.area, 4.5, 1e-14);
    ASSERT_NEAR(ct.centroid[0], 1., 1e-14);
    ASSERT_NEAR(ct.centroid[1], 1., 1e-14);

    // a cubic with an interior knot: the integration span by span is exact, against a Gauss rule on each span (10
    // points, exact up to degree 19: x y x' is of degree 8)
    gbs::BSCurve<double, 2> wave{gbs::points_vector<double, 2>{{0., 0.}, {1., 2.}, {2., -1.}, {3., 2.}, {4., 0.}},
                                 {0., 0., 0., 0., 0.3, 1., 1., 1., 1.}, 3};
    const auto w = gbs::area_properties(wave);
    double ref_a{}, ref_x{}; // the area under the wave, closed by the axis, and its moment
    for (const auto [ua, ub] : {std::pair{0., 0.3}, std::pair{0.3, 1.}})
    {
        using gauss10 = boost::math::quadrature::gauss<double, 10>;
        ref_a += gauss10::integrate([&](double u) { return wave.value(u)[1] * wave.value(u, 1)[0]; }, ua, ub);
        ref_x += gauss10::integrate([&](double u) { return wave.value(u)[0] * wave.value(u)[1] * wave.value(u, 1)[0]; }, ua, ub);
    }
    ASSERT_NEAR(w.area, -ref_a, 1e-14);
    ASSERT_NEAR(w.centroid[0], ref_x / ref_a, 1e-14);
    // the same curve of degree 6 (Kronrod exact) and 10 (Kronrod halved until it converges): the same properties
    for (size_t step : {3u, 7u})
    {
        auto raised = wave;
        raised.increaseDegree(step);
        const auto r = gbs::area_properties(raised);
        ASSERT_NEAR(r.area, w.area, 1e-13);
        ASSERT_NEAR(r.centroid[0], w.centroid[0], 1e-13);
        ASSERT_NEAR(r.centroid[1], w.centroid[1], 1e-13);
    }

    // no area
    ASSERT_THROW(gbs::area_properties(a), std::invalid_argument);
    ASSERT_THROW(gbs::area_properties(std::vector<gbs::BSCurve<double, 2>>{}), std::invalid_argument);
    // a curve and its reverse: their rounding is no area, whatever the unit and the distance to the origin; a thin
    // strip there is one
    for (double s : {1e-3, 1., 1e3})
    {
        auto at = [s](double x, double y) { return gbs::point<double, 2>{s * (x + 1e3), s * (y - 5e2)}; };
        const gbs::BSCurve<double, 2> c{gbs::points_vector<double, 2>{at(0., 0.), at(1., 2.), at(2., -1.), at(3., 2.), at(4., 0.)},
                                        {0., 0., 0., 0., 0.3, 1., 1., 1., 1.}, 3};
        const auto c_back = c.reversed();
        ASSERT_THROW(gbs::area_properties({&c, &c_back}), std::invalid_argument);
        const auto disk_edge = gbs::build_circle<double, 2>(2. * s, at(0., 0.));
        const auto disk_back = disk_edge.reversed();
        ASSERT_THROW(gbs::area_properties({&disk_edge, &disk_back}), std::invalid_argument);
        const auto p0 = at(0., 0.), p1 = at(10., 0.), p2 = at(10., 1e-5), p3 = at(0., 1e-5);
        const std::vector<gbs::BSCurve<double, 2>> strip{seg(p0, p1), seg(p1, p2), seg(p2, p3), seg(p3, p0)};
        ASSERT_NEAR(gbs::area_properties(strip).area / (s * s * 1e-4), 1., 1e-6);
    }
}

namespace
{
    // a curve that counts its evaluations, through value or derivatives
    template <typename Base>
    struct EvaluationCounter : Base
    {
        explicit EvaluationCounter(const Base &base) : Base{base} {}
        mutable size_t evaluations{};
        auto value(double u, size_t d = 0) const -> std::array<double, 2> override
        {
            ++evaluations;
            return Base::value(u, d);
        }
        auto derivatives(double u, size_t n, std::array<double, 2> *CK) const -> void override
        {
            ++evaluations;
            Base::derivatives(u, n, CK);
        }
    };
}

TEST(tests_bscanalysis, area_properties_narrow_span)
{
    // an airfoil's upper side (tbs, NACA 4412, from PR #188), its lower bound an ulp below its first knot (from a
    // projection), closed by its chord: that knot is no break (smooth_spans), so no span is an ulp wide; Boost's
    // Gauss-Kronrod bisected such a span to its depth limit (3.9 million evaluations, 90 ms; 5.9 million, 345 ms,
    // for the same curve as a NURBS)
    const gbs::points_vector<double, 2> poles{{1.0001665262873147, 0.0012489471548600977},
                                              {0.9828639016423999, 0.0060343299227571118},
                                              {0.94825865235257012, 0.015605095458551137},
                                              {0.73544649012701113, 0.066789012317935728},
                                              {0.29725476526187655, 0.12441760464154356},
                                              {0.073838852703787677, 0.061977982081785149},
                                              {0.032378423293418349, 0.038275298383154681}};
    const std::vector<double> knots{0.061462396870544334, 0.061462396870544334, 0.061462396870544334,
                                    0.061462396870544334, 0.28327254123129286,  0.50508268559204128,
                                    0.83841601892537454,  1., 1., 1., 1.};
    const double u0 = 0.061462396870544327; // an ulp below the first knot
    const gbs::BSCurve<double, 2> plain{poles, knots, 3};
    EvaluationCounter<gbs::BSCurve<double, 2>> side{plain};
    side.trim(u0, 1., false);
    EvaluationCounter<gbs::BSCurveRational<double, 2>> nurbs{gbs::BSCurveRational<double, 2>{plain}};
    nurbs.trim(u0, 1., false);
    const auto polynomial = gbs::area_properties(side);
    const auto rational = gbs::area_properties(nurbs);
    ASSERT_LT(side.evaluations, 200u);
    ASSERT_LT(nurbs.evaluations, 200u);

    // against a Gauss rule of 15 points on each knot span (exact: the integrands are of degree 8 at most), closed by
    // the chord (a segment through o: nothing for 2A)
    const auto o = plain.value(u0), e = plain.value(1.);
    double a2{}, ax{}, ay{};
    const std::array<double, 5> spans{u0, knots[4], knots[5], knots[6], 1.};
    for (size_t i = 0; i + 1 < spans.size(); ++i)
    {
        auto moment = [&](auto term) {
            return boost::math::quadrature::gauss<double, 15>::integrate(
                [&](double u) {
                    const auto p = plain.value(u), d = plain.value(u, 1);
                    return term(p[0] - o[0], p[1] - o[1], d[0], d[1]);
                },
                spans[i], spans[i + 1]);
        };
        a2 += moment([](double x, double y, double dx, double dy) { return x * dy - y * dx; });
        ax += moment([](double x, double, double, double dy) { return x * x * dy; });
        ay -= moment([](double, double y, double dx, double) { return y * y * dx; });
    }
    const double ex = e[0] - o[0], ey = e[1] - o[1];
    ax -= ey * ex * ex / 3.;
    ay += ex * ey * ey / 3.;
    for (const auto &props : {polynomial, rational})
    {
        ASSERT_NEAR(props.area, a2 / 2., 1e-14);
        ASSERT_NEAR(props.centroid[0], o[0] + ax / a2, 1e-13);
        ASSERT_NEAR(props.centroid[1], o[1] + ay / a2, 1e-13);
    }
}

TEST(tests_bscanalysis, deviation_based_params_breaks_once)
{
    // a knot on one of the n uniform parameters (0.5 with n = 3): once (twice before)
    const gbs::BSCurve<double, 2> c{gbs::points_vector<double, 2>{{0., 0.}, {1., 1.}, {2., -1.}, {3., 1.}, {4., 0.}},
                                    {0., 0., 0., 0., 0.5, 1., 1., 1., 1.}, 3};
    for (size_t n : {3u, 5u, 11u})
    {
        const auto u = gbs::deviation_based_params<double, 2>(c, n, 0.05);
        ASSERT_EQ(std::ranges::adjacent_find(u), u.end());
    }
    // a composite of two segments of the same length: its join is the middle parameter, and
    // discretize gave its point twice
    auto s1 = std::make_shared<gbs::BSCurve<double, 2>>(gbs::build_segment<double, 2>({0., 0.}, {1., 0.}));
    auto s2 = std::make_shared<gbs::BSCurve<double, 2>>(gbs::build_segment<double, 2>({1., 0.}, {1., 1.}));
    const gbs::CurveComposite<double, 2> comp{std::vector<std::shared_ptr<gbs::Curve<double, 2>>>{s1, s2}};
    ASSERT_EQ((gbs::deviation_based_params<double, 2>(comp, 3, 0.05)), (std::list<double>{0., 1., 2.}));
    ASSERT_EQ((gbs::discretize<double, 2>(comp, 3, 0.05).size()), 3u);
}
