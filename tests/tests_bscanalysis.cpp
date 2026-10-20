#include <doctest_gtest.hpp>
#include <fstream>
#include <numbers>
#include <gbs/bscanalysis.h>
#include <gbs/bscbuild.h>
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
