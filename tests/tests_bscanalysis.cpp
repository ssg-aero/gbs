#include <doctest_gtest.hpp>
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
    auto e = gbs::build_ellipse<double,2>(1.,0.2);
    auto [u, cur] = gbs::max_curvature_pos(e, 0., 1., 1e-63);
    // ASSERT_NEAR(u,0.5,1e-6);
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
