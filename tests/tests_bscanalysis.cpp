#include <gtest/gtest.h>
#include <numbers>
#include <gbs/bscanalysis.h>
#include <gbs/bscbuild.h>
#include <gbs-render/vtkcurvesrender.h>
using gbs::operator-;
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

    ASSERT_NEAR(f_u(  PI),0.5,1e-5);
    ASSERT_NEAR(f_u(2*PI),1.0,1e-5);

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
    auto [u, cur] = gbs::max_curvature_pos(e, 0., 1., 1e-3);
    // ASSERT_NEAR(u,0.5,1e-6);
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