#include <gtest/gtest.h>

import vecop;
#include <numbers>

const double tol = 1e-10;

using namespace gbs;

const std::array<double,3> a{1.,1.,1.};
const std::array<double,3> b{2.,2.,2.};
const std::array<double,3> c{3.,3.,3.};
const std::array<double,3> x{1.,0.,0.};
const std::array<double,3> y{0.,1.,0.};
const std::array<double,3> z{0.,0.,1.};

TEST(tests_vecop, sum)
{
    auto c =a+b;
    std::for_each(c.begin(),c.end(),[](const auto &v){ASSERT_DOUBLE_EQ(v,3.);});
    auto d =a-b;
    std::for_each(d.begin(),d.end(),[](const auto &v){ASSERT_DOUBLE_EQ(v,-1.);});
}

TEST(tests_vecop, norms)
{
    ASSERT_DOUBLE_EQ(sq_distance(a,b),3.);
    ASSERT_DOUBLE_EQ(distance(a,b),std::sqrt(3.));
}

TEST(tests_vecop, lentgh)
{
    std::vector< std::array<double,3> > pts =
    {
        {0.,0.,0.},
        {1.,0.,0.},
        {1.,1.,0.},
        {1.,1.,2.},
    };

    ASSERT_DOUBLE_EQ(length(pts),4.);
}

TEST(tests_vecop, cross_product)
{
    ASSERT_LT(distance(x^y,z),tol);
    ASSERT_LT(distance(y^z,x),tol);
    ASSERT_LT(distance(z^x,y),tol);
}

TEST(tests_vecop, angle)
{
    using T = double;
    ASSERT_DOUBLE_EQ(gbs::angle<T>({1,0,0}, {0,1,0}), std::numbers::pi_v<T>/2);
    ASSERT_DOUBLE_EQ(gbs::angle<T>({0,1,0}, {1,0,0}), std::numbers::pi_v<T>/2);
    ASSERT_DOUBLE_EQ(gbs::angle<T>({1,1,0}, {0,1,0}), std::numbers::pi_v<T>/4);
    ASSERT_DOUBLE_EQ(gbs::angle<T>({1,0,0}, {-1,1,0}), 3*std::numbers::pi_v<T>/4);
    ASSERT_DOUBLE_EQ(gbs::angle<T>({1,0,0}, {-1,-1,0}), 3*std::numbers::pi_v<T>/4);
    ASSERT_DOUBLE_EQ(gbs::angle<T>({0.0, 0.0, 1.0}, {0.0, -0.9001691682414991, 0.7553316170686007}), 50/180.*std::numbers::pi_v<T>);
}

TEST(tests_vecop, dot_product)
{
    ASSERT_LT(std::fabs(a*b-6.),tol);
}

TEST(tests_vecop, scale)
{
    auto as_d = a * 2.;
    auto as_f = a * 2.f;
    std::for_each(as_d.begin(),as_d.end(),[](const auto &d_){ASSERT_DOUBLE_EQ(d_,2.);});
    std::for_each(as_f.begin(),as_f.end(),[](const auto &d_){ASSERT_DOUBLE_EQ(d_,2.);});
}

TEST(tests_vecop, delta)
{
    std::vector<double> l = {0.,1.,2.,3.,4.};
    auto d = delta(l);
    ASSERT_EQ(d.size(),l.size()-1);
    std::for_each(d.begin(),d.end(),[](const auto &d_){ASSERT_LT(std::fabs(d_-1.),tol);});

    std::vector<std::array<double,3> > pts = {a,b};
    auto d_pts = delta(pts);
    ASSERT_EQ(d_pts.size(),pts.size()-1);
    std::for_each(d_pts.begin(),d_pts.end(),[](const auto &d_){ASSERT_LT(distance(d_,{1.,1.,1.}),tol);});
    
}

TEST(tests_vecop, coef_div)
{
    std::vector<std::array<double,3> > pts = {a,b};
    std::vector<double> u = {1.,2.};
    
    auto d = pts / u;
    ASSERT_EQ(d.size(),u.size());
    std::for_each(d.begin(),d.end(),[](const auto &d_){ASSERT_LT(distance(d_,{1.,1.,1.}),tol);});
}