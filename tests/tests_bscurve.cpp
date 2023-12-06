#include <gtest/gtest.h>
#include <gbs/bscurve.h>
#include <gbs/curvereparametrized.h>
#include <gbs/bscanalysis.h>
#include <gbs/bscbuild.h>
#include <gbs/knotsfunctions.h>

#include <numbers>
#include <chrono>

import vecop;

const double tol = 1e-10;
const double tol_confusion = 1e-7;


using gbs::operator-;

TEST(tests_bscurve, perf_long_double)
{
    using T = long double;
    std::vector<T> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<T,3> > poles =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;
    double u = 2.3;

    auto c1_3d_dp = gbs::BSCurve<T,3>(poles,k,p);

    const auto count_max = 10000000;

    const auto t1 = std::chrono::high_resolution_clock::now();
    auto count = count_max;
    while(count)
    {
        u = static_cast<T>(count) / static_cast<T>(count_max) * 5.;
        c1_3d_dp.value(u);
        count--;
    }
    const auto t2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<T, std::milli> ms = t2 - t1;
    std::cout << std::fixed 
                  << "gbs took " << ms.count() << " ms\n";

}

TEST(tests_bscurve, perf_double)
{
    std::vector<double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<double,3> > poles =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;
    double u = 2.3;

    auto c1_3d_dp = gbs::BSCurve<double,3>(poles,k,p);

    const auto count_max = 10000000;

    const auto t1 = std::chrono::high_resolution_clock::now();
    auto count = count_max;
    while(count)
    {
        u = double(count) / double(count_max) * 5.;
        c1_3d_dp.value(u);
        count--;
    }
    const auto t2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> ms = t2 - t1;
    std::cout << std::fixed 
                  << "gbs took " << ms.count() << " ms\n";

}

TEST(tests_bscurve, perf_float)
{
    std::vector<float> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<float,3> > poles =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;
    float u = 2.3;

    auto c1_3d_dp = gbs::BSCurve<float,3>(poles,k,p);

    const auto count_max = 10000000;

    const auto t1 = std::chrono::high_resolution_clock::now();
    auto count = count_max;
    while(count)
    {
        u = double(count) / double(count_max) * 5.;
        c1_3d_dp.value(u);
        count--;
    }
    const auto t2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> ms = t2 - t1;
    std::cout << std::fixed 
                  << "gbs float took " << ms.count() << " ms\n";

}
TEST(tests_bscurve, curve_values)
{
    std::vector<float> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<float,3> > poles =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;
    float u = 2.3;

    auto c1_3d_dp = gbs::BSCurve<float,3>(poles,k,p);

    size_t n{1000};
    auto u_lst = gbs::make_range(k.front(), k.back(), n);
    const auto t1 = std::chrono::high_resolution_clock::now();
    auto pts = c1_3d_dp.values(u_lst);
    const auto t2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> ms = t2 - t1;
    std::cout << std::fixed 
                  << "gbs float took " << ms.count() << " ms\n";
                  
    for (size_t i{}; i < n; i++)
    {
        ASSERT_LT( gbs::distance( c1_3d_dp(u_lst[i]), pts[i] ), tol);
    }

}
TEST(tests_bscurve, curve_parametrization)
{
    std::vector<std::array<double,3> > pt =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    auto k1 = gbs::curve_parametrization(pt, gbs::KnotsCalcMode::EQUALLY_SPACED);
    std::for_each(k1.begin(), k1.end(), [](auto k_) { printf("k=%f\n", k_); });

    auto k2 = gbs::curve_parametrization(pt, gbs::KnotsCalcMode::CHORD_LENGTH);
    std::for_each(k2.begin(), k2.end(), [](auto k_) { printf("k=%f\n", k_); });

    auto k3 = gbs::curve_parametrization(pt, gbs::KnotsCalcMode::CENTRIPETAL);
    std::for_each(k3.begin(), k3.end(), [](auto k_) { printf("k=%f\n", k_); });
}

TEST(tests_bscurve, curve_increase_degree)
{
    std::vector<float> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<float,3> > poles =
    {
        {0.,0.,0.},
        {0.,1.,0.},
        {1.,1.,0.},
        {1.,1.,1.},
        {1.,1.,2.},
        {3.,1.,1.},
        {0.,4.,1.},
    };
    size_t p = 2;
    auto c1_3d_dp = gbs::BSCurve<float,3>(poles,k,p);
    auto c2_3d_dp = gbs::BSCurve<float,3>(poles,k,p);
    c1_3d_dp.increaseDegree();
    auto u = gbs::make_range(k.front(),k.back(),100);
    for(auto u_ : u) 
    {
        ASSERT_LT(gbs::norm(c1_3d_dp.value(u_) - c2_3d_dp.value(u_)), 1e-6);
    }
}

TEST(tests_bscurve, curve_reparametrized)
{
    auto circle = gbs::build_circle<float,2>(1.);
    
    auto f_u = gbs::BSCfunction{ gbs::abs_curv<float,2,10>(circle,36*30) };

    gbs::CurveReparametrized<float,2> reparam(std::make_shared<gbs::BSCurveRational<float,2>>(circle),f_u);

    auto tol = 1e-6;

    ASSERT_NEAR(reparam(std::numbers::pi/2.)[1], 1.,tol);
    ASSERT_NEAR(reparam(std::numbers::pi/2.)[0], 0.,tol);
    ASSERT_NEAR(reparam(std::numbers::pi   )[1], 0.,tol);
    ASSERT_NEAR(reparam(std::numbers::pi   )[0],-1.,tol);
}

TEST(tests_bscurve, eval_except)
{
    auto crv =  gbs::build_segment<double,2>({0.,0.},{1.,0.});
    ASSERT_THROW(crv.value(2.), gbs::OutOfBoundsCurveEval<double>);
}