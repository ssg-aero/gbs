#include <gtest/gtest.h>
#include <gbs/bscurve.h>
#include <gbs/bscbuild.h>
#include <gbs/knotsfunctions.h>
#include <gbs-occt/curvesbuild.h>
#include <gbs-occt/export.h>
#include <GeomTools.hxx>
const double tol = 1e-10;
const double tol_confusion = 1e-7;

using gbs::operator-;

TEST(tests_bscurve, ctor)
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
    // double u = 2.3;

    int n = 1000;
    for (int i = 0; i < n; i++)
    {
        auto u = k.front() + (k.back() - k.front()) * i / (n - 1.);
        auto c1_3d_dp = gbs::BSCurve3d_d(poles, k, p);
        auto c1_3d_dp_pt1 = occt_utils::point(c1_3d_dp.value(u));

        auto h_c1_3d_dp_ref = occt_utils::BSplineCurve(c1_3d_dp);
        auto c1_3d_dp_pt1_ref = h_c1_3d_dp_ref->Value(u);

        ASSERT_LT(c1_3d_dp_pt1_ref.Distance(c1_3d_dp_pt1), tol);

        auto c1_3d_dp_tg1 = occt_utils::vector(c1_3d_dp.value(u, 1));
        auto c1_3d_dp_tg1_ref = h_c1_3d_dp_ref->DN(u, 1);

        ASSERT_LT((c1_3d_dp_tg1_ref - c1_3d_dp_tg1).Magnitude(), tol);
    }
}

TEST(tests_bscurve, ctor_rational)
{
    std::vector<double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<double, 3>> poles_nr =
        {
            {0., 0., 0.},
            {0., 1., 0.},
            {1., 1., 0.},
            {1., 1., 1.},
            {1., 1., 2.},
            {3., 1., 1.},
            {0., 4., 1.},
        };
    std::vector<std::array<double, 4>> poles =
        {
            {0., 0., 0., 1.},
            {0., 1., 0., 2.},
            {1., 1., 0., 3.},
            {1., 1., 1., 4.},
            {1., 1., 2., 5.},
            {3., 1., 1., 6.},
            {0., 4., 1., 7.},
        };
    size_t p = 2;
    int n = 1000;
    auto c1_3d_dp = gbs::BSCurveRational3d_d(poles, k, p);
    auto c2_3d_dp = gbs::BSCurveRational3d_d(poles_nr, k, p);
    auto c3_3d_dp = gbs::BSCurve3d_d(poles_nr, k, p);
    auto h_c1_3d_dp_ref = occt_utils::NURBSplineCurve(c1_3d_dp);
    auto h_c2_3d_dp_ref = occt_utils::BSplineCurve(c3_3d_dp);
    for (int i = 0; i < n; i++)
    {
        auto u = k.front() + (k.back() - k.front()) * i / (n - 1.);
        for(auto d =0 ; d <=2 ; d++)
        {
            auto c1_3d_dp_cu1 = occt_utils::vector(c1_3d_dp.value(u, d));
            auto c1_3d_dp_cu1_ref = h_c1_3d_dp_ref->DN(u, d);
            ASSERT_LT((c1_3d_dp_cu1_ref - c1_3d_dp_cu1).Magnitude(), tol);
            auto c2_3d_dp_cu1 = occt_utils::vector(c2_3d_dp.value(u, d));
            auto c2_3d_dp_cu1_ref = h_c2_3d_dp_ref->DN(u, d);
            ASSERT_LT((c2_3d_dp_cu1_ref - c2_3d_dp_cu1).Magnitude(), tol);
        }
    }
}

TEST(tests_bscurve, perf_occt)
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
    auto h_c1_3d_dp_ref = occt_utils::BSplineCurve(c1_3d_dp);

    const auto count_max = 10000000;

	auto count = count_max;
	const auto t1_ref = std::chrono::high_resolution_clock::now();
	while (count)
	{
		u = double(count) / double(count_max) * 5.;
		h_c1_3d_dp_ref->Value(u);
		count--;
	}
	const auto t2_ref = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> ms_ref = t2_ref - t1_ref;
	std::cout << std::fixed
		<< "occt curve took " << ms_ref.count() << " ms\n";

    // const auto t1 = std::chrono::high_resolution_clock::now();
    // count = count_max;
    // while(count)
    // {
    //     u = double(count) / double(count_max) * 5.;
    //     c1_3d_dp.value(u);
    //     count--;
    // }
    // const auto t2 = std::chrono::high_resolution_clock::now();
    // const std::chrono::duration<double, std::milli> ms = t2 - t1;
    // std::cout << std::fixed 
    //               << "gbs took " << ms.count() << " ms\n";

}

TEST(tests_bscurve, perf_long_double)
{
    std::vector<long double> k = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<std::array<long double,3> > poles =
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

    auto c1_3d_dp = gbs::BSCurve<long double,3>(poles,k,p);

    const auto count_max = 10000000;

    const auto t1 = std::chrono::high_resolution_clock::now();
    auto count = count_max;
    while(count)
    {
        u = long double(count) / long double(count_max) * 5.;
        c1_3d_dp.value(u);
        count--;
    }
    const auto t2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<long double, std::milli> ms = t2 - t1;
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
    auto h_c1_3d_dp_ref = occt_utils::BSplineCurve(c1_3d_dp);

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
    auto k1 = gbs::curve_parametrization(pt, gbs::KnotsCalcMode::EQUALY_SPACED);
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