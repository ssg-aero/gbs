#include <gtest/gtest.h>
#include <gbs-occt/curvesbuild.h>
#include <gbs-occt/export.h>
#include <array>
#include <vector>
#include <chrono>
#include <GeomAPI.hxx>

const double tol = 1e-10;
const double tol_confusion = 1e-7;


TEST(tests_curvebuild, bsc_c1)
{
    auto c1 = occt_utils::bscurve_c1<2>(
        {
            {0.,0.},
            {-.1,-.1}
        },
        {0.,5.},
        {1.,0.},
        false,
        1e-6
    );
    auto c2 = occt_utils::bscurve_c1<3>(
        {
            {0.,0.,0.},
            {-.1,-.1,-1.}
        },
        {0.,5.,0.},
        {1.,0.,0.3},
        false,
        1e-6
    );
    // auto c1 = occt_utils::bscurve_c1<4>(//won't compile
    //     {
    //         {0.,0.},
    //         {-.1,-.1}
    //     },
    //     {0.,5.},
    //     {1.,0.},
    //     false,
    //     1e-6
    // );
    // occt_utils::to_iges({c1},"tests/out/c1.igs");
    std::vector<Handle(Geom_Geometry)> geom_lst{GeomAPI::To3d(c1, gp_Pln()),c2};
    occt_utils::to_iges(geom_lst, "tests/out/bsc_c1.igs");
}

TEST(tests_curvebuild, bscurve_c2_approx)
{
    // auto c1 = occt_utils::bscurve_c2_approx<2>(
    //     {
    //         {0., 0.},
    //         {-.1, -.1}
    //     },
    //     {
    //         {0., 5.},
    //         {1., 0.}
    //     },
    //     1e-6);
    auto pt = occt_utils::col_from_vec<gp_Pnt2d>({std::array<double, 2>{0., 0.},
                                                  std::array<double, 2>{0.3, 0.1},
                                                  std::array<double, 2>{1., 1.}});
    auto tg = occt_utils::col_from_vec<gp_Vec2d>({std::array<double, 2>{1., 0.},
                                                  std::array<double, 2>{0., 0.},
                                                  std::array<double, 2>{0., 1.}});
    auto c1 = occt_utils::bscurve_c2_approx(pt, tg, 1e-6);
    occt_utils::to_iges(std::vector<Handle(Geom_Geometry)>{GeomAPI::To3d(c1, gp_Pln())}, "tests/out/bscurve_c2_approx.igs");
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

}

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


TEST(tests_basis_functions, eval_curve)
{
    std::vector<double> k1 = {0.,0.,0.,1.,1.,1.};
    auto u = 0.3;
    auto p = 2;

    auto crv = occt_utils::BSplineCurve(
        {gp_Pnt{0., 0., 0.},
         gp_Pnt{1., 0., 0.},
         gp_Pnt{1., 1., 1.}},
        {0., 1.},
        {3, 3},
        2);
    auto pt1 = crv->Value(u);
    const std::vector<std::array<double,3> > poles = {{0., 0., 0.}, {1., 0., 0.}, {1., 1., 1.}};
    auto pt2 = occt_utils::point(gbs::eval_value_decasteljau(u, k1, poles , 2));
    auto v1 = crv->DN(u,1);
    auto v2 = occt_utils::vector(gbs::eval_value_decasteljau(u, k1, poles , 2,1));

    ASSERT_LT(pt1.Distance(pt2),tol);
    ASSERT_LT((v1-v2).Magnitude(),tol);


}

TEST(tests_basis_functions, eval_curve_perf)
{

    std::vector<double> k1 = {0.,0.,0.,1.,1.,1.};
    const std::vector<std::array<double,3> > poles = {{0., 0., 0.}, {1., 0., 0.}, {1., 1., 1.}};
    auto u = 0.3;
    auto p = 2;

    auto crv = occt_utils::BSplineCurve(
        {gp_Pnt{0., 0., 0.},
         gp_Pnt{1., 0., 0.},
         gp_Pnt{1., 1., 1.}},
        {0., 1.},
        {3, 3},
        2);

	auto count = 100000;
	const auto t1_ref = std::chrono::high_resolution_clock::now();
	while (count)
	{
		u = (rand() % 1000) / 999.;
		crv->Value(u);
		count--;
	}
	const auto t2_ref = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double, std::milli> ms_ref = t2_ref - t1_ref;
	std::cout << std::fixed
		<< " took " << ms_ref.count() << " ms\n";

    const auto t1 = std::chrono::high_resolution_clock::now();
    count = 100000;
    while(count)
    {
        u = (rand() % 1000) / 999.;
        gbs::eval_value_decasteljau(u, k1, poles , 2);
        count--;
    }
    const auto t2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> ms = t2 - t1;
    std::cout << std::fixed 
                  << " took " << ms.count() << " ms\n";
}