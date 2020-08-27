#include <gtest/gtest.h>
#include <gbslib/bscurve.h>
#include <gbslib/bscbuild.h>
#include <gbslib/knotsfunctions.h>
#include <occt-utils/curvesbuild.h>
#include <occt-utils/export.h>
#include <GeomTools.hxx>
const double tol = 1e-10;

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
    double u = 2.3;

    auto c1_3d_dp = gbs::BSCurve(poles,k,p);
    auto c1_3d_dp_pt1 = occt_utils::point( c1_3d_dp.value(u) );

    auto h_c1_3d_dp_ref = occt_utils::BSplineCurve(c1_3d_dp);
    auto c1_3d_dp_pt1_ref = h_c1_3d_dp_ref->Value(u) ;

    ASSERT_LT(c1_3d_dp_pt1_ref.Distance(c1_3d_dp_pt1),tol);

    auto c1_3d_dp_tg1 = occt_utils::vector( c1_3d_dp.value(u,1) );
    auto c1_3d_dp_tg1_ref = h_c1_3d_dp_ref->DN(u,1) ;

    ASSERT_LT((c1_3d_dp_tg1_ref-c1_3d_dp_tg1).Magnitude(),tol);


}

TEST(tests_bscurve, perf)
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

    auto c1_3d_dp = gbs::BSCurve(poles,k,p);
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
		<< " took " << ms_ref.count() << " ms\n";

    const auto t1 = std::chrono::high_resolution_clock::now();
    count = count_max;
    while(count)
    {
        u = double(count) / double(count_max) * 5.;
        c1_3d_dp.value(u);
        count--;
    }
    const auto t2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> ms = t2 - t1;
    std::cout << std::fixed 
                  << " took " << ms.count() << " ms\n";

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

TEST(tests_bscurve, build_circle)
{

    auto r = 1.23;
    std::array<double,3> center{1.,2.,3.};
    auto c1_3d_dp = gbs::build_circle<double,3>(r,center);

    auto h_c1_3d_dp_ref = occt_utils::NURBSplineCurve(c1_3d_dp);
    
    
    c1_3d_dp.insertKnot(3./8.,2);
    h_c1_3d_dp_ref->InsertKnot(3./8.,2);

    
    double u;
    for(int i=0 ; i < 100 ;i++)
    {
        u = i / 99.;
        auto pt = c1_3d_dp.valueRationnal(u);
        ASSERT_NEAR(gbs::norm(pt-center),r,tol);
        ASSERT_LT( (occt_utils::point(pt).XYZ()-h_c1_3d_dp_ref->Value(u).XYZ()).Modulus(),tol);
        // std::cout << u * 2 * M_PI << " ; " << atan2(pt[1],pt[0]) << " ; " <<gbs::norm(pt) << std::endl;
    }

    std::vector<Handle_Geom_Curve> crv_lst;
    crv_lst.push_back(h_c1_3d_dp_ref);   

    occt_utils::to_iges(crv_lst, "C:/Users/sebastien/workspace2/gbslib/tests/out/build_circle.igs");

}

TEST(tests_bscurve, build_elipse)
{

    auto r1 = 1.23;
    auto r2 = 2.34;

    auto c1_3d_dp = gbs::build_ellipse<double,3>(r1,r2);
    auto h_c1_3d_dp_ref = occt_utils::NURBSplineCurve(c1_3d_dp);
    
    double u;
    for(int i=0 ; i < 100 ;i++)
    {
        u = i / 99.;
        auto pt = c1_3d_dp.valueRationnal(u);
        ASSERT_NEAR(pt[0]*pt[0]/(r1*r1)+pt[1]*pt[1]/(r2*r2),1,tol);
        ASSERT_LT( (occt_utils::point(pt).XYZ()-h_c1_3d_dp_ref->Value(u).XYZ()).Modulus(),tol);
        // std::cout << u * 2 * M_PI << " ; " << atan2(pt[1],pt[0]) << " ; " <<gbs::norm(pt) << std::endl;
    }

    std::vector<Handle_Geom_Curve> crv_lst;
    crv_lst.push_back(h_c1_3d_dp_ref);   

    occt_utils::to_iges(crv_lst, "C:/Users/sebastien/workspace2/gbslib/tests/out/build_ellipse.igs");

}