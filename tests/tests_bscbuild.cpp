#include <gtest/gtest.h>
#include <gbslib/bscurve.h>
#include <gbslib/bscbuild.h>
#include <occt-utils/curvesbuild.h>
#include <occt-utils/export.h>

const double tol = 1e-10;

using gbs::operator-;

TEST(tests_bscbuild, build_circle)
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

TEST(tests_bscbuild, build_elipse)
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

TEST(tests_bscbuild, build_derivate)
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
    
    auto C  = gbs::BSCurve(poles,k,p);
    auto C1 = gbs::derivate(C);
    int n = 1000;
    for( int i = 0 ; i < n ; i++)
    {
        auto u = k.front() + (k.back()-k.front()) * i / (n-1.);
        ASSERT_LT(fabs(gbs::norm(C.value(u,1)-C1.value(u))),tol);
    }
}

TEST(tests_bscbuild, build_integrate)
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
    
    auto C  = gbs::BSCurve(poles,k,p);
    auto C1 = gbs::derivate(C);
    auto Ci = gbs::integrate(C1,C.poles().front());
    int n = 1000;
    for( int i = 0 ; i < n ; i++)
    {
        auto u = k.front() + (k.back()-k.front()) * i / (n-1.);
        ASSERT_LT(fabs(gbs::norm(C.value(u)-Ci.value(u))),tol);
    }
}