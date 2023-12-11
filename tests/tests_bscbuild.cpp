#include <gtest/gtest.h>
#include <gbs/bscurve.h>
#include <gbs/bscbuild.h>
#include <gbs/bscapprox.h>
#include <gbs-render/vtkGbsRender.h>
import math;
const double tol = 1e-10;

using gbs::operator-;

TEST(tests_bscbuild, build_circle)
{

    auto r = 1.23;
    std::array<double,3> center{1.,2.,3.};
    auto c1_3d_dp = gbs::build_circle<double,3>(r,center);   
    
    c1_3d_dp.insertKnot(3./8.,2);
    
    double u;
    for(int i=0 ; i < 100 ;i++)
    {
        u = i / 99.;
        auto pt = c1_3d_dp.value(u);
        ASSERT_NEAR(gbs::norm(pt-center),r,tol);
    }

}

TEST(tests_bscbuild, build_elipse)
{

    auto r1 = 1.23;
    auto r2 = 2.34;

    auto c1_3d_dp = gbs::build_ellipse<double,3>(r1,r2);
    
    double u;
    for(int i=0 ; i < 100 ;i++)
    {
        u = i / 99.;
        auto pt = c1_3d_dp.value(u);
        ASSERT_NEAR(pt[0]*pt[0]/(r1*r1)+pt[1]*pt[1]/(r2*r2),1,tol);
    }


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
    
    auto C  = gbs::BSCurve<double,3>(poles,k,p);
    auto C1 = gbs::build_derivate(C);
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
    
    auto C  = gbs::BSCurve<double,3>(poles,k,p);
    auto C1 = gbs::build_derivate(C);

    auto Ci = gbs::build_integrate(C1,C.poles().front());
    int n = 1000;
    for( int i = 0 ; i < n ; i++)
    {
        auto u = k.front() + (k.back()-k.front()) * i / (n-1.);
        ASSERT_LT(fabs(gbs::norm(C.value(u,1)-C1.value(u))),tol);
        ASSERT_LT(fabs(gbs::norm(C.value(u)-Ci.value(u))),tol);
    }
}

TEST(tests_bscbuild, build_lenght)
{
    auto r = 1.;
    std::array<double,3> center{0.,0.,0.};
    auto c = gbs::build_circle<double,3>(r,center);
    {
        auto perimeter = gbs::length<double, 3, 10>(c, c.knotsFlats().front(), c.knotsFlats().back());
        ASSERT_NEAR(perimeter, 2 * std:numbers::pi, 1e-6);
    }

    auto np= 360;
    gbs::points_vector_3d_d pts(np);
    for(auto i = 0 ; i < np ; i++)
    {
        pts[i] = {center[0]+r*cos(2*std:numbers::pi/(np-1.)*i), center[1]+r*sin(2*std:numbers::pi/(np-1.)*i), center[2]};
    }

    auto crv = gbs::approx(pts,5,gbs::KnotsCalcMode::CHORD_LENGTH,true);
    {
        auto perimeter = gbs::length<double, 3, 10>(crv, crv.knotsFlats().front(), crv.knotsFlats().back());
        ASSERT_NEAR(perimeter, 2 * std:numbers::pi, 1e-5);
    }
    // gbs::plot(crv,c);
}