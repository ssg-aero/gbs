#include <gtest/gtest.h>
#include <gbslib/bsctools.h>

#include <occt-utils/curvesbuild.h>
#include <occt-utils/export.h>

TEST(tests_bsctools, trim)
{
    std::vector<double> k = {0., 0., 0., 1./4., 1./4., 1./2., 1./2., 3./4., 3./4., 1., 1., 1.};
    auto wi = sqrt(2.)/2.;
    auto r = 1.23;
    std::vector<std::array<double,4> > poles =
    {
        { r, 0.,0.,1.},
        { wi*r, wi*r,0.,wi},
        { 0., r,0.,1.},
        {-wi*r, wi*r,0.,wi},
        {-r, 0.,0.,1.},
        {-wi*r,-wi*r,0.,wi},
        { 0.,-r,0.,1.},
        { wi*r,-wi*r,0.,wi},
        { r, 0.,0.,1.},
    };

    size_t p = 2;
    
    // gbs::trim(p,k,poles,0.1,3./8.);
    gbs::trim(p,k,poles,0.,3./8.);
    gbs::trim(p,k,poles,0.1,3./8.);

    auto c1_3d_dp = gbs::BSCurveRational<double,3>(poles,k,p);

    std::vector<Handle_Geom_Curve> crv_lst;
    crv_lst.push_back(occt_utils::NURBSplineCurve(c1_3d_dp));   

    occt_utils::to_iges(crv_lst, "C:/Users/sebastien/workspace2/gbslib/tests/out/trim.igs");
}