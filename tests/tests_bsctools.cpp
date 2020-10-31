#include <gtest/gtest.h>
#include <gbs/bsctools.h>

#include <gbs-occt/curvesbuild.h>
#include <gbs-occt/export.h>

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

    occt_utils::to_iges(crv_lst, "trim.igs");
}

TEST(tests_bsctools, unify_knots)
{

    std::vector<double> k1 = {0., 0., 0., 1, 2, 3, 4, 5., 5., 5.};
    std::vector<double> k2 = {0., 0., 0.1, 1.3, 1.8 , 2.1, 2.6,  3., 3.};
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
    size_t p1 = 2;
    size_t p2 = 1;
    
    gbs::BSCurve3d_d c1(poles,k1,p1);
    gbs::BSCurve3d_d c2(poles,k2,p2);

    std::list<gbs::BSCurve3d_d> bsc_lst = {c1,c2};

    ASSERT_THROW( gbs::unify_knots(bsc_lst) , std::exception);
    
    gbs::unify_degree(bsc_lst);

    gbs::unify_knots(bsc_lst);

    auto k1_ = bsc_lst.front().knotsFlats();
    auto k2_ = bsc_lst.back().knotsFlats();

    ASSERT_TRUE(
        std::equal(
        k1_.begin(),
        k1_.end(),
        k2_.begin(),
        [](const auto &k1,const auto &k2)
        {
            return fabs( k1 - k2 ) < gbs::knot_eps;
        })
    );

    std::vector<Handle(Geom_Geometry)> crv_lst;
    crv_lst.push_back(occt_utils::BSplineCurve(c1));
    crv_lst.push_back(occt_utils::BSplineCurve(c2));
    crv_lst.push_back(occt_utils::BSplineCurve(bsc_lst.front()));
    crv_lst.push_back(occt_utils::BSplineCurve(bsc_lst.back()));

    occt_utils::to_iges(crv_lst, "unify_knots.igs");
}