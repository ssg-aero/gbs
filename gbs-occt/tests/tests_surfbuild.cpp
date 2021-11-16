#include <gtest/gtest.h>
#include <gbs-occt/surfacesbuild.h>
#include <gbs/knotsfunctions.h>
#include <gbs/bssurf.h>
#include <gbs/bssinterp.h>
#include <vector>
#include <GeomTools.hxx>



TEST(tests_surfacebuild, from_std_container)
{
    std::vector<double> ku_flat = {0.,0.,0.,5.,5.,5.};
    std::vector<double> kv_flat = {0.,0.,0.,1.,2.,3.,3.,3.};
    size_t p = 2;
    size_t q = 2;
    auto u = 2.5;
    auto v = 1.;

    std::vector<double> ku,kv;
    std::vector<int> mu,mv;
    gbs::unflat_knots(ku_flat,mu,ku);
    gbs::unflat_knots(kv_flat,mv,kv);

    //Pij avec j inner loop
    const std::vector<std::array<double,3> > poles =
    {
        {0,0,0},{1,0,1},{1,2,0},
        {0,1,0},{1,1,1},{1,2,1},
        {0,2,0},{1,2,1},{1,2,2},
        {0,3,0},{1,3,1},{1,2,3},
        {0,4,0},{1,4,1},{1,2,4},
    };

    //   =========> U >=========
    //   ||
    //   ||
    //   ||
    //   vv
    //   V
    //   vv
    //   ||
    //   ||

    // auto arr = occt_utils::col2_from_vec<gp_Pnt>(poles,3);
    // for(int i = 1 ; i <= 5; i++)
    // {
    //     for(int j= 1 ; j <= 3 ; j++)
    //     {
    //         arr.Value(i,j).Coord().DumpJson(std::cout);
    //         std::cout << std::endl;
    //     }
    // }

    auto srf = occt_utils::BSplineSurface(poles, ku, kv, mu, mv, p, q);
    // GeomTools::Dump(srf,std::cout);
    gbs::BSSurface<double,3> srf2(poles,ku_flat,kv_flat,p,q);
    auto pt_ref = srf->Value(u, v);
    auto pt     = occt_utils::point( srf2.value(u,v) );

    ASSERT_LT(pt.Distance(pt_ref),1e-6);

    srf->Value(ku.front(), kv.front()).Coord().DumpJson(std::cout);
    std::cout << std::endl;
    srf->Value(ku.back(), kv.front()).Coord().DumpJson(std::cout);
    std::cout << std::endl;
    srf->Value(ku.back(), kv.back()).Coord().DumpJson(std::cout);
    std::cout << std::endl;
    srf->Value(ku.front(), kv.back()).Coord().DumpJson(std::cout);
    std::cout << std::endl;
    std::cout << std::endl;

    occt_utils::point( srf2.value(ku.front(), kv.front())).Coord().DumpJson(std::cout);
    std::cout << std::endl;
    occt_utils::point( srf2.value(ku.back(), kv.front())).Coord().DumpJson(std::cout);
    std::cout << std::endl;
    occt_utils::point( srf2.value(ku.back(), kv.back())).Coord().DumpJson(std::cout);
    std::cout << std::endl;
    occt_utils::point( srf2.value(ku.front(), kv.back())).Coord().DumpJson(std::cout);
    std::cout << std::endl;
}
const double tol = 1e-10;
TEST(tests_bssurf,derivatives)
{
    const std::vector<std::array<double,3> > points =
    {
        {0,0,0},{1,0,0},
        {0,1,0},{1,1,1},
        {0,2,1},{2,1,0},
        {3,2,0},{3,2,0},
    };

    auto srf = gbs::interpolate(points,4,1,2,gbs::KnotsCalcMode::CHORD_LENGTH);

    auto srf_occt = occt_utils::BSplineSurface(srf);

    auto u=0.3,v=0.7;

    ASSERT_LT((srf_occt->DN(u,v,1,0) - occt_utils::vector( srf.value(u,v,1,0) )).Magnitude(),tol);
    ASSERT_LT((srf_occt->DN(u,v,0,1) - occt_utils::vector( srf.value(u,v,0,1) )).Magnitude(),tol);
    ASSERT_LT((srf_occt->DN(u,v,0,2) - occt_utils::vector( srf.value(u,v,0,2) )).Magnitude(),tol);
    ASSERT_LT((srf_occt->DN(u,v,1,1) - occt_utils::vector( srf.value(u,v,1,1) )).Magnitude(),tol);
    ASSERT_LT((srf_occt->DN(u,v,1,2) - occt_utils::vector( srf.value(u,v,1,2) )).Magnitude(),tol);
}