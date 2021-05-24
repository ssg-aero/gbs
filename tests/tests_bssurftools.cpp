#include <gtest/gtest.h>
#include <gbs/surfaces>
#include <gbs/bsctools.h>
#include <gbs/bsstools.h>
#include <gbs/bscapprox.h>
#include <gbs-io/print.h>
#include <gbs-render/vtkcurvesrender.h>

using gbs::operator-;

TEST(tests_bssurf, unify)
{

}
TEST(tests_bssurf, extention)
{

    const gbs::points_vector<double,3> points_srf =
    {
        {0,0,0},{3,2,0},{5,1,0},{7,0,0},
        {0,1,1},{3,3,1},{5,2,1},{7,1,1},
        {0,0,2},{2,0,2},{5,0,2},{7,-0.5,2},
    };
    // std::vector<double> ku = {0.,0.,0.,0.5,1.,1.,1.};
    std::vector<double> ku = {0.,0.,0.,0.,1.,1.,1.,1.};
    std::vector<double> kv = {0.,0.,0.,1.,1.,1.};
    size_t p = 3;
    size_t q = 2;
    std::vector<double> u = {0.,0.33,0.66,1.};
    std::vector<double> v = {0.,0.5,1.};
    auto poles_srf = gbs::build_poles(points_srf,ku,kv,u,v,p,q);

    gbs::BSSurface srf(poles_srf,ku,kv,p,q) ;

    const gbs::points_vector<double,3> points_crv =
    {
        {9,0,0},{9,1,1},{9,0,2},
    };
    auto poles_crv = gbs::build_poles(points_crv,kv,v,q);

    gbs::BSCurve crv1{
        poles_crv,kv,q
    };

    using namespace gbs;

    auto natural_end = true;
    auto max_cont = 2;
    auto srf_extention_to_curve = gbs::extention_to_curve(srf,crv1,natural_end,max_cont);
    auto srf_extended_to_curve = gbs::extended_to_curve(srf,crv1,natural_end,max_cont);
    auto l_ext = 0.5;
    auto srf_extended_u = gbs::extention(srf,l_ext,natural_end,max_cont);

    auto srf_extended_u_= gbs::extention(srf,l_ext,gbs::SurfaceBound::U_START,natural_end,max_cont);
    auto srf_extended_v_= gbs::extention(srf,l_ext,gbs::SurfaceBound::V_START,natural_end,max_cont);
    auto srf_extended_v = gbs::extention(srf,l_ext,gbs::SurfaceBound::V_END,natural_end,max_cont);



    std::array<double, 3> peacock{51. / 255., 161. / 255., 201. / 255.};
    std::array<double, 3> rosy_brown{188. / 255., 143. / 255., 143. / 255.};
    std::array<double, 3> DarkSeaGreen{143. / 255.,  188. / 255.,  143. / 255.};
    gbs::plot(
        // crv1,
        gbs::make_actor( srf     , peacock),
        // gbs::make_actor( srf_extention_to_curve , rosy_brown),
        // gbs::make_actor( srf_extended_to_curve , DarkSeaGreen),
        gbs::make_actor( srf_extended_u , rosy_brown),
        gbs::make_actor( srf_extended_v , DarkSeaGreen),
        gbs::make_actor( srf_extended_u_, rosy_brown),
        gbs::make_actor( srf_extended_v_, DarkSeaGreen)
    );

}

