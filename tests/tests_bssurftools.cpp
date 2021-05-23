#include <gtest/gtest.h>
#include <gbs/surfaces>
#include <gbs/bsctools.h>
#include <gbs/bsstools.h>
#include <gbs-io/print.h>
#include <gbs-render/vtkcurvesrender.h>

using gbs::operator-;
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

    gbs::BSCurve crv{
        poles_crv,kv,q
    };

    auto natural_end = true;
    auto max_cont = 2;
    auto srf_ext = gbs::extention_to_curve(srf,crv,natural_end,max_cont);
    std::array<double, 3> peacock{51. / 255., 161. / 255., 201. / 255.};
    std::array<double, 3> rosy_brown{188. / 255., 143. / 255., 143. / 255.};

    gbs::plot(
        crv,
        gbs::make_actor( srf     , peacock),
        gbs::make_actor( srf_ext , rosy_brown)
    );

}