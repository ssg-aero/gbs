#include <gtest/gtest.h>
#include <gbs/surfaces>
#include <gbs/bsctools.h>
#include <gbs/bsstools.h>
#include <gbs/bscapprox.h>
#include <gbs-io/print.h>
#include <gbs-render/vtkGbsRender.h>
#include <gbs/loft/loftBase.h>

using gbs::operator-;

#ifdef TEST_PLOT_ON
    const bool PLOT_ON = true;
#else
    const bool PLOT_ON = false;
#endif

TEST(tests_bssurf, unify)
{
    using namespace gbs;
    using T = double;
    constexpr size_t dim = 3;
    size_t p = 3;
    size_t q = 2;

    std::vector<T> u{0., 0.33, 0.66, 1.};
    std::vector<T> v{0., 0.5, 1.};
    auto ku = build_simple_mult_flat_knots<T>(u, p);
    auto kv = build_simple_mult_flat_knots<T>(v, q);
    points_vector<T, dim> Q{
        {0., 0., 0.}, {0.33, 0., 0.}, {0.66, 0., 0.}, {1., 0., 0.}, 
        {0., 0.5, 0.}, {0.33, 0.5, 0.1}, {0.66, 0.5, 0.2}, {1., 0.5, 0.},
        {0., 1., 0.}, {0.33, 1., 0.}, {0.66, 1., 0.}, {1., 1., 0.},
    };
    auto poles_T = build_poles(Q, ku, kv, u, v, p, q);
    BSSurface<T, dim> srf1{
        poles_T,
        ku, kv,
        p, q};

    auto srf2{srf1};
    auto srf3{srf1};

    srf2.insertKnotU(0.33, 2);
    srf2.insertKnotV(0.33, 1);

    srf3.insertKnotU(0.66, 1);
    srf3.insertKnotV(0.30, 1); 

    std::vector<BSSurfaceInfo<T, dim >> surfaces_info{srf1.info(), srf2.info(), srf3.info()};

    unify_knots(surfaces_info);

    // Check if surfaces are now compatibles
    const auto & [polesUV0, ku0, kv0, p0, d0] = surfaces_info[0];
    std::for_each(
        std::next(surfaces_info.begin()), surfaces_info.end(),
        [&kv0, &ku0](const auto& surface_info){
            const auto & [polesUV_, ku_, kv_, p_, d_] = surface_info;
            ASSERT_TRUE( std::ranges::equal( ku0, ku_) );
            ASSERT_TRUE( std::ranges::equal( kv0, kv_) );
        }
    );
    // Check on construction points
    BSSurface<T, dim> srf1_{surfaces_info[0]};
    BSSurface<T, dim> srf2_{surfaces_info[1]};
    BSSurface<T, dim> srf3_{surfaces_info[2]};
    for (auto u_ : u)
    {
        for (auto v_ : v)
        {
            ASSERT_LE(distance(srf1(u_, v_), srf1_(u_, v_)), 1e-6);
        }
    }

    if(PLOT_ON)
    {
        gbs::plot(srf1_, srf2_, srf3_);
    }
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

    gbs::BSSurface<double,3> srf(poles_srf,ku,kv,p,q) ;

    const gbs::points_vector<double,3> points_crv =
    {
        {9,0,0},{9,1,1},{9,0,2},
    };
    auto poles_crv = gbs::build_poles(points_crv,kv,v,q);

    gbs::BSCurve<double,3> crv1{
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
    if(PLOT_ON) gbs::plot(
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

