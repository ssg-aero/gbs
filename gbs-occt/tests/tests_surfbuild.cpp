#include <doctest_gtest.hpp>
#include <gbs-occt/surfacesbuild.h>
#ifdef GBS_USE_MODULES
    import knots_functions;
#endif
#include <gbs/bssurf.h>
#include <gbs/bssinterp.h>
#include <vector>
#include <GeomTools.hxx>
#include <Geom_OffsetSurface.hxx>
#include <gbs/offsets.h>
#include <cmath>
#include <memory>



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
TEST(tests_surfacebuild, offset_surface_derivatives)
{
    // gbs::OffsetSurface against Geom_OffsetSurface on a bicubic B-spline: the point and every
    // derivative up to the third
    std::vector<std::array<double, 3>> poles;
    for (int j = 0; j < 6; ++j)
        for (int i = 0; i < 6; ++i)
            poles.push_back({double(i), double(j), 0.3 * std::sin(1.3 * i) * std::cos(0.9 * j) + 0.05 * i * j});
    const std::vector<double> k{0., 0., 0., 0., 1., 2., 3., 3., 3., 3.};
    auto s = std::make_shared<gbs::BSSurface<double, 3>>(poles, k, k, 3, 3);
    for (double d : {0.2, -0.15})
    {
        const gbs::OffsetSurface<double> o{s, d};
        Handle(Geom_OffsetSurface) g = new Geom_OffsetSurface(occt_utils::BSplineSurface(*s), d);
        for (double u : {0.3, 1.45, 2.8})
            for (double v : {0.2, 1.7, 2.55})
            {
                gp_Pnt P;
                gp_Vec D1U, D1V, D2U, D2V, D2UV, D3U, D3V, D3UUV, D3UVV;
                g->D3(u, v, P, D1U, D1V, D2U, D2V, D2UV, D3U, D3V, D3UUV, D3UVV);
                auto near = [&](std::size_t a, std::size_t b, const gp_XYZ &x) {
                    const auto y = o.value(u, v, a, b);
                    const double e = std::hypot(y[0] - x.X(), y[1] - x.Y(), y[2] - x.Z());
                    EXPECT_LT(e, 1e-10 * (1. + x.Modulus())) << a << " " << b;
                };
                near(0, 0, P.XYZ());
                near(1, 0, D1U.XYZ());
                near(0, 1, D1V.XYZ());
                near(2, 0, D2U.XYZ());
                near(0, 2, D2V.XYZ());
                near(1, 1, D2UV.XYZ());
                near(3, 0, D3U.XYZ());
                near(0, 3, D3V.XYZ());
                near(2, 1, D3UUV.XYZ());
                near(1, 2, D3UVV.XYZ());
            }
    }
}
