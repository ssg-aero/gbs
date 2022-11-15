#include <gtest/gtest.h>
#include <gbs/bscurve.h>
#include <gbs/vecop.h>
#include <gbs/bssinterp.h>
#include <gbs/bscbuild.h>
#include <gbs/extrema.h>

#include <gbs-render/vtkcurvesrender.h>

using gbs::operator-;



using namespace gbs;



    template <typename... T, std::size_t... I>
    auto subtuple_(const std::tuple<T...>& t, std::index_sequence<I...>) {
    return std::make_tuple(std::get<I>(t)...);
    }

    template <int Trim, typename... T>
    auto subtuple(const std::tuple<T...>& t) {
    return subtuple_(t, std::make_index_sequence<sizeof...(T) - Trim>());
    }

    // template <typename... Targs>
    //     auto bracket(Targs... Fargs)
    //     {
    //         auto tuple = std::tie(Fargs...);
    //         if
    //     }
    


TEST(tests_extrema, PC)
{
    std::vector<gbs::constrType<double, 3, 1>> Q =
        {
            {{0., 0., 0.}},
            {{1., 0., 0.}},
            {{1., 1., 0.}},
            {{1., 1., 2.}},
            {{0., 1., 1.}},
            {{0., -1., 4.}},
        };

    gbs::points_vector<double,3> pts;

    auto crv = gbs::interpolate(Q,2,gbs::KnotsCalcMode::CHORD_LENGTH);

    {
        auto u = 0.3;
        pts.push_back(crv.value(u));
        auto [res_u, res_d] = gbs::extrema_curve_point(crv, pts.back(), 1e-6);
        ASSERT_NEAR(res_d, 0., 1e-6);
        ASSERT_NEAR(res_u, u, 1e-6);
    }
    {
        auto u = 0.7;
        pts.push_back(crv.value(u));
        auto [res_u, res_d] = gbs::extrema_curve_point(crv, pts.back(), 1e-6);
        ASSERT_NEAR(res_d, 0., 1e-6);
        ASSERT_NEAR(res_u, u, 1e-6);
    }

}

TEST(tests_extrema, PC_f)
{
    std::vector<gbs::constrType<float, 3, 1>> Q =
        {
            {{0., 0., 0.}},
            {{1., 0., 0.}},
            {{1., 1., 0.}},
            {{1., 1., 2.}},
            {{0., 1., 1.}},
            {{0., -1., 4.}},
        };

    gbs::points_vector<float,3> pts;

    auto crv = gbs::interpolate(Q,2,gbs::KnotsCalcMode::CHORD_LENGTH);

    {
        auto u = 0.3;
        pts.push_back(crv.value(u));
        auto [res_u, res_d] = gbs::extrema_curve_point<float,3>(crv, pts.back(), 1e-6);
        ASSERT_NEAR(res_d, 0., 1e-5);
        ASSERT_NEAR(res_u, u, 1e-5);
    }
    {
        auto u = 0.7;
        pts.push_back(crv.value(u));
        auto [res_u, res_d] = gbs::extrema_curve_point<float,3>(crv, pts.back(), 1e-6);
        ASSERT_NEAR(res_d, 0., 1e-5);
        ASSERT_NEAR(res_u, u, 1e-5);
    }

}

TEST(tests_extrema, PS)
{
    std::vector<double> ku_flat = {0.,0.,0.,5.,5.,5.};
    std::vector<double> kv_flat = {0.,0.,0.,1.,2.,3.,3.,3.};
    size_t p = 2;
    size_t q = 2;

    //Pij avec j inner loop
    const std::vector<std::array<double,3> > poles =
    {
        {0,0,0},{1,0,1},{2,0,2},
        {0,1,0},{1.2,1,1},{2,1,1},
        {0,2,1},{1.2,2,1},{2,2,2},
        {0,3,1},{1.1,3,1},{2,3,3},
        {0,4,0},{1,4,1},{2,4,0},
    };
    

    gbs::BSSurface<double,3> srf(poles,ku_flat,kv_flat,p,q);
    auto u = 2.5;
    auto v = 1.;

    auto pt = srf.value(u,v);

    // gbs::plot(
    //     srf,
    //     gbs::points_vector<double,3>{srf(res.u,res.v),pt}
    // );


    {
        auto [res_u, res_v,res_d] = gbs::extrema_surf_pnt(srf, pt, 1e-6);
        ASSERT_NEAR(res_d, 0., 1e-6);
        ASSERT_NEAR(res_u, u, 5e-6);
        ASSERT_NEAR(res_v, v, 5e-6);
    }


    std::vector<double> ku = {0.,0.,0.,1.,2.,3.,4.,4.,5.,5.,5.};
    std::vector<double> kv = {0.,0.,0.,1.,2.,3.,3.,3.};

    std::vector<std::array<double,4> > poles_t = {
        {0,2,4,1},{0,2,4,1}, {0,6,4,2},    {0,2,0,1},{0,2,0,1},
        {0,2,4,1},{0,2,4,1}, {0,6,4,2},    {0,2,0,1},{0,2,0,1},
        {0,2,4,1},{0,2,4,1}, {0,6,4,2},    {0,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,6,8,2}, {12,24,12,6}, {4,6,0,2},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1},
        {0,2,4,1},{4,2,4,1}, {8,6,4,2},    {4,2,0,1},{0,2,0,1}
                                                    };

    //Pij avec i inner loop
    std::vector<std::array<double,4> > poles_r(poles_t.size());
    int ni = 5 , nj =8;
    for (int i = 0; i < ni; i++)
    {
        for (int j = 0; j < nj; j++)
        {
            poles_r[j + nj * i] = poles_t[i + ni * j];
        }
    }

    gbs::BSSurfaceRational<double,3> srfNURBS(poles_r,ku,kv,p,q);

    pt = srfNURBS.value(u,v);
    // gbs::plot(
    //     srfNURBS,
    //     gbs::points_vector<double,3>{pt}
    // );

    {
        auto [res_u, res_v,res_d] = gbs::extrema_surf_pnt(srfNURBS, pt, 1e-6, nlopt::LN_COBYLA);
        ASSERT_NEAR(res_d, 0., 1e-6);
        ASSERT_NEAR(res_u, u, 1e-6);
        ASSERT_NEAR(res_v, v, 1e-6);
    }

}

TEST(tests_extrema, PS_f)
{
    std::vector<float> ku_flat = {0.,0.,0.,5.,5.,5.};
    std::vector<float> kv_flat = {0.,0.,0.,1.,2.,3.,3.,3.};
    size_t p = 2;
    size_t q = 2;

    //Pij avec j inner loop
    const std::vector<std::array<float,3> > poles =
    {
        {0,0,0},{1,0,1},{2,0,2},
        {0,1,0},{1.2,1,1},{2,1,1},
        {0,2,1},{1.2,2,1},{2,2,2},
        {0,3,1},{1.1,3,1},{2,3,3},
        {0,4,0},{1,4,1},{2,4,0},
    };
    

    gbs::BSSurface<float,3> srf(poles,ku_flat,kv_flat,p,q);
    auto u = 2.5;
    auto v = 1.;

    auto pt = srf.value(u,v);

    // gbs::plot(
    //     srf,
    //     gbs::points_vector<double,3>{srf(res.u,res.v),pt}
    // );


    {
        auto [res_u, res_v,res_d] = gbs::extrema_surf_pnt<float,3>(srf, pt, 1e-6);
        ASSERT_NEAR(res_d, 0., 1e-5);
        ASSERT_NEAR(res_u, u, 5e-5);
        ASSERT_NEAR(res_v, v, 5e-5);
    }
}

TEST(tests_extrema, CS)
{
    std::vector<double> ku_flat = {0.,0.,1.,1.};
    std::vector<double> kv_flat = {0.,0.,1.,1.};
    size_t p = 1;
    size_t q = 1;
    const std::vector<std::array<double,3> > polesS =
    {
        {0,0,0},{1,0,0},
        {0,1,0},{1,1,0},
    };
    const std::vector<std::array<double,3> > polesC =
    {
        {0.5,0.5,-1. },{0.5,0.5,1.},
    };
    gbs::BSSurface<double,3> srf(polesS,ku_flat,kv_flat,p,q);
    gbs::BSCurve<double,3> crv(polesC,ku_flat,p);

    // gbs::plot(srf,crv);

    {
        auto [res_u, res_v, res_uc,res_d] = gbs::extrema_surf_curve(srf, crv, 1.e-6);
        ASSERT_LT(res_d, 1e-6);
        ASSERT_NEAR(res_u, 0.5, 1e-6);
        ASSERT_NEAR(res_v, 0.5, 1e-6);
        ASSERT_NEAR(res_uc, 0.5, 1e-6);
    }

    std::vector<double> ku_flat2 = {0.,0.,0.,1.,1.,1.};
    std::vector<double> kv_flat2 = {0.,0.,0.,1.,1.,1.};
    p = 2;
    q = 2;
    const std::vector<std::array<double,4> > polesS2 =
    {
        {0,0,0  ,1},{0.5,0,0   ,1},{1,0,0  ,1},
        {0,0.5,0,1},{5*0.75,5*0.5,5*1.,5},{1,0.5,0,1},
        {0,1,0  ,1},{0.5,1.,0   ,1},{1,1,0 ,1},
    };
        const std::vector<std::array<double,3> > polesC2 =
    {
        {0.5,0.5,-.5 },{0.,0.5,0.},{0.5,0.5,0.75},
    };
    gbs::BSSurfaceRational<double,3> srf2(polesS2,ku_flat2,kv_flat2,p,q);
    gbs::BSCurve<double,3> crv2(polesC2,ku_flat2,p);

    {
        auto [res_u, res_v, res_uc,res_d] = gbs::extrema_surf_curve(srf2, crv2, 1.e-6);
        auto I = crv2.value(res_uc);

        ASSERT_LT(res_d, 2e-5);
    }

    // auto colors = vtkSmartPointer<vtkNamedColors>::New();
    // gbs::plot(
    //     srf2,
    //     crv2,
    //     gbs::make_actor(gbs::points_vector<double,3>{I},30.,true,colors->GetColor4d("Blue").GetData())
    //     );

}
TEST(tests_extrema, CS_f)
{
    std::vector<float> ku_flat = {0.,0.,1.,1.};
    std::vector<float> kv_flat = {0.,0.,1.,1.};
    size_t p = 1;
    size_t q = 1;
    const std::vector<std::array<float,3> > polesS =
    {
        {0,0,0},{1,0,0},
        {0,1,0},{1,1,0},
    };
    const std::vector<std::array<float,3> > polesC =
    {
        {0.5,0.5,-1. },{0.5,0.5,1.},
    };
    gbs::BSSurface<float,3> srf(polesS,ku_flat,kv_flat,p,q);
    gbs::BSCurve<float,3> crv(polesC,ku_flat,p);

    // gbs::plot(srf,crv);

    {
        auto [res_u, res_v, res_uc,res_d] = gbs::extrema_surf_curve<float,3>(srf, crv, 1.e-6);
        ASSERT_LT(res_d, 1e-6);
        ASSERT_NEAR(res_u, 0.5, 1e-6);
        ASSERT_NEAR(res_v, 0.5, 1e-6);
        ASSERT_NEAR(res_uc, 0.5, 1e-6);
    }

}
TEST(tests_extrema, CC)
{
    auto c1 = gbs::build_circle<double,3>(1.);
    auto c2 = gbs::build_segment<double,3>({0.,0.,0.},{1.,1.,0.});
    auto [res_u1,res_u2,res_d] = gbs::extrema_curve_curve(c1,c2,1.e-6,nlopt::LN_COBYLA);

    ASSERT_LT(res_d,5e-6);
    ASSERT_NEAR(res_u2,1.,5e-6);
}

TEST(tests_extrema, CC_f)
{
    auto c1 = gbs::build_circle<float,3>(1.);
    auto c2 = gbs::build_segment<float,3>({0.,0.,0.},{1.,1.,0.});
    auto [res_u1,res_u2,res_d] = gbs::extrema_curve_curve<float,3>(c1,c2,1.e-6,nlopt::LN_COBYLA);

    ASSERT_LT(res_d,5e-6);
    ASSERT_NEAR(res_u2,1.,5e-6);
}

TEST(tests_extrema, max_coordinate)
{
    gbs::BSCurve2d_d crv{
        {
            {0.0,0.},
            {0.1,0.},
            {0.4,1.},
            {0.6,1.},
            {0.9,0.},
            {1.0,0.},
        },
        {0., 0.4, 0.6, 1.},
        {4,1,1,4},
        3
    };
    auto [u,r] = gbs::max_coordinate(crv, 1, 0.5);
    ASSERT_NEAR(u,0.5,1e-6);
}