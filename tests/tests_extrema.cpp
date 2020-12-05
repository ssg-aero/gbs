#include <gtest/gtest.h>
#include <gbs/bscurve.h>
#include <gbs/vecop.h>
#include <gbs/bssinterp.h>
#include <gbs/bscbuild.h>
#include <gbs/extrema.h>

#include <gbs-render/vtkcurvesrender.h>

using gbs::operator-;

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

    auto crv = gbs::interpolate(Q,2,gbs::KnotsCalcMode::CHORD_LENGTH);
    auto u = 0.3;
    auto res = gbs::extrema_PC(crv,crv.value(u),1e-6);
    ASSERT_NEAR(res.u,u,1e-6);
    u = 0.7;
    res = gbs::extrema_PC(crv,crv.value(u),1e-6);
    ASSERT_NEAR(res.u,u,1e-6);
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
        {0,0,0},{1,0,1},{1,2,0},
        {0,1,0},{1,1,1},{1,2,1},
        {0,2,0},{1,2,1},{1,2,2},
        {0,3,0},{1,3,1},{1,2,3},
        {0,4,0},{1,4,1},{1,2,4},
    };

    gbs::BSSurface srf(poles,ku_flat,kv_flat,p,q);
    auto u = 2.5;
    auto v = 1.;

    auto pt = srf.value(u,v);

    auto res = gbs::extrema_PS(srf,pt,1e-6);
    ASSERT_NEAR(res.u,u,1e-6);
    ASSERT_NEAR(res.v,v,1e-6);

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

    res = gbs::extrema_PS(srfNURBS,pt,1e-6);
    ASSERT_NEAR(res.u,u,1e-6);
    ASSERT_NEAR(res.v,v,1e-6);

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
    gbs::BSSurface srf(polesS,ku_flat,kv_flat,p,q);
    gbs::BSCurve crv(polesC,ku_flat,p);

    // gbs::plot(srf,crv);

    auto res = gbs::extrema_CS(srf,crv,1.e-6);
    ASSERT_LT(res.d,1e-6);
    ASSERT_NEAR(res.u_s,0.5,1e-6);
    ASSERT_NEAR(res.v_s,0.5,1e-6);
    ASSERT_NEAR(res.u_c,0.5,1e-6);

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
    gbs::BSCurve crv2(polesC2,ku_flat2,p);

    res = gbs::extrema_CS(srf2,crv2,1.e-6);
    auto I = crv2.value(res.u_c);

    ASSERT_LT(res.d,2e-5);

    // auto colors = vtkSmartPointer<vtkNamedColors>::New();
    // gbs::plot(
    //     srf2,
    //     crv2,
    //     gbs::make_actor(gbs::points_vector<double,3>{I},30.,true,colors->GetColor4d("Blue").GetData())
    //     );

}

TEST(tests_extrema, CC)
{
    auto c1 = gbs::build_circle<double,3>(1.);
    auto c2 = gbs::build_segment<double,3>({0.,0.,0.},{1.,1.,0.});
    auto result = gbs::extrema_CC(c1,c2,1.e-6);

    ASSERT_LT(result.d,5e-6);
    ASSERT_NEAR(result.u2,1.,5e-6);
}