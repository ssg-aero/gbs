#include <gtest/gtest.h>
#include <gbslib/bscurve.h>
#include <gbslib/vecop.h>
#include <gbslib/bssinterp.h>
#include <gbslib/extrema.h>

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