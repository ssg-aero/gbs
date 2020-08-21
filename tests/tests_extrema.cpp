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