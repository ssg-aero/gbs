#include <gtest/gtest.h>
#include <gbs/bscbuild.h>
#include <gbs-mesh/mshedge.h>
// #include <gbs-mesh/tfi.h>

TEST(tests_mesh, msh_ed)
{
    auto r_cir = 1.2f;
    auto cir = gbs::build_circle<float,2>(r_cir);
    auto p_cir = new gbs::BSCurveRational2d_f{cir};
    auto np = 11;
    gbs::msh_edge<float,2> ed1{p_cir};
    ed1.set_points(np);
    ed1.compute_pnts();
    ASSERT_FLOAT_EQ( r_cir,ed1.points()[0][0]);
    ASSERT_FLOAT_EQ(-r_cir,ed1.points()[5][0]);
    ASSERT_FLOAT_EQ( r_cir,ed1.points()[10][0]);
}
