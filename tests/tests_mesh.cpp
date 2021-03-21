#include <gtest/gtest.h>
#include <gbs/bscbuild.h>
#include <gbs-mesh/mshedge.h>
#include <gbs-mesh/tfi.h>

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

TEST(tests_mesh, tfi_blend_functions)
{
    std::vector<float> ksi_i_f = {0., 0.5, 1.};

    const auto n = 2;
    auto a_n = gbs::build_tfi_blend_function_with_derivatives<float, n>(ksi_i_f);
    for (auto d = 0; d < n; d++)
    {
        for (auto d_ = 0; d_ < n; d_++)
        {
            for (auto i = 0; i < ksi_i_f.size(); i++)
            {
                for (auto i_ = 0; i_ < ksi_i_f.size(); i_++)
                {
                    ASSERT_FLOAT_EQ(
                        a_n[i][d].value(ksi_i_f[i_], d_)[0],
                        gbs::kronecker<float>(i, i_) * gbs::kronecker<float>(d, d_));
                }
            }
        }
    }
}