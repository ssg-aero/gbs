#include <gtest/gtest.h>
#include <gbs/bssurf.h>

#include <array>
#include <cmath>
#include <vector>

namespace
{
    // Paraboloid S(u, v) = (u, v, u^2 + v^2) as a bidegree-(2,2) Bezier
    // patch on [0, 1]^2. Pole grid is U-first, V-outer to match the
    // library's storage convention.
    //   S_u  = (1, 0, 2u),        S_v  = (0, 1, 2v)
    //   S_uu = (0, 0, 2),         S_vv = (0, 0, 2)
    gbs::BSSurface<double, 3> make_paraboloid()
    {
        const std::vector<double> ku{0., 0., 0., 1., 1., 1.};
        const std::vector<double> kv{0., 0., 0., 1., 1., 1.};
        const gbs::points_vector<double, 3> poles{
            {0., 0., 0.}, {0.5, 0., 0.}, {1., 0., 1.},
            {0., 0.5, 0.}, {0.5, 0.5, 0.}, {1., 0.5, 1.},
            {0., 1., 1.}, {0.5, 1., 1.}, {1., 1., 2.},
        };
        return gbs::BSSurface<double, 3>(poles, ku, kv, 2, 2);
    }
}

TEST(tests_surface_derivatives, d_dmu_is_unit_tangent_along_u)
{
    const auto srf = make_paraboloid();
    constexpr double tol = 1e-12;

    for (double u : {0.0, 0.25, 0.5, 0.75, 1.0})
    for (double v : {0.0, 0.25, 0.5, 0.75, 1.0})
    {
        const auto t = srf.d_dmu(u, v);

        // Expected: (1, 0, 2u) / sqrt(1 + 4u^2)
        const double inv_n = 1.0 / std::sqrt(1.0 + 4.0 * u * u);
        EXPECT_NEAR(t[0],        inv_n, tol) << "u=" << u << " v=" << v;
        EXPECT_NEAR(t[1],           0., tol) << "u=" << u << " v=" << v;
        EXPECT_NEAR(t[2], 2.0 * u * inv_n, tol) << "u=" << u << " v=" << v;
        EXPECT_NEAR(std::sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]), 1.0, tol);
    }
}

TEST(tests_surface_derivatives, d_dmv_is_unit_tangent_along_v)
{
    const auto srf = make_paraboloid();
    constexpr double tol = 1e-12;

    for (double u : {0.0, 0.25, 0.5, 0.75, 1.0})
    for (double v : {0.0, 0.25, 0.5, 0.75, 1.0})
    {
        const auto t = srf.d_dmv(u, v);

        // Expected: (0, 1, 2v) / sqrt(1 + 4v^2)
        const double inv_n = 1.0 / std::sqrt(1.0 + 4.0 * v * v);
        EXPECT_NEAR(t[0],           0., tol) << "u=" << u << " v=" << v;
        EXPECT_NEAR(t[1],        inv_n, tol) << "u=" << u << " v=" << v;
        EXPECT_NEAR(t[2], 2.0 * v * inv_n, tol) << "u=" << u << " v=" << v;
        EXPECT_NEAR(std::sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]), 1.0, tol);
    }
}

TEST(tests_surface_derivatives, d_dmu2_matches_closed_form_along_u)
{
    const auto srf = make_paraboloid();
    constexpr double tol = 1e-12;

    // For S(u, v) = (u, v, u^2 + v^2), along the u-iso (v fixed):
    //   S_u = (1, 0, 2u), |S_u|^2 = 1 + 4u^2, S_uu = (0, 0, 2)
    //   S_u . S_uu = 4u
    //   d2S/dmu2 = (S_uu |S_u|^2 - S_u (S_u.S_uu)) / |S_u|^4
    //            = (-4u, 0, 2) / (1 + 4u^2)^2
    for (double u : {0.0, 0.25, 0.5, 0.75, 1.0})
    for (double v : {0.0, 0.25, 0.5, 0.75, 1.0})
    {
        const auto a = srf.d_dmu2(u, v);

        const double N  = 1.0 + 4.0 * u * u;
        const double N2 = N * N;
        EXPECT_NEAR(a[0], -4.0 * u / N2, tol) << "u=" << u << " v=" << v;
        EXPECT_NEAR(a[1],            0., tol) << "u=" << u << " v=" << v;
        EXPECT_NEAR(a[2],   2.0      / N2, tol) << "u=" << u << " v=" << v;
    }
}

TEST(tests_surface_derivatives, d_dmv2_matches_closed_form_along_v)
{
    const auto srf = make_paraboloid();
    constexpr double tol = 1e-12;

    // Symmetric: d2S/dmv2 = (0, -4v, 2) / (1 + 4v^2)^2
    for (double u : {0.0, 0.25, 0.5, 0.75, 1.0})
    for (double v : {0.0, 0.25, 0.5, 0.75, 1.0})
    {
        const auto a = srf.d_dmv2(u, v);

        const double N  = 1.0 + 4.0 * v * v;
        const double N2 = N * N;
        EXPECT_NEAR(a[0],            0., tol) << "u=" << u << " v=" << v;
        EXPECT_NEAR(a[1], -4.0 * v / N2, tol) << "u=" << u << " v=" << v;
        EXPECT_NEAR(a[2],   2.0      / N2, tol) << "u=" << u << " v=" << v;
    }
}

TEST(tests_surface_derivatives, vectorized_arc_length_derivatives_match_scalar)
{
    const auto srf = make_paraboloid();
    constexpr double tol = 1e-14;

    std::vector<std::array<double, 2>> uv_lst;
    for (double u : {0.0, 0.25, 0.5, 0.75, 1.0})
    for (double v : {0.0, 0.25, 0.5, 0.75, 1.0})
        uv_lst.push_back({u, v});

    // The vectorized arc-length derivatives must match the scalar overloads
    // entry by entry, in order.
    const auto du   = srf.d_dmus(uv_lst);
    const auto dv   = srf.d_dmvs(uv_lst);
    const auto du2  = srf.d_dmu2s(uv_lst);
    const auto dv2  = srf.d_dmv2s(uv_lst);

    ASSERT_EQ(du.size(),  uv_lst.size());
    ASSERT_EQ(dv.size(),  uv_lst.size());
    ASSERT_EQ(du2.size(), uv_lst.size());
    ASSERT_EQ(dv2.size(), uv_lst.size());

    for (std::size_t i{}; i < uv_lst.size(); ++i)
    {
        const auto u = uv_lst[i][0];
        const auto v = uv_lst[i][1];
        const auto su   = srf.d_dmu(u, v);
        const auto sv   = srf.d_dmv(u, v);
        const auto su2  = srf.d_dmu2(u, v);
        const auto sv2  = srf.d_dmv2(u, v);
        for (std::size_t k{}; k < 3; ++k)
        {
            EXPECT_NEAR(du[i][k],  su[k],  tol) << "i=" << i << " k=" << k;
            EXPECT_NEAR(dv[i][k],  sv[k],  tol) << "i=" << i << " k=" << k;
            EXPECT_NEAR(du2[i][k], su2[k], tol) << "i=" << i << " k=" << k;
            EXPECT_NEAR(dv2[i][k], sv2[k], tol) << "i=" << i << " k=" << k;
        }
    }
}
