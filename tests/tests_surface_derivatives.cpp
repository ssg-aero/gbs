#include "doctest_gtest.hpp"
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

// The surface eval path span-reduces unconditionally (no d==0 guard), so before
// the find_span half-open fix its DERIVATIVES were silently wrong on interior
// knots of multiplicity > 1 (undetected because the other tests use a Bezier
// patch with no interior knots). This guards span-reduced vs full-span surface
// derivatives at C^0 joints in both u and v. See the curve counterpart
// tests_bscurve.deriv_span_reduction_matches_full_span.
TEST(tests_surface_derivatives, deriv_span_reduction_matches_full_span)
{
    using T = double;
    constexpr double tol = 1e-10;

    // Bidegree (2,2). Interior knots u=2 and v=2 both have multiplicity 2 (==p,
    // ==q) => C^0 joints with discontinuous derivatives there.
    const std::vector<T> ku{0.,0.,0., 1., 2.,2., 3., 4.,4.,4.};
    const std::vector<T> kv{0.,0.,0., 1., 2.,2., 3., 4.,4.,4.};
    const size_t p = 2, q = 2;
    const size_t nu = ku.size() - p - 1;   // 7
    const size_t nv = kv.size() - q - 1;   // 7

    gbs::points_vector<T, 3> poles(nu * nv);   // storage: U-first, V-outer
    for (size_t j = 0; j < nv; ++j)
        for (size_t i = 0; i < nu; ++i)
            poles[j * nu + i] = {T(i), T(j), 0.3 * T(i) + 0.7 * T(j) * T(j) - 0.2 * T(i) * T(j)};

    gbs::BSSurface<T, 3> srf(poles, ku, kv, p, q);

    // Full-span reference: sum over ALL basis functions (no span reduction).
    auto full = [&](T u, T v, size_t du, size_t dv) {
        std::array<T, 3> pt{0., 0., 0.};
        for (size_t j = 0; j < nv; ++j)
        {
            T Nv = gbs::basis_function(v, j, q, dv, kv);
            for (size_t i = 0; i < nu; ++i)
            {
                T Nu = gbs::basis_function(u, i, p, du, ku);
                for (size_t c = 0; c < 3; ++c) pt[c] += Nu * Nv * poles[j * nu + i][c];
            }
        }
        return pt;
    };

    // Sample the grid plus every distinct knot value (hits the multiple knots).
    std::vector<T> us, vs;
    for (size_t i = 0; i <= 8; ++i) { us.push_back(4.0 * T(i) / 8.0); vs.push_back(4.0 * T(i) / 8.0); }
    for (T kvv : ku) us.push_back(kvv);
    for (T kvv : kv) vs.push_back(kvv);

    for (size_t du : {size_t{0}, size_t{1}, size_t{2}})
    for (size_t dv : {size_t{0}, size_t{1}, size_t{2}})
        for (T u : us)
        for (T v : vs)
        {
            const auto got = srf.value(u, v, du, dv);
            const auto ref = full(u, v, du, dv);
            T err2 = 0;
            for (size_t c = 0; c < 3; ++c) { T e = got[c] - ref[c]; err2 += e * e; }
            ASSERT_LT(std::sqrt(err2), tol)
                << "du=" << du << " dv=" << dv << " u=" << u << " v=" << v;
        }
}

// The one-pass multi-order derivatives(u, v, nu, nv, SKL) must agree with
// calling value(u, v, ku, kv) per mixed order, sampled on a grid plus every
// distinct knot, with interior knots of multiplicity p/q (C^0 joints).
TEST(tests_surface_derivatives, derivatives_match_per_order_value)
{
    using T = double;
    constexpr double tol = 1e-9;

    // Bidegree (2,2) with interior knots u=2 (mult 2) and v=2 (mult 2).
    const std::vector<T> ku{0.,0.,0., 1., 2.,2., 3., 4.,4.,4.};
    const std::vector<T> kv{0.,0.,0., 1., 2.,2., 3., 4.,4.,4.};
    const size_t p = 2, q = 2;
    const size_t nu = ku.size() - p - 1;
    const size_t nv = kv.size() - q - 1;

    gbs::points_vector<T, 3> poles(nu * nv);
    for (size_t j = 0; j < nv; ++j)
        for (size_t i = 0; i < nu; ++i)
            poles[j * nu + i] = {T(i), T(j), 0.3 * T(i) + 0.7 * T(j) * T(j) - 0.2 * T(i) * T(j)};

    gbs::BSSurface<T, 3> srf(poles, ku, kv, p, q);

    std::vector<T> us, vs;
    for (size_t i = 0; i <= 6; ++i) { us.push_back(4.0 * T(i) / 6.0); vs.push_back(4.0 * T(i) / 6.0); }
    for (T x : ku) us.push_back(x);
    for (T x : kv) vs.push_back(x);

    const size_t Nu = 3, Nv = 3; // request orders 0..3 in each direction
    std::array<std::array<T, 3>, (Nu + 1) * (Nv + 1)> SKL;
    for (T u : us)
        for (T v : vs)
        {
            srf.derivatives(u, v, Nu, Nv, SKL.data());
            for (size_t ku_ = 0; ku_ <= Nu; ++ku_)
                for (size_t kv_ = 0; kv_ <= Nv; ++kv_)
                {
                    const auto ref = srf.value(u, v, ku_, kv_);
                    const auto &got = SKL[ku_ * (Nv + 1) + kv_];
                    T err2 = 0;
                    for (size_t c = 0; c < 3; ++c) { T e = got[c] - ref[c]; err2 += e * e; }
                    ASSERT_LT(std::sqrt(err2), tol)
                        << "ku=" << ku_ << " kv=" << kv_ << " u=" << u << " v=" << v;
                }
        }
}
