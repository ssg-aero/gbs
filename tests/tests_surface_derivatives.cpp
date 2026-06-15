#include <doctest_gtest.hpp>
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
    const std::vector<double> paraboloid_ku{0., 0., 0., 1., 1., 1.};
    const std::vector<double> paraboloid_kv{0., 0., 0., 1., 1., 1.};
    const gbs::points_vector<double, 3> paraboloid_poles{
        {0., 0., 0.}, {0.5, 0., 0.}, {1., 0., 1.},
        {0., 0.5, 0.}, {0.5, 0.5, 0.}, {1., 0.5, 1.},
        {0., 1., 1.}, {0.5, 1., 1.}, {1., 1., 2.},
    };

    gbs::BSSurface<double, 3> make_paraboloid()
    {
        return gbs::BSSurface<double, 3>(paraboloid_poles, paraboloid_ku, paraboloid_kv, 2, 2);
    }

    // Same control net as the paraboloid, promoted to rational with all weights
    // equal to 1: the projected surface and all its derivatives must reproduce
    // the non-rational ones to machine precision.
    gbs::BSSurfaceRational<double, 3> make_unit_weight_paraboloid()
    {
        const std::vector<double> w(paraboloid_poles.size(), 1.0);
        return gbs::BSSurfaceRational<double, 3>(paraboloid_poles, w, paraboloid_ku, paraboloid_kv, 2, 2);
    }

    // A genuinely weighted bidegree-(2,2) patch: a quarter-circle profile in the
    // (x, z) plane (NURBS unit arc, weights 1, 1/sqrt(2), 1) swept linearly in y.
    // The non-trivial middle weight exercises the A4.4 quotient rule.
    gbs::BSSurfaceRational<double, 3> make_weighted_patch()
    {
        const std::vector<double> ku{0., 0., 0., 1., 1., 1.};
        const std::vector<double> kv{0., 0., 0., 1., 1., 1.};
        const double s = 1.0 / std::sqrt(2.0);
        // Arc control points (x, _, z); y filled per v-row below.
        const gbs::points_vector<double, 3> poles{
            {1., 0., 0.}, {1., 0., 1.}, {0., 0., 1.},
            {1., 1., 0.}, {1., 1., 1.}, {0., 1., 1.},
            {1., 2., 0.}, {1., 2., 1.}, {0., 2., 1.},
        };
        const std::vector<double> weights{
            1., s, 1.,
            1., s, 1.,
            1., s, 1.,
        };
        return gbs::BSSurfaceRational<double, 3>(poles, weights, ku, kv, 2, 2);
    }

    // Exact NURBS quarter-cylinder of radius R: a rational quadratic unit-circle
    // arc (control points (R,0),(R,R),(0,R), weights 1, 1/sqrt(2), 1) in the
    // (x, z) plane, swept linearly (degree 1) along y. The u-iso curves are exact
    // circles of radius R, so the arc-length second derivative d_dmu2 has the
    // closed form -(x, 0, z)/R^2 (curvature vector, |.| = 1/R); the straight
    // v-rulings give d_dmv2 == 0.
    constexpr double cylinder_R = 2.0;
    gbs::BSSurfaceRational<double, 3> make_quarter_cylinder()
    {
        const std::vector<double> ku{0., 0., 0., 1., 1., 1.}; // degU = 2
        const std::vector<double> kv{0., 0., 1., 1.};         // degV = 1
        const double R = cylinder_R;
        const gbs::points_vector<double, 3> poles{
            {R, 0., 0.}, {R, 0., R}, {0., 0., R}, // y = 0 row
            {R, R, 0.}, {R, R, R}, {0., R, R},    // y = R row
        };
        const double s = 1.0 / std::sqrt(2.0);
        const std::vector<double> weights{
            1., s, 1.,
            1., s, 1.,
        };
        return gbs::BSSurfaceRational<double, 3>(poles, weights, ku, kv, 2, 1);
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

// A rational surface with all weights == 1 must reproduce the non-rational
// mixed derivatives to machine precision (the quotient rule collapses to the
// homogeneous numerator). This catches the bulk of A4.4 transcription bugs.
TEST(tests_surface_derivatives, rational_all_weights_one_matches_non_rational)
{
    const auto srf  = make_paraboloid();
    const auto rsrf = make_unit_weight_paraboloid();
    constexpr double tol = 1e-12;

    for (double u : {0.0, 0.2, 0.5, 0.8, 1.0})
    for (double v : {0.0, 0.2, 0.5, 0.8, 1.0})
        for (size_t du : {size_t{0}, size_t{1}, size_t{2}})
        for (size_t dv : {size_t{0}, size_t{1}, size_t{2}})
        {
            CAPTURE(u);
            CAPTURE(v);
            CAPTURE(du);
            CAPTURE(dv);
            const auto ref = srf.value(u, v, du, dv);
            const auto got = rsrf.value(u, v, du, dv);
            for (size_t c = 0; c < 3; ++c)
                EXPECT_NEAR(got[c], ref[c], tol);
        }

    // The arc-length helpers route through derivatives(); they must now work
    // (no throw) on rational surfaces and match the non-rational values.
    for (double u : {0.2, 0.5, 0.8})
    for (double v : {0.2, 0.5, 0.8})
    {
        CAPTURE(u);
        CAPTURE(v);
        const auto u2  = rsrf.d_dmu2(u, v);
        const auto v2  = rsrf.d_dmv2(u, v);
        const auto u2r = srf.d_dmu2(u, v);
        const auto v2r = srf.d_dmv2(u, v);
        for (size_t c = 0; c < 3; ++c)
        {
            EXPECT_NEAR(u2[c], u2r[c], tol);
            EXPECT_NEAR(v2[c], v2r[c], tol);
        }
    }
}

// Weighted (genuine NURBS) derivatives must match central finite differences
// of value(u, v, 0, 0). Looser tolerance because the reference is itself an
// O(h^2) approximation.
TEST(tests_surface_derivatives, rational_weighted_matches_finite_differences)
{
    const auto srf = make_weighted_patch();
    const double h = 1e-4;
    const double tol = 1e-6;

    auto S = [&](double u, double v) { return srf.value(u, v, 0, 0); };

    for (double u : {0.25, 0.5, 0.75})
    for (double v : {0.25, 0.5, 0.75})
    {
        CAPTURE(u);
        CAPTURE(v);

        const auto Su  = srf.value(u, v, 1, 0);
        const auto Sv  = srf.value(u, v, 0, 1);
        const auto Suu = srf.value(u, v, 2, 0);
        const auto Svv = srf.value(u, v, 0, 2);
        const auto Suv = srf.value(u, v, 1, 1);

        const auto Sp0 = S(u + h, v), Sm0 = S(u - h, v);
        const auto S0p = S(u, v + h), S0m = S(u, v - h);
        const auto S00 = S(u, v);
        const auto Spp = S(u + h, v + h), Spm = S(u + h, v - h);
        const auto Smp = S(u - h, v + h), Smm = S(u - h, v - h);

        for (size_t c = 0; c < 3; ++c)
        {
            const double fd_u  = (Sp0[c] - Sm0[c]) / (2 * h);
            const double fd_v  = (S0p[c] - S0m[c]) / (2 * h);
            const double fd_uu = (Sp0[c] - 2 * S00[c] + Sm0[c]) / (h * h);
            const double fd_vv = (S0p[c] - 2 * S00[c] + S0m[c]) / (h * h);
            const double fd_uv = (Spp[c] - Spm[c] - Smp[c] + Smm[c]) / (4 * h * h);
            EXPECT_NEAR(Su[c],  fd_u,  tol);
            EXPECT_NEAR(Sv[c],  fd_v,  tol);
            EXPECT_NEAR(Suu[c], fd_uu, 1e-4);
            EXPECT_NEAR(Svv[c], fd_vv, 1e-4);
            EXPECT_NEAR(Suv[c], fd_uv, 1e-4);
        }
    }
}

// Arc-length second derivative on an exact NURBS circle. For the radius-R
// quarter-cylinder, the u-iso curves are exact circles, so d_dmu2 is the
// curvature vector: it points from the surface point to the axis with
// magnitude 1/R, i.e. d_dmu2 == -(x, 0, z)/R^2. The straight v-rulings give
// d_dmv2 == 0. Both hold to machine precision (exact rational geometry,
// closed-form A4.4 derivatives -- no finite differencing).
TEST(tests_surface_derivatives, d_dmu2_on_nurbs_circle_is_curvature_vector)
{
    const auto srf = make_quarter_cylinder();
    const double R = cylinder_R;
    constexpr double tol = 1e-11;

    for (double u : {0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0})
    for (double v : {0.0, 0.3, 0.7, 1.0})
    {
        CAPTURE(u);
        CAPTURE(v);
        const auto P = srf.value(u, v, 0, 0);

        // Sanity: the (x, z) section lies exactly on the circle of radius R.
        EXPECT_NEAR(P[0] * P[0] + P[2] * P[2], R * R, tol);

        // Curvature vector of the u-iso circle: toward the axis, magnitude 1/R.
        const auto a = srf.d_dmu2(u, v);
        EXPECT_NEAR(a[0], -P[0] / (R * R), tol);
        EXPECT_NEAR(a[1], 0.0,             tol);
        EXPECT_NEAR(a[2], -P[2] / (R * R), tol);
        EXPECT_NEAR(std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]), 1.0 / R, tol);

        // The v-rulings are straight lines: zero arc-length curvature.
        const auto b = srf.d_dmv2(u, v);
        EXPECT_NEAR(b[0], 0.0, tol);
        EXPECT_NEAR(b[1], 0.0, tol);
        EXPECT_NEAR(b[2], 0.0, tol);
    }
}
