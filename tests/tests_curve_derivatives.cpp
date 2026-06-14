#include "doctest_gtest.hpp"
#include <gbs/bscurve.h>

#include <array>
#include <cmath>
#include <vector>

namespace
{
    // Parabola C(u) = (u, u^2) expressed as a degree-2 Bezier on u in [0, 1].
    // Control points (0,0), (0.5,0), (1,1) reproduce the parabola exactly,
    // which gives closed-form derivatives to validate d_dm / d_dm2 against:
    //   C'(u)  = (1, 2u)
    //   C''(u) = (0, 2)
    gbs::BSCurve<double, 2> make_parabola()
    {
        const std::vector<double> knots{0., 0., 0., 1., 1., 1.};
        const std::vector<std::array<double, 2>> poles{{0., 0.}, {0.5, 0.}, {1., 1.}};
        return gbs::BSCurve<double, 2>(poles, knots, 2);
    }
}

TEST(tests_curve_derivatives, d_dm_is_unit_tangent_on_parabola)
{
    const auto crv = make_parabola();
    constexpr double tol = 1e-12;

    for (double u : {0.0, 0.25, 0.5, 0.75, 1.0})
    {
        const auto t = crv.d_dm(u);

        // Expected unit tangent: (1, 2u) / sqrt(1 + 4u^2)
        const double inv_n = 1.0 / std::sqrt(1.0 + 4.0 * u * u);
        const double tx = 1.0 * inv_n;
        const double ty = 2.0 * u * inv_n;

        EXPECT_NEAR(t[0], tx, tol) << "u=" << u;
        EXPECT_NEAR(t[1], ty, tol) << "u=" << u;
        EXPECT_NEAR(std::sqrt(t[0] * t[0] + t[1] * t[1]), 1.0, tol) << "u=" << u;
    }
}

TEST(tests_curve_derivatives, d_dm2_matches_closed_form_on_parabola)
{
    const auto crv = make_parabola();
    constexpr double tol = 1e-12;

    // For C(u) = (u, u^2):
    //   d2C/dm2 = (C'' |C'|^2 - C' (C'.C'')) / |C'|^4
    //           = (-4u, 2) / (1 + 4u^2)^2
    for (double u : {0.0, 0.25, 0.5, 0.75, 1.0})
    {
        const auto a = crv.d_dm2(u);

        const double N  = 1.0 + 4.0 * u * u;
        const double N2 = N * N;
        const double ax = -4.0 * u / N2;
        const double ay =  2.0      / N2;

        EXPECT_NEAR(a[0], ax, tol) << "u=" << u;
        EXPECT_NEAR(a[1], ay, tol) << "u=" << u;
    }
}

namespace
{
    template <size_t dim>
    double diff_norm(const std::array<double, dim> &a, const std::array<double, dim> &b)
    {
        double s = 0.;
        for (size_t c = 0; c < dim; ++c) { double e = a[c] - b[c]; s += e * e; }
        return std::sqrt(s);
    }

    // Clamped knot vector of degree p with one interior knot of multiplicity p
    // (a C^0 joint, where one-sided derivatives are discontinuous).
    std::vector<double> c0_knots(size_t p)
    {
        std::vector<double> k;
        for (size_t i = 0; i <= p; ++i) k.push_back(0.);
        k.push_back(1.);
        for (size_t i = 0; i < p; ++i) k.push_back(2.); // interior mult == p
        k.push_back(3.);
        for (size_t i = 0; i <= p; ++i) k.push_back(4.);
        return k;
    }
}

// The one-pass multi-order derivatives(u, n, CK) must agree with calling
// value(u, k) per order, for every order including beyond the degree, sampled
// on a grid plus every distinct knot (incl. a multiplicity-p C^0 joint).
TEST(tests_curve_derivatives, derivatives_match_per_order_value_nonrational)
{
    using T = double;
    constexpr size_t nmax = 4;
    for (size_t p : {size_t{2}, size_t{3}, size_t{5}})
    {
        const auto k = c0_knots(p);
        const size_t n_poles = k.size() - p - 1;
        std::vector<std::array<T, 3>> poles(n_poles);
        for (size_t i = 0; i < n_poles; ++i)
            poles[i] = {T(i), T((i * 5) % 7), T((i * 2) % 3)};
        gbs::BSCurve<T, 3> crv(poles, k, p);

        auto us = gbs::make_range(k.front(), k.back(), 40);
        for (auto kv : k) us.push_back(kv);

        std::array<std::array<T, 3>, nmax + 1> CK;
        for (auto u : us)
        {
            crv.derivatives(u, nmax, CK.data());
            for (size_t d = 0; d <= nmax; ++d)
                ASSERT_LT(diff_norm<3>(CK[d], crv.value(u, d)), 1e-9)
                    << "p=" << p << " d=" << d << " u=" << u;
        }
    }
}

// Same for the rational path: the one-pass A4.2 derivatives must match the
// per-order rational value() (eval_rational_value_simple).
TEST(tests_curve_derivatives, derivatives_match_per_order_value_rational)
{
    using T = double;
    constexpr size_t nmax = 3;
    for (size_t p : {size_t{2}, size_t{3}})
    {
        const auto k = c0_knots(p);
        const size_t n_poles = k.size() - p - 1;
        std::vector<std::array<T, 3>> poles(n_poles);
        std::vector<T> weights(n_poles);
        for (size_t i = 0; i < n_poles; ++i)
        {
            poles[i] = {T(i), T((i * 5) % 7), T((i * 2) % 3)};
            weights[i] = T(1) + T((i * 3) % 5) * T(0.5); // non-trivial weights
        }
        gbs::BSCurveRational<T, 3> crv(poles, k, weights, p);

        auto us = gbs::make_range(k.front(), k.back(), 40);
        for (auto kv : k) us.push_back(kv);

        std::array<std::array<T, 3>, nmax + 1> CK;
        for (auto u : us)
        {
            crv.derivatives(u, nmax, CK.data());
            for (size_t d = 0; d <= nmax; ++d)
                ASSERT_LT(diff_norm<3>(CK[d], crv.value(u, d)), 1e-9)
                    << "p=" << p << " d=" << d << " u=" << u;
        }
    }
}
