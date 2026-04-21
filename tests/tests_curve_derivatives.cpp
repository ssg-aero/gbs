#include <gtest/gtest.h>
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
