#pragma once
/**
 * @file quadrature.h
 * @brief Adaptive Gauss-Kronrod quadrature over an interval of any width.
 */
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/tools/precision.hpp>

#include <cmath>
#include <cstddef>

namespace gbs
{
    /**
     * @brief The integral of f over [a, b] by Boost's adaptive Gauss-Kronrod rule of N points, taken on the variable
     * x of [-1, 1] (u = (a + b) / 2 + x (b - a) / 2), so that its stopping test does not depend on the width of
     * [a, b].
     *
     * Boost compares the error of its rule on [-1, 1], not scaled, with a tolerance scaled by the half-width of the
     * interval: an interval narrower than about 2 eps / tol (3e-8 by default) is bisected down to max_depth whatever
     * f (983 025 evaluations of a constant for N = 15), and a wide one stops before reaching tol. On [-1, 1] both
     * are at the same scale: tol is relative to the integral, whatever the width. Infinite bounds are left to Boost
     * (its own changes of variable).
     *
     * @param f the integrand, of u
     * @param error if not null, the estimated error of the integral
     * @param L1 if not null, the integral of |f| over [a, b]
     */
    template <typename T, std::size_t N = 15, typename F>
    auto integrate_kronrod(F f, T a, T b, unsigned max_depth = 15, T tol = boost::math::tools::root_epsilon<T>(),
                           T *error = nullptr, T *L1 = nullptr)
    {
        using rule = boost::math::quadrature::gauss_kronrod<T, N>;
        if (!std::isfinite(a) || !std::isfinite(b))
            return rule::integrate(f, a, b, max_depth, tol, error, L1);
        const T mid = a / 2 + b / 2, half = b / 2 - a / 2; // no overflow for the widest finite bounds
        return rule::integrate([&](T x) { return half * f(mid + half * x); }, T(-1), T(1), max_depth, tol, error, L1);
    }
} // namespace gbs
