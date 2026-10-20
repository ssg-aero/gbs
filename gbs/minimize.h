#pragma once
/**
 * @file minimize.h
 * @brief Derivative-free minimization of a function of a few variables in a box: the
 * Nelder-Mead simplex method.
 *
 * For objectives without derivatives: shape parameters fitted to angles, distances measured
 * in a transformed space, non-smooth sums of absolute values. Where the problem is a closest
 * point between geometries, the Newton projections of `gbs/extrema.h` are faster and exact;
 * a function of one variable is better served by Brent (`boost::math::tools::brent_find_minima`).
 *
 * The bounds are held by a change of variables, the simplex moving in z: x = lo + (hi - lo)
 * (1 + sin z) / 2 between two bounds, x = lo + z² or hi - z² with one bound, x = z without.
 * Clamping the trial points on the bounds instead lays the simplex flat on a face of the box,
 * from which it cannot leave. Trial points are shifted by whole periods to the period of the
 * best vertex, so that the centroid of the simplex is not a mirror point of the box. Once
 * converged, the method starts again from the best point with a fresh simplex, as long as this
 * lowers f (at most 5 restarts): a simplex collapsed onto a line stops short of the minimum,
 * the restart recovers it. Design:
 * docs/sources/design/nelder_mead.md.
 */
#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <limits>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <type_traits>

namespace gbs
{
    /// Result of `nelder_mead`.
    template <std::floating_point T, std::size_t n>
    struct MinimizeResult
    {
        std::array<T, n> x{};      ///< best point found, inside the box
        T f{};                     ///< f(x), +inf when f was defined nowhere the simplex went
        std::size_t evaluations{}; ///< evaluations of f, restarts included
        bool converged{};          ///< false when stopped by the maximum number of evaluations, or f(x) is not finite
    };

    namespace minimize_detail
    {
        /// The change of variables x(z) holding x in [lo, hi].
        template <std::floating_point T, std::size_t n>
        struct Box
        {
            std::array<T, n> lo, hi;

            auto x(const std::array<T, n> &z) const -> std::array<T, n>
            {
                std::array<T, n> x;
                for (std::size_t i{}; i < n; ++i)
                {
                    const bool l = std::isfinite(lo[i]), h = std::isfinite(hi[i]);
                    if (l && h)
                        x[i] = std::clamp(lo[i] + (hi[i] - lo[i]) * (T(1) + std::sin(z[i])) / T(2), lo[i], hi[i]);
                    else if (l)
                        x[i] = lo[i] + z[i] * z[i];
                    else if (h)
                        x[i] = hi[i] - z[i] * z[i];
                    else
                        x[i] = z[i];
                }
                return x;
            }

            /// z shifted by whole periods (2 pi, between two bounds) to the period of ref, the
            /// same x: vertices spread over several periods have their centroid at a mirror point
            /// of the box. The symmetric z (pi - z, or -z with one bound) are left alone: a vertex
            /// beyond a bound is how the simplex straddles it, a minimum on a bound being a smooth
            /// one in z (folding them back lays the simplex against the bound).
            auto unwrap(std::array<T, n> z, const std::array<T, n> &ref) const -> std::array<T, n>
            {
                constexpr T two_pi = T(2) * std::numbers::pi_v<T>;
                for (std::size_t i{}; i < n; ++i)
                    if (std::isfinite(lo[i]) && std::isfinite(hi[i]))
                        z[i] -= two_pi * std::round((z[i] - ref[i]) / two_pi);
                return z;
            }

            /// z of one variable, x in [lo, hi].
            auto z(std::size_t i, T x) const -> T
            {
                const bool l = std::isfinite(lo[i]), h = std::isfinite(hi[i]);
                if (l && h)
                    return hi[i] > lo[i] ? std::asin(std::clamp(T(2) * (x - lo[i]) / (hi[i] - lo[i]) - T(1), T(-1), T(1))) : T(0);
                if (l)
                    return std::sqrt(x - lo[i]);
                if (h)
                    return std::sqrt(hi[i] - x);
                return x;
            }
        };

        /// f(x) as T, NaN taken as +inf: a point where the objective is not defined is the worst one.
        template <std::floating_point T, std::size_t n, typename F>
        auto evaluate(F &f, const std::array<T, n> &x) -> T
        {
            const T y = static_cast<T>(std::invoke(f, x));
            return std::isnan(y) ? std::numeric_limits<T>::infinity() : y;
        }

        /**
         * One Nelder-Mead run in z from x0, the initial simplex stepping by h on each x_i (into
         * the box). Stops when, on every coordinate, the vertices are within tol_x (1 + |x_i|) of
         * the best one in x, or when `budget` evaluations are spent (`evaluations` counts them).
         */
        template <std::floating_point T, std::size_t n, typename F>
        auto run(F &f, const Box<T, n> &box, const std::array<T, n> &x0, T f0, const std::array<T, n> &h, T tol_x,
                 std::size_t budget, std::size_t &evaluations) -> MinimizeResult<T, n>
        {
            // reflection, expansion, contraction and shrink coefficients of Nelder and Mead
            constexpr T alpha = T(1), beta = T(2), gamma = T(0.5), delta = T(0.5);
            constexpr T inf = std::numeric_limits<T>::infinity();

            // vertices in z, their points in x and f there
            std::array<std::array<T, n>, n + 1> z, x;
            std::array<T, n + 1> fv;
            std::size_t best{}; // the best vertex, reference of the unwrapping
            auto set = [&](std::size_t k, const std::array<T, n> &zk) {
                z[k] = box.unwrap(zk, z[best]);
                x[k] = box.x(zk);
                fv[k] = evaluations < budget ? (++evaluations, evaluate<T, n>(f, x[k])) : inf;
            };
            // c + t (y - c)
            auto along = [](const std::array<T, n> &c, const std::array<T, n> &y, T t) {
                std::array<T, n> r;
                for (std::size_t i{}; i < n; ++i)
                    r[i] = c[i] + t * (y[i] - c[i]);
                return r;
            };

            for (std::size_t i{}; i < n; ++i)
                z[0][i] = box.z(i, x0[i]);
            x[0] = x0;
            fv[0] = f0;
            for (std::size_t j{}; j < n; ++j)
            {
                auto zj = z[0];
                // step into the box: away from the upper bound when x0 is on it
                const T xj = x0[j] + h[j] <= box.hi[j] ? x0[j] + h[j] : x0[j] - h[j];
                zj[j] = box.z(j, std::clamp(xj, box.lo[j], box.hi[j]));
                set(j + 1, zj);
            }

            std::array<std::size_t, n + 1> order;
            for (;;)
            {
                std::iota(order.begin(), order.end(), std::size_t{});
                std::ranges::stable_sort(order, {}, [&](std::size_t k) { return fv[k]; });
                best = order[0];
                const std::size_t worst = order[n], second = order[n - 1];

                bool small = true;
                for (std::size_t i{}; i < n && small; ++i)
                    for (std::size_t k{}; k <= n && small; ++k)
                        small = std::abs(x[k][i] - x[best][i]) <= tol_x * (T(1) + std::abs(x[best][i]));
                if (small || evaluations >= budget)
                    return {x[best], fv[best], evaluations, small && std::isfinite(fv[best])};

                std::array<T, n> c{};
                for (std::size_t k{}; k <= n; ++k)
                    if (k != worst)
                        for (std::size_t i{}; i < n; ++i)
                            c[i] += z[k][i] / T(n);

                // trial vertices are evaluated in the slot of the worst one, which is kept aside
                const auto zw = z[worst];
                const auto xw = x[worst];
                const T fw = fv[worst];
                auto restore = [&] {
                    z[worst] = zw;
                    x[worst] = xw;
                    fv[worst] = fw;
                };

                set(worst, along(c, zw, -alpha));
                const auto zr = z[worst];
                const auto xr = x[worst];
                const T fr = fv[worst];
                if (fr < fv[best])
                {
                    set(worst, along(c, zw, -alpha * beta));
                    if (!(fv[worst] < fr))
                    {
                        z[worst] = zr;
                        x[worst] = xr;
                        fv[worst] = fr;
                    }
                    continue;
                }
                if (fr < fv[second])
                    continue;
                // contraction: outside (towards the reflected point) when it improves on the worst point, inside otherwise
                const bool outside = fr < fw;
                set(worst, outside ? along(c, zr, gamma) : along(c, zw, gamma));
                if (fv[worst] < (outside ? fr : fw))
                    continue;
                restore();
                // shrink towards the best point
                for (std::size_t k{}; k <= n; ++k)
                    if (k != best)
                        set(k, along(z[best], z[k], delta));
            }
        }
    }

    /**
     * @brief Minimizes f over the box [lo, hi] by the Nelder-Mead simplex method, without derivatives.
     *
     * @param f callable taking a `const std::array<T, n> &` and returning a value convertible to T,
     *        called non-const (a `mutable` lambda or a stateful functor counts its calls); NaN is
     *        taken as +inf (the point is rejected)
     * @param x0 start, clamped into the box
     * @param lo, hi bounds, lo[i] <= hi[i]; infinite bounds leave the variable free, lo[i] == hi[i] holds it
     * @param tol_x stop when every vertex of the simplex is within tol_x (1 + |x_i|) of the best one
     * @param max_evaluations maximum number of evaluations of f, restarts included
     * @param initial_step the steps of the initial simplex (and of the restarts), relative to the
     *        box on each variable (to max(1, |x0_i|) when the variable is unbounded)
     * @return the best point, f there, the number of evaluations and whether tol_x was reached
     * @throw std::invalid_argument if tol_x <= 0, max_evaluations == 0, initial_step is not
     *        positive and finite, lo[i] > hi[i] or x0 is not finite
     *
     * The minimum found is local: start close to the expected one. The initial step chooses
     * among the basins around x0: a larger one looks further (and may jump to another basin), a
     * smaller one stays closer to x0.
     */
    template <std::floating_point T, std::size_t n, typename F>
        requires(n > 0) && std::invocable<F &, const std::array<T, n> &> &&
                std::convertible_to<std::invoke_result_t<F &, const std::array<T, n> &>, T>
    auto nelder_mead(F &&f, std::array<T, n> x0, const std::array<T, n> &lo, const std::array<T, n> &hi,
                     std::type_identity_t<T> tol_x, std::size_t max_evaluations = 1000 * n,
                     std::type_identity_t<T> initial_step = T(0.1)) -> MinimizeResult<T, n>
    {
        if (!(tol_x > T(0)))
            throw std::invalid_argument("nelder_mead: tol_x must be positive");
        if (max_evaluations == 0)
            throw std::invalid_argument("nelder_mead: max_evaluations must be positive");
        if (!(initial_step > T(0)) || !std::isfinite(initial_step))
            throw std::invalid_argument("nelder_mead: initial_step must be positive and finite");
        std::array<T, n> h;
        for (std::size_t i{}; i < n; ++i)
        {
            if (!(lo[i] <= hi[i]))
                throw std::invalid_argument("nelder_mead: lo must not be above hi");
            if (!std::isfinite(x0[i]))
                throw std::invalid_argument("nelder_mead: x0 must be finite");
            x0[i] = std::clamp(x0[i], lo[i], hi[i]);
            const T width = hi[i] - lo[i];
            h[i] = initial_step * (std::isfinite(width) ? width : std::max(T(1), std::abs(x0[i])));
        }

        const minimize_detail::Box<T, n> box{lo, hi};
        std::size_t evaluations{1};
        auto best = minimize_detail::run<T, n>(f, box, x0, minimize_detail::evaluate<T, n>(f, x0), h, tol_x,
                                               max_evaluations, evaluations);
        // restart from the best point while it lowers f (a collapsed simplex stops short)
        for (int restart{}; restart < 5 && best.converged && evaluations < max_evaluations; ++restart)
        {
            auto next = minimize_detail::run<T, n>(f, box, best.x, best.f, h, tol_x, max_evaluations, evaluations);
            if (!(next.f < best.f))
                break;
            best = next;
        }
        best.evaluations = evaluations;
        return best;
    }
}
