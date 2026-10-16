#pragma once
/**
 * @file projection.h
 * @brief Closest points by Newton's method: point-curve, point-surface, curve-curve and
 * curve-surface, in any dimension, on bounded parameter ranges.
 *
 * Each problem minimizes f(x) = |r(x)|² / 2 for a residual r (the difference of the two
 * points) of 1 to 3 parameters x kept in their box. A step solves H dx = -g with the full
 * Hessian H = JᵀJ + Σ r_k ∇²r_k, so that the convergence stays quadratic when the distance
 * is not zero; Gauss-Newton (JᵀJ) is used where H is not positive definite, a diagonal step
 * where JᵀJ is singular. Variables on a bound with the gradient pushing out of the box are
 * held for the step, and the step is halved until f does not increase. Stopped on a bound
 * where the opposite bound gives the same distance (a seam), Newton starts again across it.
 *
 * The local functions (`newton_*`) refine from a start; the global ones sample the ranges and
 * run Newton from the two nearest samples (on a closed curve or surface the nearest one may
 * sit at the wrong end of the seam). Design: docs/sources/design/nlopt_removal.md.
 */
#include <gbs/bssurf.h>
#include <gbs/execution.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs
{
    namespace projection_detail
    {
        template <std::floating_point T, std::size_t n>
        using vec = std::array<T, n>;
        template <std::floating_point T, std::size_t n>
        using mat = std::array<std::array<T, n>, n>;

        /// Gradient, Hessian and Gauss-Newton matrix JᵀJ of f = |r|² / 2.
        template <std::floating_point T, std::size_t n>
        struct Model
        {
            vec<T, n> g{};
            mat<T, n> H{}; ///< full Hessian (JᵀJ alone when the second derivatives are not available)
            mat<T, n> G{}; ///< JᵀJ
        };

        /// Solves A x = b by Cholesky; false when A is not (numerically) positive definite.
        template <std::floating_point T, std::size_t n>
        auto cholesky_solve(const mat<T, n> &A, const vec<T, n> &b, vec<T, n> &x) -> bool
        {
            T scale{};
            for (std::size_t i{}; i < n; ++i)
                scale = std::max(scale, std::abs(A[i][i]));
            if (!(scale > T(0)))
                return false;
            mat<T, n> L{};
            for (std::size_t j{}; j < n; ++j)
            {
                T d = A[j][j];
                for (std::size_t k{}; k < j; ++k)
                    d -= L[j][k] * L[j][k];
                if (!(d > T(16) * std::numeric_limits<T>::epsilon() * scale))
                    return false;
                L[j][j] = std::sqrt(d);
                for (std::size_t i = j + 1; i < n; ++i)
                {
                    T s = A[i][j];
                    for (std::size_t k{}; k < j; ++k)
                        s -= L[i][k] * L[j][k];
                    L[i][j] = s / L[j][j];
                }
            }
            vec<T, n> y{};
            for (std::size_t i{}; i < n; ++i)
            {
                T s = b[i];
                for (std::size_t k{}; k < i; ++k)
                    s -= L[i][k] * y[k];
                y[i] = s / L[i][i];
            }
            for (std::size_t i = n; i-- > 0;)
            {
                T s = y[i];
                for (std::size_t k = i + 1; k < n; ++k)
                    s -= L[k][i] * x[k];
                x[i] = s / L[i][i];
            }
            return true;
        }

        /**
         * Damped Newton on f in the box [lo, hi]. `value(x)` gives f, `model(x)` its derivatives.
         * Stops when every step |dx_i| <= max(tol_x², 4 eps) (1 + |x_i|), when f cannot decrease
         * any more, or after 100 iterations. Then, at most `seams` times, starts again across a
         * seam (see below). Returns x and f(x).
         */
        template <std::floating_point T, std::size_t n, typename FValue, typename FModel>
        auto minimize(const FValue &value, const FModel &model, vec<T, n> x, const vec<T, n> &lo, const vec<T, n> &hi, T tol_x,
                      int seams = int(n)) -> std::pair<vec<T, n>, T>
        {
            const T tau = std::max(tol_x * tol_x, T(4) * std::numeric_limits<T>::epsilon());
            for (std::size_t i{}; i < n; ++i)
                x[i] = std::clamp(x[i], lo[i], hi[i]);
            T f = value(x);
            for (int it{}; it < 100 && f > T(0); ++it)
            {
                const Model<T, n> m = model(x);
                // variables on a bound with the gradient pushing out of the box are held
                std::array<bool, n> held{};
                bool all_held = true;
                vec<T, n> b{};
                for (std::size_t i{}; i < n; ++i)
                {
                    held[i] = (x[i] <= lo[i] && m.g[i] > T(0)) || (x[i] >= hi[i] && m.g[i] < T(0));
                    all_held = all_held && held[i];
                    b[i] = held[i] ? T(0) : -m.g[i];
                }
                if (all_held)
                    break;
                auto reduced = [&](mat<T, n> A) {
                    for (std::size_t i{}; i < n; ++i)
                        if (held[i])
                            for (std::size_t j{}; j < n; ++j)
                                A[i][j] = A[j][i] = T(i == j);
                    return A;
                };
                vec<T, n> dx{};
                if (!cholesky_solve(reduced(m.H), b, dx) && !cholesky_solve(reduced(m.G), b, dx))
                    for (std::size_t i{}; i < n; ++i) // singular JᵀJ (parallel tangents, null derivative)
                        dx[i] = !held[i] && m.G[i][i] > T(0) ? b[i] / m.G[i][i] : T(0);

                // halve the step until f does not increase
                vec<T, n> x_new{};
                T f_new{};
                bool moved{};
                T alpha{1};
                for (int ls{}; ls < 30; ++ls, alpha /= T(2))
                {
                    for (std::size_t i{}; i < n; ++i)
                        x_new[i] = std::clamp(x[i] + alpha * dx[i], lo[i], hi[i]);
                    f_new = value(x_new);
                    if (f_new <= f)
                    {
                        moved = true;
                        break;
                    }
                }
                if (!moved)
                    break; // minimum reached at the machine precision
                bool small = true;
                for (std::size_t i{}; i < n; ++i)
                    small = small && std::abs(x_new[i] - x[i]) <= tau * (T(1) + std::abs(x[i]));
                x = x_new;
                f = f_new;
                if (small)
                    break;
            }
            // stopped on a bound where the opposite bound gives the same distance (the seam of a
            // closed curve or surface): the minimum may lie on the other side, start again there
            for (std::size_t i{}; i < n && seams > 0 && f > T(0); ++i)
            {
                if (x[i] != lo[i] && x[i] != hi[i])
                    continue;
                auto y = x;
                y[i] = x[i] == lo[i] ? hi[i] : lo[i];
                if (!(std::abs(value(y) - f) <= T(16) * std::sqrt(std::numeric_limits<T>::epsilon()) * f))
                    continue;
                if (const auto other = minimize<T, n>(value, model, y, lo, hi, tol_x, seams - 1); other.second < f)
                {
                    x = other.first;
                    f = other.second;
                }
            }
            return {x, f};
        }

        /// Second derivatives of a geometry that does not provide them (such as an offset
        /// surface) are given up after the first failure: Gauss-Newton is used instead.
        template <typename F>
        auto second_order(bool &enabled, const F &f) -> void
        {
            if (!enabled)
                return;
            try
            {
                f();
            }
            catch (const std::exception &)
            {
                enabled = false;
            }
        }

        /// The two samples nearest to the target, as starts for Newton.
        template <std::floating_point T, std::size_t n>
        struct TwoBest
        {
            std::array<std::pair<T, vec<T, n>>, 2> best{std::pair{std::numeric_limits<T>::max(), vec<T, n>{}},
                                                        std::pair{std::numeric_limits<T>::max(), vec<T, n>{}}};
            auto push(T d, const vec<T, n> &x) -> void
            {
                if (d < best[0].first)
                {
                    best[1] = best[0];
                    best[0] = {d, x};
                }
                else if (d < best[1].first)
                    best[1] = {d, x};
            }
        };

        /// Distances of every sample (in parallel above gbs::parallel_min_size), then the two nearest.
        template <std::floating_point T, std::size_t n, typename F>
        auto two_best(const std::vector<vec<T, n>> &samples, const F &dist) -> TwoBest<T, n>
        {
            std::vector<T> d(samples.size());
            transform_threshold(samples.begin(), samples.end(), d.begin(), dist);
            TwoBest<T, n> tb;
            for (std::size_t k{}; k < samples.size(); ++k)
                tb.push(d[k], samples[k]);
            return tb;
        }

        /// n samples of [a, b], ends included (n >= 2).
        template <std::floating_point T>
        auto samples(T a, T b, std::size_t n) -> std::vector<T>
        {
            n = std::max<std::size_t>(n, 2);
            std::vector<T> s(n);
            for (std::size_t k{}; k < n; ++k)
                s[k] = k + 1 == n ? b : a + (b - a) * T(k) / T(n - 1);
            return s;
        }

        /// Runs `local` from the two starts and keeps the closest result (the last entry is the distance).
        template <typename R, std::floating_point T, std::size_t n, typename F>
        auto from_two_best(const TwoBest<T, n> &tb, const F &local) -> R
        {
            R best = local(tb.best[0].second);
            if (tb.best[1].first < std::numeric_limits<T>::max())
                if (const R other = local(tb.best[1].second); other.back() < best.back())
                    best = other;
            return best;
        }
    } // namespace projection_detail

    /**
     * @brief Projection of p on c restricted to range, by Newton from t0.
     * @return {t, distance}
     */
    template <std::floating_point T, std::size_t dim>
    auto newton_point_curve(const Curve<T, dim> &c, const point<T, dim> &p, T t0, const std::array<T, 2> &range, T tol_x)
        -> std::array<T, 2>
    {
        using namespace projection_detail;
        bool second{true};
        auto value = [&](const vec<T, 1> &x) { return sq_norm(c.value(x[0]) - p) / T(2); };
        auto model = [&](const vec<T, 1> &x) {
            const auto r = c.value(x[0]) - p;
            const auto d1 = c.value(x[0], 1);
            Model<T, 1> m;
            m.g = {d1 * r};
            m.G = {{{d1 * d1}}};
            m.H = m.G;
            second_order(second, [&] { m.H[0][0] += c.value(x[0], 2) * r; });
            return m;
        };
        const auto [x, f] = minimize<T, 1>(value, model, {t0}, {range[0]}, {range[1]}, tol_x);
        return {x[0], std::sqrt(T(2) * f)};
    }

    /**
     * @brief Projection of p on s restricted to box {u1, u2, v1, v2}, by Newton from (u0, v0).
     * @return {u, v, distance}
     */
    template <std::floating_point T, std::size_t dim>
    auto newton_point_surface(const Surface<T, dim> &s, const point<T, dim> &p, T u0, T v0, const std::array<T, 4> &box, T tol_x)
        -> std::array<T, 3>
    {
        using namespace projection_detail;
        bool second{true};
        auto value = [&](const vec<T, 2> &x) { return sq_norm(s.value(x[0], x[1]) - p) / T(2); };
        auto model = [&](const vec<T, 2> &x) {
            const auto r = s.value(x[0], x[1]) - p;
            const auto su = s.value(x[0], x[1], 1, 0), sv = s.value(x[0], x[1], 0, 1);
            Model<T, 2> m;
            m.g = {su * r, sv * r};
            m.G = {{{su * su, su * sv}, {su * sv, sv * sv}}};
            m.H = m.G;
            second_order(second, [&] {
                const T huu = s.value(x[0], x[1], 2, 0) * r, huv = s.value(x[0], x[1], 1, 1) * r, hvv = s.value(x[0], x[1], 0, 2) * r;
                m.H[0][0] += huu;
                m.H[0][1] += huv;
                m.H[1][0] += huv;
                m.H[1][1] += hvv;
            });
            return m;
        };
        const auto [x, f] = minimize<T, 2>(value, model, {u0, v0}, {box[0], box[2]}, {box[1], box[3]}, tol_x);
        return {x[0], x[1], std::sqrt(T(2) * f)};
    }

    /**
     * @brief Closest points of c1 on range1 and c2 on range2, by Newton from (s0, t0).
     * @return {s, t, distance}
     */
    template <std::floating_point T, std::size_t dim>
    auto newton_curve_curve(const Curve<T, dim> &c1, const Curve<T, dim> &c2, T s0, T t0, const std::array<T, 2> &range1,
                            const std::array<T, 2> &range2, T tol_x) -> std::array<T, 3>
    {
        using namespace projection_detail;
        bool second{true};
        auto value = [&](const vec<T, 2> &x) { return sq_norm(c1.value(x[0]) - c2.value(x[1])) / T(2); };
        auto model = [&](const vec<T, 2> &x) {
            const auto r = c1.value(x[0]) - c2.value(x[1]);
            const auto d1 = c1.value(x[0], 1), d2 = c2.value(x[1], 1);
            Model<T, 2> m;
            m.g = {d1 * r, -(d2 * r)};
            m.G = {{{d1 * d1, -(d1 * d2)}, {-(d1 * d2), d2 * d2}}};
            m.H = m.G;
            second_order(second, [&] {
                const T h1 = c1.value(x[0], 2) * r, h2 = -(c2.value(x[1], 2) * r);
                m.H[0][0] += h1;
                m.H[1][1] += h2;
            });
            return m;
        };
        const auto [x, f] = minimize<T, 2>(value, model, {s0, t0}, {range1[0], range2[0]}, {range1[1], range2[1]}, tol_x);
        return {x[0], x[1], std::sqrt(T(2) * f)};
    }

    /**
     * @brief Closest points of c on range and s on box, by Newton from (t0, u0, v0).
     * @return {t, u, v, distance}
     */
    template <std::floating_point T, std::size_t dim>
    auto newton_curve_surface(const Curve<T, dim> &c, const Surface<T, dim> &s, T t0, T u0, T v0, const std::array<T, 2> &range,
                              const std::array<T, 4> &box, T tol_x) -> std::array<T, 4>
    {
        using namespace projection_detail;
        bool second{true};
        auto value = [&](const vec<T, 3> &x) { return sq_norm(c.value(x[0]) - s.value(x[1], x[2])) / T(2); };
        auto model = [&](const vec<T, 3> &x) {
            const auto r = c.value(x[0]) - s.value(x[1], x[2]);
            const std::array<point<T, dim>, 3> J{c.value(x[0], 1), T(-1) * s.value(x[1], x[2], 1, 0), T(-1) * s.value(x[1], x[2], 0, 1)};
            Model<T, 3> m;
            for (std::size_t i{}; i < 3; ++i)
            {
                m.g[i] = J[i] * r;
                for (std::size_t j{}; j < 3; ++j)
                    m.G[i][j] = J[i] * J[j];
            }
            m.H = m.G;
            second_order(second, [&] {
                const T htt = c.value(x[0], 2) * r, huu = -(s.value(x[1], x[2], 2, 0) * r), huv = -(s.value(x[1], x[2], 1, 1) * r),
                        hvv = -(s.value(x[1], x[2], 0, 2) * r);
                m.H[0][0] += htt;
                m.H[1][1] += huu;
                m.H[1][2] += huv;
                m.H[2][1] += huv;
                m.H[2][2] += hvv;
            });
            return m;
        };
        const auto [x, f] = minimize<T, 3>(value, model, {t0, u0, v0}, {range[0], box[0], box[2]}, {range[1], box[1], box[3]}, tol_x);
        return {x[0], x[1], x[2], std::sqrt(T(2) * f)};
    }

    /**
     * @brief Projection of p on c restricted to range: n samples, Newton from the two nearest.
     * @return {t, distance}
     */
    template <std::floating_point T, std::size_t dim>
    auto project_point_curve(const Curve<T, dim> &c, const point<T, dim> &p, const std::array<T, 2> &range, std::size_t n, T tol_x)
        -> std::array<T, 2>
    {
        using namespace projection_detail;
        std::vector<vec<T, 1>> x;
        for (const auto t : samples(range[0], range[1], n))
            x.push_back({t});
        const auto tb = two_best<T, 1>(x, [&](const vec<T, 1> &t) { return sq_norm(c.value(t[0]) - p); });
        return from_two_best<std::array<T, 2>>(tb, [&](const vec<T, 1> &t) { return newton_point_curve(c, p, t[0], range, tol_x); });
    }

    /**
     * @brief Projection of p on s restricted to box: nu x nv samples, Newton from the two nearest.
     * @return {u, v, distance}
     */
    template <std::floating_point T, std::size_t dim>
    auto project_point_surface(const Surface<T, dim> &s, const point<T, dim> &p, const std::array<T, 4> &box, std::size_t nu,
                               std::size_t nv, T tol_x) -> std::array<T, 3>
    {
        using namespace projection_detail;
        std::vector<vec<T, 2>> x;
        const auto us = samples(box[0], box[1], nu);
        for (const auto v : samples(box[2], box[3], nv))
            for (const auto u : us)
                x.push_back({u, v});
        const auto tb = two_best<T, 2>(x, [&](const vec<T, 2> &uv) { return sq_norm(s.value(uv[0], uv[1]) - p); });
        return from_two_best<std::array<T, 3>>(tb, [&](const vec<T, 2> &uv) { return newton_point_surface(s, p, uv[0], uv[1], box, tol_x); });
    }

    /**
     * @brief Closest points of c1 on range1 and c2 on range2: n1 x n2 samples, Newton from the
     * two nearest pairs.
     * @return {s, t, distance}
     */
    template <std::floating_point T, std::size_t dim>
    auto closest_curve_curve(const Curve<T, dim> &c1, const Curve<T, dim> &c2, const std::array<T, 2> &range1,
                             const std::array<T, 2> &range2, std::size_t n1, std::size_t n2, T tol_x) -> std::array<T, 3>
    {
        using namespace projection_detail;
        const auto s1 = samples(range1[0], range1[1], n1), s2 = samples(range2[0], range2[1], n2);
        points_vector<T, dim> p1(s1.size()), p2(s2.size());
        std::transform(s1.begin(), s1.end(), p1.begin(), [&](T s) { return c1.value(s); });
        std::transform(s2.begin(), s2.end(), p2.begin(), [&](T t) { return c2.value(t); });
        TwoBest<T, 2> tb;
        for (std::size_t i{}; i < s1.size(); ++i)
            for (std::size_t j{}; j < s2.size(); ++j)
                tb.push(sq_norm(p1[i] - p2[j]), {s1[i], s2[j]});
        return from_two_best<std::array<T, 3>>(
            tb, [&](const vec<T, 2> &st) { return newton_curve_curve(c1, c2, st[0], st[1], range1, range2, tol_x); });
    }

    /**
     * @brief Closest points of c on range and s on box: n x nu x nv samples, Newton from the two
     * nearest triples.
     * @return {t, u, v, distance}
     */
    template <std::floating_point T, std::size_t dim>
    auto closest_curve_surface(const Curve<T, dim> &c, const Surface<T, dim> &s, const std::array<T, 2> &range, const std::array<T, 4> &box,
                               std::size_t n, std::size_t nu, std::size_t nv, T tol_x) -> std::array<T, 4>
    {
        using namespace projection_detail;
        const auto ts = samples(range[0], range[1], n), us = samples(box[0], box[1], nu), vs = samples(box[2], box[3], nv);
        points_vector<T, dim> pc(ts.size());
        std::transform(ts.begin(), ts.end(), pc.begin(), [&](T t) { return c.value(t); });
        std::vector<vec<T, 2>> uv;
        for (const auto v : vs)
            for (const auto u : us)
                uv.push_back({u, v});
        points_vector<T, dim> ps(uv.size());
        transform_threshold(uv.begin(), uv.end(), ps.begin(), [&](const vec<T, 2> &x) { return s.value(x[0], x[1]); });
        TwoBest<T, 3> tb;
        for (std::size_t i{}; i < ts.size(); ++i)
            for (std::size_t k{}; k < uv.size(); ++k)
                tb.push(sq_norm(pc[i] - ps[k]), {ts[i], uv[k][0], uv[k][1]});
        return from_two_best<std::array<T, 4>>(
            tb, [&](const vec<T, 3> &x) { return newton_curve_surface(c, s, x[0], x[1], x[2], range, box, tol_x); });
    }
} // namespace gbs
