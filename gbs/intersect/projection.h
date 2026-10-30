#pragma once
/**
 * @file projection.h
 * @brief Closest points by Newton's method: point-curve, point-surface, curve-curve and
 * curve-surface, in any dimension, on bounded parameter ranges.
 *
 * Each problem minimizes f(x) = |r(x)|² / 2 for a residual r (the difference of the two
 * points) of 1 to 3 parameters x kept in their box. A step solves H dx = -g with the full
 * Hessian H = JᵀJ + Σ r_k ∇²r_k, so that the convergence stays quadratic when the distance
 * is not zero; Gauss-Newton (JᵀJ) is used where H is not positive definite (or the geometry
 * has no second derivatives), Levenberg-Marquardt where JᵀJ is singular. Variables on a bound with the gradient pushing out of the box are
 * held for the step, and the step is halved until f does not increase. Stopped on a bound
 * where the opposite bound gives the same distance (a seam), Newton starts again across it.
 *
 * The local functions (`newton_*`) refine from a start; the global ones sample the ranges and
 * run Newton from the two nearest samples (on a closed curve or surface the nearest one may
 * sit at the wrong end of the seam). Design: docs/sources/design/nlopt_removal.md.
 */
#include <gbs/intersect/surface_pieces.h>
#include <gbs/execution.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <numeric>
#include <queue>
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
         * Variables on a bound are held when the gradient or the Newton step points out of the
         * box; when the step cannot move x, a scaled steepest descent on the projected gradient
         * is tried before giving up. Stops when every step |dx_i| <= max(tol_x², 4 eps) (1 + |x_i|), when f cannot decrease
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
                auto on_lo = [&](std::size_t i) { return x[i] <= lo[i]; };
                auto on_hi = [&](std::size_t i) { return x[i] >= hi[i]; };
                // variables on a bound with the gradient pushing out of the box are held
                std::array<bool, n> held{};
                for (std::size_t i{}; i < n; ++i)
                    held[i] = (on_lo(i) && m.g[i] > T(0)) || (on_hi(i) && m.g[i] < T(0));
                const auto pushed_out = held;

                // Newton step on the free variables
                auto newton_step = [&]() {
                    vec<T, n> b{}, dx{};
                    for (std::size_t i{}; i < n; ++i)
                        b[i] = held[i] ? T(0) : -m.g[i];
                    auto reduced = [&](mat<T, n> A) {
                        for (std::size_t i{}; i < n; ++i)
                            if (held[i])
                                for (std::size_t j{}; j < n; ++j)
                                    A[i][j] = A[j][i] = T(i == j);
                        return A;
                    };
                    if (cholesky_solve(reduced(m.H), b, dx) || cholesky_solve(reduced(m.G), b, dx))
                        return dx;
                    // singular JᵀJ (parallel tangents): Levenberg-Marquardt, JᵀJ + mu diag(JᵀJ)
                    auto lm = m.G;
                    for (std::size_t i{}; i < n; ++i)
                        lm[i][i] += T(1e-6) * (m.G[i][i] > T(0) ? m.G[i][i] : T(1));
                    if (cholesky_solve(reduced(lm), b, dx))
                        return dx;
                    for (std::size_t i{}; i < n; ++i) // null derivatives
                        dx[i] = !held[i] && m.G[i][i] > T(0) ? b[i] / m.G[i][i] : T(0);
                    return dx;
                };
                // a free variable on a bound whose step points out of the box (through the coupling
                // of H) would be cancelled by the clamp: it is held too and the step recomputed
                vec<T, n> dx{};
                bool any_free{};
                for (std::size_t pass{}; pass <= n; ++pass)
                {
                    any_free = std::find(held.begin(), held.end(), false) != held.end();
                    if (!any_free)
                        break;
                    dx = newton_step();
                    bool again{};
                    for (std::size_t i{}; i < n; ++i)
                        if (!held[i] && ((on_lo(i) && dx[i] < T(0)) || (on_hi(i) && dx[i] > T(0))))
                            held[i] = again = true;
                    if (!again)
                        break;
                }

                auto small_step = [&](const vec<T, n> &y) {
                    for (std::size_t i{}; i < n; ++i)
                        if (std::abs(y[i] - x[i]) > tau * (T(1) + std::abs(x[i])))
                            return false;
                    return true;
                };
                // halve the step until f decreases (or stays equal on a small step, at the machine
                // precision; a large step at the same f would swap between equal minima forever);
                // false if x cannot move
                vec<T, n> x_new{};
                T f_new{};
                auto line_search = [&](const vec<T, n> &d) {
                    T alpha{1};
                    for (int ls{}; ls < 30; ++ls, alpha /= T(2))
                    {
                        for (std::size_t i{}; i < n; ++i)
                            x_new[i] = std::clamp(x[i] + alpha * d[i], lo[i], hi[i]);
                        f_new = value(x_new);
                        if (f_new < f || (f_new <= f && small_step(x_new)))
                            return x_new != x;
                    }
                    return false;
                };
                bool moved = any_free && line_search(dx);
                if (!moved)
                {
                    // no progress while the projected gradient is not null (variables held by the
                    // clamp with the gradient pointing into the box): scaled steepest descent
                    vec<T, n> dg{};
                    bool descent{};
                    for (std::size_t i{}; i < n; ++i)
                    {
                        dg[i] = !pushed_out[i] && m.G[i][i] > T(0) ? -m.g[i] / m.G[i][i] : T(0);
                        descent = descent || dg[i] != T(0);
                    }
                    moved = descent && line_search(dg);
                }
                if (!moved)
                    break; // minimum reached at the machine precision
                const bool small = small_step(x_new);
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

        /// The two samples nearest to the target, as starts for Newton. Samples at the same
        /// distance as the nearest one are taken for the same point and not kept as the second
        /// start: the samples of a degenerate edge (the pole of a sphere) are all the same point,
        /// where Newton stays (null derivative); the copies of a point across a seam are left to
        /// the restart of `minimize`.
        template <std::floating_point T, std::size_t n>
        struct TwoBest
        {
            std::array<std::pair<T, vec<T, n>>, 2> best{std::pair{std::numeric_limits<T>::max(), vec<T, n>{}},
                                                        std::pair{std::numeric_limits<T>::max(), vec<T, n>{}}};
            static auto same(T a, T b) -> bool { return std::abs(a - b) <= T(64) * std::numeric_limits<T>::epsilon() * std::max(a, b); }
            auto push(T d, const vec<T, n> &x) -> void
            {
                if (d < best[0].first)
                {
                    if (!same(d, best[0].first))
                        best[1] = best[0];
                    best[0] = {d, x};
                }
                else if (d < best[1].first && !same(d, best[0].first))
                    best[1] = {d, x};
            }
        };

        /// The two nearest of several partial results, in their order (deterministic).
        template <std::floating_point T, std::size_t n>
        auto merge(const std::vector<TwoBest<T, n>> &parts) -> TwoBest<T, n>
        {
            TwoBest<T, n> tb;
            for (const auto &p : parts)
                for (const auto &[d, x] : p.best)
                    if (d < std::numeric_limits<T>::max())
                        tb.push(d, x);
            return tb;
        }

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

        /**
         * Best-first branch and bound on pieces of the geometries. `lower(node)` bounds the
         * distance from below on the node. `leaf(node)` lowers `best` with the exact points it
         * knows (ends, corners) and, when the node is flat (or cannot be halved any more), with
         * Newton from it, and then returns true: the node is done. Otherwise `split(node)` returns
         * its two halves. A node is dropped when its lower bound cannot beat `best` beyond the
         * rounding errors. Returns false when `max_nodes` were examined before the end (`best` is
         * then only the best found).
         *
         * Certifying the minimum to a margin delta instead (halving the flat nodes on until their
         * lower bound reaches best - delta) was measured: near a minimum at a non-zero distance the
         * lower bounds converge like the flatness, so the number of pieces of surface grows like
         * 1/delta (53 000 for 1e-6 on a sphere), and near-ties (a segment almost parallel to a
         * circle) exhaust any budget. Newton from a flat piece finds its minimum unless the piece,
         * flat within `flat_ratio`, holds two distinct minima.
         */
        template <std::floating_point T, typename Node, typename FLower, typename FLeaf, typename FSplit>
        auto branch_and_bound(std::vector<Node> nodes, const FLower &lower, const FLeaf &leaf, const FSplit &split, const T &best,
                              std::size_t max_nodes = 100000) -> bool
        {
            auto worth = [&](T lb) { return lb < best - T(16) * std::numeric_limits<T>::epsilon() * (T(1) + best); };
            using Item = std::pair<T, std::size_t>;
            std::priority_queue<Item, std::vector<Item>, std::greater<>> queue;
            for (std::size_t i{}; i < nodes.size(); ++i)
                queue.push({lower(nodes[i]), i});
            std::size_t examined{};
            while (!queue.empty())
            {
                const auto [lb, i] = queue.top();
                queue.pop();
                if (!worth(lb))
                    break; // the others are farther
                if (++examined > max_nodes)
                    return false;
                Node node = std::move(nodes[i]);
                if (leaf(node) || !worth(lb))
                    continue;
                auto halves = split(node);
                for (auto *h : {&halves.first, &halves.second})
                    if (const T l = lower(*h); worth(l))
                    {
                        nodes.push_back(std::move(*h));
                        queue.push({l, nodes.size() - 1});
                    }
            }
            return true;
        }

        /// A piece flat enough to start Newton from its chord or its centre, after which it is done
        /// (the pieces are halved until they are). 0.05 and 0.2 gave the same results on the
        /// bench (bench/bench_extrema.cpp); 0.1 halves the time of curve-surface against 0.05.
        inline constexpr double flat_ratio = 0.1;

        template <std::floating_point T, std::size_t dim>
        auto flat(const CurvePiece<T, dim> &pc, T small) -> bool
        {
            return pc.flatness <= T(flat_ratio) * pc.chord() || pc.size() <= small;
        }

        template <std::floating_point T>
        auto flat(const SurfacePiece<T> &pc, T small) -> bool
        {
            return pc.flatness <= T(flat_ratio) * pc.size() || pc.size() <= small;
        }

        // Lower bounds of the distance. A piece of curve lies within `flatness` of its chord, a
        // patch within `flatness` of the plane of its corners (its slab): these bounds are much
        // tighter than the boxes once the pieces are flat, and are combined with them.

        /// Signed distance of x to the plane of a patch (0 when its normal is undefined).
        template <std::floating_point T>
        auto to_plane(const SurfacePiece<T> &q, const point<T, 3> &x) -> T
        {
            return (x - q.centre) * q.normal;
        }

        template <std::floating_point T, std::size_t dim>
        auto lower_bound(const CurvePiece<T, dim> &pc, const point<T, dim> &p) -> T
        {
            return std::max(pc.box.distance(p), intersect_detail::distance_to_segment(p, pc.first, pc.last) - pc.flatness);
        }

        template <std::floating_point T>
        auto lower_bound(const SurfacePiece<T> &q, const point<T, 3> &p) -> T
        {
            const T slab = sq_norm(q.normal) > T(0) ? std::abs(to_plane(q, p)) : norm(p - q.centre);
            return std::max(q.box.distance(p), slab - q.flatness);
        }

        template <std::floating_point T, std::size_t dim>
        auto lower_bound(const CurvePiece<T, dim> &a, const CurvePiece<T, dim> &b) -> T
        {
            const T chords = intersect_detail::closest_on_chords(a.first, a.last, b.first, b.last)[2];
            return std::max(a.box.distance(b.box), chords - a.flatness - b.flatness);
        }

        template <std::floating_point T>
        auto lower_bound(const CurvePiece<T, 3> &a, const SurfacePiece<T> &q) -> T
        {
            T slab{};
            if (sq_norm(q.normal) > T(0))
            {
                // the chord on one side of the slab: its nearest end
                const T d0 = to_plane(q, a.first), d1 = to_plane(q, a.last);
                if ((d0 > T(0)) == (d1 > T(0)))
                    slab = std::min(std::abs(d0), std::abs(d1)) - a.flatness - q.flatness;
            }
            return std::max(a.box.distance(q.box), slab);
        }

        /// Parameter of the point of the chord of a piece nearest to p.
        template <std::floating_point T, std::size_t dim>
        auto on_chord(const CurvePiece<T, dim> &pc, const point<T, dim> &p) -> T
        {
            const auto ch = pc.last - pc.first;
            const T l2 = sq_norm(ch);
            const T s = l2 > T(0) ? std::clamp((p - pc.first) * ch / l2, T(0), T(1)) : T(0.5);
            return pc.a + s * (pc.b - pc.a);
        }

        /// Parameters of the point of a patch nearest to p, estimated on the bilinear patch of its
        /// corners (Gauss-Newton, no evaluation of the surface): with the centre, a start of
        /// Newton on a flat patch. Each of them alone fails on a degenerate parametrization, where
        /// Newton stops on a stationary point: from the centre, on a line where Su = 0 (two equal
        /// poles) between the centre and the foot; from the bilinear estimate, at a pole whose
        /// corners all coincide.
        template <std::floating_point T>
        auto on_patch(const SurfacePiece<T> &q, const point<T, 3> &p) -> std::array<T, 2>
        {
            const auto p00 = q.pole(0, 0), p10 = q.pole(q.nu - 1, 0), p01 = q.pole(0, q.nv - 1), p11 = q.pole(q.nu - 1, q.nv - 1);
            T s{T(0.5)}, t{T(0.5)};
            for (int it{}; it < 10; ++it)
            {
                const auto r = (T(1) - s) * (T(1) - t) * p00 + s * (T(1) - t) * p10 + (T(1) - s) * t * p01 + s * t * p11 - p;
                const auto bs = (T(1) - t) * (p10 - p00) + t * (p11 - p01), bt = (T(1) - s) * (p01 - p00) + s * (p11 - p10);
                const T a = bs * bs, b = bs * bt, c = bt * bt, g1 = bs * r, g2 = bt * r, det = a * c - b * b;
                if (!(det > T(0)))
                    break;
                s = std::clamp(s - (c * g1 - b * g2) / det, T(0), T(1));
                t = std::clamp(t - (a * g2 - b * g1) / det, T(0), T(1));
            }
            return {q.uv[0] + s * (q.uv[1] - q.uv[0]), q.uv[2] + t * (q.uv[3] - q.uv[2])};
        }

        /// Size under which a piece is not halved any more: a fraction of the size of the geometry.
        template <std::floating_point T, typename Pieces>
        auto small_size(const Pieces &pieces) -> T
        {
            T size{};
            for (const auto &pc : pieces)
                size = std::max(size, pc.size());
            return T(1e-7) * size;
        }

        /// Keeps the closest result (the last entry is the distance).
        template <typename R>
        auto keep(R &best, const R &r) -> void
        {
            if (r.back() < best.back())
                best = r;
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
        const bool second = c.has_second_derivatives();
        auto value = [&](const vec<T, 1> &x) { return sq_norm(c.value(x[0]) - p) / T(2); };
        auto model = [&](const vec<T, 1> &x) {
            const auto r = c.value(x[0]) - p;
            const auto d1 = c.value(x[0], 1);
            Model<T, 1> m;
            m.g = {d1 * r};
            m.G = {{{d1 * d1}}};
            m.H = m.G;
            if (second)
                m.H[0][0] += c.value(x[0], 2) * r;
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
        const bool second = s.has_second_derivatives();
        auto value = [&](const vec<T, 2> &x) { return sq_norm(s.value(x[0], x[1]) - p) / T(2); };
        auto model = [&](const vec<T, 2> &x) {
            const auto r = s.value(x[0], x[1]) - p;
            const auto su = s.value(x[0], x[1], 1, 0), sv = s.value(x[0], x[1], 0, 1);
            Model<T, 2> m;
            m.g = {su * r, sv * r};
            m.G = {{{su * su, su * sv}, {su * sv, sv * sv}}};
            m.H = m.G;
            if (second)
            {
                const T huu = s.value(x[0], x[1], 2, 0) * r, huv = s.value(x[0], x[1], 1, 1) * r, hvv = s.value(x[0], x[1], 0, 2) * r;
                m.H[0][0] += huu;
                m.H[0][1] += huv;
                m.H[1][0] += huv;
                m.H[1][1] += hvv;
            }
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
        const bool second = c1.has_second_derivatives() && c2.has_second_derivatives();
        auto value = [&](const vec<T, 2> &x) { return sq_norm(c1.value(x[0]) - c2.value(x[1])) / T(2); };
        auto model = [&](const vec<T, 2> &x) {
            const auto r = c1.value(x[0]) - c2.value(x[1]);
            const auto d1 = c1.value(x[0], 1), d2 = c2.value(x[1], 1);
            Model<T, 2> m;
            m.g = {d1 * r, -(d2 * r)};
            m.G = {{{d1 * d1, -(d1 * d2)}, {-(d1 * d2), d2 * d2}}};
            m.H = m.G;
            if (second)
            {
                const T h1 = c1.value(x[0], 2) * r, h2 = -(c2.value(x[1], 2) * r);
                m.H[0][0] += h1;
                m.H[1][1] += h2;
            }
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
        const bool second = c.has_second_derivatives() && s.has_second_derivatives();
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
            if (second)
            {
                const T htt = c.value(x[0], 2) * r, huu = -(s.value(x[1], x[2], 2, 0) * r), huv = -(s.value(x[1], x[2], 1, 1) * r),
                        hvv = -(s.value(x[1], x[2], 0, 2) * r);
                m.H[0][0] += htt;
                m.H[1][1] += huu;
                m.H[1][2] += huv;
                m.H[2][1] += huv;
                m.H[2][2] += hvv;
            }
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
        // the two nearest pairs of each row (in parallel above gbs::parallel_min_size), then of all rows
        std::vector<std::size_t> rows(s1.size());
        std::iota(rows.begin(), rows.end(), std::size_t{});
        std::vector<TwoBest<T, 2>> row_best(rows.size());
        transform_threshold(rows.begin(), rows.end(), row_best.begin(), [&](std::size_t i) {
            TwoBest<T, 2> tb;
            for (std::size_t j{}; j < s2.size(); ++j)
                tb.push(sq_norm(p1[i] - p2[j]), {s1[i], s2[j]});
            return tb;
        });
        const auto tb = merge(row_best);
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
        // the two nearest triples of each sample of the surface (in parallel above
        // gbs::parallel_min_size), then of all of them
        std::vector<std::size_t> rows(uv.size());
        std::iota(rows.begin(), rows.end(), std::size_t{});
        std::vector<TwoBest<T, 3>> row_best(rows.size());
        transform_threshold(rows.begin(), rows.end(), row_best.begin(), [&](std::size_t k) {
            TwoBest<T, 3> tb;
            for (std::size_t i{}; i < ts.size(); ++i)
                tb.push(sq_norm(pc[i] - ps[k]), {ts[i], uv[k][0], uv[k][1]});
            return tb;
        });
        const auto tb = merge(row_best);
        return from_two_best<std::array<T, 4>>(
            tb, [&](const vec<T, 3> &x) { return newton_curve_surface(c, s, x[0], x[1], x[2], range, box, tol_x); });
    }
    // Global closest points by subdivision. The geometries are cut into pieces enclosed by boxes
    // (gbs/intersect/curve_pieces.h, surface_pieces.h); the pieces are examined by increasing
    // distance of their boxes, those that cannot be closer than the best point found are
    // dropped, the others are halved until flat, and Newton starts from each flat piece. For a
    // B-spline (rational or not) the box of the poles of a Bézier piece encloses it: the result
    // is the global minimum. Other geometries are sampled, their boxes are not guaranteed.
    // Surfaces are cut in 3D only: in another dimension the sampled versions above are used.

    /**
     * @brief Global projection of p on c restricted to range, by subdivision.
     * @return {t, distance}
     */
    template <std::floating_point T, std::size_t dim>
    auto global_point_curve(const Curve<T, dim> &c, const point<T, dim> &p, const std::array<T, 2> &range, T tol_x) -> std::array<T, 2>
    {
        using namespace projection_detail;
        CurvePieces<T, dim> pieces{c};
        auto roots = pieces.pieces(range[0], range[1]);
        const T small = small_size<T>(roots);
        std::array<T, 2> best{range[0], norm(c.value(range[0]) - p)};
        keep(best, {range[1], norm(c.value(range[1]) - p)});
        branch_and_bound<T>(
            std::move(roots), [&](const CurvePiece<T, dim> &pc) { return lower_bound(pc, p); },
            [&](const CurvePiece<T, dim> &pc) {
                keep(best, {pc.a, norm(pc.first - p)});
                keep(best, {pc.b, norm(pc.last - p)});
                if (flat(pc, small))
                {
                    keep(best, newton_point_curve(c, p, on_chord(pc, p), range, tol_x));
                    return true;
                }
                return false;
            },
            [&](const CurvePiece<T, dim> &pc) { return pieces.split(pc); }, best.back());
        return best;
    }

    /**
     * @brief Global projection of p on s restricted to box, by subdivision (3D; sampled with
     * nu x nv samples in another dimension).
     * @return {u, v, distance}
     */
    template <std::floating_point T, std::size_t dim>
    auto global_point_surface(const Surface<T, dim> &s, const point<T, dim> &p, const std::array<T, 4> &box, T tol_x,
                              std::size_t nu = 30, std::size_t nv = 30) -> std::array<T, 3>
    {
        using namespace projection_detail;
        if constexpr (dim != 3)
            return project_point_surface(s, p, box, nu, nv, tol_x);
        else
        {
            SurfacePieces<T> pieces{s};
            auto roots = pieces.pieces(box);
            const T small = small_size<T>(roots);
            std::array<T, 3> best{box[0], box[2], norm(s.value(box[0], box[2]) - p)};
            auto corners = [&](const SurfacePiece<T> &q) {
                for (const auto &[i, j] : {std::pair{0, 0}, std::pair{1, 0}, std::pair{0, 1}, std::pair{1, 1}})
                    keep(best, {q.uv[i], q.uv[2 + j], norm(q.pole(i * (q.nu - 1), j * (q.nv - 1)) - p)});
            };
            branch_and_bound<T>(
                std::move(roots), [&](const SurfacePiece<T> &q) { return lower_bound(q, p); },
                [&](const SurfacePiece<T> &q) {
                    corners(q);
                    if (flat(q, small))
                    {
                        for (const auto &uv : {q.center(), on_patch(q, p)})
                            keep(best, newton_point_surface(s, p, uv[0], uv[1], box, tol_x));
                        return true;
                    }
                    return false;
                },
                [&](const SurfacePiece<T> &q) { return pieces.split(q, SurfacePieces<T>::longer_along_u(q)); }, best.back());
            return best;
        }
    }

    /**
     * @brief Global closest points of c1 on range1 and c2 on range2, by subdivision.
     * @return {s, t, distance}
     */
    template <std::floating_point T, std::size_t dim>
    auto global_curve_curve(const Curve<T, dim> &c1, const Curve<T, dim> &c2, const std::array<T, 2> &range1, const std::array<T, 2> &range2,
                            T tol_x) -> std::array<T, 3>
    {
        using namespace projection_detail;
        using Piece = CurvePiece<T, dim>;
        using Pair = std::pair<Piece, Piece>;
        CurvePieces<T, dim> P1{c1}, P2{c2};
        const auto r1 = P1.pieces(range1[0], range1[1]), r2 = P2.pieces(range2[0], range2[1]);
        const T small1 = small_size<T>(r1), small2 = small_size<T>(r2);
        std::vector<Pair> roots;
        for (const auto &a : r1)
            for (const auto &b : r2)
                roots.emplace_back(a, b);
        std::array<T, 3> best{range1[0], range2[0], norm(c1.value(range1[0]) - c2.value(range2[0]))};
        branch_and_bound<T>(
            std::move(roots), [&](const Pair &pq) { return lower_bound(pq.first, pq.second); },
            [&](const Pair &pq) {
                const auto &[a, b] = pq;
                for (const auto &[s, x] : {std::pair{a.a, a.first}, std::pair{a.b, a.last}})
                    for (const auto &[t, y] : {std::pair{b.a, b.first}, std::pair{b.b, b.last}})
                        keep(best, {s, t, norm(x - y)});
                if (flat(a, small1) && flat(b, small2) || (a.size() <= small1 && b.size() <= small2))
                {
                    const auto [ls, lt, ld] = intersect_detail::closest_on_chords(a.first, a.last, b.first, b.last);
                    keep(best, newton_curve_curve(c1, c2, a.a + ls * (a.b - a.a), b.a + lt * (b.b - b.a), range1, range2, tol_x));
                    return true;
                }
                return false;
            },
            [&](const Pair &pq) {
                const auto &[a, b] = pq;
                // halve the piece that is not flat, the larger one if both are not
                if (!flat(a, small1) && (flat(b, small2) || a.size() >= b.size()))
                {
                    auto [l, r] = P1.split(a);
                    return std::pair{Pair{std::move(l), b}, Pair{std::move(r), b}};
                }
                auto [l, r] = P2.split(b);
                return std::pair{Pair{a, std::move(l)}, Pair{a, std::move(r)}};
            },
            best.back());
        return best;
    }

    /**
     * @brief Global closest points of c on range and s on box, by subdivision (3D; sampled with
     * n x nu x nv samples in another dimension).
     * @return {t, u, v, distance}
     */
    template <std::floating_point T, std::size_t dim>
    auto global_curve_surface(const Curve<T, dim> &c, const Surface<T, dim> &s, const std::array<T, 2> &range, const std::array<T, 4> &box,
                              T tol_x, std::size_t n = 30, std::size_t nu = 30, std::size_t nv = 30) -> std::array<T, 4>
    {
        using namespace projection_detail;
        if constexpr (dim != 3)
            return closest_curve_surface(c, s, range, box, n, nu, nv, tol_x);
        else
        {
            using CPiece = CurvePiece<T, 3>;
            using SPiece = SurfacePiece<T>;
            using Pair = std::pair<CPiece, SPiece>;
            CurvePieces<T, 3> CP{c};
            SurfacePieces<T> SP{s};
            const auto rc = CP.pieces(range[0], range[1]);
            const auto rs = SP.pieces(box);
            const T small_c = small_size<T>(rc), small_s = small_size<T>(rs);
            std::vector<Pair> roots;
            for (const auto &a : rc)
                for (const auto &b : rs)
                    roots.emplace_back(a, b);
            std::array<T, 4> best{range[0], box[0], box[2], norm(c.value(range[0]) - s.value(box[0], box[2]))};
            branch_and_bound<T>(
                std::move(roots), [&](const Pair &pq) { return lower_bound(pq.first, pq.second); },
                [&](const Pair &pq) {
                    const auto &[a, q] = pq;
                    for (const auto &[t, x] : {std::pair{a.a, a.first}, std::pair{a.b, a.last}})
                        for (const auto &[i, j] : {std::pair{0, 0}, std::pair{1, 0}, std::pair{0, 1}, std::pair{1, 1}})
                            keep(best, {t, q.uv[i], q.uv[2 + j], norm(x - q.pole(i * (q.nu - 1), j * (q.nv - 1)))});
                    if (flat(a, small_c) && flat(q, small_s) || (a.size() <= small_c && q.size() <= small_s))
                    {
                        const auto t0 = on_chord(a, s.value(q.center()[0], q.center()[1]));
                        for (const auto &uv : {q.center(), on_patch(q, c.value(t0))})
                            keep(best, newton_curve_surface(c, s, t0, uv[0], uv[1], range, box, tol_x));
                        return true;
                    }
                    return false;
                },
                [&](const Pair &pq) {
                    const auto &[a, q] = pq;
                    if (!flat(a, small_c) && (flat(q, small_s) || a.size() >= q.size()))
                    {
                        auto [l, r] = CP.split(a);
                        return std::pair{Pair{std::move(l), q}, Pair{std::move(r), q}};
                    }
                    auto [l, r] = SP.split(q, SurfacePieces<T>::longer_along_u(q));
                    return std::pair{Pair{a, std::move(l)}, Pair{a, std::move(r)}};
                },
                best.back());
            return best;
        }
    }
} // namespace gbs
