#pragma once
/**
 * @file curves.h
 * @brief Intersection of two curves: every crossing and every common portion on bounded
 * parameter ranges.
 *
 * The curves are cut into pieces (`CurvePieces`): pairs of pieces whose boxes, or fat lines, are
 * farther than the tolerance are dropped, the others are halved until both pieces are flat. Two
 * flat pieces give a crossing candidate from the closest points of their chords, refined by
 * Newton on the whole curves. A piece lying within the tolerance of the other curve, and parallel
 * to it, seeds a common portion and is not subdivided further; the seeds that follow each other
 * are grown to where the curves part, then told apart: an overlap, or a tangency. Written for any
 * dimension: in 2D it crosses pcurves in the parameter space of a face, in 3D it finds the points
 * where two curves meet within the tolerance.
 */
#include <gbs/intersect/curve_pieces.h>
#include <gbs/intersect/projection.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs
{
    /// Options of the intersection algorithms.
    template <std::floating_point T>
    struct IntersectOptions
    {
        T tol{extrema_tol<T>};         ///< distance under which two points are the same
        T angular_tol{T(1e-3)};        ///< sine of the angle under which a crossing is tangent
        std::size_t max_pairs{200000}; ///< cap on the pairs of pieces examined
        bool closed_forms{true};       ///< surface-surface: use the closed forms of the pairs of primitives
        T max_turn{T(0.05)};           ///< surface-surface marching: largest turn of the tangent per step (radians)
    };

    /// A crossing of two curves: c1(t1) = c2(t2) = p within the tolerance.
    template <std::floating_point T, std::size_t dim>
    struct CurveCurveHit
    {
        T t1{}, t2{};
        point<T, dim> p{};
        bool tangent{}; ///< the curves are tangent at p
    };

    /// A common portion of two curves: c1 on [t1a, t1b] (t1a < t1b) lies on c2 from t2a to t2b;
    /// t2a > t2b when the curves run in opposite directions.
    template <std::floating_point T>
    struct CurveOverlap
    {
        T t1a{}, t1b{}, t2a{}, t2b{};
    };

    template <std::floating_point T, std::size_t dim>
    struct CurveCurveResult
    {
        std::vector<CurveCurveHit<T, dim>> hits; ///< sorted by t1, none inside an overlap
        std::vector<CurveOverlap<T>> overlaps;   ///< sorted by t1a
        bool complete{true};                     ///< false when `max_pairs` was reached: the result is partial
    };

    namespace intersect_detail
    {
        /**
         * Gauss-Newton on |c1(s) - c2(t)|² from (s, t), parameters kept in their ranges, each step
         * halved until the distance decreases: a local minimum of the distance (a crossing, the
         * contact of a tangency or of a near miss), not the first point found within the
         * tolerance. The 2 x 2 system is singular relative to its own diagonal only (the speeds of
         * the curves may differ by orders of magnitude, a grazing crossing is well posed); with
         * parallel tangents, or when its step does not decrease the distance (against a bound),
         * each curve is projected on the other. Stopped at a saddle (the top of the hump between two
         * crossings), it goes on along the direction of negative curvature, towards `side`, and sets
         * `escaped`. Returns the parameters and the distance.
         */
        template <std::floating_point T, std::size_t dim>
        auto refine(const Curve<T, dim> &c1, T a1, T b1, const Curve<T, dim> &c2, T a2, T b2, T s, T t, T side = T(1), bool *escaped = nullptr)
            -> std::array<T, 3>
        {
            constexpr T eps = std::numeric_limits<T>::epsilon();
            T d = norm(c1.value(s) - c2.value(t));
            bool at_rounding = false;
            // halves the step (ds, dt) until the distance decreases; false if it does not
            auto step = [&](T ds, T dt) {
                for (int k{}; k < 30; ++k, ds /= T(2), dt /= T(2))
                {
                    const T s_n = std::clamp(s + ds, a1, b1), t_n = std::clamp(t + dt, a2, b2);
                    if (s_n == s && t_n == t)
                        return false;
                    if (const T d_n = norm(c1.value(s_n) - c2.value(t_n)); d_n < d)
                    {
                        at_rounding = std::abs(s_n - s) <= T(4) * eps * (T(1) + std::abs(s)) && std::abs(t_n - t) <= T(4) * eps * (T(1) + std::abs(t));
                        s = s_n;
                        t = t_n;
                        d = d_n;
                        return true;
                    }
                }
                return false;
            };
            for (int escape{}; escape < 3; ++escape)
            {
                for (int it{}; it < 100 && d > T(0) && !at_rounding; ++it)
                {
                    const auto r = c1.value(s) - c2.value(t);
                    const auto d1 = c1.value(s, 1), d2 = c2.value(t, 1);
                    const auto A = d1 * d1, B = -(d1 * d2), C = d2 * d2;
                    const auto g1 = d1 * r, g2 = -(d2 * r);
                    const auto det = A * C - B * B;
                    const bool regular = det > T(1e-14) * A * C && det > T(0);
                    if (!(regular && step(-(C * g1 - B * g2) / det, -(A * g2 - B * g1) / det)) &&
                        !step(A > T(0) ? -g1 / A : T(0), C > T(0) ? -g2 / C : T(0)))
                        break;
                }
                // stopped above the rounding of the points: a saddle when the Hessian of |r|² / 2 has a
                // negative eigenvalue (the top of the hump between two crossings, a line grazing a
                // circle); along its eigenvector, the quadratic model vanishes at d / sqrt(-lambda)
                const auto r = c1.value(s) - c2.value(t);
                if (d <= T(64) * eps * (T(1) + norm(c1.value(s))))
                    break;
                const auto d1 = c1.value(s, 1), d2 = c2.value(t, 1);
                const T h11 = d1 * d1 + r * c1.value(s, 2), h22 = d2 * d2 - r * c2.value(t, 2), h12 = -(d1 * d2);
                const T mid = (h11 + h22) / T(2), rad = std::hypot((h11 - h22) / T(2), h12);
                const T lambda = mid > T(0) ? (h11 * h22 - h12 * h12) / (mid + rad) : mid - rad; // the smaller eigenvalue, without cancellation
                if (!(lambda < T(0)))
                    break;
                std::array<T, 2> v = std::abs(h11 - lambda) >= std::abs(h22 - lambda) ? std::array{-h12, h11 - lambda} : std::array{h22 - lambda, -h12};
                const T nv = std::hypot(v[0], v[1]);
                if (!(nv > T(0)))
                    break;
                const T a = side * d / std::sqrt(-lambda) / nv;
                at_rounding = false;
                if (!step(a * v[0], a * v[1]) && !step(-a * v[0], -a * v[1]))
                    break;
                if (escaped)
                    *escaped = true;
            }
            return {s, t, d};
        }

        /**
         * Projection of `p` on `c` restricted to [a, b], by Newton from the two nearest of a few
         * samples (on a closed curve the nearest may sit at the wrong end, against its bound).
         */
        template <std::floating_point T, std::size_t dim>
        auto project(const Curve<T, dim> &c, T a, T b, const point<T, dim> &p) -> std::pair<T, T>
        {
            const auto [t, d] = project_point_curve(c, p, {a, b}, 17, T(0));
            return {t, d};
        }

        /// Foot of p on c on [a, b] near the parameter t: Gauss-Newton from t, a few steps (p near
        /// c(t), or c nearly straight). Returns the parameter and the distance.
        template <std::floating_point T, std::size_t dim>
        auto local_foot(const Curve<T, dim> &c, const point<T, dim> &p, T t, T a, T b) -> std::pair<T, T>
        {
            for (int it{}; it < 6; ++it)
            {
                const auto d1 = c.value(t, 1);
                const T l2 = d1 * d1;
                if (!(l2 > T(0)))
                    break;
                const T t_n = std::clamp(t - d1 * (c.value(t) - p) / l2, a, b);
                if (std::abs(t_n - t) <= T(4) * std::numeric_limits<T>::epsilon() * (T(1) + std::abs(t)))
                    break;
                t = t_n;
            }
            return {t, norm(c.value(t) - p)};
        }

        /// Sine of the angle between two vectors (0 if one is null).
        template <std::floating_point T, std::size_t dim>
        auto sine(const point<T, dim> &u, const point<T, dim> &v) -> T
        {
            const auto nu = norm(u), nv = norm(v);
            if (nu <= T(0) || nv <= T(0))
                return T(0);
            const auto c = std::clamp(u * v / (nu * nv), T(-1), T(1));
            return std::sqrt(std::max(T(0), T(1) - c * c));
        }

        template <std::floating_point T>
        struct Seed
        {
            T s0, s1, t0, t1;
        };

        /// Extent of a piece along the unit vector n: between its poles (convex hull of a Bézier
        /// piece), or within its flatness of its chord (sampled piece).
        template <std::floating_point T, std::size_t dim>
        auto extent(const CurvePiece<T, dim> &pc, const point<T, dim> &n) -> std::array<T, 2>
        {
            if (!pc.bezier)
            {
                const T u = pc.first * n, v = pc.last * n;
                return {std::min(u, v) - pc.flatness, std::max(u, v) + pc.flatness};
            }
            std::array<T, 2> r{std::numeric_limits<T>::max(), std::numeric_limits<T>::lowest()};
            for (const auto &h : pc.poles)
            {
                const T x = cartesian<T, dim>(h) * n;
                r = {std::min(r[0], x), std::max(r[1], x)};
            }
            return r;
        }

        /// Fat line: is q farther than `gap` from p across the chord of p? Along the normal to the
        /// chord of p that points to the middle of the chord of q (in 2D, the normal of the chord).
        template <std::floating_point T, std::size_t dim>
        auto apart(const CurvePiece<T, dim> &p, const CurvePiece<T, dim> &q, T gap) -> bool
        {
            const auto c = p.last - p.first;
            const T l2 = sq_norm(c);
            auto n = (q.first + q.last) / T(2) - p.first;
            if (l2 > T(0))
                n = n - (n * c / l2) * c;
            const T ln = norm(n);
            if (!(ln > T(0)))
                return false;
            n = n / ln;
            return extent(q, n)[0] > extent(p, n)[1] + gap;
        }

        /// Parameter at the fraction f of the way from x to y on [a, b]: round the seam of a closed
        /// curve when that way is shorter.
        template <std::floating_point T>
        auto between(T x, T y, T a, T b, bool closed, T f) -> T
        {
            T d = y - x;
            if (closed && std::abs(d) > (b - a) / T(2))
                d -= std::copysign(b - a, d);
            T u = x + f * d;
            if (u < a)
                u += b - a;
            else if (u > b)
                u -= b - a;
            return std::clamp(u, a, b);
        }

        /// Projection of p on c on the way from x to y (see `between`): parameter and distance.
        template <std::floating_point T, std::size_t dim>
        auto project_on_way(const Curve<T, dim> &c, T x, T y, T a, T b, bool closed, const point<T, dim> &p) -> std::pair<T, T>
        {
            const T lo = std::min(x, y), hi = std::max(x, y);
            if (closed && hi - lo > (b - a) / T(2))
                return std::min(project(c, hi, b, p), project(c, a, lo, p), [](const auto &l, const auto &r) { return l.second < r.second; });
            return lo < hi ? project(c, lo, hi, p) : std::pair{lo, norm(c.value(lo) - p)};
        }
    }

    /**
     * @brief Every crossing and common portion of c1 on [a1, b1] and c2 on [a2, b2].
     *
     * A crossing is reported once, at the parameters where the curves are the closest (both
     * within `tol`); the crossings at the ends of the ranges are included. Crossings between
     * which the curves stay within `tol` of each other are one contact, reported once and tangent
     * (a line grazing a circle within the tolerance). A common portion is an overlap when the
     * curves stay together over more than ten times the tolerance besides the parts where they
     * leave each other tangentially, and the crossings on it are dropped. `complete` is false when
     * `max_pairs` pairs of pieces were not enough: the result is partial.
     */
    template <std::floating_point T, std::size_t dim>
    [[nodiscard]] auto intersect(const Curve<T, dim> &c1, std::type_identity_t<T> a1, std::type_identity_t<T> b1,
                                 const Curve<T, dim> &c2, std::type_identity_t<T> a2, std::type_identity_t<T> b2,
                                 IntersectOptions<T> o = {}) -> CurveCurveResult<T, dim>
    {
        using namespace intersect_detail;
        using Piece = CurvePiece<T, dim>;
        const T tol = o.tol;
        const T small = T(10) * tol; // pieces under this size are not halved any more
        const T gap = T(5) * small;  // pieces near a patch border stay unseeded up to twice `small`

        CurveCurveResult<T, dim> res;
        CurvePieces<T, dim> P1{c1}, P2{c2};
        // a stack per pair of initial pieces, examined in turn: curves running side by side just
        // beyond the tolerance take many pairs, they do not use up the budget of the other pairs
        std::vector<std::vector<std::pair<Piece, Piece>>> stacks;
        const auto pieces2 = P2.pieces(a2, b2);
        for (auto &p : P1.pieces(a1, b1))
            for (auto &q : pieces2)
                if (p.box.overlaps(q.box, tol))
                    stacks.emplace_back(1, std::pair{p, q});

        std::vector<std::array<T, 3>> candidates; // s, t, distance
        std::vector<Seed<T>> seeds;

        // Does the smaller piece lie on the other curve? Then it is a common portion.
        auto try_overlap = [&](const Piece &p, const Piece &q) -> bool {
            const bool p_small = p.size() <= q.size();
            const Piece &sm = p_small ? p : q, &ot = p_small ? q : p;
            if (sm.size() <= small || !ot.box.contains(sm.box, tol))
                return false;
            const auto &cs = p_small ? c1 : c2, &co = p_small ? c2 : c1;
            constexpr int n = 4;
            // curves within tol of each other over a step h are parallel within about 8 tol / h:
            // curves crossing at the samples are not a common portion
            const T max_sine = T(8 * n) * tol / sm.size();
            std::array<T, n + 1> proj{};
            for (int k = 0; k <= n; ++k)
            {
                const int kk = k == 0 ? 0 : (k == 1 ? n : k - 1); // the ends first: they fail most often
                const auto u = kk == n ? sm.b : sm.a + (sm.b - sm.a) * T(kk) / T(n);
                const auto [t, d] = project(co, ot.a, ot.b, cs.value(u));
                if (d > tol || sine(cs.value(u, 1), co.value(t, 1)) > max_sine)
                    return false;
                proj[kk] = t;
            }
            // on a closed curve, an end may project on either end of the range: take the one next
            // to the projection of its neighbour
            for (int k : {0, n})
            {
                const auto nb = proj[k == 0 ? 1 : n - 1];
                const auto alt = proj[k] == ot.a ? ot.b : (proj[k] == ot.b ? ot.a : proj[k]);
                if (alt != proj[k] && norm(co.value(alt) - co.value(proj[k])) <= tol && std::abs(alt - nb) < std::abs(proj[k] - nb))
                    proj[k] = alt;
            }
            if (p_small)
                seeds.push_back({sm.a, sm.b, proj[0], proj[n]});
            else if (proj[0] <= proj[n])
                seeds.push_back({proj[0], proj[n], sm.a, sm.b});
            else // opposite directions: ordered along c1
                seeds.push_back({proj[n], proj[0], sm.b, sm.a});
            return true;
        };

        // a crossing candidate by Newton from (s, t); from a saddle, on both sides. A start next to a
        // crossing just found (the pieces around it come in a row) gives it again: skipped, within
        // `small` of it on both curves, or within (ds, dt) of its parameters (a sample of a profile)
        std::array<std::pair<std::array<T, 2>, point<T, dim>>, 8> recent{};
        std::size_t n_recent{};
        auto follow = [&](T s, T t, T ds = T(0), T dt = T(0)) {
            const auto p1 = c1.value(s), p2 = c2.value(t);
            for (std::size_t k{}; k < std::min(n_recent, recent.size()); ++k)
            {
                const auto &[st, x] = recent[k];
                if (std::abs(s - st[0]) <= ds && std::abs(t - st[1]) <= dt)
                    return;
                if (norm(p1 - x) <= small && norm(p2 - x) <= small && norm(c2.value((t + st[1]) / T(2)) - x) <= small)
                    return;
            }
            bool escaped = false;
            for (const T side : {T(1), T(-1)})
            {
                const auto r = refine(c1, a1, b1, c2, a2, b2, s, t, side, &escaped);
                candidates.push_back(r);
                if (r[2] <= tol)
                    recent[n_recent++ % recent.size()] = {{r[0], r[1]}, c1.value(r[0])};
                if (!escaped)
                    break;
            }
        };

        // Flat, nearly parallel pieces: the distance from p to q at a few points of p, and Newton from
        // its local minima that may hide a crossing, and from the ends of q within the tolerance of
        // p. False, for the halves, when the pieces stay within the tolerance over two samples in a
        // row (a contact or a common portion: its seeds).
        auto side_by_side = [&](const Piece &p, const Piece &q) -> bool {
            constexpr int n = 16;
            std::array<T, n + 1> x{}, t{}, d{};
            const auto qc = q.last - q.first;
            const T q2 = sq_norm(qc);
            int run{};
            for (int k = 0; k <= n; ++k)
            {
                x[k] = k == n ? p.b : p.a + (p.b - p.a) * T(k) / T(n);
                const auto pt = c1.value(x[k]);
                const T l = std::clamp((pt - q.first) * qc / q2, T(0), T(1));
                const auto [tk, dk] = local_foot(c2, pt, q.a + l * (q.b - q.a), q.a, q.b);
                t[k] = tk;
                d[k] = dk;
                run = dk <= tol ? run + 1 : 0;
                if (run >= 2)
                    return false;
            }
            const auto pc = p.last - p.first;
            for (const T tq : {q.a, q.b})
            {
                const auto qt = c2.value(tq);
                const T l = std::clamp((qt - p.first) * pc / sq_norm(pc), T(0), T(1));
                if (const auto [sp, dq] = local_foot(c1, qt, p.a + l * (p.b - p.a), p.a, p.b); dq <= tol)
                    follow(sp, tq);
            }
            // between two samples, the distance dips below the lower one by less than its change over
            // a step (a crossing between them) plus its second difference (the curves bending apart)
            for (int k = 0; k <= n; ++k)
            {
                if (!((k == 0 || d[k] < d[k - 1]) && (k == n || d[k] <= d[k + 1])))
                    continue;
                const T dl = k > 0 ? d[k - 1] : d[k + 1], dr = k < n ? d[k + 1] : d[k - 1];
                if (d[k] <= tol + std::max(dl, dr) - d[k] + (dl - T(2) * d[k] + dr))
                    follow(x[k], t[k], (p.b - p.a) / T(n), std::abs(q.b - q.a) / T(n));
            }
            return true;
        };

        auto examine = [&](std::vector<std::pair<Piece, Piece>> &stack) {
            auto [p, q] = std::move(stack.back());
            stack.pop_back();
            if (!p.box.overlaps(q.box, tol) || apart(p, q, tol) || apart(q, p, tol))
                return;
            // a flat piece is a thin tube around its chord: test it against the other box
            if (p.flatness <= T(1e-3) * p.chord() && !segment_overlaps(p.first, p.last, q.box, p.flatness + tol))
                return;
            if (q.flatness <= T(1e-3) * q.chord() && !segment_overlaps(q.first, q.last, p.box, q.flatness + tol))
                return;
            if (try_overlap(p, q))
                return;
            const bool p_flat = p.flatness <= T(1e-3) * p.chord() || p.size() <= small;
            const bool q_flat = q.flatness <= T(1e-3) * q.chord() || q.size() <= small;
            const bool tiny = p.size() <= small && q.size() <= small;
            if (p_flat && q_flat)
            {
                const auto [s, t, d] = closest_on_chords(p.first, p.last, q.first, q.last);
                if (d > p.flatness + q.flatness + tol)
                    return;
                const T s0 = p.a + s * (p.b - p.a), t0 = q.a + t * (q.b - q.a);
                const T sin_pq = sine(p.last - p.first, q.last - q.first);
                if (tiny || sin_pq > T(0.1) || !(p.chord() > T(0)) || !(q.chord() > T(0)))
                {
                    follow(s0, t0);
                    return;
                }
                // nearly parallel pieces may cross twice, or run side by side just beyond the tolerance
                // (halved, their number would grow like sqrt(R / tol)): followed along p
                if (side_by_side(p, q))
                    return;
            }
            if (tiny)
            {
                follow((p.a + p.b) / T(2), (q.a + q.b) / T(2));
                return;
            }
            // halve the larger piece that is not flat
            const bool split_p = q_flat && !p_flat ? true : (p_flat && !q_flat ? false : p.size() >= q.size());
            if (split_p)
            {
                auto [l, r] = P1.split(p);
                stack.emplace_back(std::move(l), q);
                stack.emplace_back(std::move(r), std::move(q));
            }
            else
            {
                auto [l, r] = P2.split(q);
                stack.emplace_back(p, std::move(l));
                stack.emplace_back(std::move(p), std::move(r));
            }
        };

        constexpr std::size_t quantum = 256;
        std::size_t pairs{};
        for (bool pending = true; pending && res.complete;)
        {
            pending = false;
            for (auto &stack : stacks)
            {
                for (std::size_t k{}; k < quantum && !stack.empty() && res.complete; ++k)
                {
                    if (pairs++ >= o.max_pairs)
                        res.complete = false;
                    else
                        examine(stack);
                }
                pending = pending || !stack.empty();
            }
        }

        // ---- common portions: merge the seeds that follow each other on both curves. Across the
        // seam of a closed c2 the points meet but t2 jumps from one end of its range to the other:
        // two portions, each grown up to the seam
        auto inside = [](T x, T a, T b) { return (x - a) * (x - b) <= T(0); };
        auto along2 = [&](T t, T u) { // c2 goes from t to u within `gap`, not round its range
            const auto p = c2.value(t);
            return norm(c2.value(u) - p) <= gap && norm(c2.value((t + u) / T(2)) - p) <= gap;
        };
        std::ranges::sort(seeds, {}, &Seed<T>::s0);
        std::vector<Seed<T>> merged;
        for (const auto &sd : seeds)
        {
            if (!merged.empty())
            {
                auto &m = merged.back();
                const bool along1 = sd.s0 <= m.s1 || norm(c1.value(sd.s0) - c1.value(m.s1)) <= gap;
                if (along1 && (inside(sd.t0, m.t0, m.t1) || along2(m.t1, sd.t0)))
                {
                    if (sd.s1 > m.s1)
                    {
                        m.s1 = sd.s1;
                        m.t1 = sd.t1;
                    }
                    continue;
                }
            }
            merged.push_back(sd);
        }

        // distance from c1(s) to c2 on [lo, hi], by Newton from t
        auto foot = [&](T s, T t, T lo, T hi) { return newton_point_curve(c2, c1.value(s), t, std::array{lo, hi}, T(0)); };
        struct End
        {
            T s, t;
            bool departure; ///< the curves part there (rather than a range ending)
        };
        // From (s, t), within the tolerance, along c1 towards s_lim and along c2 towards t_lim: the
        // end of the common portion, where the curves part (d = tol) or where a range ends. The
        // steps double up to h_max (eight lengths of the seed): longer, they would hop over the gaps
        // between two curves crossing again and again
        auto grow = [&](T s, T t, T s_lim, T t_lim, T h_max) -> End {
            const T lo = std::min(t, t_lim), hi = std::max(t, t_lim);
            T h = small / std::max(norm(c1.value(s, 1)), std::numeric_limits<T>::min());
            h_max = std::max(h, h_max);
            while (s != s_lim)
            {
                const T s_n = s_lim > s ? std::min(s + h, s_lim) : std::max(s - h, s_lim);
                const auto [t_n, d_n] = local_foot(c2, c1.value(s_n), t, lo, hi);
                if (t_n == t_lim)
                {
                    // c2 ends first: the portion ends where c1 passes the end of c2
                    const auto [s_e, d_e] = newton_point_curve(c1, c2.value(t_lim), s_n, std::array{std::min(s, s_lim), std::max(s, s_lim)}, T(0));
                    if (d_e <= tol)
                        return {s_e, t_lim, false};
                }
                if (d_n > tol)
                {
                    T s_out = s_n; // bisection on d = tol
                    for (int it{}; it < 60 && norm(c1.value(s_out) - c1.value(s)) > T(1e-3) * tol; ++it)
                    {
                        const T s_m = (s + s_out) / T(2);
                        if (const auto [t_m, d_m] = local_foot(c2, c1.value(s_m), t, lo, hi); d_m <= tol)
                        {
                            s = s_m;
                            t = t_m;
                        }
                        else
                            s_out = s_m;
                    }
                    return {s, t, true};
                }
                s = s_n;
                t = t_n;
                h = std::min(T(2) * h, h_max);
            }
            return {s, t, false};
        };
        struct Stretch
        {
            T s0, s1, t0, t1;
            bool dep0, dep1;
        };
        std::vector<Stretch> stretches;
        for (const auto &m : merged)
        {
            const bool up = m.t1 != m.t0 ? m.t1 > m.t0 : c1.value(m.s0, 1) * c2.value(m.t0, 1) >= T(0);
            const auto e0 = grow(m.s0, m.t0, a1, up ? a2 : b2, T(8) * (m.s1 - m.s0)), e1 = grow(m.s1, m.t1, b1, up ? b2 : a2, T(8) * (m.s1 - m.s0));
            stretches.push_back({e0.s, e1.s, e0.t, e1.t, e0.departure, e1.departure});
        }
        // grown, two of them may cover each other
        std::ranges::sort(stretches, {}, &Stretch::s0);
        std::vector<Stretch> portions;
        for (const auto &st : stretches)
        {
            if (!portions.empty())
            {
                auto &m = portions.back();
                if (st.s0 <= m.s1 && (inside(st.t0, m.t0, m.t1) || along2(m.t1, st.t0)))
                {
                    if (st.s1 > m.s1)
                    {
                        m.s1 = st.s1;
                        m.t1 = st.t1;
                        m.dep1 = st.dep1;
                    }
                    continue;
                }
            }
            portions.push_back(st);
        }

        // A tangency also keeps the curves within the tolerance, over a length ~ sqrt(tol / relative
        // curvature). Along it, as along a crossing at a small angle or a line grazing a circle, the
        // gap c1 - c2 is a quadratic of the parameter; from an end where the curves part, the
        // distance grows as a parabola from its vertex, d'/d'' away. A common portion is not a
        // quadratic (an exact part, two curves waving within the tolerance of each other), or has a
        // core where the curves stay together besides the departures at its ends. Such a portion
        // longer than ten tolerances is an overlap; otherwise the local minima of the distance along
        // it are crossing candidates (one contact, for the deduplication below). Each portion is
        // measured on its own range of c2.
        std::vector<Stretch> contacts; // the curves within the tolerance all along: one crossing
        for (const auto &st : portions)
        {
            const T w = st.s1 - st.s0;
            if (!(w > T(0)))
            {
                follow(st.s0, st.t0);
                continue;
            }
            constexpr int n = 16;
            std::array<T, n + 1> x{}, t{}, d{};
            std::array<point<T, dim>, n + 1> r{};
            const T lo = std::min(st.t0, st.t1), hi = std::max(st.t0, st.t1);
            T len{}, d_max{};
            for (int k = 0; k <= n; ++k)
            {
                x[k] = k == n ? st.s1 : st.s0 + w * T(k) / T(n);
                const auto [tk, dk] = foot(x[k], st.t0 + (st.t1 - st.t0) * T(k) / T(n), lo, hi);
                t[k] = tk;
                d[k] = dk;
                r[k] = c1.value(x[k]) - c2.value(tk);
                d_max = std::max(d_max, dk);
                if (k > 0)
                    len += norm(c1.value(x[k]) - c1.value(x[k - 1]));
            }
            // least squares quadratic of the gap, on the basis 1, xi, xi² - mean (orthogonal on the
            // samples xi = 2k / n - 1): its largest residual
            auto residual = [&] {
                std::array<std::array<T, n + 1>, 3> P{};
                T m2{};
                for (int k = 0; k <= n; ++k)
                    m2 += std::pow(T(2 * k) / T(n) - T(1), 2) / T(n + 1);
                for (int k = 0; k <= n; ++k)
                {
                    const T xi = T(2 * k) / T(n) - T(1);
                    P[0][k] = T(1);
                    P[1][k] = xi;
                    P[2][k] = xi * xi - m2;
                }
                auto fit = r;
                for (auto &f : fit)
                    f = point<T, dim>{};
                for (const auto &Pi : P)
                {
                    point<T, dim> num{};
                    T den{};
                    for (int k = 0; k <= n; ++k)
                    {
                        num = num + Pi[k] * r[k];
                        den += Pi[k] * Pi[k];
                    }
                    for (int k = 0; k <= n; ++k)
                        fit[k] = fit[k] + (Pi[k] / den) * num;
                }
                T e{};
                for (int k = 0; k <= n; ++k)
                    e = std::max(e, norm(r[k] - fit[k]));
                return e;
            };
            // from the end e inwards (dir), to the vertex of the parabola; not beyond two samples in a
            // row where the curves coincide (a departure growing slower than a parabola, d'' small at
            // the end, then an exact part: a tangency touches at one point)
            auto to_vertex = [&](int e, int dir) {
                const T h = w / T(n), d1 = d[e + dir], d2 = d[e + 2 * dir];
                const T slope = (T(3) * d[e] - T(4) * d1 + d2) / (T(2) * h), curv = (d[e] - T(2) * d1 + d2) / (h * h);
                const T v = slope > T(0) && curv > T(0) ? slope / curv : std::numeric_limits<T>::infinity(); // no parabola
                for (int k = e; k + dir >= 0 && k + dir <= n; k += dir)
                    if (std::max(d[k], d[k + dir]) <= T(1e-3) * tol)
                        return std::min(v, std::abs(x[k] - x[e]));
                return v;
            };
            // an end of a range where the curves already part as a parabola, from a vertex inside the
            // portion or at its other end, is a departure too, cut short (a carrier clipped where it
            // leaves the tolerance: a circle touching a border, restricted to the box of a face, maybe
            // cut at its seam); a linear gap (two lines at a small angle) has its vertex far away
            auto departure = [&](bool dep, int e, int dir) {
                const T v = to_vertex(e, dir);
                return dep || v <= T(2) * w ? v : T(0);
            };
            const T core = w - departure(st.dep0, 0, 1) - departure(st.dep1, n, -1);
            if (d_max <= tol && len > small && ((core > w / T(20) && core / w * len > small) || residual() > tol / T(8)))
            {
                res.overlaps.push_back({st.s0, st.s1, st.t0, st.t1});
                continue;
            }
            if (d_max <= tol)
                contacts.push_back(st);
            for (int k = 0; k <= n; ++k)
                if ((k == 0 || d[k] < d[k - 1]) && (k == n || d[k] <= d[k + 1]))
                    follow(x[k], t[k]);
        }

        // ---- crossings: those within the tolerance and off the overlaps, sorted, one per contact
        struct Found
        {
            CurveCurveHit<T, dim> h;
            T d;
            bool spread; ///< several crossings, more than `small` apart, between which the curves stay together
        };
        auto on_overlap = [&](const CurveCurveHit<T, dim> &h) {
            return std::ranges::any_of(res.overlaps, [&](const CurveOverlap<T> &ov) {
                return (h.t1 >= ov.t1a && h.t1 <= ov.t1b) || norm(h.p - c1.value(ov.t1a)) <= gap || norm(h.p - c1.value(ov.t1b)) <= gap;
            });
        };
        std::vector<Found> found;
        for (const auto &[s, t, d] : candidates)
            if (d <= tol)
                if (const CurveCurveHit<T, dim> h{s, t, (c1.value(s) + c2.value(t)) / T(2), false}; !on_overlap(h))
                    found.push_back({h, d, false});
        std::ranges::sort(found, {}, [](const Found &f) { return f.h.t1; });

        // Two crossings are one when the curves stay within the tolerance of each other between
        // them, on the ways between them on both curves, round the seam of a closed curve if
        // shorter (a crossing at the seam is found at both ends of the range): close to each other
        // and to the middles of the ways, on a contact portion (where that was checked), or within
        // the tolerance and parallel at a few points of the ways. Those points are off any regular
        // grid (the crossings of a sine wave with its axis, of two interpolations of a circle, fall
        // on one) and the parallelism rejects curves crossing there. The closest points are kept.
        const bool closed1 = norm(c1.value(a1) - c1.value(b1)) <= tol, closed2 = norm(c2.value(a2) - c2.value(b2)) <= tol;
        auto on_contact = [&](const Stretch &st, const CurveCurveHit<T, dim> &h) {
            const bool on1 = inside(h.t1, st.s0, st.s1) || norm(h.p - c1.value(st.s0)) <= small || norm(h.p - c1.value(st.s1)) <= small;
            return on1 && (inside(h.t2, st.t0, st.t1) || norm(h.p - c2.value(st.t0)) <= small || norm(h.p - c2.value(st.t1)) <= small);
        };
        auto near = [&](const CurveCurveHit<T, dim> &u, const CurveCurveHit<T, dim> &h) {
            return norm(u.p - h.p) <= small && norm(c1.value(between(u.t1, h.t1, a1, b1, closed1, T(0.5))) - h.p) <= small &&
                   norm(c2.value(between(u.t2, h.t2, a2, b2, closed2, T(0.5))) - h.p) <= small;
        };
        auto same_contact = [&](const CurveCurveHit<T, dim> &u, const CurveCurveHit<T, dim> &h) {
            if (near(u, h))
                return true;
            // an overlap between them separates them
            auto on_way = [&](T x) {
                const T lo = std::min(u.t1, h.t1), hi = std::max(u.t1, h.t1);
                return closed1 && hi - lo > (b1 - a1) / T(2) ? (x > hi || x < lo) : (x > lo && x < hi);
            };
            if (std::ranges::any_of(res.overlaps, [&](const CurveOverlap<T> &ov) { return on_way(ov.t1a) || on_way(ov.t1b); }))
                return false;
            if (std::ranges::any_of(contacts, [&](const Stretch &st) { return on_contact(st, u) && on_contact(st, h); }))
                return true;
            const auto m1 = c1.value(between(u.t1, h.t1, a1, b1, closed1, T(0.5)));
            const T max_sine = T(32) * tol / std::max(norm(m1 - u.p) + norm(h.p - m1), tol);
            for (const auto &e : {u, h}) // first at the crossings themselves
                if (sine(c1.value(e.t1, 1), c2.value(e.t2, 1)) > max_sine)
                    return false;
            for (T f : {T(0.47), T(0.29), T(0.71), T(0.13), T(0.87), T(0.38), T(0.61)})
            {
                const T s = between(u.t1, h.t1, a1, b1, closed1, f), t = between(u.t2, h.t2, a2, b2, closed2, f);
                const auto [t_s, d_s] = project_on_way(c2, u.t2, h.t2, a2, b2, closed2, c1.value(s));
                if (d_s > tol || sine(c1.value(s, 1), c2.value(t_s, 1)) > max_sine)
                    return false;
                const auto [s_t, d_t] = project_on_way(c1, u.t1, h.t1, a1, b1, closed1, c2.value(t));
                if (d_t > tol || sine(c1.value(s_t, 1), c2.value(t, 1)) > max_sine)
                    return false;
            }
            return true;
        };
        auto absorb = [&](Found &u, const Found &f) {
            u.spread = u.spread || f.spread || norm(u.h.p - f.h.p) > small;
            if (f.d < u.d)
            {
                u.h = f.h;
                u.d = f.d;
            }
        };
        // first, in one pass along c1, the copies of a point (from the pieces around it) and the
        // points of a contact portion, by their parameters
        auto portion_of = [&](const CurveCurveHit<T, dim> &h) {
            const auto it = std::ranges::find_if(contacts, [&](const Stretch &st) { return inside(h.t1, st.s0, st.s1) && inside(h.t2, st.t0, st.t1); });
            return it - contacts.begin();
        };
        std::vector<Found> points;
        auto last_portion = contacts.end() - contacts.begin();
        for (const auto &f : found)
        {
            const auto portion = portion_of(f.h);
            if (!points.empty() && ((portion == last_portion && portion < std::ssize(contacts)) || near(points.back().h, f.h)))
                absorb(points.back(), f);
            else
                points.push_back(f);
            last_portion = portion;
        }
        // then the crossings of a contact, which follow each other along c1 (round its seam if closed)
        std::vector<Found> unique;
        for (const auto &f : points)
        {
            if (!unique.empty() && same_contact(unique.back().h, f.h))
                absorb(unique.back(), f);
            else
                unique.push_back(f);
        }
        if (closed1 && unique.size() > 1 && same_contact(unique.back().h, unique.front().h))
        {
            absorb(unique.front(), unique.back());
            unique.pop_back();
        }
        std::ranges::sort(unique, {}, [](const Found &f) { return f.h.t1; });
        for (auto &f : unique)
        {
            f.h.tangent = f.spread || sine(c1.value(f.h.t1, 1), c2.value(f.h.t2, 1)) < o.angular_tol;
            res.hits.push_back(f.h);
        }
        return res;
    }

    /// Every crossing and common portion of two curves on their whole bounds.
    template <std::floating_point T, std::size_t dim>
    [[nodiscard]] auto intersect(const Curve<T, dim> &c1, const Curve<T, dim> &c2, IntersectOptions<T> o = {}) -> CurveCurveResult<T, dim>
    {
        const auto b1 = c1.bounds(), b2 = c2.bounds();
        return intersect(c1, b1[0], b1[1], c2, b2[0], b2[1], o);
    }
}
