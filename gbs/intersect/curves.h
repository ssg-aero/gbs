#pragma once
/**
 * @file curves.h
 * @brief Intersection of two curves: every crossing and every common portion on bounded
 * parameter ranges.
 *
 * The curves are cut into pieces (`CurvePieces`): pairs of pieces whose boxes are farther
 * than the tolerance are dropped, the others are halved until both pieces are flat. Two flat
 * pieces give a crossing candidate from the closest points of their chords, refined by
 * Gauss-Newton on the whole curves. A piece lying within the tolerance of the other curve is a
 * common portion (overlap) and is not subdivided further. Written for any dimension: in 2D it
 * crosses pcurves in the parameter space of a face, in 3D it finds the points where two curves
 * meet within the tolerance.
 */
#include <gbs/intersect/curve_pieces.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
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
        bool complete{true};                     ///< false when `max_pairs` was reached
    };

    namespace intersect_detail
    {
        /// Gauss-Newton on |c1(s) - c2(t)|², parameters kept in their ranges.
        /// Returns the parameters and the distance reached.
        template <std::floating_point T, std::size_t dim>
        auto refine(const Curve<T, dim> &c1, T a1, T b1, const Curve<T, dim> &c2, T a2, T b2, T s, T t, T tol) -> std::array<T, 3>
        {
            auto d = norm(c1.value(s) - c2.value(t));
            for (int it{}; it < 100; ++it)
            {
                const auto r = c1.value(s) - c2.value(t);
                const auto d1 = c1.value(s, 1), d2 = c2.value(t, 1);
                const auto A = d1 * d1, B = -(d1 * d2), C = d2 * d2;
                const auto g1 = d1 * r, g2 = -(d2 * r);
                const auto det = A * C - B * B;
                T ds, dt;
                if (det > T(1e-14) * A * C && det > T(0))
                {
                    ds = -(C * g1 - B * g2) / det;
                    dt = -(A * g2 - B * g1) / det;
                }
                else // parallel tangents: steepest projection of each curve on the other
                {
                    ds = A > T(0) ? -g1 / A : T(0);
                    dt = C > T(0) ? -g2 / C : T(0);
                }
                const auto s_new = std::clamp(s + ds, a1, b1), t_new = std::clamp(t + dt, a2, b2);
                const auto d_new = norm(c1.value(s_new) - c2.value(t_new));
                const auto step = std::abs(s_new - s) * std::sqrt(A) + std::abs(t_new - t) * std::sqrt(C);
                if (d_new > d && d <= tol)
                    break; // no further progress
                s = s_new;
                t = t_new;
                d = d_new;
                if (step <= tol * T(1e-6) || d <= tol * T(1e-9))
                    break;
            }
            return {s, t, d};
        }

        /// Projection of `p` on `c` restricted to [a, b], from the nearest of a few samples.
        template <std::floating_point T, std::size_t dim>
        auto project(const Curve<T, dim> &c, T a, T b, const point<T, dim> &p) -> std::pair<T, T>
        {
            constexpr int n = 8;
            T t = a, best = norm(c.value(a) - p);
            for (int k = 1; k <= n; ++k)
            {
                const auto u = k == n ? b : a + (b - a) * T(k) / T(n);
                const auto d = norm(c.value(u) - p);
                if (d < best)
                {
                    best = d;
                    t = u;
                }
            }
            for (int it{}; it < 50; ++it)
            {
                const auto r = c.value(t) - p;
                const auto d1 = c.value(t, 1), d2 = c.value(t, 2);
                const auto f = d1 * r, df = d2 * r + d1 * d1;
                if (df <= T(0))
                    break;
                const auto t_new = std::clamp(t - f / df, a, b);
                const auto moved = std::abs(t_new - t);
                t = t_new;
                if (moved * std::sqrt(d1 * d1) <= T(1e-14) * (T(1) + norm(p)))
                    break;
            }
            return {t, norm(c.value(t) - p)};
        }

        /// Parameters of the closest points of the chords [p1, q1] and [p2, q2], in [0, 1].
        template <std::floating_point T, std::size_t dim>
        auto closest_on_chords(const point<T, dim> &p1, const point<T, dim> &q1, const point<T, dim> &p2, const point<T, dim> &q2)
            -> std::array<T, 3>
        {
            const auto d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
            const auto a = d1 * d1, e = d2 * d2, f = d2 * r;
            T s{}, t{};
            if (a <= T(0) && e <= T(0))
                return {T(0), T(0), norm(r)};
            if (a <= T(0))
                t = std::clamp(f / e, T(0), T(1));
            else
            {
                const auto c = d1 * r;
                if (e <= T(0))
                    s = std::clamp(-c / a, T(0), T(1));
                else
                {
                    const auto b = d1 * d2, den = a * e - b * b;
                    s = den > T(0) ? std::clamp((b * f - c * e) / den, T(0), T(1)) : T(0);
                    t = (b * s + f) / e;
                    if (t < T(0))
                    {
                        t = T(0);
                        s = std::clamp(-c / a, T(0), T(1));
                    }
                    else if (t > T(1))
                    {
                        t = T(1);
                        s = std::clamp((b - c) / a, T(0), T(1));
                    }
                }
            }
            return {s, t, norm(p1 + s * d1 - (p2 + t * d2))};
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
    }

    /**
     * @brief Every crossing and common portion of c1 on [a1, b1] and c2 on [a2, b2].
     *
     * A crossing is reported once, at the parameters where the curves are the closest (both
     * within `tol`); the crossings at the ends of the ranges are included. A common portion
     * longer than about ten times the tolerance is an overlap, and the crossings on it are
     * dropped.
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

        CurveCurveResult<T, dim> res;
        CurvePieces<T, dim> P1{c1}, P2{c2};
        std::vector<std::pair<Piece, Piece>> stack;
        for (auto &p : P1.pieces(a1, b1))
            for (auto &q : P2.pieces(a2, b2))
                if (p.box.overlaps(q.box, tol))
                    stack.emplace_back(p, q);

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
            std::array<T, n + 1> proj{};
            for (int k = 0; k <= n; ++k)
            {
                const int kk = k == 0 ? 0 : (k == 1 ? n : k - 1); // the ends first: they fail most often
                const auto u = kk == n ? sm.b : sm.a + (sm.b - sm.a) * T(kk) / T(n);
                const auto [t, d] = project(co, ot.a, ot.b, cs.value(u));
                if (d > tol)
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

        std::size_t pairs{};
        while (!stack.empty())
        {
            auto [p, q] = std::move(stack.back());
            stack.pop_back();
            if (++pairs > o.max_pairs)
            {
                res.complete = false;
                break;
            }
            if (!p.box.overlaps(q.box, tol))
                continue;
            // a flat piece is a thin tube around its chord: test it against the other box
            if (p.flatness <= T(1e-3) * p.chord() && !segment_overlaps(p.first, p.last, q.box, p.flatness + tol))
                continue;
            if (q.flatness <= T(1e-3) * q.chord() && !segment_overlaps(q.first, q.last, p.box, q.flatness + tol))
                continue;
            if (try_overlap(p, q))
                continue;
            const bool p_flat = p.flatness <= T(1e-3) * p.chord() || p.size() <= small;
            const bool q_flat = q.flatness <= T(1e-3) * q.chord() || q.size() <= small;
            const bool tiny = p.size() <= small && q.size() <= small;
            if (p_flat && q_flat)
            {
                const auto [s, t, d] = closest_on_chords(p.first, p.last, q.first, q.last);
                if (d > p.flatness + q.flatness + tol)
                    continue;
                // nearly parallel pieces may cross twice: halve them until they are tiny
                if (tiny || sine(p.last - p.first, q.last - q.first) > T(0.1))
                {
                    candidates.push_back(refine(c1, a1, b1, c2, a2, b2, p.a + s * (p.b - p.a), q.a + t * (q.b - q.a), tol));
                    continue;
                }
            }
            if (tiny)
            {
                candidates.push_back(refine(c1, a1, b1, c2, a2, b2, (p.a + p.b) / T(2), (q.a + q.b) / T(2), tol));
                continue;
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
        }

        // ---- overlaps: merge the seeds that follow each other on both curves
        std::ranges::sort(seeds, {}, &Seed<T>::s0);
        const T gap = T(5) * small; // pieces near a patch border stay unseeded up to twice `small`
        std::vector<Seed<T>> merged;
        for (const auto &sd : seeds)
        {
            if (!merged.empty())
            {
                auto &m = merged.back();
                const bool along1 = sd.s0 <= m.s1 || norm(c1.value(sd.s0) - c1.value(m.s1)) <= gap;
                const bool along2 = norm(c2.value(sd.t0) - c2.value(m.t1)) <= gap;
                if (along1 && along2)
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

        // A tangency also keeps the curves within the tolerance over a length ~ sqrt(tol / curvature
        // gap): along a true overlap the distance hardly varies, along a tangency it grows from about
        // zero at the contact to the tolerance at the ends. Tangencies become crossing candidates.
        std::erase_if(merged, [&](const Seed<T> &m) {
            constexpr int n = 8;
            T d_min = std::numeric_limits<T>::max(), d_max{}, s_min = m.s0, t_min = m.t0;
            for (int k = 1; k < n; ++k)
            {
                const auto s = m.s0 + (m.s1 - m.s0) * T(k) / T(n);
                const auto [t, d] = project(c2, std::min(m.t0, m.t1), std::max(m.t0, m.t1), c1.value(s));
                if (d < d_min)
                {
                    d_min = d;
                    s_min = s;
                    t_min = t;
                }
                d_max = std::max(d_max, d);
            }
            if (d_max - d_min <= tol / T(4))
                return false;
            candidates.push_back(refine(c1, a1, b1, c2, a2, b2, s_min, t_min, tol));
            return true;
        });

        // ---- crossings: keep those within the tolerance, sorted, without duplicates
        std::vector<CurveCurveHit<T, dim>> hits;
        for (const auto &[s, t, d] : candidates)
            if (d <= tol)
                hits.push_back({s, t, (c1.value(s) + c2.value(t)) / T(2), false});
        std::ranges::sort(hits, {}, &CurveCurveHit<T, dim>::t1);
        std::vector<CurveCurveHit<T, dim>> unique;
        for (const auto &h : hits)
        {
            const bool dup = std::ranges::any_of(unique, [&](const auto &u) {
                return norm(u.p - h.p) <= small && norm(c1.value((u.t1 + h.t1) / T(2)) - h.p) <= small &&
                       norm(c2.value((u.t2 + h.t2) / T(2)) - h.p) <= small;
            });
            if (!dup)
                unique.push_back(h);
        }

        // crossings at the ends of an overlap extend it; those inside it are dropped
        for (auto &m : merged)
        {
            const auto e0 = c1.value(m.s0), e1 = c1.value(m.s1);
            std::erase_if(unique, [&](const auto &h) {
                const bool inside = h.t1 >= m.s0 && h.t1 <= m.s1;
                const bool at_start = norm(h.p - e0) <= gap, at_end = norm(h.p - e1) <= gap;
                if (at_start && h.t1 < m.s0)
                {
                    m.s0 = h.t1;
                    m.t0 = h.t2;
                }
                if (at_end && h.t1 > m.s1)
                {
                    m.s1 = h.t1;
                    m.t1 = h.t2;
                }
                return inside || at_start || at_end;
            });
        }
        for (const auto &m : merged)
            if (norm(c1.value(m.s1) - c1.value(m.s0)) > small || norm(c1.value((m.s0 + m.s1) / T(2)) - c1.value(m.s0)) > small)
                res.overlaps.push_back({m.s0, m.s1, m.t0, m.t1});

        // a closed curve meets the other curve at both ends of its range: once
        auto closed_dup = [&](auto key, const Curve<T, dim> &c, T a, T b) {
            if (unique.size() < 2 || norm(c.value(a) - c.value(b)) > tol)
                return;
            auto lo = std::ranges::min_element(unique, {}, key), hi = std::ranges::max_element(unique, {}, key);
            if (lo != hi && std::invoke(key, *lo) == a && std::invoke(key, *hi) == b && norm(lo->p - hi->p) <= small)
                unique.erase(hi);
        };
        closed_dup(&CurveCurveHit<T, dim>::t1, c1, a1, b1);
        closed_dup(&CurveCurveHit<T, dim>::t2, c2, a2, b2);

        for (auto &h : unique)
            h.tangent = sine(c1.value(h.t1, 1), c2.value(h.t2, 1)) < o.angular_tol;
        res.hits = std::move(unique);
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
