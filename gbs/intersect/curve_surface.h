#pragma once
/**
 * @file curve_surface.h
 * @brief Intersection of a curve and a surface: every crossing and every portion of the curve
 * lying on the surface, on bounded parameter ranges.
 *
 * Closed forms: a line (`LineCurve`, or a B-spline of degree 1 with two poles) against a plane,
 * a cylinder, a cone or a sphere is a quadratic equation, against a torus a quartic; a circle or
 * an ellipse against a plane is a trigonometric equation. Any other pair, and the configurations
 * where the equation vanishes (the line lies on the surface), go through subdivision: pieces of
 * the curve against patches of the surface (`CurvePieces`, `SurfacePieces`), then Gauss-Newton
 * on (t, u, v).
 */
#include <gbs/intersect/curves.h>
#include <gbs/intersect/polynomial.h>
#include <gbs/intersect/surface_pieces.h>
#include <gbs/primitives.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs
{
    /// A crossing of a curve and a surface: c(t) = s(u, v) = p within the tolerance.
    template <std::floating_point T>
    struct CurveSurfaceHit
    {
        T t{}, u{}, v{};
        point<T, 3> p{};
        bool tangent{}; ///< the curve is tangent to the surface at p
    };

    /// A portion of the curve lying on the surface: c on [ta, tb] (ta < tb), from s(ua, va) to s(ub, vb).
    template <std::floating_point T>
    struct CurveSurfacePortion
    {
        T ta{}, tb{}, ua{}, va{}, ub{}, vb{};
    };

    template <std::floating_point T>
    struct CurveSurfaceResult
    {
        std::vector<CurveSurfaceHit<T>> hits;         ///< sorted by t, none inside a portion
        std::vector<CurveSurfacePortion<T>> portions; ///< sorted by ta
        bool complete{true};                          ///< false when `max_pairs` was reached
    };

    namespace intersect_detail
    {
        /// Is (u, v) in the box, up to a relative rounding margin?
        template <std::floating_point T>
        bool in_box(const std::array<T, 4> &box, T u, T v)
        {
            const T eu = T(1e-10) * (T(1) + std::abs(box[0]) + std::abs(box[1]));
            const T ev = T(1e-10) * (T(1) + std::abs(box[2]) + std::abs(box[3]));
            return u >= box[0] - eu && u <= box[1] + eu && v >= box[2] - ev && v <= box[3] + ev;
        }

        /// Gauss-Newton on |c(t) - s(u, v)|², parameters kept in their ranges. Returns t, u, v, distance.
        template <std::floating_point T>
        auto refine_cs(const Curve<T, 3> &c, T a, T b, const Surface<T, 3> &s, const std::array<T, 4> &box, T t, T u, T v, T tol)
            -> std::array<T, 4>
        {
            auto d = norm(c.value(t) - s.value(u, v));
            for (int it{}; it < 100; ++it)
            {
                const auto r = c.value(t) - s.value(u, v);
                const auto ct = c.value(t, 1), su = s.value(u, v, 1, 0), sv = s.value(u, v, 0, 1);
                // normal equations of J = [ct, -su, -sv]
                const std::array<point<T, 3>, 3> J{ct, T(-1) * su, T(-1) * sv};
                T A[3][3], g[3];
                for (int i = 0; i < 3; ++i)
                {
                    g[i] = J[i] * r;
                    for (int j = 0; j < 3; ++j)
                        A[i][j] = J[i] * J[j];
                }
                const T det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                              A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
                std::array<T, 3> dx{};
                if (std::abs(det) > T(1e-14) * A[0][0] * A[1][1] * A[2][2] && det != T(0))
                {
                    for (int k = 0; k < 3; ++k)
                    {
                        T M[3][3];
                        for (int i = 0; i < 3; ++i)
                            for (int j = 0; j < 3; ++j)
                                M[i][j] = j == k ? -g[i] : A[i][j];
                        dx[k] = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
                                 M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])) /
                                det;
                    }
                }
                else // tangent: project the curve point on the surface, and the surface point on the curve
                {
                    const T d2 = A[1][1] * A[2][2] - A[1][2] * A[1][2];
                    if (d2 > T(0))
                    {
                        dx[1] = -(A[2][2] * g[1] - A[1][2] * g[2]) / d2;
                        dx[2] = -(A[1][1] * g[2] - A[1][2] * g[1]) / d2;
                    }
                    dx[0] = A[0][0] > T(0) ? -g[0] / A[0][0] : T(0);
                }
                const auto t_new = std::clamp(t + dx[0], a, b), u_new = std::clamp(u + dx[1], box[0], box[1]),
                           v_new = std::clamp(v + dx[2], box[2], box[3]);
                const auto d_new = norm(c.value(t_new) - s.value(u_new, v_new));
                const auto step = std::abs(t_new - t) * norm(ct) + std::abs(u_new - u) * norm(su) + std::abs(v_new - v) * norm(sv);
                if (d_new > d && d <= tol)
                    break;
                t = t_new;
                u = u_new;
                v = v_new;
                d = d_new;
                if (step <= tol * T(1e-6) || d <= tol * T(1e-9))
                    break;
            }
            return {t, u, v, d};
        }

        /**
         * Projection of `p` on the surface restricted to `box`, by Newton from the two
         * nearest points of a grid (on a closed surface the nearest may sit on the wrong side of
         * the seam, against the bound of the box). Returns u, v, distance.
         */
        template <std::floating_point T>
        auto project_on_surface(const Surface<T, 3> &s, const std::array<T, 4> &box, const point<T, 3> &p) -> std::array<T, 3>
        {
            return project_point_surface(s, p, box, 5, 5, T(0));
        }

        /// The line carried by a curve, if it is one: c(t) = o + t d.
        template <std::floating_point T>
        auto as_line(const Curve<T, 3> &c) -> std::optional<std::array<point<T, 3>, 2>>
        {
            if (const auto *l = dynamic_cast<const LineCurve<T> *>(&c))
                return std::array<point<T, 3>, 2>{l->origin(), l->vector()};
            if (const auto *bs = dynamic_cast<const BSCurve<T, 3> *>(&c); bs && bs->degree() == 1 && bs->poles().size() == 2)
            {
                const auto &k = bs->knotsFlats();
                const T a = k.front(), b = k.back();
                const auto d = (bs->poles()[1] - bs->poles()[0]) / (b - a);
                return std::array<point<T, 3>, 2>{bs->poles()[0] - a * d, d};
            }
            return std::nullopt;
        }

        /// u (and v for a torus) of an angle, in the turn of the box.
        template <std::floating_point T>
        T box_turn(T x, T lo, T hi) { return primitives_detail::wrap(x, T(0.5) * (lo + hi) - std::numbers::pi_v<T>); }

        /**
         * Candidate parameters of the curve from a closed form, or nullopt when there is none for
         * the pair or the equation vanishes (curve on the surface).
         */
        template <std::floating_point T>
        auto closed_form(const Curve<T, 3> &c, T a, T b, const Surface<T, 3> &s) -> std::optional<std::vector<T>>
        {
            const auto *es = dynamic_cast<const ElementarySurface<T> *>(&s);
            if (!es)
                return std::nullopt;
            const auto &f = es->frame();
            auto finish = [&](const std::vector<T> &coef) -> std::optional<std::vector<T>> {
                T scale{};
                for (T x : coef)
                    scale = std::max(scale, std::abs(x));
                bool vanishes = true;
                for (T x : coef)
                    vanishes = vanishes && std::abs(x) <= T(1e-12) * std::max(scale, T(1));
                if (vanishes)
                    return std::nullopt;
                const auto r = polynomial_roots(coef, a, b);
                std::vector<T> ts = r.roots;
                ts.insert(ts.end(), r.extremes.begin(), r.extremes.end());
                return ts;
            };
            if (const auto line = as_line(c))
            {
                const auto q = (*line)[0] - f.o, d = (*line)[1];
                // in the frame of the surface: q + t d = (x, y, z)
                const T qx = q * f.x, qy = q * f.y, qz = q * f.z, dx = d * f.x, dy = d * f.y, dz = d * f.z;
                if (dynamic_cast<const Plane<T> *>(es))
                    return finish({qz, dz});
                // radial part: rho² = A t² + B t + C
                const T A = dx * dx + dy * dy, B = T(2) * (qx * dx + qy * dy), C = qx * qx + qy * qy;
                if (const auto *cy = dynamic_cast<const CylindricalSurface<T> *>(es))
                    return finish({C - cy->radius() * cy->radius(), B, A});
                if (const auto *co = dynamic_cast<const ConicalSurface<T> *>(es))
                {
                    // rho² = (R + z tan a)², z = qz + t dz
                    const T k = std::tan(co->semi_angle()), r0 = co->radius() + k * qz, r1 = k * dz;
                    return finish({C - r0 * r0, B - T(2) * r0 * r1, A - r1 * r1});
                }
                if (const auto *sp = dynamic_cast<const SphericalSurface<T> *>(es))
                    return finish({C + qz * qz - sp->radius() * sp->radius(), B + T(2) * qz * dz, A + dz * dz});
                if (const auto *to = dynamic_cast<const ToroidalSurface<T> *>(es))
                {
                    // (|x|² + R² - r²)² = 4 R² rho²
                    const T R = to->major_radius(), r = to->minor_radius();
                    const T a2 = A + dz * dz, a1 = B + T(2) * qz * dz, a0 = C + qz * qz + R * R - r * r;
                    const T R4 = T(4) * R * R;
                    return finish({a0 * a0 - R4 * C, T(2) * a0 * a1 - R4 * B, a1 * a1 + T(2) * a0 * a2 - R4 * A, T(2) * a1 * a2, a2 * a2});
                }
                return std::nullopt;
            }
            if (const auto *el = dynamic_cast<const EllipseCurve<T> *>(&c); el && dynamic_cast<const Plane<T> *>(es))
            {
                // (C(t) - O).z = d0 + d1 cos t + d2 sin t = 0
                const auto &g = el->frame();
                const T d0 = (g.o - f.o) * f.z, d1 = el->radius1() * (g.x * f.z), d2 = el->radius2() * (g.y * f.z);
                const T rho = std::hypot(d1, d2);
                const T scale = std::max({std::abs(d0), rho, T(1)});
                if (rho <= T(1e-12) * scale)
                    return std::abs(d0) <= T(1e-12) * scale ? std::nullopt : std::optional<std::vector<T>>{std::vector<T>{}};
                const T phi = std::atan2(d2, d1);
                std::vector<T> ts;
                const T x = -d0 / rho;
                if (x > T(1) + T(1e-12) || x < T(-1) - T(1e-12))
                    return ts;
                const T w = std::acos(std::clamp(x, T(-1), T(1)));
                for (T base : {phi + w, phi - w})
                    for (int k = -2; k <= 2; ++k)
                    {
                        const T t = base + T(k) * T(2) * std::numbers::pi_v<T>;
                        if (t >= a - T(1e-12) && t <= b + T(1e-12))
                            ts.push_back(std::clamp(t, a, b));
                    }
                return ts;
            }
            return std::nullopt;
        }
    }

    /**
     * @brief Every crossing of the curve c on [a, b] with the surface s on the parameter box
     * `box` = {u1, u2, v1, v2}, and every portion of the curve lying on the surface.
     *
     * A crossing is reported once, at the parameters where the curve and the surface are the
     * closest (within `tol`); a portion longer than about ten times the tolerance along which
     * the distance hardly varies lies on the surface, and the crossings on it are dropped.
     */
    template <std::floating_point T>
    [[nodiscard]] auto intersect(const Curve<T, 3> &c, std::type_identity_t<T> a, std::type_identity_t<T> b, const Surface<T, 3> &s,
                                 const std::array<T, 4> &box, IntersectOptions<T> o = {}) -> CurveSurfaceResult<T>
    {
        using namespace intersect_detail;
        const T tol = o.tol;
        const T small = T(10) * tol;
        CurveSurfaceResult<T> res;
        std::vector<std::array<T, 4>> candidates; // t, u, v, distance
        struct Seed
        {
            T t0, t1;
            std::array<T, 2> uv0, uv1;
        };
        std::vector<Seed> seeds;

        const auto *es = dynamic_cast<const ElementarySurface<T> *>(&s);
        auto uv_of = [&](const point<T, 3> &p) -> std::array<T, 2> {
            auto uv = es->parameters(p);
            if (!dynamic_cast<const Plane<T> *>(es))
                uv[0] = box_turn(uv[0], box[0], box[1]);
            if (dynamic_cast<const ToroidalSurface<T> *>(es))
                uv[1] = box_turn(uv[1], box[2], box[3]);
            return uv;
        };

        if (const auto ts = closed_form(c, a, b, s))
        {
            for (T t : *ts)
            {
                const auto p = c.value(t);
                const auto uv = uv_of(p);
                if (!in_box(box, uv[0], uv[1]))
                    continue;
                candidates.push_back(refine_cs(c, a, b, s, box, t, std::clamp(uv[0], box[0], box[1]), std::clamp(uv[1], box[2], box[3]), tol));
            }
        }
        else
        {
            using CPiece = CurvePiece<T, 3>;
            using SPiece = SurfacePiece<T>;
            CurvePieces<T, 3> CP{c};
            SurfacePieces<T> SP{s};
            std::vector<std::pair<CPiece, SPiece>> stack;
            const auto spieces = SP.pieces(box);
            for (auto &p : CP.pieces(a, b))
                for (const auto &q : spieces)
                    if (p.box.overlaps(q.box, tol))
                        stack.emplace_back(p, q);

            // a piece of the curve lying on the surface patch: a seed of portion
            auto try_portion = [&](const CPiece &p, const SPiece &q) -> bool {
                if (p.size() <= small || !q.box.contains(p.box, tol))
                    return false;
                constexpr int n = 4;
                std::array<std::array<T, 3>, n + 1> uvd{};
                for (int k = 0; k <= n; ++k)
                {
                    const int kk = k == 0 ? 0 : (k == 1 ? n : k - 1);
                    const auto t = kk == n ? p.b : p.a + (p.b - p.a) * T(kk) / T(n);
                    uvd[kk] = project_on_surface(s, q.uv, c.value(t));
                    if (uvd[kk][2] > tol)
                        return false;
                }
                seeds.push_back({p.a, p.b, {uvd[0][0], uvd[0][1]}, {uvd[n][0], uvd[n][1]}});
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
                // a flat piece of curve is a thin tube around its chord
                if (p.flatness <= T(1e-3) * p.chord() && !segment_overlaps(p.first, p.last, q.box, p.flatness + tol))
                    continue;
                if (try_portion(p, q))
                    continue;
                const bool p_flat = p.flatness <= T(1e-3) * p.chord() || p.size() <= small;
                const bool q_flat = q.flatness <= T(1e-3) * q.size() || q.size() <= small;
                const bool tiny = p.size() <= small && q.size() <= small;
                const auto quv = q.center();
                if (p_flat && q_flat)
                {
                    // the chord against the plane of the patch
                    const auto c0 = s.value(quv[0], quv[1]);
                    auto n = cross(s.value(quv[0], quv[1], 1, 0), s.value(quv[0], quv[1], 0, 1));
                    const auto nn = norm(n);
                    const auto ch = p.last - p.first;
                    const T sn = nn > T(0) && norm(ch) > T(0) ? std::abs(ch * n) / (nn * norm(ch)) : T(0);
                    if (tiny || sn > T(0.1))
                    {
                        T lambda = T(0.5);
                        if (nn > T(0) && std::abs(ch * n) > T(0))
                            lambda = std::clamp(((c0 - p.first) * n) / (ch * n), T(0), T(1));
                        candidates.push_back(refine_cs(c, a, b, s, box, p.a + lambda * (p.b - p.a), quv[0], quv[1], tol));
                        continue;
                    }
                }
                if (tiny)
                {
                    candidates.push_back(refine_cs(c, a, b, s, box, (p.a + p.b) / T(2), quv[0], quv[1], tol));
                    continue;
                }
                const bool split_p = q_flat && !p_flat ? true : (p_flat && !q_flat ? false : p.size() >= q.size());
                if (split_p)
                {
                    auto [l, r] = CP.split(p);
                    stack.emplace_back(std::move(l), q);
                    stack.emplace_back(std::move(r), std::move(q));
                }
                else
                {
                    auto [l, r] = SP.split(q, SurfacePieces<T>::longer_along_u(q));
                    stack.emplace_back(p, std::move(l));
                    stack.emplace_back(std::move(p), std::move(r));
                }
            }
        }

        // ---- portions: merge the seeds that follow each other along the curve
        std::ranges::sort(seeds, {}, &Seed::t0);
        const T gap = T(5) * small; // pieces near a patch border stay unseeded up to twice `small`
        std::vector<Seed> merged;
        for (const auto &sd : seeds)
        {
            if (!merged.empty() && (sd.t0 <= merged.back().t1 || norm(c.value(sd.t0) - c.value(merged.back().t1)) <= gap))
            {
                if (sd.t1 > merged.back().t1)
                {
                    merged.back().t1 = sd.t1;
                    merged.back().uv1 = sd.uv1;
                }
                continue;
            }
            merged.push_back(sd);
        }
        // a tangency keeps the curve within the tolerance over a short length: the distance varies
        // along it, from about zero at the contact to the tolerance at the ends (see curves.h)
        std::erase_if(merged, [&](const Seed &m) {
            constexpr int n = 8;
            T d_min = std::numeric_limits<T>::max(), d_max{};
            std::array<T, 3> best{m.t0, m.uv0[0], m.uv0[1]};
            for (int k = 1; k < n; ++k)
            {
                const T t = m.t0 + (m.t1 - m.t0) * T(k) / T(n);
                const auto [u, v, d] = project_on_surface(s, box, c.value(t));
                if (d < d_min)
                {
                    d_min = d;
                    best = {t, u, v};
                }
                d_max = std::max(d_max, d);
            }
            if (d_max - d_min <= tol / T(4))
                return false;
            candidates.push_back(refine_cs(c, a, b, s, box, best[0], best[1], best[2], tol));
            return true;
        });

        // ---- crossings: within the tolerance, sorted, without duplicates
        std::vector<CurveSurfaceHit<T>> hits;
        for (const auto &[t, u, v, d] : candidates)
            if (d <= tol && in_box(box, u, v))
                hits.push_back({t, u, v, (c.value(t) + s.value(u, v)) / T(2), false});
        std::ranges::sort(hits, {}, &CurveSurfaceHit<T>::t);
        std::vector<CurveSurfaceHit<T>> unique;
        for (const auto &h : hits)
        {
            const bool dup = std::ranges::any_of(unique, [&](const auto &x) {
                return norm(x.p - h.p) <= small && norm(c.value((x.t + h.t) / T(2)) - h.p) <= small;
            });
            if (!dup)
                unique.push_back(h);
        }
        for (auto &m : merged)
        {
            const auto e0 = c.value(m.t0), e1 = c.value(m.t1);
            std::erase_if(unique, [&](const auto &h) {
                const bool inside = h.t >= m.t0 && h.t <= m.t1;
                const bool at_start = norm(h.p - e0) <= gap, at_end = norm(h.p - e1) <= gap;
                if (at_start && h.t < m.t0)
                {
                    m.t0 = h.t;
                    m.uv0 = {h.u, h.v};
                }
                if (at_end && h.t > m.t1)
                {
                    m.t1 = h.t;
                    m.uv1 = {h.u, h.v};
                }
                return inside || at_start || at_end;
            });
        }
        for (const auto &m : merged)
            if (norm(c.value(m.t1) - c.value(m.t0)) > small || norm(c.value((m.t0 + m.t1) / T(2)) - c.value(m.t0)) > small)
                res.portions.push_back({m.t0, m.t1, m.uv0[0], m.uv0[1], m.uv1[0], m.uv1[1]});
        // a closed curve crossing the surface at its start point: once
        if (unique.size() >= 2 && norm(c.value(a) - c.value(b)) <= tol && unique.front().t == a && unique.back().t == b &&
            norm(unique.front().p - unique.back().p) <= small)
            unique.pop_back();

        for (auto &h : unique)
        {
            const auto n = cross(s.value(h.u, h.v, 1, 0), s.value(h.u, h.v, 0, 1));
            const auto d = c.value(h.t, 1);
            const T nn = norm(n), nd = norm(d);
            h.tangent = nn > T(0) && nd > T(0) && std::abs(n * d) / (nn * nd) < o.angular_tol;
        }
        res.hits = std::move(unique);
        return res;
    }

    /// Every crossing of a curve with a surface, on their whole bounds.
    template <std::floating_point T>
    [[nodiscard]] auto intersect(const Curve<T, 3> &c, const Surface<T, 3> &s, IntersectOptions<T> o = {}) -> CurveSurfaceResult<T>
    {
        const auto b = c.bounds();
        return intersect(c, b[0], b[1], s, s.bounds(), o);
    }
}
