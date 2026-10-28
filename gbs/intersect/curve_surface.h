#pragma once
/**
 * @file curve_surface.h
 * @brief Intersection of a curve and a surface: every crossing and every portion of the curve
 * lying on the surface, on bounded parameter ranges.
 *
 * Closed forms: a line (`LineCurve`, or a B-spline of degree 1 with two poles) against a plane,
 * a cylinder, a cone or a sphere is a quadratic equation, against a torus a quartic; a circle or
 * an ellipse against a plane is a trigonometric equation. Any other pair, and the curves that
 * stay within the tolerance of the surface (lying on it), go through subdivision: pieces of the
 * curve against patches of the surface (`CurvePieces`, `SurfacePieces`), then Gauss-Newton on
 * (t, u, v).
 *
 * A stretch of the curve within the tolerance of the surface is a portion, or a contact (a
 * tangency, a crossing at a small angle) when the distance grows from its minimum to the
 * tolerance towards the ends where the curve leaves it: one crossing, at the minimum.
 */
#include <gbs/intersect/curves.h>
#include <gbs/intersect/polynomial.h>
#include <gbs/intersect/surface_pieces.h>
#include <gbs/primitives.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <numbers>
#include <optional>
#include <tuple>
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
        bool tangent{}; ///< the curve is tangent to the surface at p, or grazes it
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

        /// Is the parametrization singular at (u, v): a pole or an apex, where Su or Sv vanishes?
        template <std::floating_point T>
        bool singular_at(const Surface<T, 3> &s, T u, T v)
        {
            const auto su = s.value(u, v, 1, 0), sv = s.value(u, v, 0, 1);
            return norm(cross(su, sv)) <= T(1e-6) * (sq_norm(su) + sq_norm(sv));
        }

        /**
         * `refine_cs`, started again around a singular point of the parametrization where it
         * stops away from the surface: at a pole Su = 0, the steps cannot turn u towards the
         * crossing (nor can a projection). The parameter that vanishes there is sampled over the
         * box, the other one moved off the singular line on both sides, by about the distance left
         * (the crossing is about as near; farther, the steps fall back on the line).
         */
        template <std::floating_point T>
        auto refine_cs_restart(const Curve<T, 3> &c, T a, T b, const Surface<T, 3> &s, const std::array<T, 4> &box, T t, T u, T v, T tol)
            -> std::array<T, 4>
        {
            auto r = refine_cs(c, a, b, s, box, t, u, v, tol);
            if (r[3] <= tol || !singular_at(s, r[1], r[2]))
                return r;
            const auto su = s.value(r[1], r[2], 1, 0), sv = s.value(r[1], r[2], 0, 1);
            const bool along_u = sq_norm(su) <= sq_norm(sv); // Su vanishes
            const int i = along_u ? 0 : 2, j = along_u ? 2 : 0; // sampled, moved
            const T speed = norm(along_u ? sv : su);
            const T off = std::min(speed > T(0) ? T(2) * r[3] / speed : T(0), T(1e-3) * (box[j + 1] - box[j]));
            constexpr int n = 8;
            for (int k = 0; k < n && r[3] > tol; ++k)
                for (T side : {T(-1), T(1)})
                {
                    std::array<T, 2> x{};
                    x[i / 2] = box[i] + (box[i + 1] - box[i]) * (T(k) + T(0.5)) / T(n);
                    x[j / 2] = r[1 + j / 2] + side * off;
                    if (x[j / 2] < box[j] || x[j / 2] > box[j + 1])
                        continue;
                    if (const auto r2 = refine_cs(c, a, b, s, box, r[0], x[0], x[1], tol); r2[3] < r[3])
                        r = r2;
                }
            return r;
        }

        /// The line carried by a curve, if it is one: c(t) = o + t d.
        template <std::floating_point T>
        auto as_line(const Curve<T, 3> &c) -> std::optional<std::array<point<T, 3>, 2>>
        {
            const auto *bs = dynamic_cast<const BSCurve<T, 3> *>(&c);
            if (!dynamic_cast<const LineCurve<T> *>(&c) && !(bs && bs->degree() == 1 && bs->poles().size() == 2))
                return std::nullopt;
            const T a = c.bounds()[0];
            const auto d = c.value(a, 1);
            return std::array<point<T, 3>, 2>{c.value(a) - a * d, d};
        }

        /// u (and v for a torus) of an angle, in the turn of the box.
        template <std::floating_point T>
        T box_turn(T x, T lo, T hi) { return primitives_detail::wrap(x, T(0.5) * (lo + hi) - std::numbers::pi_v<T>); }

        /// Parameters of the angles x + 2k pi in [a, b] (ends pulled in from a rounding margin).
        template <std::floating_point T>
        void push_turns(std::vector<T> &ts, T x, T a, T b)
        {
            constexpr T two_pi = T(2) * std::numbers::pi_v<T>;
            const T e = T(1e-12) * (T(1) + std::abs(a) + std::abs(b));
            for (T k = std::floor((a - x) / two_pi); x + k * two_pi <= b + e; ++k)
                if (const T t = x + k * two_pi; t >= a - e)
                    ts.push_back(std::clamp(t, a, b));
        }

        /**
         * Candidate parameters of the curve from a closed form: the roots, the extremes of the
         * equation (tangencies, near misses) and the ends of the range within the tolerance (a
         * crossing there may be missed by the rounding). nullopt when there is none for the pair,
         * or the curve stays within the tolerance of the surface at a few points (it lies on it):
         * the subdivision gives the portions and their clipping by the box.
         */
        template <std::floating_point T>
        auto closed_form(const Curve<T, 3> &c, T a, T b, const Surface<T, 3> &s, T tol) -> std::optional<std::vector<T>>
        {
            const auto *es = dynamic_cast<const ElementarySurface<T> *>(&s);
            if (!es || !(b > a))
                return std::nullopt;
            const auto &f = es->frame();
            if (const auto line = as_line(c))
            {
                const auto d = (*line)[1];
                const T nd = norm(d);
                if (!(nd > T(0)))
                    return std::nullopt;
                // by the arc length s from the foot of the frame origin, c = O + q + s e with q
                // orthogonal to e: the coefficients do not depend on the parameter origin of the
                // line (far away, or knots at 1e7) and do not cancel
                const auto e = d / nd;
                const auto pa = c.value(a);
                const T sf = (f.o - pa) * e;
                const auto q = pa + sf * e - f.o;
                const T s0 = -sf, s1 = (b - a) * nd - sf;
                const T qx = q * f.x, qy = q * f.y, qz = q * f.z, ex = e * f.x, ey = e * f.y, ez = e * f.z;
                // radial part: rho² = A s² + B s + C
                const T A = ex * ex + ey * ey, B = T(2) * (qx * ex + qy * ey), C = qx * qx + qy * qy;
                auto rho = [&](T x) { return std::hypot(qx + x * ex, qy + x * ey); };
                auto z = [&](T x) { return qz + x * ez; };
                std::vector<T> coef;
                std::function<T(T)> dist; // distance of c to the surface at s
                if (dynamic_cast<const Plane<T> *>(es))
                {
                    coef = {qz, ez};
                    dist = [&](T x) { return std::abs(z(x)); };
                }
                else if (const auto *cy = dynamic_cast<const CylindricalSurface<T> *>(es))
                {
                    const T R = cy->radius();
                    coef = {C - R * R, B, A};
                    dist = [&, R](T x) { return std::abs(rho(x) - R); };
                }
                else if (const auto *co = dynamic_cast<const ConicalSurface<T> *>(es))
                {
                    // rho² = (R + z tan a)², z = qz + s ez: both nappes
                    const T R = co->radius(), k = std::tan(co->semi_angle()), r0 = R + k * qz, r1 = k * ez;
                    coef = {C - r0 * r0, B - T(2) * r0 * r1, A - r1 * r1};
                    dist = [&, R, k](T x) { return std::abs(rho(x) - std::abs(R + k * z(x))) / std::sqrt(T(1) + k * k); };
                }
                else if (const auto *sp = dynamic_cast<const SphericalSurface<T> *>(es))
                {
                    const T R = sp->radius();
                    coef = {C + qz * qz - R * R, B + T(2) * qz * ez, A + ez * ez};
                    dist = [&, R](T x) { return std::abs(std::hypot(rho(x), z(x)) - R); };
                }
                else if (const auto *to = dynamic_cast<const ToroidalSurface<T> *>(es))
                {
                    // (|x|² + R² - r²)² = 4 R² rho²
                    const T R = to->major_radius(), r = to->minor_radius();
                    const T a2 = A + ez * ez, a1 = B + T(2) * qz * ez, a0 = C + qz * qz + R * R - r * r;
                    const T R4 = T(4) * R * R;
                    coef = {a0 * a0 - R4 * C, T(2) * a0 * a1 - R4 * B, a1 * a1 + T(2) * a0 * a2 - R4 * A, T(2) * a1 * a2, a2 * a2};
                    dist = [&, R, r](T x) { return std::abs(std::hypot(rho(x) - R, z(x)) - r); };
                }
                else
                    return std::nullopt;
                constexpr int n = 4;
                bool on = true;
                for (int k = 0; k <= n && on; ++k)
                    on = dist(s0 + (s1 - s0) * T(k) / T(n)) <= tol;
                if (on)
                    return std::nullopt;
                const auto r = polynomial_roots(coef, s0, s1);
                std::vector<T> ts;
                for (const auto &xs : {r.roots, r.extremes})
                    for (T x : xs)
                        ts.push_back(std::clamp(a + (x - s0) / nd, a, b));
                if (dist(s0) <= tol)
                    ts.push_back(a);
                if (dist(s1) <= tol)
                    ts.push_back(b);
                return ts;
            }
            if (const auto *el = dynamic_cast<const EllipseCurve<T> *>(&c); el && dynamic_cast<const Plane<T> *>(es))
            {
                // (C(t) - O).z = d0 + d1 cos t + d2 sin t = d0 + rho cos(t - phi)
                const auto &g = el->frame();
                const T d0 = (g.o - f.o) * f.z, d1 = el->radius1() * (g.x * f.z), d2 = el->radius2() * (g.y * f.z);
                const T rho = std::hypot(d1, d2);
                if (std::abs(d0) + rho <= tol)
                    return std::nullopt; // in the plane within the tolerance
                std::vector<T> ts;
                if (!(rho > T(0)))
                    return ts; // parallel to the plane, farther than the tolerance
                const T phi = std::atan2(d2, d1);
                // the extremes phi and phi + pi: a tangency, or a near miss within the tolerance
                std::vector<T> xs{phi, phi + std::numbers::pi_v<T>};
                if (const T x = -d0 / rho; std::abs(x) <= T(1))
                {
                    const T w = std::acos(x);
                    xs.insert(xs.end(), {phi + w, phi - w});
                }
                for (T x : xs)
                    push_turns(ts, x, a, b);
                auto dist = [&](T t) { return std::abs(d0 + d1 * std::cos(t) + d2 * std::sin(t)); };
                if (dist(a) <= tol)
                    ts.push_back(a);
                if (dist(b) <= tol)
                    ts.push_back(b);
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
     * closest (within `tol`); a tangency and the crossings next to it, the curve within the
     * tolerance between them, are one grazing contact, reported tangent. A stretch longer than
     * about ten times the tolerance along which the curve stays within it is a portion, and the
     * crossings on it are dropped; unless the distance grows from its minimum to the tolerance
     * towards the ends where the curve leaves it (a tangency, a crossing at a small angle): one
     * crossing then, at the minimum.
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
        auto refine = [&](T t, T u, T v) { return refine_cs_restart(c, a, b, s, box, t, u, v, tol); };
        // projection of c(t) on the surface on the box, by Newton from (u, v): u, v, distance
        auto near = [&](T t, const std::array<T, 2> &uv) { return newton_point_surface(s, c.value(t), uv[0], uv[1], box, T(0)); };
        // does the curve stay within the tolerance between t0 and t1 (at a few points)?
        auto within = [&](T t0, T t1, std::array<T, 2> uv) {
            constexpr int n = 4;
            for (int k = 1; k < n; ++k)
            {
                const auto [u, v, d] = near(t0 + (t1 - t0) * T(k) / T(n), uv);
                if (d > tol)
                    return false;
                uv = {u, v};
            }
            return true;
        };

        if (const auto ts = closed_form(c, a, b, s, tol))
        {
            // parameters of the point, clamped into the box: at a pole or an apex u is arbitrary,
            // the refinement and the tolerance decide
            const auto &es = dynamic_cast<const ElementarySurface<T> &>(s);
            const bool angle_u = !dynamic_cast<const Plane<T> *>(&es), torus = dynamic_cast<const ToroidalSurface<T> *>(&es) != nullptr,
                       cone = dynamic_cast<const ConicalSurface<T> *>(&es) != nullptr;
            for (T t : *ts)
            {
                const auto p = c.value(t);
                auto uv = es.parameters(p);
                if (angle_u)
                    uv[0] = box_turn(uv[0], box[0], box[1]);
                if (torus)
                    uv[1] = box_turn(uv[1], box[2], box[3]);
                if (cone) // on the other nappe the radial direction is reversed: of u and u + pi, the nearer
                    if (const T w = box_turn(uv[0] + std::numbers::pi_v<T>, box[0], box[1]); norm(es.value(w, uv[1]) - p) < norm(es.value(uv[0], uv[1]) - p))
                        uv[0] = w;
                candidates.push_back(refine(t, std::clamp(uv[0], box[0], box[1]), std::clamp(uv[1], box[2], box[3])));
            }
        }
        else
        {
            using CPiece = CurvePiece<T, 3>;
            using SPiece = SurfacePiece<T>;
            CurvePieces<T, 3> CP{c};
            SurfacePieces<T> SP{s};
            // pieces apart: their boxes, the tube of a flat piece of curve against the box of the
            // patch, the chord against the slab of the plane of the patch (axis-aligned boxes do not
            // follow a curve lying on a tilted or curved surface)
            auto apart = [&](const CPiece &p, const SPiece &q) {
                return !p.box.overlaps(q.box, tol) ||
                       (p.flatness <= T(1e-3) * p.chord() && !segment_overlaps(p.first, p.last, q.box, p.flatness + tol)) ||
                       projection_detail::lower_bound(p, q) > tol;
            };
            std::vector<std::pair<CPiece, SPiece>> stack;
            const auto spieces = SP.pieces(box);
            for (auto &p : CP.pieces(a, b))
                for (const auto &q : spieces)
                    if (!apart(p, q))
                        stack.emplace_back(p, q);

            // a piece of the curve lying on the surface: a seed of portion. Its points are projected
            // on the whole box from the patch: the piece may run over the next patches.
            auto try_portion = [&](const CPiece &p, const SPiece &q) -> bool {
                if (p.size() <= small)
                    return false;
                constexpr int n = 4;
                std::array<std::array<T, 3>, n + 1> uvd{};
                for (int k = 0; k <= n; ++k)
                {
                    const int kk = k == 0 ? 0 : (k == 1 ? n : k - 1); // the ends first: they fail most often
                    const auto t = kk == n ? p.b : p.a + (p.b - p.a) * T(kk) / T(n);
                    const auto from = k <= 1 ? projection_detail::on_patch(q, c.value(t)) : std::array<T, 2>{uvd[kk - 1][0], uvd[kk - 1][1]};
                    uvd[kk] = near(t, from);
                    if (uvd[kk][2] > tol)
                        return false;
                }
                seeds.push_back({p.a, p.b, {uvd[0][0], uvd[0][1]}, {uvd[n][0], uvd[n][1]}});
                return true;
            };
            auto seeded = [&](const CPiece &p) { return std::ranges::any_of(seeds, [&](const Seed &sd) { return sd.t0 <= p.a && p.b <= sd.t1; }); };

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
                if (seeded(p) || apart(p, q) || try_portion(p, q))
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
                        candidates.push_back(refine(p.a + lambda * (p.b - p.a), quv[0], quv[1]));
                        continue;
                    }
                }
                if (tiny)
                {
                    candidates.push_back(refine((p.a + p.b) / T(2), quv[0], quv[1]));
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

        // ---- stretches within the tolerance: the seeds that overlap, or with the curve within the
        // tolerance between them (fragments of one stretch, left by the borders of the patches)
        std::ranges::sort(seeds, {}, &Seed::t0);
        std::vector<Seed> merged;
        for (const auto &sd : seeds)
        {
            if (!merged.empty() && (sd.t0 <= merged.back().t1 || within(merged.back().t1, sd.t0, merged.back().uv1)))
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
        // the end of a stretch, from t towards `end`: where the curve leaves the tolerance (steps
        // doubled from about `small`, then bisection), or `end` itself (open: the curve stays in it)
        auto extend = [&](T t, std::array<T, 2> uv, T end) -> std::tuple<T, std::array<T, 2>, bool> {
            const T speed = norm(c.value(t, 1));
            T h = speed > T(0) ? small / speed : std::abs(end - t) / T(16), lo = t, hi = end;
            for (;;)
            {
                const bool last = std::abs(end - lo) <= h;
                const T tn = last ? end : lo + std::copysign(h, end - lo);
                const auto [u, v, d] = near(tn, uv);
                if (d > tol)
                {
                    hi = tn;
                    break;
                }
                lo = tn;
                uv = {u, v};
                if (last)
                    return {end, uv, true};
                h *= T(2);
            }
            for (int it{}; it < 64 && norm(c.value(hi) - c.value(lo)) > tol; ++it)
            {
                const T m = (lo + hi) / T(2);
                if (m == lo || m == hi)
                    break;
                const auto [u, v, d] = near(m, uv);
                if (d > tol)
                    hi = m;
                else
                {
                    lo = m;
                    uv = {u, v};
                }
            }
            return {lo, uv, false};
        };
        // A stretch is a contact when the distance stays near its minimum (within an eighth of the
        // way to the tolerance) on less than half of it: along a tangency it grows like a parabola,
        // faster along a crossing at a small angle, and in both ways when the curve dips under the
        // surface by less than the tolerance. Along a portion the curve stays close to the surface,
        // and a stretch that the curve does not leave (open at both ends) is a portion.
        std::vector<std::array<T, 2>> portion_spans, contact_spans; // t0, t1
        std::vector<std::array<T, 4>> contacts;                     // t, u, v, distance
        for (const auto &m : merged)
        {
            const auto [t0, uv0, open0] = extend(m.t0, m.uv0, a);
            const auto [t1, uv1, open1] = extend(m.t1, m.uv1, b);
            constexpr int n = 16;
            std::array<T, n + 1> d{};
            std::array<T, 4> best{t0, uv0[0], uv0[1], std::numeric_limits<T>::max()};
            std::array<T, 2> uv = uv0;
            for (int k = 0; k <= n; ++k)
            {
                const T t = k == n ? t1 : t0 + (t1 - t0) * T(k) / T(n);
                const auto r = near(t, uv);
                uv = {r[0], r[1]};
                d[k] = r[2];
                if (r[2] < best[3])
                    best = {t, r[0], r[1], r[2]};
            }
            // the minimum: Newton from the nearest sample
            if (const auto mn = newton_curve_surface(c, s, best[0], best[1], best[2], {t0, t1}, box, T(0)); mn[3] < best[3])
                best = mn;
            const auto close = std::ranges::count_if(d, [&](T x) { return x - best[3] <= (tol - best[3]) / T(8); });
            const bool longer = norm(c.value(t1) - c.value(t0)) > small || norm(c.value((t0 + t1) / T(2)) - c.value(t0)) > small;
            if (!longer || ((!open0 || !open1) && 2 * close < n + 1))
            {
                contact_spans.push_back({t0, t1});
                contacts.push_back(best);
            }
            else
            {
                portion_spans.push_back({t0, t1});
                res.portions.push_back({t0, t1, uv0[0], uv0[1], uv1[0], uv1[1]});
            }
        }

        // ---- crossings: within the tolerance and outside the portions, sorted; on a contact, its
        // minimum and the crossings found exactly (two steep crossings near an apex stay apart)
        auto in = [&](const auto &spans, T t) { return std::ranges::any_of(spans, [&](const auto &sp) { return t >= sp[0] && t <= sp[1]; }); };
        std::vector<std::pair<CurveSurfaceHit<T>, T>> hits; // with their distance
        auto add = [&](T t, T u, T v, T d) {
            const auto n = cross(s.value(u, v, 1, 0), s.value(u, v, 0, 1));
            const auto dc = c.value(t, 1);
            const T nn = norm(n), nd = norm(dc);
            const bool tangent = nn > T(0) && nd > T(0) && std::abs(n * dc) / (nn * nd) < o.angular_tol;
            hits.push_back({{t, u, v, (c.value(t) + s.value(u, v)) / T(2), tangent}, d});
        };
        for (const auto &[t, u, v, d] : candidates)
            if (d <= tol && in_box(box, u, v) && !in(portion_spans, t) && (d <= T(1e-3) * tol || !in(contact_spans, t)))
                add(t, u, v, d);
        for (const auto &[t, u, v, d] : contacts)
            if (d <= tol)
                add(t, u, v, d);
        std::ranges::sort(hits, {}, [](const auto &h) { return h.first.t; });
        // the same crossing found several times, or a tangency with the crossings next to it, the
        // curve within the tolerance between them (a grazing contact): once, tangent, at the nearest
        std::vector<std::pair<CurveSurfaceHit<T>, T>> unique;
        for (const auto &h : hits)
        {
            if (!unique.empty())
            {
                const auto &[x, dx] = unique.back();
                const bool same = norm(x.p - h.first.p) <= tol;
                if (same || ((x.tangent || h.first.tangent) && within(x.t, h.first.t, {x.u, x.v})))
                {
                    const bool tangent = x.tangent || h.first.tangent;
                    if (h.second < dx)
                        unique.back() = h;
                    unique.back().first.tangent = tangent || !same;
                    continue;
                }
            }
            unique.push_back(h);
        }
        // a closed curve crossing the surface at its start point: once
        if (unique.size() >= 2 && norm(c.value(a) - c.value(b)) <= tol && norm(unique.front().first.p - unique.back().first.p) <= small &&
            norm(unique.front().first.p - c.value(a)) <= small)
            unique.pop_back();
        for (const auto &h : unique)
            res.hits.push_back(h.first);
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
