#pragma once
/**
 * @file surface_surface.h
 * @brief Intersection of two surfaces on bounded parameter boxes: curves with a 3D curve and a
 * pcurve on each surface on the same parameter, contact points, coincidence.
 *
 * Pairs of primitives with a closed form (`analytic_surfaces.h`) give exact lines, circles and
 * ellipses, restricted to the parts lying in both boxes and cut where a periodic parameter
 * crosses the seam of a box. A pcurve is exact (a segment) when the parameters of the surface
 * vary linearly along the curve, as for a circle on a coaxial surface or a generatrix, and is
 * interpolated otherwise, with its deviation measured. The other pairs go through marching
 * (`marching.h`).
 */
#include <gbs/intersect/analytic_surfaces.h>
#include <gbs/intersect/intersection_curve.h>
#include <gbs/intersect/marching.h>
#include <gbs/intersect/curve_surface.h>
#include <gbs/bscinterp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <tuple>
#include <numbers>
#include <optional>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs
{
    namespace intersect_detail
    {
        /// Parameters of a point of a primitive, angles in the turn of the box.
        template <std::floating_point T>
        auto primitive_uv(const ElementarySurface<T> &s, const std::array<T, 4> &box, const point<T, 3> &p) -> std::array<T, 2>
        {
            auto uv = s.parameters(p);
            if (!dynamic_cast<const Plane<T> *>(&s))
                uv[0] = box_turn(uv[0], box[0], box[1]);
            if (dynamic_cast<const ToroidalSurface<T> *>(&s))
                uv[1] = box_turn(uv[1], box[2], box[3]);
            return uv;
        }

        template <std::floating_point T>
        auto periodic(const ElementarySurface<T> &s) -> std::array<bool, 2>
        {
            return {dynamic_cast<const Plane<T> *>(&s) == nullptr, dynamic_cast<const ToroidalSurface<T> *>(&s) != nullptr};
        }

        /// Parameter range of a line inside a box, or nullopt.
        template <std::floating_point T>
        auto clip_line(const LineCurve<T> &l, const Box<T, 3> &b) -> std::optional<std::array<T, 2>>
        {
            T lo = std::numeric_limits<T>::lowest(), hi = std::numeric_limits<T>::max();
            const auto &o = l.origin(), &d = l.vector();
            for (int i = 0; i < 3; ++i)
            {
                if (d[i] == T(0))
                {
                    if (o[i] < b.min[i] || o[i] > b.max[i])
                        return std::nullopt;
                    continue;
                }
                T t1 = (b.min[i] - o[i]) / d[i], t2 = (b.max[i] - o[i]) / d[i];
                if (t1 > t2)
                    std::swap(t1, t2);
                lo = std::max(lo, t1);
                hi = std::min(hi, t2);
            }
            if (!(lo < hi))
                return std::nullopt;
            return std::array<T, 2>{lo, hi};
        }

        /// Pcurve of c on [a, b] on a primitive: a segment when its parameters vary linearly, interpolated otherwise.
        template <std::floating_point T>
        auto primitive_pcurve(const Curve<T, 3> &c, T a, T b, const ElementarySurface<T> &s, const std::array<T, 4> &box, T tol)
            -> std::tuple<std::shared_ptr<Curve<T, 2>>, T, bool>
        {
            const auto per = periodic(s);
            const T pi = std::numbers::pi_v<T>;
            // parameters along the curve, unwrapped from the middle outwards (ends may sit on the seam)
            auto sample = [&](std::size_t n) {
                std::vector<T> t(n + 1);
                std::vector<std::array<T, 2>> uv(n + 1);
                for (std::size_t k = 0; k <= n; ++k)
                    t[k] = k == n ? b : a + (b - a) * T(k) / T(n);
                const std::size_t mid = n / 2;
                uv[mid] = primitive_uv(s, box, c.value(t[mid]));
                auto follow = [&](std::size_t k, std::size_t ref) {
                    auto x = primitive_uv(s, box, c.value(t[k]));
                    for (int i = 0; i < 2; ++i)
                        if (per[i])
                            x[i] = uv[ref][i] + std::remainder(x[i] - uv[ref][i], T(2) * pi);
                    uv[k] = x;
                };
                for (std::size_t k = mid + 1; k <= n; ++k)
                    follow(k, k - 1);
                for (std::size_t k = mid; k-- > 0;)
                    follow(k, k + 1);
                // at a pole or an apex the angle is undefined: it is extrapolated from the
                // neighbouring samples (a generatrix keeps its angle up to the apex)
                auto singular = [&](std::size_t k) {
                    return norm(s.value(uv[k][0], uv[k][1], 1, 0)) <= T(1e-8) * norm(s.value(uv[k][0], uv[k][1], 0, 1));
                };
                for (std::size_t k = 0; k <= n; ++k)
                {
                    if (!singular(k))
                        continue;
                    if (k + 2 <= n && !singular(k + 1) && !singular(k + 2))
                        uv[k][0] = T(2) * uv[k + 1][0] - uv[k + 2][0];
                    else if (k >= 2 && !singular(k - 1) && !singular(k - 2))
                        uv[k][0] = T(2) * uv[k - 1][0] - uv[k - 2][0];
                }
                return std::pair{t, uv};
            };
            auto deviation = [&](const Curve<T, 2> &pc) {
                T d{};
                constexpr int m = 1000;
                for (int k = 0; k < m; ++k)
                {
                    const T t = a + (b - a) * (T(k) + T(0.5)) / T(m);
                    const auto q = pc.value(t);
                    d = std::max(d, norm(s.value(q[0], q[1]) - c.value(t)));
                }
                return d;
            };
            {
                const auto [t, uv] = sample(8);
                T err{}, scale = T(1);
                for (std::size_t k = 0; k <= 8; ++k)
                {
                    const T l = (t[k] - a) / (b - a);
                    for (int i = 0; i < 2; ++i)
                    {
                        err = std::max(err, std::abs(uv[k][i] - (uv[0][i] + l * (uv[8][i] - uv[0][i]))));
                        scale = std::max(scale, std::abs(uv[k][i]));
                    }
                }
                if (err <= T(1e-10) * scale)
                {
                    auto pc = std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{uv[0], uv[8]}, std::vector<T>{a, a, b, b}, 1);
                    return {pc, deviation(*pc), true};
                }
            }
            std::shared_ptr<Curve<T, 2>> pc;
            T dev{};
            for (std::size_t n = 16; n <= 2048; n *= 2)
            {
                const auto [t, uv] = sample(n);
                std::vector<constrType<T, 2, 1>> q(uv.size());
                for (std::size_t k = 0; k < uv.size(); ++k)
                    q[k] = {uv[k]};
                pc = std::make_shared<BSCurve<T, 2>>(interpolate<T, 2>(q, t, std::size_t{3}));
                dev = deviation(*pc);
                if (dev <= tol / T(10))
                    break;
            }
            return {pc, dev, false};
        }

        /**
         * Parts of a carrier inside both boxes, cut where a periodic parameter crosses the seam
         * of its box. The pcurve of the carrier on each surface (periodic parameters unwrapped) is
         * crossed with the lines of the bounds of the box, repeated every turn for a periodic
         * parameter (curve-curve intersection: all the crossings, however short the part between
         * them); the pieces between the crossings are kept when their middle is in both boxes,
         * and joined again where no periodic parameter jumps.
         */
        template <std::floating_point T>
        auto inside_intervals(const Curve<T, 3> &c, T ta, T tb, const ElementarySurface<T> &s1, const std::array<T, 4> &b1,
                              const ElementarySurface<T> &s2, const std::array<T, 4> &b2, bool full_turn, T tol)
            -> std::vector<std::array<T, 2>>
        {
            const T two_pi = T(2) * std::numbers::pi_v<T>;
            const auto per1 = periodic(s1), per2 = periodic(s2);
            auto uv = [&](T t) {
                const auto p = c.value(t);
                const auto a = primitive_uv(s1, b1, p), b = primitive_uv(s2, b2, p);
                return std::array<T, 4>{a[0], a[1], b[0], b[1]};
            };
            auto inside = [&](T t) {
                const auto x = uv(t);
                return in_box(b1, x[0], x[1]) && in_box(b2, x[2], x[3]);
            };
            const std::array<bool, 4> per{per1[0], per1[1], per2[0], per2[1]};
            const std::array<std::array<T, 2>, 4> lim{{{b1[0], b1[1]}, {b1[2], b1[3]}, {b2[0], b2[1]}, {b2[2], b2[3]}}};
            auto on_end = [&](T x, int i) {
                const T eps = T(1e-9) * (T(1) + std::abs(lim[i][0]) + std::abs(lim[i][1]));
                return std::abs(x - lim[i][0]) <= eps || std::abs(x - lim[i][1]) <= eps;
            };
            // a periodic parameter jumping by more than half a turn between t and t': the seam
            // (a value on an end of the box is on the seam: it belongs to both sides)
            auto jump = [&](const std::array<T, 4> &x, const std::array<T, 4> &y) {
                for (int i = 0; i < 4; ++i)
                    if (per[i] && std::abs(x[i] - y[i]) > std::numbers::pi_v<T> && !on_end(x[i], i) && !on_end(y[i], i))
                        return true;
                return false;
            };

            // ---- crossings of the pcurves with the bounds of the boxes
            std::vector<T> cuts{ta, tb};
            for (int side = 0; side < 2; ++side)
            {
                const auto &s = side == 0 ? s1 : s2;
                const auto &b = side == 0 ? b1 : b2;
                const auto [pc, dev, exact] = primitive_pcurve(c, ta, tb, s, b, tol);
                if (!pc)
                    continue;
                Box<T, 2> ext;
                for (int k = 0; k <= 64; ++k)
                    ext.add(pc->value(ta + (tb - ta) * T(k) / T(64)));
                const T margin = T(1e-3) * (T(1) + ext.diagonal());
                IntersectOptions<T> o;
                o.tol = T(1e-10) * (T(1) + ext.diagonal());
                for (int i = 0; i < 2; ++i)
                {
                    const T lo = ext.min[i] - margin, hi = ext.max[i] + margin;
                    const T olo = ext.min[1 - i] - margin, ohi = ext.max[1 - i] + margin;
                    for (T bound : {b[2 * i], b[2 * i + 1]})
                    {
                        // every turn of a periodic parameter over the unwrapped range
                        const bool turns = per[2 * side + i];
                        const int k0 = turns ? static_cast<int>(std::floor((lo - bound) / two_pi)) : 0;
                        const int k1 = turns ? static_cast<int>(std::ceil((hi - bound) / two_pi)) : 0;
                        for (int k = k0; k <= k1; ++k)
                        {
                            const T x = bound + T(k) * two_pi;
                            if (x < lo || x > hi)
                                continue;
                            point<T, 2> p0, p1;
                            p0[i] = p1[i] = x;
                            p0[1 - i] = olo;
                            p1[1 - i] = ohi;
                            const BSCurve<T, 2> line{points_vector<T, 2>{p0, p1}, std::vector<T>{0, 0, 1, 1}, 1};
                            const auto r = intersect(*pc, ta, tb, line, T(0), T(1), o);
                            for (const auto &h : r.hits)
                                cuts.push_back(h.t1);
                            for (const auto &ov : r.overlaps)
                            {
                                cuts.push_back(ov.t1a);
                                cuts.push_back(ov.t1b);
                            }
                        }
                    }
                }
            }
            std::ranges::sort(cuts);
            const T eps_t = T(1e-12) * (T(1) + std::abs(ta) + std::abs(tb));
            cuts.erase(std::unique(cuts.begin(), cuts.end(), [&](T x, T y) { return y - x <= eps_t; }), cuts.end());

            // ---- the pieces inside both boxes, joined where nothing jumps
            std::vector<std::array<T, 2>> r;
            for (std::size_t k = 0; k + 1 < cuts.size(); ++k)
            {
                const T a = cuts[k], b = cuts[k + 1];
                if (!inside((a + b) / T(2)))
                    continue;
                const T d = T(1e-6) * (b - a);
                if (!r.empty() && r.back()[1] == a)
                {
                    const T dl = T(1e-6) * (r.back()[1] - r.back()[0]);
                    if (!jump(uv(a - dl), uv(a + d)))
                    {
                        r.back()[1] = b;
                        continue;
                    }
                }
                r.push_back({a, b});
            }
            // parts reduced to a point (the curve touching a corner of a box)
            std::erase_if(r, [&](const std::array<T, 2> &i) {
                const auto p = c.value(i[0]);
                return norm(c.value(i[1]) - p) <= T(10) * tol && norm(c.value((i[0] + i[1]) / T(2)) - p) <= T(10) * tol;
            });
            // a full turn cut by a box in its middle: the parts at both ends are one curve
            if (full_turn && r.size() >= 2 && r.front()[0] == ta && r.back()[1] == tb)
            {
                const T d = T(1e-6) * (tb - ta);
                if (!jump(uv(tb - d), uv(ta + d)))
                {
                    r.front() = {r.back()[0], r.front()[1] + (tb - ta)};
                    r.pop_back();
                }
            }
            return r;
        }

        /// The carrier restricted to [a, b]: a new primitive curve on that range.
        template <std::floating_point T>
        auto restricted(const ElementaryCurve<T> &c, T a, T b) -> std::shared_ptr<Curve<T, 3>>
        {
            if (const auto *l = dynamic_cast<const LineCurve<T> *>(&c))
                return std::make_shared<LineCurve<T>>(l->origin(), l->vector(), a, b);
            const auto &e = dynamic_cast<const EllipseCurve<T> &>(c);
            return std::make_shared<EllipseCurve<T>>(e.radius1(), e.radius2(), e.frame(), a, b);
        }
    }

    /**
     * @brief Intersection of s1 on `box1` and s2 on `box2` (each {u1, u2, v1, v2}).
     *
     * Pairs of primitives with a closed form give exact curves; the other pairs go through
     * marching (B-spline curves and pcurves within the tolerance). The curves are cut at the seams
     * of the boxes, so that each pcurve stays continuous inside its box.
     */
    template <std::floating_point T>
    [[nodiscard]] auto intersect(const Surface<T, 3> &s1, const std::array<T, 4> &box1, const Surface<T, 3> &s2, const std::array<T, 4> &box2,
                                 IntersectOptions<T> o = {}) -> SurfaceSurfaceResult<T>
    {
        using namespace intersect_detail;
        SurfaceSurfaceResult<T> res;
        const auto *e1 = dynamic_cast<const ElementarySurface<T> *>(&s1);
        const auto *e2 = dynamic_cast<const ElementarySurface<T> *>(&s2);
        const auto ai = e1 && e2 && o.closed_forms ? analytic_intersection(*e1, *e2, o.tol) : std::nullopt;
        if (!ai)
            return marching(s1, box1, s2, box2, o);
        if (ai->coincident)
        {
            res.coincident = true;
            return res;
        }
        // the region where both surfaces are: lines are clipped to it
        auto region = SurfacePieces<T>{s1}.pieces(box1).front().box;
        {
            Box<T, 3> b2;
            for (const auto &p : SurfacePieces<T>{s2}.pieces(box2))
                b2.add(p.box);
            for (int i = 0; i < 3; ++i)
            {
                region.min[i] = std::max(region.min[i], b2.min[i]);
                region.max[i] = std::min(region.max[i], b2.max[i]);
            }
            region.inflate(o.tol);
        }
        const T small = T(10) * o.tol;
        for (const auto &[carrier, tangent] : ai->carriers)
        {
            const auto [ta, tb] = carrier->bounds();
            std::array<T, 2> range{ta, tb};
            if (const auto *l = dynamic_cast<const LineCurve<T> *>(carrier.get()))
            {
                if (region.empty())
                    continue;
                const auto r = clip_line(*l, region);
                if (!r)
                    continue;
                range = *r;
            }
            const bool full_turn = !dynamic_cast<const LineCurve<T> *>(carrier.get());
            for (const auto &[a, b] : inside_intervals(*carrier, range[0], range[1], *e1, box1, *e2, box2, full_turn, o.tol))
            {
                if (norm(carrier->value(b) - carrier->value(a)) <= small && norm(carrier->value((a + b) / T(2)) - carrier->value(a)) <= small)
                    continue;
                IntersectionCurve<T> ic;
                ic.curve = restricted(*carrier, a, b);
                ic.a = a;
                ic.b = b;
                const auto [p1, d1, x1] = primitive_pcurve(*ic.curve, a, b, *e1, box1, o.tol);
                const auto [p2, d2, x2] = primitive_pcurve(*ic.curve, a, b, *e2, box2, o.tol);
                ic.pc1 = p1;
                ic.pc2 = p2;
                ic.deviation = std::max(d1, d2);
                ic.exact = x1 && x2;
                ic.closed = full_turn && b - a >= T(2) * std::numbers::pi_v<T> * (T(1) - T(1e-12));
                ic.tangent = tangent;
                res.curves.push_back(std::move(ic));
            }
        }
        for (const auto &p : ai->points)
        {
            const auto a = primitive_uv(*e1, box1, p), b = primitive_uv(*e2, box2, p);
            if (in_box(box1, a[0], a[1]) && in_box(box2, b[0], b[1]))
                res.points.push_back({a[0], a[1], b[0], b[1], p});
        }
        return res;
    }

    /// Intersection of two surfaces on their whole bounds.
    template <std::floating_point T>
    [[nodiscard]] auto intersect(const Surface<T, 3> &s1, const Surface<T, 3> &s2, IntersectOptions<T> o = {}) -> SurfaceSurfaceResult<T>
    {
        return intersect(s1, s1.bounds(), s2, s2.bounds(), o);
    }
}
