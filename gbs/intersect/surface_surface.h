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
#include <gbs/intersect/curve_pieces.h>
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
            -> PCurveFit<T>
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
         * Parameters in [a, b] where the coordinate i of the 2D B-spline pc equals x: Bézier
         * pieces of that coordinate (`CurvePieces` in one dimension) halved until they no longer
         * reach x or are monotone, then bisection. A span where the coordinate stays at x within
         * eps gives its two ends.
         */
        template <std::floating_point T>
        void coordinate_roots(const BSCurve<T, 2> &pc, int i, T x, T a, T b, T eps, std::vector<T> &roots)
        {
            points_vector<T, 1> poles;
            for (const auto &q : pc.poles())
                poles.push_back({q[static_cast<std::size_t>(i)] - x});
            const BSCurve<T, 1> f{poles, pc.knotsFlats(), pc.degree()};
            CurvePieces<T, 1> pieces{f};
            auto stack = pieces.pieces(a, b);
            const T t_min = T(1e-13) * (T(1) + std::abs(a) + std::abs(b));
            auto value = [&](T t) { return f.value(t)[0]; };
            while (!stack.empty())
            {
                auto pce = std::move(stack.back());
                stack.pop_back();
                if (pce.box.min[0] > eps || pce.box.max[0] < -eps)
                    continue; // does not reach x
                if (pce.box.max[0] - pce.box.min[0] <= eps)
                {
                    roots.push_back(pce.a); // along x over the span
                    roots.push_back(pce.b);
                    continue;
                }
                bool up = true, down = true;
                for (std::size_t k = 1; k < pce.poles.size(); ++k)
                {
                    const T d = pce.poles[k][0] / pce.poles[k][1] - pce.poles[k - 1][0] / pce.poles[k - 1][1];
                    up = up && d >= T(0);
                    down = down && d <= T(0);
                }
                if (up || down)
                {
                    // monotone: one root at most, by bisection on the sign
                    T lo = pce.a, hi = pce.b, flo = value(lo), fhi = value(hi);
                    if (std::abs(flo) <= eps)
                        roots.push_back(lo);
                    if (std::abs(fhi) <= eps)
                        roots.push_back(hi);
                    if ((flo < T(0)) == (fhi < T(0)))
                        continue;
                    for (int it = 0; it < 200 && hi - lo > t_min; ++it)
                    {
                        const T m = (lo + hi) / T(2), fm = value(m);
                        if ((fm < T(0)) == (flo < T(0)))
                        {
                            lo = m;
                            flo = fm;
                        }
                        else
                        {
                            hi = m;
                            fhi = fm;
                        }
                    }
                    // a last secant step: exact where the coordinate is linear (a segment)
                    roots.push_back(fhi != flo ? std::clamp(lo - flo * (hi - lo) / (fhi - flo), lo, hi) : (lo + hi) / T(2));
                    continue;
                }
                // no pole beyond x on one side: the coordinate only comes back to x, a touch (a
                // double root) at its extremum when it reaches x
                bool below = false, above = false;
                for (const auto &q : pce.poles)
                {
                    below = below || q[0] / q[1] < -eps;
                    above = above || q[0] / q[1] > eps;
                }
                if (!(below && above))
                {
                    const T sg = above ? T(1) : T(-1); // the side of the coordinate: its extremum is a minimum of sg f
                    T lo = pce.a, hi = pce.b;
                    for (int it = 0; it < 200 && hi - lo > t_min; ++it)
                    {
                        const T m1 = lo + (hi - lo) / T(3), m2 = hi - (hi - lo) / T(3);
                        (sg * value(m1) < sg * value(m2) ? hi : lo) = (sg * value(m1) < sg * value(m2) ? m2 : m1);
                    }
                    const T t = (lo + hi) / T(2);
                    if (std::abs(value(t)) <= eps)
                        roots.push_back(t);
                    continue;
                }
                if (pce.b - pce.a <= t_min)
                {
                    roots.push_back((pce.a + pce.b) / T(2));
                    continue;
                }
                auto [l, r] = pieces.split(pce);
                stack.push_back(std::move(l));
                stack.push_back(std::move(r));
            }
        }

        /// Points of a primitive where its parametrization degenerates: apex of a cone, poles of a sphere.
        template <std::floating_point T>
        auto singular_points(const ElementarySurface<T> &s) -> std::vector<point<T, 3>>
        {
            const auto &f = s.frame();
            if (const auto *k = dynamic_cast<const ConicalSurface<T> *>(&s))
                return {f.o - (k->radius() / std::tan(k->semi_angle())) * f.z};
            if (const auto *sp = dynamic_cast<const SphericalSurface<T> *>(&s))
                return {f.o + sp->radius() * f.z, f.o - sp->radius() * f.z};
            return {};
        }

        /// A span of a carrier between its singular points, and its pcurves on both surfaces there (for the crossings).
        template <std::floating_point T>
        struct RestrictionSpan
        {
            T ta{}, tb{};
            std::array<PCurveFit<T>, 2> pcurves; ///< periodic parameters unwrapped
        };

        /// The restriction of a carrier to two boxes: its parts, and the pcurves of its spans.
        template <std::floating_point T>
        struct Restriction
        {
            std::vector<std::array<T, 2>> parts;
            std::vector<RestrictionSpan<T>> spans;
            bool accurate{true}; ///< the pcurves of the spans holding the parts reached their tolerance
        };

        /**
         * Parts of a carrier inside both boxes, cut where a periodic parameter crosses the seam
         * of its box. The carrier is first cut where it goes through an apex or a pole (its
         * pcurves jump there). On each span, the pcurve of the carrier on each surface (periodic
         * parameters unwrapped) is crossed with the lines of the bounds of the box, repeated every
         * turn for a periodic parameter: every root of the coordinate minus the bound, however
         * short the part between them. The pieces between the crossings are kept when their middle
         * is in both boxes; a full turn is joined across the origin of the carrier unless a seam or
         * a singular point is there.
         */
        template <std::floating_point T>
        auto inside_intervals(const Curve<T, 3> &c, T ta, T tb, const ElementarySurface<T> &s1, const std::array<T, 4> &b1,
                              const ElementarySurface<T> &s2, const std::array<T, 4> &b2, bool full_turn, T tol) -> Restriction<T>
        {
            const T two_pi = T(2) * std::numbers::pi_v<T>;
            auto inside = [&](T t) {
                const auto p = c.value(t);
                const auto a = primitive_uv(s1, b1, p), b = primitive_uv(s2, b2, p);
                return in_box(b1, a[0], a[1]) && in_box(b2, b[0], b[1]);
            };
            Restriction<T> res;
            struct Cut
            {
                T t;
                bool seam; ///< on a line of a periodic parameter, or at a singular point: the parts on both sides stay apart
            };
            std::vector<Cut> cuts{{ta, false}, {tb, false}};

            // ---- spans between the singular points the carrier goes through
            std::vector<T> ends{ta, tb};
            for (const auto *s : {&s1, &s2})
                for (const auto &q : singular_points(*s))
                {
                    // a carrier reduced to the singular point (a circle of radius 0 at an apex): a contact point
                    if (norm(c.value(ta) - q) <= tol && norm(c.value(tb) - q) <= tol && norm(c.value((ta + tb) / T(2)) - q) <= tol)
                        return res;
                    const auto [t, d] = project(c, ta, tb, q);
                    if (d <= tol && t > ta + T(1e-9) * (tb - ta) && t < tb - T(1e-9) * (tb - ta))
                    {
                        ends.push_back(t);
                        cuts.push_back({t, true});
                    }
                }
            std::ranges::sort(ends);

            // ---- crossings of the pcurves with the bounds of the boxes, span by span
            for (std::size_t sp = 0; sp + 1 < ends.size(); ++sp)
            {
                RestrictionSpan<T> span{ends[sp], ends[sp + 1]};
                if (!(span.tb > span.ta))
                    continue;
            for (int side = 0; side < 2; ++side)
            {
                const auto &s = side == 0 ? s1 : s2;
                const auto &b = side == 0 ? b1 : b2;
                const auto per = periodic(s);
                auto fit = primitive_pcurve(c, span.ta, span.tb, s, b, tol);
                const auto pc = std::dynamic_pointer_cast<BSCurve<T, 2>>(fit.pcurve);
                span.pcurves[static_cast<std::size_t>(side)] = std::move(fit);
                // the exact extent of the pcurve: the box of its poles
                Box<T, 2> ext;
                for (const auto &q : pc->poles())
                    ext.add(q);
                const T eps = T(1e-10) * (T(1) + ext.diagonal());
                for (int i = 0; i < 2; ++i)
                    for (T bound : {b[2 * i], b[2 * i + 1]})
                    {
                        // every turn of a periodic parameter over the unwrapped range
                        const bool turns = per[static_cast<std::size_t>(i)];
                        const int k0 = turns ? static_cast<int>(std::floor((ext.min[i] - bound) / two_pi)) : 0;
                        const int k1 = turns ? static_cast<int>(std::ceil((ext.max[i] - bound) / two_pi)) : 0;
                        for (int k = k0; k <= k1; ++k)
                        {
                            const T x = bound + T(k) * two_pi;
                            if (x < ext.min[i] - eps || x > ext.max[i] + eps)
                                continue;
                            std::vector<T> roots;
                            coordinate_roots(*pc, i, x, span.ta, span.tb, eps, roots);
                            for (T t : roots)
                                cuts.push_back({t, turns});
                        }
                    }
            }
                res.spans.push_back(std::move(span));
            }
            // sorted, and merged when closer than eps_t to the first of a group (keeping a seam)
            std::ranges::sort(cuts, {}, &Cut::t);
            const T eps_t = T(1e-12) * (T(1) + std::abs(ta) + std::abs(tb));
            std::vector<Cut> merged;
            for (const auto &x : cuts)
                if (!merged.empty() && x.t - merged.back().t <= eps_t)
                    merged.back().seam = merged.back().seam || x.seam;
                else
                    merged.push_back(x);

            // ---- the pieces inside both boxes; two pieces meeting at a cut stay apart: a seam, or a
            //      point where the curve touches a bound (a contact, which gets its vertex)
            auto &r = res.parts;
            for (std::size_t k = 0; k + 1 < merged.size(); ++k)
            {
                const T a = merged[k].t, b = merged[k + 1].t;
                if (!inside((a + b) / T(2)))
                    continue;
                if (!r.empty() && r.back()[1] == a && !merged[k].seam)
                    r.back()[1] = b;
                else
                    r.push_back({a, b});
            }
            // parts reduced to a point (the curve touching a corner of a box)
            std::erase_if(r, [&](const std::array<T, 2> &i) {
                const auto p = c.value(i[0]);
                return norm(c.value(i[1]) - p) <= T(10) * tol && norm(c.value((i[0] + i[1]) / T(2)) - p) <= T(10) * tol;
            });
            // a full turn cut by a box in its middle: the parts at both ends are one curve, unless
            // a seam passes through the origin of the carrier
            if (full_turn && r.size() >= 2 && r.front()[0] == ta && r.back()[1] == tb && !merged.front().seam && !merged.back().seam)
            {
                r.front() = {r.back()[0], r.front()[1] + (tb - ta)};
                r.pop_back();
            }
            // the crossings are as accurate as the pcurves of the spans holding the parts (beyond
            // an apex the carrier may run on the other nappe, outside the box, where no pcurve fits)
            for (const auto &sp : res.spans)
                for (const auto &part : r)
                    if (part[0] < sp.tb && part[1] > sp.ta)
                        for (const auto &f : sp.pcurves)
                            res.accurate = res.accurate && f.deviation <= tol / T(10);
            return res;
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
            const auto rs = inside_intervals(*carrier, range[0], range[1], *e1, box1, *e2, box2, full_turn, o.tol);
            res.complete = res.complete && rs.accurate;
            for (const auto &[a, b] : rs.parts)
            {
                if (norm(carrier->value(b) - carrier->value(a)) <= small && norm(carrier->value((a + b) / T(2)) - carrier->value(a)) <= small)
                    continue;
                IntersectionCurve<T> ic;
                ic.curve = restricted(*carrier, a, b);
                ic.a = a;
                ic.b = b;
                const auto f1 = primitive_pcurve(*ic.curve, a, b, *e1, box1, o.tol);
                const auto f2 = primitive_pcurve(*ic.curve, a, b, *e2, box2, o.tol);
                ic.pc1 = f1.pcurve;
                ic.pc2 = f2.pcurve;
                ic.deviation = std::max(f1.deviation, f2.deviation);
                ic.exact = f1.exact && f2.exact;
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
