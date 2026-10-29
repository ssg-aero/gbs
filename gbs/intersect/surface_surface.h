#pragma once
/**
 * @file surface_surface.h
 * @brief Intersection of two surfaces on bounded parameter boxes: curves with a 3D curve and a
 * pcurve on each surface on the same parameter, contact points, coincidence.
 *
 * Pairs of primitives with a closed form (`analytic_surfaces.h`) give exact lines, circles and
 * ellipses, clipped to the common region of the 3D boxes of the surfaces, then restricted to the
 * parts lying in both boxes and cut where a periodic parameter crosses the seam of a box. A
 * pcurve is exact (a segment) when the parameters of the surface vary linearly along the curve,
 * as for a circle on a coaxial surface or a generatrix, and is interpolated otherwise (samples
 * refined where it misses tol / 10), with its deviation measured. The other pairs go through
 * marching (`marching.h`).
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
#include <ranges>
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

        /// Parameter range of a line inside a box, or nullopt (also for a line or a box not finite).
        template <std::floating_point T>
        auto clip_line(const LineCurve<T> &l, const Box<T, 3> &b) -> std::optional<std::array<T, 2>>
        {
            T lo = std::numeric_limits<T>::lowest(), hi = std::numeric_limits<T>::max();
            const auto &o = l.origin(), &d = l.vector();
            for (int i = 0; i < 3; ++i)
            {
                // (std::max and std::min would keep lo and hi against a NaN)
                if (!std::isfinite(o[i]) || !std::isfinite(d[i]) || !std::isfinite(b.min[i]) || !std::isfinite(b.max[i]))
                    return std::nullopt;
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
            if (!(lo < hi) || !std::isfinite(lo) || !std::isfinite(hi) || lo == std::numeric_limits<T>::lowest())
                return std::nullopt;
            return std::array<T, 2>{lo, hi};
        }

        /**
         * Parameter ranges of the full turn of an ellipse inside a box: each coordinate is
         * c + a cos t + b sin t, which crosses a bound at phi -+ acos((bound - c) / hypot(a, b)). An
         * arc across the origin of the turn is one range, ending beyond the turn.
         */
        template <std::floating_point T>
        auto clip_ellipse(const EllipseCurve<T> &e, const Box<T, 3> &box) -> std::vector<std::array<T, 2>>
        {
            const auto [t0, t1] = e.bounds();
            const auto &f = e.frame();
            std::vector<T> cuts{t0, t1};
            for (std::size_t i = 0; i < 3; ++i)
            {
                const T a = e.radius1() * f.x[i], b = e.radius2() * f.y[i], r = std::hypot(a, b);
                if (!(r > T(0)))
                    continue;
                const T phi = std::atan2(b, a);
                for (T bound : {box.min[i], box.max[i]})
                    if (const T k = (bound - f.o[i]) / r; std::abs(k) < T(1))
                        for (T s : {T(-1), T(1)})
                            cuts.push_back(primitives_detail::wrap(phi + s * std::acos(k), t0));
            }
            std::ranges::sort(cuts);
            std::vector<std::array<T, 2>> r;
            for (std::size_t k = 0; k + 1 < cuts.size(); ++k)
            {
                const T a = cuts[k], b = cuts[k + 1];
                if (!(b > a) || !box.contains(e.value((a + b) / T(2))))
                    continue;
                if (!r.empty() && r.back()[1] == a)
                    r.back()[1] = b;
                else
                    r.push_back({a, b});
            }
            if (r.size() >= 2 && r.front()[0] == t0 && r.back()[1] == t1)
            {
                r.front() = {r.back()[0], r.front()[1] + (t1 - t0)};
                r.pop_back();
            }
            return r;
        }

        /// A pcurve and its fit error: the largest distance between S(pcurve(t)) and the point of the
        /// surface at the parameters of C(t), whatever the distance of the curve to the surface (a
        /// tangency or a near-parallel pair is accepted within tol by the closed forms).
        template <std::floating_point T>
        struct PrimitiveFit
        {
            PCurveFit<T> fit;
            T error{};
        };

        /**
         * Pcurve of c on [a, b] on a primitive: a segment when its parameters vary linearly,
         * interpolated otherwise. The interpolation samples are halved where the fit misses tol / 10
         * in the middle of an interval or at one of 1000 check points: they gather where the pcurve
         * turns fast, as on a circle passing near a pole, where u turns by pi over a length of the
         * order of the distance to the pole.
         */
        template <std::floating_point T>
        auto fit_primitive_pcurve(const Curve<T, 3> &c, T a, T b, const ElementarySurface<T> &s, const std::array<T, 4> &box, T tol)
            -> PrimitiveFit<T>
        {
            const auto per = periodic(s);
            const T two_pi = T(2) * std::numbers::pi_v<T>;
            struct Sample
            {
                T t;
                point<T, 3> p, foot;   ///< the curve, the surface at its parameters
                std::array<T, 2> uv;   ///< in the turn of the box
                bool singular;         ///< at a pole or an apex: no angle
            };
            auto sample = [&](T t, bool interpolated = true) {
                Sample x{.t = t, .p = c.value(t), .foot = {}, .uv = {}, .singular = false};
                x.uv = primitive_uv(s, box, x.p);
                x.foot = s.value(x.uv[0], x.uv[1]);
                x.singular = interpolated && norm(s.value(x.uv[0], x.uv[1], 1, 0)) <= T(1e-8) * norm(s.value(x.uv[0], x.uv[1], 0, 1));
                return x;
            };
            constexpr std::size_t n0 = 16;
            const T t_mid = a + (b - a) * T(n0 / 2) / T(n0); // a sample, the reference of the turns
            // parameters along the curve, unwrapped from the middle outwards (ends may sit on the
            // seam); at a pole or an apex the angle is undefined: it is taken linear in t from the
            // nearest samples (a generatrix keeps its angle up to the apex)
            auto unwrap = [&](const std::vector<Sample> &x) {
                const std::size_t n = x.size();
                std::vector<std::array<T, 2>> uv(n);
                const auto mid = static_cast<std::size_t>(std::ranges::lower_bound(x, t_mid, {}, &Sample::t) - x.begin());
                uv[mid] = x[mid].uv;
                auto follow = [&](std::size_t k, std::size_t ref) {
                    uv[k] = x[k].uv;
                    for (std::size_t i = 0; i < 2; ++i)
                        if (per[i])
                            uv[k][i] = uv[ref][i] + std::remainder(uv[k][i] - uv[ref][i], two_pi);
                };
                for (std::size_t k = mid + 1; k < n; ++k)
                    follow(k, k - 1);
                for (std::size_t k = mid; k-- > 0;)
                    follow(k, k + 1);
                std::vector<std::size_t> regular;
                for (std::size_t k = 0; k < n; ++k)
                    if (!x[k].singular)
                        regular.push_back(k);
                if (regular.size() >= 2)
                    for (std::size_t k = 0; k < n; ++k)
                        if (x[k].singular)
                        {
                            const auto it = std::ranges::upper_bound(regular, k);
                            const std::size_t i = it == regular.begin() ? regular[0] : it == regular.end() ? regular[regular.size() - 2] : *(it - 1);
                            const std::size_t j = it == regular.begin() ? regular[1] : it == regular.end() ? regular.back() : *it;
                            uv[k][0] = uv[i][0] + (uv[j][0] - uv[i][0]) * (x[k].t - x[i].t) / (x[j].t - x[i].t);
                        }
                return uv;
            };
            std::vector<Sample> x;
            for (std::size_t k = 0; k <= n0; ++k)
                x.push_back(sample(k == n0 ? b : a + (b - a) * T(k) / T(n0)));
            {
                const auto uv = unwrap(x);
                T err{}, scale = T(1);
                for (std::size_t k = 0; k <= n0; ++k)
                {
                    const T l = (x[k].t - a) / (b - a);
                    for (std::size_t i = 0; i < 2; ++i)
                    {
                        err = std::max(err, std::abs(uv[k][i] - (uv[0][i] + l * (uv[n0][i] - uv[0][i]))));
                        scale = std::max(scale, std::abs(uv[k][i]));
                    }
                }
                if (err <= T(1e-10) * scale)
                {
                    auto pc = std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{uv[0], uv[n0]}, std::vector<T>{a, a, b, b}, 1);
                    T dev{};
                    constexpr int m = 1000;
                    // (the ends included: a gap varying along the curve, as to a slightly tilted or
                    // offset surface, often peaks there)
                    for (int k = 0; k <= m; ++k)
                    {
                        const T t = a + (b - a) * T(k) / T(m);
                        const auto q = pc->value(t);
                        dev = std::max(dev, norm(s.value(q[0], q[1]) - c.value(t)));
                    }
                    return {{pc, dev, true}, T(0)};
                }
            }
            // the samples at the middle of the intervals, and the check points (once, when needed)
            std::vector<Sample> mids, checks;
            for (std::size_t k = 0; k < n0; ++k)
                mids.push_back(sample((x[k].t + x[k + 1].t) / T(2)));
            const T target = tol / T(10), t_min = T(1e-13) * (T(1) + std::abs(a) + std::abs(b));
            constexpr std::size_t max_samples = 4096, max_rounds = 64;
            std::shared_ptr<BSCurve<T, 2>> pc;
            T dev{}, err{};
            for (std::size_t round = 0;; ++round)
            {
                const auto uv = unwrap(x);
                std::vector<constrType<T, 2, 1>> q(x.size());
                std::vector<T> t(x.size());
                for (std::size_t k = 0; k < x.size(); ++k)
                {
                    q[k] = {uv[k]};
                    t[k] = x[k].t;
                }
                pc = std::make_shared<BSCurve<T, 2>>(interpolate<T, 2>(q, t, std::size_t{3}));
                dev = err = T(0);
                auto missed = [&](const Sample &y) {
                    const auto r = pc->value(y.t);
                    const auto p = s.value(r[0], r[1]);
                    dev = std::max(dev, norm(p - y.p));
                    const T e = norm(p - y.foot);
                    err = std::max(err, e);
                    return e > target;
                };
                std::vector<bool> split(mids.size());
                std::size_t n_split{};
                auto mark = [&](std::size_t k) {
                    if (!split[k] && x[k + 1].t - x[k].t > t_min)
                    {
                        split[k] = true;
                        ++n_split;
                    }
                };
                for (std::size_t k = 0; k < mids.size(); ++k)
                    if (missed(mids[k]))
                        mark(k);
                const bool stop = round + 1 >= max_rounds || x.size() + n_split > max_samples;
                if (n_split == 0 || stop)
                {
                    if (checks.empty())
                    {
                        constexpr int m = 1000;
                        for (int k = 0; k < m; ++k)
                            checks.push_back(sample(a + (b - a) * (T(k) + T(0.5)) / T(m), false));
                    }
                    for (const auto &y : checks)
                        if (missed(y))
                            mark(static_cast<std::size_t>(std::ranges::upper_bound(x, y.t, {}, &Sample::t) - x.begin()) - 1);
                }
                if (n_split == 0 || stop || x.size() + n_split > max_samples)
                    break;
                // the middles of the split intervals become samples, each half gets its middle
                std::vector<Sample> nx, nm;
                for (std::size_t k = 0; k < x.size(); ++k)
                {
                    nx.push_back(x[k]);
                    if (k + 1 == x.size())
                        break;
                    if (!split[k])
                    {
                        nm.push_back(mids[k]);
                        continue;
                    }
                    nm.push_back(sample((x[k].t + mids[k].t) / T(2)));
                    nm.push_back(sample((mids[k].t + x[k + 1].t) / T(2)));
                    nx.push_back(mids[k]);
                }
                x = std::move(nx);
                mids = std::move(nm);
            }
            // (at the samples, which it interpolates, the gap of the curve to the surface: largest
            // at an end for a curve on a slightly offset surface)
            for (const auto &y : x)
            {
                const auto r = pc->value(y.t);
                dev = std::max(dev, norm(s.value(r[0], r[1]) - y.p));
            }
            return {{pc, dev, false}, err};
        }

        /// Pcurve of c on [a, b] on a primitive: a segment when its parameters vary linearly, interpolated otherwise.
        template <std::floating_point T>
        auto primitive_pcurve(const Curve<T, 3> &c, T a, T b, const ElementarySurface<T> &s, const std::array<T, 4> &box, T tol)
            -> PCurveFit<T>
        {
            return fit_primitive_pcurve(c, a, b, s, box, tol).fit;
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
            std::array<T, 2> errors{};           ///< their fit errors (`PrimitiveFit`)
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
            const T eps_t = T(1e-12) * (T(1) + std::abs(ta) + std::abs(tb)); // cuts closer than this are one
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
                auto [fit, error] = fit_primitive_pcurve(c, span.ta, span.tb, s, b, tol);
                const auto pc = std::dynamic_pointer_cast<BSCurve<T, 2>>(fit.pcurve);
                span.pcurves[static_cast<std::size_t>(side)] = std::move(fit);
                span.errors[static_cast<std::size_t>(side)] = error;
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
                            // a root of the pcurve, moved to where the parameter of the curve itself
                            // crosses x when it does so close by: at a small angle, the error of the
                            // pcurve would slide the crossing along the curve (by error / angle)
                            auto g = [&](T t) {
                                const T r = pc->value(t)[static_cast<std::size_t>(i)];
                                const T v = primitive_uv(s, b, c.value(t))[static_cast<std::size_t>(i)];
                                return (turns ? r + std::remainder(v - r, two_pi) : v) - x;
                            };
                            auto polish = [&](T t) {
                                const T g0 = g(t), d = pc->value(t, 1)[static_cast<std::size_t>(i)];
                                if (g0 == T(0) || !(std::abs(d) > T(0)))
                                    return t;
                                const T e = std::clamp(t - T(2) * g0 / d, span.ta, span.tb); // beyond the Newton step
                                if (e == t || !(std::abs(e - t) <= T(1e-3) * (span.tb - span.ta)) || (g(e) < T(0)) == (g0 < T(0)))
                                    return t; // a touch, or no crossing of the curve close by
                                T lo = t, hi = e;
                                for (int it = 0; it < 100 && std::abs(hi - lo) > eps_t; ++it)
                                {
                                    const T m = (lo + hi) / T(2);
                                    ((g(m) < T(0)) == (g0 < T(0)) ? lo : hi) = m;
                                }
                                return (lo + hi) / T(2);
                            };
                            std::vector<T> roots;
                            coordinate_roots(*pc, i, x, span.ta, span.tb, eps, roots);
                            for (T t : roots)
                                cuts.push_back({polish(t), turns});
                        }
                    }
            }
                res.spans.push_back(std::move(span));
            }
            // sorted, and merged when closer than eps_t to the first of a group (keeping a seam)
            std::ranges::sort(cuts, {}, &Cut::t);
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
            // the crossings are as accurate as the pcurves of the spans holding the parts (not a span
            // touching a part at an end, the cuts being merged within eps_t): fitted within tol / 10,
            // on a carrier within tol of the surfaces (the closed forms accept a tangency, an apex or
            // a near-parallel pair within tol)
            for (const auto &sp : res.spans)
                for (const auto &part : r)
                    if (std::min(part[1], sp.tb) - std::max(part[0], sp.ta) > eps_t)
                        for (std::size_t i = 0; i < 2; ++i)
                            res.accurate = res.accurate && sp.errors[i] <= tol / T(10) && sp.pcurves[i].deviation <= tol + tol / T(10);
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

    namespace intersect_detail
    {
        /// The common region of the 3D boxes of s1 on `box1` and s2 on `box2`, grown by tol;
        /// nullopt when the boxes are more than 2 tol apart.
        template <std::floating_point T>
        auto common_region(const Surface<T, 3> &s1, const std::array<T, 4> &box1, const Surface<T, 3> &s2, const std::array<T, 4> &box2, T tol)
            -> std::optional<Box<T, 3>>
        {
            Box<T, 3> region, b1, b2;
            for (const auto &p : SurfacePieces<T>{s1}.pieces(box1))
                b1.add(p.box);
            for (const auto &p : SurfacePieces<T>{s2}.pieces(box2))
                b2.add(p.box);
            for (int i = 0; i < 3; ++i)
            {
                if (std::max(b1.min[i], b2.min[i]) > std::min(b1.max[i], b2.max[i]) + T(2) * tol)
                    return std::nullopt;
                region.min[i] = std::max(b1.min[i], b2.min[i]);
                region.max[i] = std::min(b1.max[i], b2.max[i]);
            }
            region.inflate(tol);
            return region;
        }

        /**
         * Largest distance to the plane `to` of the points of the rectangle of the plane `from` on
         * `box_from` whose foot on `to` lies in its rectangle `box_to` (grown by `margin`): the part
         * of the first rectangle facing the second, clipped by its four sides, a convex polygon at
         * whose corners the distance, affine, is largest. nullopt when no point faces it (or the
         * rectangle is not finite).
         */
        template <std::floating_point T>
        auto facing_gap(const Plane<T> &from, const std::array<T, 4> &box_from, const Plane<T> &to, const std::array<T, 4> &box_to, T margin)
            -> std::optional<T>
        {
            if (!std::ranges::all_of(box_from, [](T x) { return std::isfinite(x); }))
                return std::nullopt;
            const auto &f = to.frame();
            std::vector<point<T, 3>> poly{from.value(box_from[0], box_from[2]), from.value(box_from[1], box_from[2]),
                                          from.value(box_from[1], box_from[3]), from.value(box_from[0], box_from[3])};
            // the sides of the second rectangle on the parameters of the feet, s ((p - o) . a - bound) >= 0
            for (const auto &[a, s, bound] : {std::tuple{f.x, T(1), box_to[0] - margin}, std::tuple{f.x, T(-1), box_to[1] + margin},
                                              std::tuple{f.y, T(1), box_to[2] - margin}, std::tuple{f.y, T(-1), box_to[3] + margin}})
            {
                auto inside = [&](const point<T, 3> &p) { return s * ((p - f.o) * a - bound); };
                std::vector<point<T, 3>> clipped;
                for (std::size_t k = 0; k < poly.size(); ++k)
                {
                    const auto &p = poly[k], &q = poly[(k + 1) % poly.size()];
                    const T gp = inside(p), gq = inside(q);
                    if (gp >= T(0))
                        clipped.push_back(p);
                    if ((gp >= T(0)) != (gq >= T(0)))
                        clipped.push_back(p + (gp / (gp - gq)) * (q - p));
                }
                poly = std::move(clipped);
                if (poly.empty())
                    return std::nullopt;
            }
            return std::ranges::max(poly | std::views::transform([&](const point<T, 3> &p) { return std::abs((p - f.o) * f.z); }));
        }
    }

    /**
     * @brief The primitives s1 on `box1` and s2 on `box2` lie on each other within `tol`: a
     * coincidence decided at another distance than the tolerance of an intersection, such as a
     * confusion radius. Two planes by their largest gap where their rectangles face each other
     * (not across the common region of their 3D boxes: for two flat faces aligned with the axes,
     * that region is a thin tube along their line of intersection, within tol of both however
     * they cross); the other pairs by the closed forms across that region (coaxial within tol, the
     * same meridian within tol), to be confirmed where the region is thin.
     */
    template <std::floating_point T>
    [[nodiscard]] bool coincident(const ElementarySurface<T> &s1, const std::array<T, 4> &box1, const ElementarySurface<T> &s2,
                                  const std::array<T, 4> &box2, T tol)
    {
        const auto *p1 = dynamic_cast<const Plane<T> *>(&s1), *p2 = dynamic_cast<const Plane<T> *>(&s2);
        if (p1 && p2)
        {
            const auto g12 = intersect_detail::facing_gap(*p1, box1, *p2, box2, tol), g21 = intersect_detail::facing_gap(*p2, box2, *p1, box1, tol);
            return (g12 || g21) && g12.value_or(T(0)) <= tol && g21.value_or(T(0)) <= tol;
        }
        const auto region = intersect_detail::common_region(s1, box1, s2, box2, tol);
        if (!region)
            return false;
        const auto ai = intersect_detail::analytic_intersection(s1, s2, tol, *region);
        return ai && ai->coincident;
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
        if (!e1 || !e2 || !o.closed_forms)
            return marching(s1, box1, s2, box2, o);
        // the region where both surfaces are: the carriers are clipped to it, and a near-parallel
        // pair is decided across it
        const auto common = common_region(s1, box1, s2, box2, o.tol);
        if (!common)
            return res; // the boxes are apart: nothing in common
        const auto &region = *common;
        const auto ai = analytic_intersection(*e1, *e2, o.tol, region);
        if (!ai)
            return marching(s1, box1, s2, box2, o);
        if (ai->coincident)
        {
            res.coincident = true;
            return res;
        }
        // (an ellipse clipped with a margin: its cuts are found by acos, ill-conditioned at a bound
        // it barely reaches)
        auto clip = region;
        clip.inflate(o.tol);
        const T small = T(10) * o.tol;
        for (const auto &[carrier, tangent] : ai->carriers)
        {
            const auto [ta, tb] = carrier->bounds();
            const auto *l = dynamic_cast<const LineCurve<T> *>(carrier.get());
            std::vector<std::array<T, 2>> ranges;
            if (l)
            {
                if (const auto r = clip_line(*l, region))
                    ranges.push_back(*r);
            }
            else
                ranges = clip_ellipse(dynamic_cast<const EllipseCurve<T> &>(*carrier), clip);
            for (const auto &range : ranges)
            {
                const bool full_turn = !l && range[0] == ta && range[1] == tb;
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
