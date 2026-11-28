#pragma once

/**
 * @file closure.h
 * @brief Native BREP core, stage 1: geometric queries used by the face
 * builders: closure / degeneracy of a surface's boundary, signed area of a
 * wire in a face's parametric space.
 *
 * Design: docs/sources/design/brep_core.md, sections 3.3 and 6.1;
 * architecture note docs/sources/design/brep_pr04_face_natural.md.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <ranges>
#include <span>
#include <vector>

#include <boost/math/quadrature/gauss.hpp>

#include <gbs/curvebreaks.h>
#include <gbs/offsets.h>
#include <gbs-brep/model.h>

namespace gbs::brep
{
    /**
     * @brief How the parametric rectangle of a surface closes on itself.
     *
     * `closed_u`: S(u1, v) == S(u2, v) for every v (seam along u = u1).
     * `degenerate_u1`: the iso u = u1 is a single point (pole, apex); same for
     * the three other sides.
     */
    struct SurfaceClosure
    {
        bool closed_u{false}, closed_v{false};
        bool degenerate_u1{false}, degenerate_u2{false};
        bool degenerate_v1{false}, degenerate_v2{false};

        [[nodiscard]] constexpr bool degenerate() const noexcept
        {
            return degenerate_u1 || degenerate_u2 || degenerate_v1 || degenerate_v2;
        }
    };

    /**
     * @brief Detects closure and degenerate sides numerically, by sampling
     * the four boundary isos at `n_samples` points and comparing within `tol`.
     * Works for any `Surface<T,3>` (NURBS, revolution, offset…) without
     * relying on the surface type; a side detected as degenerate is never
     * reported as a seam.
     */
    template <std::floating_point T>
    [[nodiscard]] auto surface_closure(const Surface<T, 3> &srf, T tol, std::size_t n_samples = 9) -> SurfaceClosure
    {
        n_samples = std::max<std::size_t>(n_samples, 3);
        const auto [u1, u2, v1, v2] = srf.bounds();
        auto param = [n_samples](T a, T b, std::size_t i) { return a + (b - a) * T(i) / T(n_samples - 1); };

        T d_closed_u{}, d_closed_v{}, d_u1{}, d_u2{}, d_v1{}, d_v2{};
        const auto p_u1v1 = srf.value(u1, v1), p_u2v1 = srf.value(u2, v1), p_u1v2 = srf.value(u1, v2);
        for (std::size_t i{}; i < n_samples; ++i)
        {
            const T v = param(v1, v2, i);
            const T u = param(u1, u2, i);
            const auto a = srf.value(u1, v), b = srf.value(u2, v);
            const auto c = srf.value(u, v1), d = srf.value(u, v2);
            d_closed_u = std::max(d_closed_u, distance(a, b));
            d_closed_v = std::max(d_closed_v, distance(c, d));
            d_u1 = std::max(d_u1, distance(a, p_u1v1));
            d_u2 = std::max(d_u2, distance(b, p_u2v1));
            d_v1 = std::max(d_v1, distance(c, p_u1v1));
            d_v2 = std::max(d_v2, distance(d, p_u1v2));
        }
        SurfaceClosure c;
        c.degenerate_u1 = d_u1 <= tol;
        c.degenerate_u2 = d_u2 <= tol;
        c.degenerate_v1 = d_v1 <= tol;
        c.degenerate_v2 = d_v2 <= tol;
        c.closed_u = d_closed_u <= tol && !c.degenerate_u1 && !c.degenerate_u2;
        c.closed_v = d_closed_v <= tol && !c.degenerate_v1 && !c.degenerate_v2;
        return c;
    }

    /**
     * @brief Points of a chain of co-edges in the (u,v) space of their face,
     * `n_per_coedge - 1` per co-edge (the last point of a co-edge is the first
     * of the next), co-edge senses applied. Every co-edge must carry a pcurve.
     */
    template <std::floating_point T>
    [[nodiscard]] auto uv_polygon(const Model<T> &m, std::span<const CoEdge<T>> coedges, std::size_t n_per_coedge = 16)
        -> std::vector<point<T, 2>>
    {
        n_per_coedge = std::max<std::size_t>(n_per_coedge, 2);
        std::vector<point<T, 2>> poly;
        for (const auto &ce : coedges)
        {
            const auto &e = m.edge(ce.edge);
            for (std::size_t i{}; i + 1 < n_per_coedge; ++i)
                poly.push_back(coedge_uv(m, ce, e.u1 + (e.u2 - e.u1) * T(i) / T(n_per_coedge - 1)));
        }
        return poly;
    }

    /// Polygon of a chain of co-edges in (u,v) with the largest chord sag left (`uv_polygon_adaptive`).
    template <std::floating_point T>
    struct UVPolygon
    {
        std::vector<point<T, 2>> points;
        T sag{}; ///< largest distance from the pcurve point at the middle parameter of a segment to its middle
    };

    /**
     * @brief Points of a chain of co-edges in the (u,v) space of their face,
     * as `uv_polygon`, refined where the pcurves bend: each co-edge is cut
     * into `n_segments` segments, each halved, at most `max_depth` times,
     * while its chord sag is larger than `sag`. A polygon inscribed in a
     * curved boundary cuts its convex arcs and overhangs its concave ones by
     * about the sag left.
     */
    template <std::floating_point T>
    [[nodiscard]] auto uv_polygon_adaptive(const Model<T> &m, std::span<const CoEdge<T>> coedges, T sag, std::size_t n_segments = 15,
                                           std::size_t max_depth = 10) -> UVPolygon<T>
    {
        n_segments = std::max<std::size_t>(n_segments, 1);
        struct Piece
        {
            T t0, t1;
            point<T, 2> p0, p1;
            std::size_t depth;
        };
        UVPolygon<T> poly;
        std::vector<Piece> stack;
        for (const auto &ce : coedges)
        {
            const auto &e = m.edge(ce.edge);
            auto at = [&](T t) { return coedge_uv(m, ce, t); };
            for (std::size_t i{}; i < n_segments; ++i)
            {
                const T t0 = e.u1 + (e.u2 - e.u1) * T(i) / T(n_segments), t1 = e.u1 + (e.u2 - e.u1) * T(i + 1) / T(n_segments);
                stack.push_back({t0, t1, at(t0), at(t1), 0});
                while (!stack.empty())
                {
                    const auto s = stack.back();
                    stack.pop_back();
                    const T tm = (s.t0 + s.t1) / 2;
                    const auto pm = at(tm);
                    const T d = distance(pm, point<T, 2>{(s.p0[0] + s.p1[0]) / 2, (s.p0[1] + s.p1[1]) / 2});
                    if (s.depth < max_depth && d > sag)
                    {
                        stack.push_back({tm, s.t1, pm, s.p1, s.depth + 1}); // the first half is taken first
                        stack.push_back({s.t0, tm, s.p0, pm, s.depth + 1});
                        continue;
                    }
                    poly.points.push_back(s.p0);
                    poly.sag = std::max(poly.sag, d);
                }
            }
        }
        return poly;
    }

    /// Distance from `p` to the boundary of the closed polygon `poly`.
    template <std::floating_point T>
    [[nodiscard]] auto distance_to_polygon(const std::vector<point<T, 2>> &poly, const point<T, 2> &p) -> T
    {
        T d = std::numeric_limits<T>::max();
        for (std::size_t i{}, j = poly.size() - 1; i < poly.size(); j = i++)
        {
            const auto &a = poly[j], &b = poly[i];
            const T ex = b[0] - a[0], ey = b[1] - a[1], l2 = ex * ex + ey * ey;
            const T s = l2 > T(0) ? std::clamp(((p[0] - a[0]) * ex + (p[1] - a[1]) * ey) / l2, T(0), T(1)) : T(0);
            d = std::min(d, std::hypot(p[0] - a[0] - s * ex, p[1] - a[1] - s * ey));
        }
        return d;
    }

    /// Signed area of a closed polygon (shoelace formula), positive counter-clockwise.
    template <std::floating_point T>
    [[nodiscard]] auto signed_area(const std::vector<point<T, 2>> &poly) -> T
    {
        T a{};
        for (std::size_t i{}; i < poly.size(); ++i)
        {
            const auto &p = poly[i];
            const auto &q = poly[(i + 1) % poly.size()];
            a += p[0] * q[1] - q[0] * p[1];
        }
        return a / T(2);
    }

    /// Even-odd test: is `p` inside the closed polygon?
    template <std::floating_point T>
    [[nodiscard]] bool inside(const std::vector<point<T, 2>> &poly, const point<T, 2> &p)
    {
        bool in = false;
        for (std::size_t i{}, j = poly.size() - 1; i < poly.size(); j = i++)
        {
            const auto &a = poly[i];
            const auto &b = poly[j];
            if ((a[1] > p[1]) != (b[1] > p[1]) && p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0])
                in = !in;
        }
        return in;
    }

    /**
     * @brief Region of (u,v) bounded by a closed chain of co-edges, for testing many points: its
     * polygon of 15 segments per co-edge, and where a point comes within twice the sag of it (the
     * inscribed polygon cuts the convex arcs and overhangs the concave ones), a polygon refined to a
     * sag of 1e-7 of its box, made once when first needed. The co-edges must outlive it.
     */
    template <std::floating_point T>
    class UVRegion
    {
        const Model<T> *m_;
        std::span<const CoEdge<T>> coedges_;
        UVPolygon<T> coarse_;
        mutable std::vector<point<T, 2>> fine_;

    public:
        UVRegion(const Model<T> &m, std::span<const CoEdge<T>> coedges)
            : m_{&m}, coedges_{coedges}, coarse_{uv_polygon_adaptive(m, coedges, T(0), 15, 0)}
        {
        }

        /// Even-odd test of p, on the refined polygon near the boundary.
        [[nodiscard]] bool contains(const point<T, 2> &p) const
        {
            if (distance_to_polygon(coarse_.points, p) > 2 * coarse_.sag)
                return inside(coarse_.points, p);
            if (fine_.empty())
            {
                point<T, 2> a = coarse_.points.front(), b = a;
                for (const auto &q : coarse_.points)
                    for (std::size_t k{}; k < 2; ++k)
                        a[k] = std::min(a[k], q[k]), b[k] = std::max(b[k], q[k]);
                fine_ = uv_polygon_adaptive(*m_, coedges_, T(1e-7) * distance(a, b)).points;
            }
            return inside(fine_, p);
        }
    };

    /**
     * @brief Signed area enclosed by a chain of co-edges in the (u,v) space
     * of their face, by the shoelace formula on `n_per_coedge` points per
     * co-edge. Positive for a counter-clockwise (outer) boundary, negative for
     * a clockwise one (hole).
     */
    template <std::floating_point T>
    [[nodiscard]] auto uv_signed_area(const Model<T> &m, std::span<const CoEdge<T>> coedges, std::size_t n_per_coedge = 16) -> T
    {
        return signed_area(uv_polygon(m, coedges, n_per_coedge));
    }

    /// Signed area of a wire whose co-edges all carry a pcurve.
    template <std::floating_point T>
    [[nodiscard]] auto uv_signed_area(const Model<T> &m, WireId wid, std::size_t n_per_coedge = 16) -> T
    {
        return uv_signed_area(m, std::span<const CoEdge<T>>{m.wire(wid).coedges}, n_per_coedge);
    }


    namespace detail
    {
        /// The 10-point Gauss-Legendre rule (Boost's, tabulated: exact up to degree 19).
        template <std::floating_point T>
        using gauss_10 = boost::math::quadrature::gauss<T, 10>;

        /// Calls g(t, w) for each node t of the 10-point Gauss-Legendre rule on [a, b], w its weight (negative if b < a).
        template <std::floating_point T, typename G>
        void gauss_10_nodes(T a, T b, G &&g)
        {
            const T mid = (a + b) / 2, half = (b - a) / 2;
            const auto &x = gauss_10<T>::abscissa();
            const auto &w = gauss_10<T>::weights();
            for (std::size_t i{}; i < x.size(); ++i)
            {
                g(mid - half * x[i], half * w[i]);
                if (x[i] != T(0))
                    g(mid + half * x[i], half * w[i]);
            }
        }

        /// Number of Chebyshev points of the tables of `face_volume_integral`: interpolation of degree 15.
        inline constexpr std::size_t chebyshev_count = 16;

        /// The Chebyshev points of the second kind on [-1, 1], cos(j pi / 15), from 1 down to -1.
        template <std::floating_point T>
        [[nodiscard]] auto chebyshev_points() -> std::array<T, chebyshev_count>
        {
            std::array<T, chebyshev_count> x{};
            for (std::size_t j{}; j < chebyshev_count; ++j)
                x[j] = std::cos(std::numbers::pi_v<T> * T(j) / T(chebyshev_count - 1));
            return x;
        }

        /**
         * The polynomial of degree 15 through the values y[j] at the Chebyshev points x[j], at t (barycentric
         * formula: weights of alternate signs, halved at the ends).
         */
        template <std::floating_point T>
        [[nodiscard]] auto chebyshev_interpolate(const std::array<T, chebyshev_count> &x, const T *y, T t) -> T
        {
            T num{}, den{};
            for (std::size_t j{}; j < chebyshev_count; ++j)
            {
                const T d = t - x[j];
                if (d == T(0))
                    return y[j];
                const T w = (j % 2 == 0 ? T(1) : T(-1)) * (j == 0 || j + 1 == chebyshev_count ? T(0.5) : T(1)) / d;
                num += w * y[j];
                den += w;
            }
            return num / den;
        }

        /// The distinct values of the flat knots k, in increasing order (none without knots).
        template <std::floating_point T>
        [[nodiscard]] auto distinct_knots(const std::optional<std::vector<T>> &k) -> std::vector<T>
        {
            std::vector<T> r;
            if (k)
                for (T x : *k)
                    if (r.empty() || x != r.back())
                        r.push_back(x);
            return r;
        }

        /// x brought onto [a, b] when it lies beyond them by at most 1e-5 of b - a; as it is otherwise.
        template <std::floating_point T>
        [[nodiscard]] auto onto_bounds(T x, T a, T b) -> T
        {
            const T band = std::isfinite(b - a) ? T(1e-5) * (b - a) : T(0);
            if (x < a && !(x < a - band))
                return a;
            if (x > b && !(x > b + band))
                return b;
            return x;
        }

        /**
         * The point of the surface at (u, v), a parameter beyond its bounds by at most 1e-5 of their range
         * taken on them (a pcurve within its tolerance off a boundary, where a B-spline refuses it).
         */
        template <std::floating_point T>
        [[nodiscard]] auto value_onto_bounds(const Surface<T, 3> &srf, const point<T, 2> &uv) -> point<T, 3>
        {
            const auto [u1, u2, v1, v2] = srf.bounds();
            return srf.value(onto_bounds(uv[0], u1, u2), onto_bounds(uv[1], v1, v2));
        }

        /**
         * The parameter in [t0, t1] where the coordinate d of the curve c equals x, by bisection; x lies strictly
         * between its values at t0 and t1, increasing from t0 to t1 if `up`.
         */
        template <std::floating_point T>
        [[nodiscard]] auto crossing(const Curve<T, 2> &c, std::size_t d, T t0, T t1, bool up, T x) -> T
        {
            for (int i{}; i < 64 && t1 - t0 > std::numeric_limits<T>::epsilon() * (std::abs(t0) + std::abs(t1)); ++i)
            {
                const T t = (t0 + t1) / 2;
                ((c.value(t)[d] < x) == up ? t0 : t1) = t;
            }
            return (t0 + t1) / 2;
        }
    } // namespace detail

    /**
     * @brief Volume integral of a face, (1/3) ∬ S · (Su × Sv) du dv over the
     * region of (u,v) bounded by its wires (outer minus holes), turned by
     * Green's theorem into the integral over its wires of F(u, v) dv, where
     * F(u, v) = (1/3) ∫ from u0 to u of S · (Su × Sv)(s, v) ds. Each wire
     * counts with the sense expected of it, the outer one counter-clockwise in
     * (u,v) and the holes clockwise, whatever its own (the region is the same).
     *
     * Quadratures, by 10-point Gauss-Legendre rules (exact up to degree 19),
     * k = max(2, n / 16):
     * - along each co-edge, on the smooth spans of its pcurve, cut where it
     *   crosses a knot of the surface (of a B-spline, or of the basis of an
     *   offset: `surface_flat_knots`), each piece into up to k by its (u,v)
     *   length;
     * - along u, from u0 the smallest u of the nodes of these rules, on k
     *   panels cut at the knots in u. The integrals from u0 to each panel are
     *   tabulated at 16 Chebyshev points of each of k intervals of v (cut at
     *   the knots in v) where a node falls, and interpolated (degree 15): a
     *   node costs one rule along u, not one per panel.
     *
     * Exact, up to rounding, on a planar face with straight pcurves and on a
     * polynomial B-spline face of degree up to 3 with straight pcurves,
     * whatever their holes; on a curved one, the error falls as a high power
     * of 1 / n. A parameter beyond the bounds of the surface by at most 1e-5
     * of their range (a pcurve within its tolerance off a boundary) is taken
     * on them, where a B-spline refuses it. Every co-edge of the face must
     * carry a pcurve. Contribution of the face to the volume of a closed shell
     * when its normal Su × Sv points outwards.
     */
    template <std::floating_point T>
    [[nodiscard]] auto face_volume_integral(const Model<T> &m, FaceId fid, std::size_t n = 64) -> T
    {
        const auto &f = m.face(fid);
        if (!f.surface || f.wires.empty())
            return T(0);
        const auto &srf = *f.surface;
        const std::size_t k = std::max<std::size_t>(n / 16, 2);
        // the polygon of each wire (16 points per co-edge): its sense, and the size of the outer one in (u,v)
        std::vector<T> areas;
        point<T, 2> lo{}, hi{};
        for (std::size_t w{}; w < f.wires.size(); ++w)
        {
            const auto poly = uv_polygon(m, std::span<const CoEdge<T>>{m.wire(f.wires[w]).coedges}, 16);
            areas.push_back(signed_area(poly));
            if (w == 0)
            {
                lo = hi = poly.front();
                for (const auto &p : poly)
                    for (std::size_t c{}; c < 2; ++c)
                        lo[c] = std::min(lo[c], p[c]), hi[c] = std::max(hi[c], p[c]);
            }
        }
        const T size = std::max(hi[0] - lo[0], hi[1] - lo[1]);
        if (!(size > T(0)))
            return T(0);
        const auto flat = surface_flat_knots(srf);
        const std::array<std::vector<T>, 2> knots{detail::distinct_knots(flat[0]), detail::distinct_knots(flat[1])};

        // the nodes along the wires: (u, v) and the weight of F there (that of the rule times dv/dt, signed by
        // the senses of the co-edge and of its wire); an iso v (dv = 0, a degenerate edge at a pole) gives none
        struct Node
        {
            T u, v, w;
        };
        std::vector<Node> nodes;
        for (std::size_t w{}; w < f.wires.size(); ++w)
        {
            const T wire_sign = (w == 0 ? areas[w] < T(0) : areas[w] > T(0)) ? T(-1) : T(1);
            for (const auto &ce : m.wire(f.wires[w]).coedges)
            {
                const auto &e = m.edge(ce.edge);
                const auto &pc = *ce.pcurve; // checked by uv_polygon
                const T sign = ce.orient == Orientation::Forward ? wire_sign : -wire_sign; // a Reversed co-edge runs backwards
                const auto spans = smooth_spans(pc, e.u1, e.u2);
                for (const auto [a, b] : std::views::zip(spans, spans | std::views::drop(1)))
                {
                    // a chord polygon of 8 sides: the (u,v) length of the span, and where it crosses the knots
                    std::array<T, 9> t{};
                    std::array<point<T, 2>, 9> p{};
                    T len{};
                    for (std::size_t j{}; j < 9; ++j)
                    {
                        t[j] = a + (b - a) * T(j) / T(8);
                        p[j] = pc.value(t[j]);
                        if (j > 0)
                            len += distance(p[j - 1], p[j]);
                    }
                    std::vector<T> cuts{a, b};
                    for (std::size_t d{}; d < 2; ++d)
                        for (std::size_t j{}; j < 8; ++j)
                        {
                            const T x0 = p[j][d], x1 = p[j + 1][d];
                            const auto first = std::ranges::upper_bound(knots[d], std::min(x0, x1));
                            const auto last = std::ranges::lower_bound(knots[d], std::max(x0, x1));
                            for (auto it = first; it < last; ++it)
                                cuts.push_back(detail::crossing(pc, d, t[j], t[j + 1], x0 < x1, *it));
                        }
                    std::ranges::sort(cuts);
                    for (const auto [c0, c1] : std::views::zip(cuts, cuts | std::views::drop(1)))
                    {
                        if (!(c1 > c0))
                            continue;
                        // pieces in proportion to the (u,v) length, 1 to k
                        const T q = T(k) * len * (c1 - c0) / ((b - a) * size);
                        const std::size_t pieces = q >= T(k) ? k : q > T(1) ? static_cast<std::size_t>(std::ceil(q)) : 1;
                        for (std::size_t j{}; j < pieces; ++j)
                            detail::gauss_10_nodes(c0 + (c1 - c0) * T(j) / T(pieces), c0 + (c1 - c0) * T(j + 1) / T(pieces), [&](T s, T ws) {
                                const T dv = pc.value(s, 1)[1];
                                if (dv == T(0))
                                    return;
                                const auto uv = pc.value(s);
                                nodes.push_back({uv[0], uv[1], sign * ws * dv});
                            });
                    }
                }
            }
        }
        if (nodes.empty())
            return T(0);

        // the box of the nodes: the panels along u from u0, its smallest u, and the intervals of v of the tables,
        // each a regular grid of k and the knots
        T ulo = nodes.front().u, uhi = ulo, vlo = nodes.front().v, vhi = vlo;
        for (const auto &q : nodes)
            ulo = std::min(ulo, q.u), uhi = std::max(uhi, q.u), vlo = std::min(vlo, q.v), vhi = std::max(vhi, q.v);
        auto grid = [k](T a, T b, const std::vector<T> &kn) {
            std::vector<T> r;
            for (std::size_t i{}; i <= k; ++i)
                r.push_back(a + (b - a) * T(i) / T(k));
            for (T x : kn)
                if (a < x && x < b)
                    r.push_back(x);
            std::ranges::sort(r);
            const auto [first, last] = std::ranges::unique(r);
            r.erase(first, last);
            if (r.size() == 1)
                r.push_back(r.front()); // a single value: one interval of no width
            return r;
        };
        const auto ub = grid(ulo, uhi, knots[0]), vb = grid(vlo, vhi, knots[1]);
        const auto [bu1, bu2, bv1, bv2] = srf.bounds();
        auto integrand = [&](T u, T v) {
            u = detail::onto_bounds(u, bu1, bu2), v = detail::onto_bounds(v, bv1, bv2);
            return srf.value(u, v) * (srf.value(u, v, 1, 0) ^ srf.value(u, v, 0, 1));
        };
        auto along_u = [&](T v, T a, T b) {
            T s{};
            detail::gauss_10_nodes(a, b, [&](T x, T w) { s += w * integrand(x, v); });
            return s;
        };
        // in the interval iv of v, at its Chebyshev points, the integral along u from u0 to each break (made
        // when a node first needs it)
        const auto xc = detail::chebyshev_points<T>();
        constexpr std::size_t nc = detail::chebyshev_count;
        std::vector<std::vector<T>> tables(vb.size() - 1);
        auto table = [&](std::size_t iv) -> const std::vector<T> & {
            auto &tab = tables[iv];
            if (tab.empty())
            {
                tab.assign(ub.size() * nc, T(0));
                for (std::size_t c{}; c < nc; ++c)
                {
                    const T v = (vb[iv] + vb[iv + 1]) / 2 + (vb[iv + 1] - vb[iv]) / 2 * xc[c];
                    T s{};
                    for (std::size_t i = 1; i < ub.size(); ++i)
                        tab[i * nc + c] = s += along_u(v, ub[i - 1], ub[i]);
                }
            }
            return tab;
        };
        // F(u, v) but for the factor 1/3: the table up to the last break at or before u, then the rest of the panel
        auto F = [&](T u, T v) {
            const auto j = static_cast<std::size_t>(std::max<std::ptrdiff_t>(std::ranges::upper_bound(ub, u) - ub.begin() - 1, 0));
            const auto iv = static_cast<std::size_t>(std::clamp<std::ptrdiff_t>(std::ranges::upper_bound(vb, v) - vb.begin() - 1, 0, std::ssize(vb) - 2));
            const T va = vb[iv], vz = vb[iv + 1];
            const T full = j == 0 ? T(0) : detail::chebyshev_interpolate(xc, table(iv).data() + j * nc, vz > va ? (2 * v - va - vz) / (vz - va) : T(0));
            return full + (ub[j] == u ? T(0) : along_u(v, ub[j], u));
        };
        T sum{};
        for (const auto &q : nodes)
            sum += q.w * F(q.u, q.v);
        return sum / T(3);
    }

    /**
     * @brief Signed volume enclosed by a closed shell (divergence theorem):
     * positive when the effective normals (face-use sense applied) point
     * outwards. The sum of `face_volume_integral` over its faces: exact to
     * the rounding on planar faces and polynomial B-spline faces of degree up
     * to 3 with straight pcurves, holes included; on curved and trimmed faces
     * the error falls as a high power of 1 / n.
     */
    template <std::floating_point T>
    [[nodiscard]] auto signed_volume(const Model<T> &m, ShellId sid, std::size_t n = 64) -> T
    {
        T v{};
        for (const auto &fu : m.shell(sid).faces)
            v += (fu.orient == Orientation::Forward ? T(1) : T(-1)) * face_volume_integral(m, fu.face, n);
        return v;
    }

} // namespace gbs::brep
