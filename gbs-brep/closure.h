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
#include <span>
#include <vector>

#include <gbs/curvebreaks.h>
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
        /// Gauss-Legendre rule of 8 points on [-1, 1]: the positive abscissas and their weights (exact up to degree 15).
        template <std::floating_point T>
        inline constexpr std::array<std::array<T, 2>, 4> gauss_legendre_8{{
            {T(0.18343464249564980493947614236018398), T(0.36268378337836198296515044927719561)},
            {T(0.52553240991632898581773904918924635), T(0.31370664587788728733796220198660131)},
            {T(0.79666647741362673959155393647583044), T(0.22238103445337447054435599442624088)},
            {T(0.96028985649753623168356086856947299), T(0.10122853629037625915253135430996219)},
        }};

        /// The integral of f over [a, b] by the 8-point Gauss-Legendre rule.
        template <std::floating_point T, typename F>
        [[nodiscard]] auto gauss_8(F &&f, T a, T b) -> T
        {
            const T mid = (a + b) / 2, half = (b - a) / 2;
            T s{};
            for (const auto &[x, w] : gauss_legendre_8<T>)
                s += w * (f(mid - half * x) + f(mid + half * x));
            return s * half;
        }

        /// The knots of a B-spline surface (rational or not) in u, empty for any other surface.
        template <std::floating_point T>
        [[nodiscard]] auto surface_knots_u(const Surface<T, 3> &srf) -> std::vector<T>
        {
            if (const auto *b = dynamic_cast<const BSSurfaceGeneral<T, 3, false> *>(&srf))
                return b->knotsU();
            if (const auto *b = dynamic_cast<const BSSurfaceGeneral<T, 3, true> *>(&srf))
                return b->knotsU();
            return {};
        }
    } // namespace detail

    /**
     * @brief Volume integral of a face, (1/3) ∬ S · (Su × Sv) du dv over the
     * region of (u,v) bounded by its wires (outer minus holes), turned by
     * Green's theorem into the integral over its wires (holes included, by
     * their own clockwise sense) of F(u, v) dv, where
     * F(u, v) = (1/3) ∫ from u0 to u of S · (Su × Sv)(s, v) ds and u0 is the
     * smallest u of the outer boundary. Both integrals are 8-point
     * Gauss-Legendre rules: along each co-edge on the smooth spans of its
     * pcurve, each cut into up to max(2, n / 16) pieces by its (u,v) length;
     * along u on panels of 1 / max(2, n / 16) of the face's width in u, cut
     * at the knots of a B-spline surface. Exact, up to rounding, on a planar
     * face with straight pcurves, whatever its holes; on a curved one, the
     * error falls as a high power of 1 / n. Every co-edge of the face must
     * carry a pcurve. Contribution of the face to the volume of a closed
     * shell when its normal Su × Sv points outwards.
     */
    template <std::floating_point T>
    [[nodiscard]] auto face_volume_integral(const Model<T> &m, FaceId fid, std::size_t n = 64) -> T
    {
        const auto &f = m.face(fid);
        if (!f.surface || f.wires.empty())
            return T(0);
        const auto &srf = *f.surface;
        const std::size_t k = std::max<std::size_t>(n / 16, 2);
        // the box of the outer boundary in (u,v)
        const auto outer = uv_polygon(m, std::span<const CoEdge<T>>{m.wire(f.wires.front()).coedges}, 16);
        point<T, 2> lo = outer[0], hi = outer[0];
        for (const auto &p : outer)
            for (std::size_t c{}; c < 2; ++c)
                lo[c] = std::min(lo[c], p[c]), hi[c] = std::max(hi[c], p[c]);
        const T width = hi[0] - lo[0], size = std::max(width, hi[1] - lo[1]);
        if (!(size > T(0)))
            return T(0);
        // the panels along u, from u0 = lo[0]: a regular grid and the knots of the surface
        std::vector<T> breaks;
        for (std::size_t i{}; i <= k; ++i)
            breaks.push_back(lo[0] + width * T(i) / T(k));
        for (T u : detail::surface_knots_u(srf))
            if (lo[0] < u && u < hi[0])
                breaks.push_back(u);
        std::ranges::sort(breaks);
        const T u0 = breaks.front();
        auto integrand = [&](T u, T v) {
            return srf.value(u, v) * (srf.value(u, v, 1, 0) ^ srf.value(u, v, 0, 1));
        };
        // F(u, v) but for the factor 1/3: the panels between u0 and u, then the part of the last one
        auto F = [&](T u, T v) {
            auto line = [&](T s) { return integrand(s, v); };
            T s{}, a = u0;
            for (std::size_t i = 1; i < breaks.size() && breaks[i] < u; ++i)
            {
                s += detail::gauss_8(line, a, breaks[i]);
                a = breaks[i];
            }
            return s + (a == u ? T(0) : detail::gauss_8(line, a, u)); // u < u0: integrated backwards
        };
        T sum{};
        for (auto wid : f.wires)
            for (const auto &ce : m.wire(wid).coedges)
            {
                if (!ce.pcurve)
                    throw BuildFailure({BuildErrc::MissingPCurve, "co-edge without pcurve", {ce.edge}});
                const auto &e = m.edge(ce.edge);
                const auto &pc = *ce.pcurve;
                // along the pcurve, from e.u1 to e.u2; a Reversed co-edge runs backwards
                auto term = [&](T t) {
                    const auto dv = pc.value(t, 1)[1];
                    if (dv == T(0))
                        return T(0); // an iso v: no contribution (and a degenerate edge at a pole)
                    const auto uv = pc.value(t);
                    return F(uv[0], uv[1]) * dv;
                };
                const auto spans = smooth_spans(pc, e.u1, e.u2);
                T s{};
                for (std::size_t i = 1; i < spans.size(); ++i)
                {
                    const T a = spans[i - 1], b = spans[i];
                    // pieces in proportion to the (u,v) length of the span, by a chord polygon of 8 sides
                    T len{};
                    for (std::size_t j{}; j < 8; ++j)
                        len += distance(pc.value(a + (b - a) * T(j) / T(8)), pc.value(a + (b - a) * T(j + 1) / T(8)));
                    const auto pieces = std::clamp<std::size_t>(static_cast<std::size_t>(std::ceil(T(k) * len / size)), 2, k);
                    for (std::size_t j{}; j < pieces; ++j)
                        s += detail::gauss_8(term, a + (b - a) * T(j) / T(pieces), a + (b - a) * T(j + 1) / T(pieces));
                }
                sum += ce.orient == Orientation::Forward ? s : -s;
            }
        return sum / T(3);
    }

    /**
     * @brief Signed volume enclosed by a closed shell (divergence theorem):
     * positive when the effective normals (face-use sense applied) point
     * outwards. The sum of `face_volume_integral` over its faces: exact to
     * the rounding on planar faces with straight pcurves, holes included;
     * on curved and trimmed faces the error falls as a high power of 1 / n.
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
