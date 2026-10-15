#pragma once
/**
 * @file classify.h
 * @brief Classification of a point against a face (in its parameter space) and against a solid,
 * and an interior point of a face.
 *
 * - Point in a face: a ray from (u, v) in the parameter space, crossed with the pcurves of the
 *   wires (`gbs/intersect/curves.h`): an odd number of crossings is inside. A ray through a
 *   vertex, tangent to a pcurve or along one is ambiguous: another direction is tried.
 * - Point in a solid: a ray from the point, crossed with the surface of every face
 *   (`gbs/intersect/curve_surface.h`), each hit being classified in its face; same parity rule,
 *   same retries when the ray meets an edge, grazes a face or runs along it.
 * The directions follow a fixed sequence: the result does not depend on chance.
 */
#include <gbs-brep/explore.h>
#include <gbs-brep/model.h>
#include <gbs/intersect/curve_surface.h>
#include <gbs/intersect/surface_surface.h>
#include <gbs/intersect/curves.h>
#include <gbs/curveline.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs::brep
{
    /// Position of a point against a face or a solid.
    enum class Location
    {
        In,
        Out,
        On
    };

    namespace classify_detail
    {
        /// Box of the pcurves of a face in its parameter space (the surface bounds for a natural face).
        template <std::floating_point T>
        auto face_uv_box(const Model<T> &m, FaceId f) -> std::array<T, 4>
        {
            const auto &face = m.face(f);
            Box<T, 2> b;
            for (auto w : face.wires)
                for (const auto &ce : m.wire(w).coedges)
                    if (ce.pcurve)
                    {
                        const auto &e = m.edge(ce.edge);
                        for (int k = 0; k <= 16; ++k)
                            b.add(ce.pcurve->value(e.u1 + (e.u2 - e.u1) * T(k) / T(16)));
                    }
            const auto sb = face.surface->bounds();
            if (b.empty())
                return sb;
            // a little margin, inside the bounds of the surface
            const T mu = T(1e-9) * (T(1) + b.max[0] - b.min[0]), mv = T(1e-9) * (T(1) + b.max[1] - b.min[1]);
            return {std::max(sb[0], b.min[0] - mu), std::min(sb[1], b.max[0] + mu), std::max(sb[2], b.min[1] - mv), std::min(sb[3], b.max[1] + mv)};
        }

        /// Directions of the rays: fixed, irregular angles.
        template <std::floating_point T>
        auto angle(int k) -> T
        {
            constexpr std::array<double, 8> a{0.1234, 1.9876, 3.5791, 5.0123, 0.7777, 2.6181, 4.3333, 5.8888};
            return T(a[static_cast<std::size_t>(k) % a.size()]);
        }

        template <std::floating_point T>
        auto directions_3d(int k) -> point<T, 3>
        {
            // points spread on the sphere (golden spiral), none along an axis
            const T z = T(1) - T(2) * (T(k) + T(0.5)) / T(16);
            const T r = std::sqrt(std::max(T(0), T(1) - z * z));
            const T phi = T(k) * T(2.39996322972865332) + T(0.31);
            return {r * std::cos(phi), r * std::sin(phi), z};
        }
    }

    /**
     * @brief Position of the parameters `uv` against the face `f`: `On` when the point of the
     * surface is within `tol` (and the tolerances of the edges) of the boundary.
     */
    template <std::floating_point T>
    [[nodiscard]] auto classify(const Model<T> &m, FaceId f, const point<T, 2> &uv, std::type_identity_t<T> tol = brep_default_tolerance<T>)
        -> Location
    {
        const auto &face = m.face(f);
        const auto box = classify_detail::face_uv_box(m, f);
        const auto sb = face.surface->bounds();
        // (a margin for the rounding of the bounds, as 2.4e-16 for 0 in a file)
        const T su = T(1e-12) * (sb[1] - sb[0]), sv = T(1e-12) * (sb[3] - sb[2]);
        if (uv[0] < sb[0] - su || uv[0] > sb[1] + su || uv[1] < sb[2] - sv || uv[1] > sb[3] + sv)
            return Location::Out;
        // on the boundary: the point of the surface near an edge
        const auto p = face.surface->value(uv[0], uv[1]);
        for (auto w : face.wires)
            for (const auto &ce : m.wire(w).coedges)
            {
                const auto &e = m.edge(ce.edge);
                const T t = std::max(tol, e.tol);
                if (e.degenerate)
                {
                    if (norm(m.vertex(e.v1).pnt - p) <= t)
                        return Location::On;
                    continue;
                }
                if (intersect_detail::project(*e.curve, e.u1, e.u2, p).second <= t)
                    return Location::On;
            }
        if (uv[0] < box[0] || uv[0] > box[1] || uv[1] < box[2] || uv[1] > box[3])
            return Location::Out;
        // parity of the crossings of a ray with the pcurves
        const T diag = std::hypot(box[1] - box[0], box[3] - box[2]);
        const T L = T(2) * diag + std::hypot(uv[0] - box[0], uv[1] - box[2]) + T(1);
        IntersectOptions<T> o;
        o.tol = T(1e-10) * (T(1) + diag);
        std::size_t crossings{};
        for (int k = 0; k < 8; ++k)
        {
            const T a = classify_detail::angle<T>(k);
            const BSCurve<T, 2> ray{points_vector<T, 2>{uv, point<T, 2>{uv[0] + L * std::cos(a), uv[1] + L * std::sin(a)}}, std::vector<T>{0, 0, 1, 1}, 1};
            bool ambiguous = false;
            crossings = 0;
            for (auto w : face.wires)
            {
                for (const auto &ce : m.wire(w).coedges)
                {
                    if (!ce.pcurve)
                        continue;
                    const auto &e = m.edge(ce.edge);
                    const auto r = intersect(*ce.pcurve, e.u1, e.u2, ray, T(0), T(1), o);
                    if (!r.overlaps.empty())
                        ambiguous = true;
                    for (const auto &h : r.hits)
                    {
                        // through a vertex (counted by both co-edges) or tangent: try another ray
                        const T span = e.u2 - e.u1;
                        if (h.tangent || h.t1 - e.u1 <= T(1e-9) * span || e.u2 - h.t1 <= T(1e-9) * span)
                            ambiguous = true;
                        ++crossings;
                    }
                    if (ambiguous)
                        break;
                }
                if (ambiguous)
                    break;
            }
            if (!ambiguous)
                break;
        }
        return crossings % 2 == 1 ? Location::In : Location::Out;
    }

    /**
     * @brief Parameters of points inside the face, away from its boundary: on a few lines of
     * constant v, the middle of the longest inside chord; the longest chords first.
     */
    template <std::floating_point T>
    [[nodiscard]] auto interior_points(const Model<T> &m, FaceId f, std::type_identity_t<T> tol = brep_default_tolerance<T>)
        -> std::vector<point<T, 2>>
    {
        const auto &face = m.face(f);
        const auto box = classify_detail::face_uv_box(m, f);
        const T du = box[1] - box[0];
        IntersectOptions<T> o;
        o.tol = T(1e-10) * (T(1) + du + (box[3] - box[2]));
        std::vector<std::pair<T, point<T, 2>>> found;
        for (T fr : {T(0.5), T(0.31), T(0.69), T(0.17), T(0.83), T(0.43), T(0.57)})
        {
            std::optional<point<T, 2>> best;
            T best_len{};
            const T v = box[2] + fr * (box[3] - box[2]);
            const BSCurve<T, 2> line{points_vector<T, 2>{{box[0] - du, v}, {box[1] + du, v}}, std::vector<T>{0, 0, 1, 1}, 1};
            std::vector<T> us;
            for (auto w : face.wires)
                for (const auto &ce : m.wire(w).coedges)
                    if (ce.pcurve)
                    {
                        const auto &e = m.edge(ce.edge);
                        for (const auto &h : intersect(*ce.pcurve, e.u1, e.u2, line, T(0), T(1), o).hits)
                            us.push_back(h.p[0]);
                    }
            std::ranges::sort(us);
            for (std::size_t k = 0; k + 1 < us.size(); ++k)
            {
                const T len = us[k + 1] - us[k];
                if (len <= best_len)
                    continue;
                const point<T, 2> c{(us[k] + us[k + 1]) / T(2), v};
                if (classify(m, f, c, tol) == Location::In)
                {
                    best = c;
                    best_len = len;
                }
            }
            if (best)
                found.emplace_back(best_len, *best);
        }
        std::ranges::stable_sort(found, [](const auto &a, const auto &b) { return a.first > b.first; });
        std::vector<point<T, 2>> r;
        for (const auto &[len, uv] : found)
            r.push_back(uv);
        return r;
    }

    /**
     * @brief Parameters of a point inside the face, away from its boundary: the middle of the
     * longest inside chord along a few lines of constant v; nullopt for a face without inside.
     */
    template <std::floating_point T>
    [[nodiscard]] auto interior_point(const Model<T> &m, FaceId f, std::type_identity_t<T> tol = brep_default_tolerance<T>)
        -> std::optional<point<T, 2>>
    {
        const auto r = interior_points(m, f, tol);
        if (r.empty())
            return std::nullopt;
        return r.front();
    }

    /**
     * @brief Position of the point `p` against the solid `s` (outer shell and cavities).
     */
    template <std::floating_point T>
    [[nodiscard]] auto classify(const Model<T> &m, SolidId s, const point<T, 3> &p, std::type_identity_t<T> tol = brep_default_tolerance<T>)
        -> Location
    {
        // (boxes sampled on the edges: grown by a part of their size, they are only a filter)
        auto bb = bounding_box(m, ShapeId{s});
        bb.inflate(T(0.05) * norm(bb.max - bb.min) + tol);
        if (!bb.contains(p))
            return Location::Out;
        const auto faces = explore<FaceId>(m, ShapeId{s});
        // on a face: a ray starting on a surface may miss it by rounding
        for (auto f : faces)
        {
            auto fb = bounding_box(m, ShapeId{f});
            fb.inflate(T(0.05) * norm(fb.max - fb.min) + T(10) * tol);
            if (!fb.contains(p))
                continue;
            const auto &srf = *m.face(f).surface;
            const auto box = classify_detail::face_uv_box(m, f);
            std::array<T, 3> q;
            if (const auto *es = dynamic_cast<const ElementarySurface<T> *>(&srf))
            {
                const auto uv = intersect_detail::primitive_uv(*es, box, p);
                q = {uv[0], uv[1], norm(es->value(uv[0], uv[1]) - p)};
            }
            else
                q = intersect_detail::project_on_surface(srf, box, p);
            if (q[2] <= tol && classify(m, f, point<T, 2>{q[0], q[1]}, tol) != Location::Out)
                return Location::On;
        }
        const T L = norm(bb.max - bb.min) + norm(p - bb.min) + T(1);
        IntersectOptions<T> o;
        o.tol = tol;
        std::size_t crossings{};
        for (int k = 0; k < 16; ++k)
        {
            const auto d = classify_detail::directions_3d<T>(k);
            const LineCurve<T> ray{p, d, T(0), L};
            bool ambiguous = false;
            crossings = 0;
            for (auto f : faces)
            {
                const auto &face = m.face(f);
                const auto box = classify_detail::face_uv_box(m, f);
                const auto r = intersect(ray, T(0), L, *face.surface, box, o);
                for (const auto &h : r.hits)
                {
                    const auto where = classify(m, f, point<T, 2>{h.u, h.v}, tol);
                    if (where == Location::Out)
                        continue;
                    if (h.t <= tol)
                        return Location::On; // the point lies on the face
                    if (where == Location::On || h.tangent)
                        ambiguous = true; // through an edge, or grazing the face
                    ++crossings;
                }
                for (const auto &pr : r.portions)
                {
                    // the ray runs along the surface: on the face if the point is on it
                    if (pr.ta <= tol && classify(m, f, point<T, 2>{pr.ua, pr.va}, tol) != Location::Out)
                        return Location::On;
                    ambiguous = true;
                }
                if (ambiguous)
                    break;
            }
            if (!ambiguous)
                break;
        }
        return crossings % 2 == 1 ? Location::In : Location::Out;
    }

    /**
     * @brief Position of the face `f` against the solid `s`, from interior points of the face:
     * the first one off the boundary of the solid (a face may touch the solid at a point or along
     * a line), `On` when all of them are on it.
     */
    template <std::floating_point T>
    [[nodiscard]] auto classify(const Model<T> &m, FaceId f, SolidId s, std::type_identity_t<T> tol = brep_default_tolerance<T>)
        -> std::optional<Location>
    {
        const auto uvs = interior_points(m, f, tol);
        if (uvs.empty())
            return std::nullopt;
        for (const auto &uv : uvs)
            if (const auto l = classify(m, s, m.face(f).surface->value(uv[0], uv[1]), tol); l != Location::On)
                return l;
        return Location::On;
    }
}
