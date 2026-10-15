#pragma once
/**
 * @file section.h
 * @brief Section of two shapes, or of a shape by a surface: the edges where their faces meet,
 * each with its pcurve on both faces, chained into wires.
 *
 * 1. Pairs of faces whose bounding boxes overlap.
 * 2. Intersection of their surfaces on the boxes of their pcurves (`gbs/intersect/surface_surface.h`).
 * 3. Restriction to the faces: each curve is cut where its pcurves cross the boundaries of the
 *    faces (curve-curve in the parameter spaces), and the parts whose middle is inside both faces
 *    (`classify.h`) become edges.
 * 4. Vertices: the ends of the edges are merged within the tolerance with the vertices of the
 *    shapes and with each other: a section edge ending on the boundary of a face shares the
 *    vertex of its neighbours.
 * 5. Wires: the edges are chained by their vertices; a branching component stays as edges.
 */
#include <gbs-brep/builders.h>
#include <gbs-brep/classify.h>
#include <gbs-brep/explore.h>
#include <gbs/intersect/surface_surface.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs::brep
{
    template <std::floating_point T>
    struct SectionOptions
    {
        T tol{brep_default_tolerance<T>};  ///< confusion distance; edge tolerances grow to the deviation of their curves
        IntersectOptions<T> intersect{};   ///< options of the surface intersections (tol replaced by `tol`)
    };

    /// An edge of a section and the faces that carry it.
    template <std::floating_point T>
    struct SectionEdge
    {
        EdgeId edge;
        FaceId face_a;                        ///< face of the first shape
        FaceId face_b;                        ///< face of the second shape (invalid for a section by a surface)
        std::shared_ptr<Curve<T, 2>> pcurve_a; ///< on face_a, same parameter as the edge
        std::shared_ptr<Curve<T, 2>> pcurve_b; ///< on face_b, or on the cutting surface
        bool tangent{};                        ///< the faces are tangent along the edge
    };

    template <std::floating_point T>
    struct SectionResult
    {
        CompoundId compound;                              ///< wires and edges of the section
        std::vector<WireId> wires;
        std::vector<SectionEdge<T>> edges;                ///< new edges, inside both faces
        std::vector<EdgeId> boundary_edges;               ///< edges of the shapes lying in the section (along a face boundary)
        std::vector<std::pair<EdgeId, FaceId>> edges_on_faces; ///< an edge of one shape and a face of the other it lies on
        std::vector<VertexId> points;                     ///< isolated contact points
        std::vector<std::pair<FaceId, FaceId>> coincident; ///< pairs of faces on the same surface (face_b invalid for a surface)
        bool singular{};                                  ///< an intersection met a singular point
        bool complete{true};                              ///< false when an intersection was capped
    };

    namespace section_detail
    {
        /// Vertices reused by the section: those of the shapes, then the new ones.
        template <std::floating_point T>
        struct VertexPool
        {
            Model<T> &m;
            std::vector<VertexId> ids;
            T tol;

            /// The vertex within the confusion radius of p (10 tol: the ends of interpolated
            /// intersection curves), or a new one; and its distance to p.
            std::pair<VertexId, T> at(const point<T, 3> &p)
            {
                VertexId best;
                T dmin = std::numeric_limits<T>::max();
                for (auto v : ids)
                {
                    const auto &x = m.vertex(v);
                    const T d = norm(x.pnt - p);
                    if (d <= std::max(T(10) * tol, x.tol) && d < dmin)
                    {
                        dmin = d;
                        best = v;
                    }
                }
                if (best.valid())
                    return {best, dmin};
                const auto v = m.add(Vertex<T>{p, tol});
                ids.push_back(v);
                return {v, T(0)};
            }
        };

        /**
         * Parameters of a curve where its pcurve crosses the boundary of a face, and the
         * intervals where it runs along it.
         */
        template <std::floating_point T>
        void boundary_crossings(const Model<T> &m, FaceId f, const Curve<T, 2> &pc, T a, T b, std::vector<T> &ts,
                                std::vector<std::array<T, 2>> &along)
        {
            const auto box = classify_detail::face_uv_box(m, f);
            IntersectOptions<T> o;
            o.tol = T(1e-9) * (T(1) + std::hypot(box[1] - box[0], box[3] - box[2]));
            for (auto w : m.face(f).wires)
                for (const auto &ce : m.wire(w).coedges)
                {
                    if (!ce.pcurve)
                        continue;
                    const auto &e = m.edge(ce.edge);
                    const auto r = intersect(pc, a, b, *ce.pcurve, e.u1, e.u2, o);
                    for (const auto &h : r.hits)
                        ts.push_back(h.t1);
                    for (const auto &ov : r.overlaps)
                    {
                        ts.push_back(ov.t1a);
                        ts.push_back(ov.t1b);
                        along.push_back({std::min(ov.t1a, ov.t1b), std::max(ov.t1a, ov.t1b)});
                    }
                }
        }

        /// The edge of the face passing through p.
        template <std::floating_point T>
        auto boundary_edge(const Model<T> &m, FaceId f, const point<T, 3> &p, T tol) -> std::optional<EdgeId>
        {
            for (auto w : m.face(f).wires)
                for (const auto &ce : m.wire(w).coedges)
                {
                    const auto &e = m.edge(ce.edge);
                    if (!e.degenerate && intersect_detail::project(*e.curve, e.u1, e.u2, p).second <= std::max(tol, e.tol))
                        return ce.edge;
                }
            return std::nullopt;
        }

        template <std::floating_point T>
        auto face_box(const Model<T> &m, FaceId f, T tol) -> BoundingBox<T>
        {
            // sampled on the edges: grown by a part of its size, it is only a filter
            auto b = bounding_box(m, ShapeId{f});
            b.inflate(T(0.05) * norm(b.max - b.min) + tol);
            return b;
        }

        template <std::floating_point T>
        bool boxes_overlap(const BoundingBox<T> &a, const BoundingBox<T> &b)
        {
            for (int i = 0; i < 3; ++i)
                if (a.min[i] > b.max[i] || b.min[i] > a.max[i])
                    return false;
            return true;
        }

        /**
         * Section of the faces `fa` (of the first shape) and `fb` (of the second, invalid for a
         * surface `tool` on `tool_box`).
         */
        template <std::floating_point T>
        void section_faces(Model<T> &m, FaceId fa, FaceId fb, const Surface<T, 3> *tool, const std::array<T, 4> &tool_box,
                           const SectionOptions<T> &opts, VertexPool<T> &pool, SectionResult<T> &res)
        {
            const T tol = opts.tol;
            const auto &sa = *m.face(fa).surface;
            const auto &sb = tool ? *tool : *m.face(fb).surface;
            const auto box_a = classify_detail::face_uv_box(m, fa);
            const auto box_b = tool ? tool_box : classify_detail::face_uv_box(m, fb);
            auto o = opts.intersect;
            o.tol = tol;
            const auto ss = intersect(sa, box_a, sb, box_b, o);
            res.singular = res.singular || ss.singular;
            res.complete = res.complete && ss.complete;
            if (ss.coincident)
            {
                res.coincident.emplace_back(fa, fb);
                return;
            }
            auto inside_b = [&](const point<T, 2> &uv) { return tool ? Location::In : classify(m, fb, uv, tol); };
            for (const auto &ic : ss.curves)
            {
                std::vector<T> ts;
                std::vector<std::array<T, 2>> along;
                boundary_crossings(m, fa, *ic.pc1, ic.a, ic.b, ts, along);
                if (!tool)
                    boundary_crossings(m, fb, *ic.pc2, ic.a, ic.b, ts, along);
                // parameters closer than the tolerance along the curve are one cut (both ends of a
                // closed curve are the same point, not the same cut); a crossing of a face boundary
                // wins over an end of the curve, which may lie on the margin of a parameter box
                auto same = [&](T x, T y) {
                    return norm(ic.curve->value(x) - ic.curve->value(y)) <= tol && norm(ic.curve->value((x + y) / T(2)) - ic.curve->value(y)) <= tol;
                };
                for (T end : {ic.a, ic.b})
                    if (std::ranges::none_of(ts, [&](T x) { return same(x, end); }))
                        ts.push_back(end);
                std::ranges::sort(ts);
                std::vector<T> cuts;
                for (T t : ts)
                    if (cuts.empty() || !same(t, cuts.back()))
                        cuts.push_back(t);
                for (std::size_t k = 0; k + 1 < cuts.size(); ++k)
                {
                    const T t0 = cuts[k], t1 = cuts[k + 1], tm = (t0 + t1) / T(2);
                    const auto p0 = ic.curve->value(t0), p1 = ic.curve->value(t1);
                    if (norm(p1 - p0) <= T(10) * tol && norm(ic.curve->value(tm) - p0) <= T(10) * tol)
                        continue;
                    const auto la = classify(m, fa, ic.pc1->value(tm), tol), lb = inside_b(ic.pc2->value(tm));
                    if (la == Location::Out || lb == Location::Out)
                        continue;
                    if (la == Location::On || lb == Location::On)
                    {
                        // along the boundary of a face: the section is an edge of the shapes there;
                        // a piece near the boundary but not along it is a contact at a point
                        if (std::ranges::none_of(along, [&](const auto &iv) { return iv[0] <= tm && tm <= iv[1]; }))
                            continue;
                        const auto pm = ic.curve->value(tm);
                        if (auto e = boundary_edge(m, la == Location::On ? fa : fb, pm, tol))
                            if (std::ranges::find(res.boundary_edges, *e) == res.boundary_edges.end())
                                res.boundary_edges.push_back(*e);
                        // the edges and the faces of the other shape they lie on (both faces when
                        // their boundaries run along each other)
                        auto on_face = [&](FaceId own, FaceId other) {
                            if (!other.valid())
                                return;
                            if (auto e = boundary_edge(m, own, pm, tol))
                                if (std::ranges::find(res.edges_on_faces, std::pair{*e, other}) == res.edges_on_faces.end())
                                    res.edges_on_faces.emplace_back(*e, other);
                        };
                        if (la == Location::On)
                            on_face(fa, fb);
                        if (lb == Location::On && !tool)
                            on_face(fb, fa);
                        continue;
                    }
                    const auto [v0, d0] = pool.at(p0);
                    const auto [v1, d1] = norm(p1 - p0) <= tol ? std::pair{v0, norm(p1 - m.vertex(v0).pnt)} : pool.at(p1);
                    // the tolerance of the edge covers its deviation and the distance to its vertices
                    const T etol = std::max({tol, ic.deviation, d0, d1});
                    auto e = make_edge(m, ic.curve, t0, t1, v0, v1, etol);
                    if (!e)
                        continue; // shorter than the tolerance
                    res.edges.push_back({*e, fa, fb, ic.pc1, ic.pc2, ic.tangent});
                }
            }
            for (const auto &h : ss.points)
                if (classify(m, fa, point<T, 2>{h.u1, h.v1}, tol) != Location::Out &&
                    inside_b(point<T, 2>{h.u2, h.v2}) != Location::Out)
                {
                    const auto v = pool.at(h.p).first;
                    if (std::ranges::find(res.points, v) == res.points.end())
                        res.points.push_back(v);
                }
        }

        /// Chains the section edges into wires, by connected components of their vertices.
        template <std::floating_point T>
        void chain(Model<T> &m, SectionResult<T> &res, T tol)
        {
            std::map<std::uint32_t, std::size_t> parent_of;
            std::vector<std::size_t> parent;
            auto id = [&](VertexId v) {
                auto [it, fresh] = parent_of.try_emplace(v.index, parent.size());
                if (fresh)
                    parent.push_back(parent.size());
                return it->second;
            };
            std::function<std::size_t(std::size_t)> root = [&](std::size_t x) { return parent[x] == x ? x : parent[x] = root(parent[x]); };
            std::vector<EdgeId> all;
            for (const auto &se : res.edges)
                all.push_back(se.edge);
            all.insert(all.end(), res.boundary_edges.begin(), res.boundary_edges.end());
            for (auto e : all)
                parent[root(id(m.edge(e).v1))] = root(id(m.edge(e).v2));
            std::map<std::size_t, std::vector<EdgeId>> comps;
            for (auto e : all)
                comps[root(id(m.edge(e).v1))].push_back(e);
            std::vector<ShapeId> shapes;
            for (auto &[r, edges] : comps)
            {
                auto w = make_wire(m, edges, tol);
                if (w)
                {
                    res.wires.push_back(*w);
                    shapes.push_back(*w);
                }
                else
                    for (auto e : edges)
                        shapes.push_back(e);
            }
            for (auto v : res.points)
                shapes.push_back(v);
            res.compound = unwrap(make_compound(m, shapes));
        }
    }

    /**
     * @brief Section of the shapes `a` and `b`: the edges where their faces meet, with their
     * pcurves on both faces, chained into wires (one compound). The shapes are not modified,
     * except that the section edges reuse their vertices.
     */
    template <std::floating_point T>
    [[nodiscard]] auto section(Model<T> &m, const ShapeId &a, const ShapeId &b, std::type_identity_t<SectionOptions<T>> opts = {})
        -> BuildResult<SectionResult<T>>
    {
        if (!m.alive(a) || !m.alive(b))
            return detail::fail(BuildErrc::InvalidId, "section: dead or invalid shape");
        SectionResult<T> res;
        const auto fas = explore<FaceId>(m, a), fbs = explore<FaceId>(m, b);
        section_detail::VertexPool<T> pool{m, {}, opts.tol};
        for (const auto &s : {a, b})
            for (auto v : explore<VertexId>(m, s))
                pool.ids.push_back(v);
        std::vector<BoundingBox<T>> bb_b;
        for (auto f : fbs)
            bb_b.push_back(section_detail::face_box(m, f, opts.tol));
        for (auto fa : fas)
        {
            const auto bb_a = section_detail::face_box(m, fa, opts.tol);
            for (std::size_t k = 0; k < fbs.size(); ++k)
                if (fbs[k] != fa && section_detail::boxes_overlap(bb_a, bb_b[k]))
                    section_detail::section_faces(m, fa, fbs[k], static_cast<const Surface<T, 3> *>(nullptr), {}, opts, pool, res);
        }
        section_detail::chain(m, res, opts.tol);
        return res;
    }

    /**
     * @brief Section of the shape `a` by the surface `tool` on its parameter box `tool_box`
     * (its bounds by default): a plane of cut, a contour of machining…
     */
    template <std::floating_point T>
    [[nodiscard]] auto section(Model<T> &m, const ShapeId &a, std::type_identity_t<std::shared_ptr<Surface<T, 3>>> tool,
                               std::optional<std::array<T, 4>> tool_box = std::nullopt, std::type_identity_t<SectionOptions<T>> opts = {})
        -> BuildResult<SectionResult<T>>
    {
        if (!m.alive(a))
            return detail::fail(BuildErrc::InvalidId, "section: dead or invalid shape");
        if (!tool)
            return detail::fail(BuildErrc::NullSurface, "section: no cutting surface");
        SectionResult<T> res;
        section_detail::VertexPool<T> pool{m, {}, opts.tol};
        for (auto v : explore<VertexId>(m, a))
            pool.ids.push_back(v);
        const auto box = tool_box.value_or(tool->bounds());
        for (auto fa : explore<FaceId>(m, a))
            section_detail::section_faces(m, fa, FaceId{}, tool.get(), box, opts, pool, res);
        section_detail::chain(m, res, opts.tol);
        return res;
    }
}
