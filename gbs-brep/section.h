#pragma once
/**
 * @file section.h
 * @brief Section of two shapes, or of a shape by a surface: the edges where their faces meet,
 * each with its pcurve on both faces, chained into wires.
 *
 * 1. Pairs of faces whose bounding boxes overlap.
 * 2. Intersection of their surfaces on the boxes of their pcurves (`gbs/intersect/surface_surface.h`);
 *    faces on the same surface are reported as coincident (closed forms for the primitives, points
 *    of each face projected on the surface of the other for the other pairs).
 * 3. Restriction to the faces: each curve is cut where its pcurves cross the boundaries of the
 *    faces (curve-curve in the parameter spaces), and where curves cross at a branch point (the
 *    spans of the marching stepping over it, crossed with each other); the parts whose middle is
 *    inside both faces (`classify.h`) become edges.
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
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <ranges>
#include <span>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs::brep
{
    /// The face `face` is intersected with `other` on these boxes of its parameters only.
    template <std::floating_point T>
    struct SectionRestriction
    {
        FaceId face, other;
        std::vector<std::array<T, 4>> boxes;
    };

    template <std::floating_point T>
    struct SectionOptions
    {
        T tol{brep_default_tolerance<T>};  ///< confusion distance; edge tolerances grow to the deviation of their curves
        IntersectOptions<T> intersect{};   ///< options of the surface intersections (tol replaced by `tol`)
        /// pairs of faces intersected on parts of the parameters of the first only (no part: not
        /// at all), the pieces of a partial coincidence (`partial.h`)
        std::vector<SectionRestriction<T>> restricted{};
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
            /// intersection curves), or a new one; and its distance to p. (Vertices merged away
            /// since they were pooled, by an imprint, are skipped.)
            std::pair<VertexId, T> at(const point<T, 3> &p)
            {
                VertexId best;
                T dmin = std::numeric_limits<T>::max();
                for (auto v : ids)
                {
                    if (!m.alive(v))
                        continue;
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
         * intervals where it runs along it. False when an intersection was capped (a crossing
         * may be missing).
         */
        template <std::floating_point T>
        [[nodiscard]] bool boundary_crossings(const Model<T> &m, FaceId f, const Curve<T, 2> &pc, T a, T b, std::vector<T> &ts,
                                              std::vector<std::array<T, 2>> &along)
        {
            const auto box = classify_detail::face_uv_box(m, f);
            IntersectOptions<T> o;
            o.tol = T(1e-9) * (T(1) + std::hypot(box[1] - box[0], box[3] - box[2]));
            bool complete = true;
            for (auto w : m.face(f).wires)
                for (const auto &ce : m.wire(w).coedges)
                {
                    if (!ce.pcurve)
                        continue;
                    const auto &e = m.edge(ce.edge);
                    const auto r = intersect(pc, a, b, *ce.pcurve, e.u1, e.u2, o);
                    complete = complete && r.complete;
                    for (const auto &h : r.hits)
                        ts.push_back(h.t1);
                    for (const auto &ov : r.overlaps)
                    {
                        ts.push_back(ov.t1a);
                        ts.push_back(ov.t1b);
                        along.push_back({std::min(ov.t1a, ov.t1b), std::max(ov.t1a, ov.t1b)});
                    }
                }
            return complete;
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

        /**
         * The boundary edges of the face f that the curve c runs along on [t0, t1], from a vertex
         * of the boundary to another: its ends within the confusion radius of these vertices (10
         * tol, or their tolerance), its inner samples within that of the edges (10 tol, or the
         * tolerance of the edge). Empty otherwise: a curve touching the boundary, or running along
         * a part of an edge, is not that boundary.
         */
        template <std::floating_point T>
        auto edges_along(const Model<T> &m, FaceId f, const Curve<T, 3> &c, T t0, T t1, T tol) -> std::vector<EdgeId>
        {
            const auto &wires = m.face(f).wires;
            auto at_vertex = [&](const point<T, 3> &p) {
                for (auto w : wires)
                    for (const auto &ce : m.wire(w).coedges)
                        for (auto v : {m.edge(ce.edge).v1, m.edge(ce.edge).v2})
                            if (norm(m.vertex(v).pnt - p) <= std::max(T(10) * tol, m.vertex(v).tol))
                                return true;
                return false;
            };
            if (!at_vertex(c.value(t0)) || !at_vertex(c.value(t1)))
                return {};
            constexpr int n = 6;
            std::vector<EdgeId> r;
            for (int k = 1; k < n; ++k)
            {
                const auto p = c.value(t0 + (t1 - t0) * T(k) / T(n));
                std::optional<EdgeId> near;
                for (auto w : wires)
                    for (const auto &ce : m.wire(w).coedges)
                    {
                        const auto &e = m.edge(ce.edge);
                        if (!near && !e.degenerate && intersect_detail::project(*e.curve, e.u1, e.u2, p).second <= std::max(T(10) * tol, e.tol))
                            near = ce.edge;
                    }
                if (!near)
                    return {};
                if (!std::ranges::contains(r, *near))
                    r.push_back(*near);
            }
            return r;
        }

        /**
         * Points of the face `from` projected on the surface `to` (on `to_box`): nullopt when one
         * of them is off the surface inside its box (the surfaces differ), otherwise the number of
         * them on the surface. Points projecting onto the border of the box are beyond `to`, and
         * ignored without being classified in `from`.
         */
        template <std::floating_point T>
        auto on_surface(const Model<T> &m, FaceId from, const Surface<T, 3> &to, const std::array<T, 4> &to_box, T tol) -> std::optional<int>
        {
            const auto &s = *m.face(from).surface;
            const auto box = classify_detail::face_uv_box(m, from);
            const T eu = T(1e-7) * (to_box[1] - to_box[0]), ev = T(1e-7) * (to_box[3] - to_box[2]);
            int on{};
            constexpr int n = 6;
            for (int j = 0; j < n; ++j)
                for (int i = 0; i < n; ++i)
                {
                    const point<T, 2> uv{box[0] + (box[1] - box[0]) * (T(i) + T(0.5)) / T(n), box[2] + (box[3] - box[2]) * (T(j) + T(0.5)) / T(n)};
                    const auto q = intersect_detail::project_on_surface(to, to_box, s.value(uv[0], uv[1]));
                    if (q[0] <= to_box[0] + eu || q[0] >= to_box[1] - eu || q[1] <= to_box[2] + ev || q[1] >= to_box[3] - ev)
                        continue;
                    if (classify(m, from, uv, tol) != Location::In)
                        continue;
                    if (q[2] > tol)
                        return std::nullopt;
                    ++on;
                }
            return on;
        }

        /// Verdict of the coincidence test of two faces.
        enum class Coincidence
        {
            No,       ///< a point of one face is off the surface of the other
            Yes,      ///< the faces face each other, and lie on each other's surface
            Undecided ///< no point of either face faces the other (faces touching along their borders, or apart)
        };

        /**
         * Faces on the same surface, for the pairs the closed forms do not cover (B-splines,
         * a primitive and a B-spline of the same shape): the points of each face lie on the
         * surface of the other wherever they face it.
         */
        template <std::floating_point T>
        auto coincident_faces(const Model<T> &m, FaceId fa, FaceId fb, T tol) -> Coincidence
        {
            const auto &sa = *m.face(fa).surface, &sb = *m.face(fb).surface;
            const auto ab = on_surface(m, fa, sb, classify_detail::face_uv_box(m, fb), tol);
            if (!ab)
                return Coincidence::No;
            const auto ba = on_surface(m, fb, sa, classify_detail::face_uv_box(m, fa), tol);
            if (!ba)
                return Coincidence::No;
            return *ab + *ba >= 3 ? Coincidence::Yes : Coincidence::Undecided;
        }

        /**
         * A boundary edge of the face `f` runs along the surface `to` (on `to_box`, its border
         * included), over a part at least: 3 of 9 samples of the edge on it (shifted faces touch
         * along a part of an edge).
         */
        template <std::floating_point T>
        bool border_on(const Model<T> &m, FaceId f, const Surface<T, 3> &to, const std::array<T, 4> &to_box, T tol)
        {
            for (auto w : m.face(f).wires)
                for (const auto &ce : m.wire(w).coedges)
                {
                    const auto &e = m.edge(ce.edge);
                    if (e.degenerate)
                        continue;
                    int on{};
                    for (int k = 0; k <= 8; ++k)
                    {
                        const auto p = e.curve->value(e.u1 + (e.u2 - e.u1) * T(k) / T(8));
                        on += intersect_detail::project_on_surface(to, to_box, p)[2] <= std::max(tol, e.tol);
                    }
                    if (on >= 3)
                        return true;
                }
            return false;
        }

        /**
         * Faces of the same surface touching along their borders (the coplanar sides of stacked
         * boxes), whose intersection is singular: a border of one lies on the surface of the
         * other, and every curve found runs along a border. Faces of different surfaces tangent
         * along a curve inside one of them do not pass.
         */
        template <std::floating_point T>
        bool touching_along_borders(const Model<T> &m, FaceId fa, FaceId fb, const std::array<T, 4> &box_a, const std::array<T, 4> &box_b,
                                    const SurfaceSurfaceResult<T> &ss, T tol)
        {
            const auto &sa = *m.face(fa).surface, &sb = *m.face(fb).surface;
            if (!border_on(m, fa, sb, box_b, tol) && !border_on(m, fb, sa, box_a, tol))
                return false;
            return std::ranges::all_of(ss.curves, [&](const IntersectionCurve<T> &ic) {
                const T tm = (ic.a + ic.b) / T(2);
                return classify(m, fa, ic.pc1->value(tm), tol) == Location::On || classify(m, fb, ic.pc2->value(tm), tol) == Location::On;
            });
        }

        template <std::floating_point T>
        auto face_box(const Model<T> &m, FaceId f, T tol) -> BoundingBox<T>
        {
            // sampled on the faces: grown by a part of its size, it is only a filter
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
            const auto restriction = std::ranges::find_if(opts.restricted, [&](const SectionRestriction<T> &r) {
                return !tool && ((r.face == fa && r.other == fb) || (r.face == fb && r.other == fa));
            });
            if (restriction != opts.restricted.end() && restriction->boxes.empty())
                return;
            const auto &sa = *m.face(fa).surface;
            const auto &sb = tool ? *tool : *m.face(fb).surface;
            const auto box_a = classify_detail::face_uv_box(m, fa);
            const auto box_b = tool ? tool_box : classify_detail::face_uv_box(m, fb);
            auto o = opts.intersect;
            o.tol = tol;
            const auto same = !tool && !(dynamic_cast<const ElementarySurface<T> *>(&sa) && dynamic_cast<const ElementarySurface<T> *>(&sb))
                                  ? coincident_faces(m, fa, fb, tol)
                                  : Coincidence::No;
            if (same == Coincidence::Yes)
            {
                res.coincident.emplace_back(fa, fb);
                return;
            }
            // (restricted: on each box of the face, within its parameters, the curves gathered)
            auto intersect_pair = [&] {
                if (restriction == opts.restricted.end())
                    return intersect(sa, box_a, sb, box_b, o);
                SurfaceSurfaceResult<T> all;
                const bool on_a = restriction->face == fa;
                for (auto r : restriction->boxes)
                {
                    const auto &whole = on_a ? box_a : box_b;
                    r = {std::max(r[0], whole[0]), std::min(r[1], whole[1]), std::max(r[2], whole[2]), std::min(r[3], whole[3])};
                    if (!(r[0] < r[1]) || !(r[2] < r[3]))
                        continue;
                    auto part = on_a ? intersect(sa, r, sb, box_b, o) : intersect(sa, box_a, sb, r, o);
                    all.singular = all.singular || part.singular;
                    all.complete = all.complete && part.complete;
                    all.coincident = all.coincident || part.coincident;
                    all.curves.append_range(part.curves | std::views::as_rvalue);
                    all.points.append_range(part.points | std::views::as_rvalue);
                }
                return all;
            };
            const auto ss = intersect_pair();
            if (same == Coincidence::Undecided && ss.singular && touching_along_borders(m, fa, fb, box_a, box_b, ss, tol))
            {
                res.coincident.emplace_back(fa, fb);
                return;
            }
            res.singular = res.singular || ss.singular;
            res.complete = res.complete && ss.complete;
            if (ss.coincident)
            {
                res.coincident.emplace_back(fa, fb);
                return;
            }
            auto inside_b = [&](const point<T, 2> &uv) { return tool ? Location::In : classify(m, fb, uv, T(2) * tol); };
            // the curves crossing at a branch point are cut there: each span stepping over a branch
            // point is crossed with the others (of the same curve or of another one)
            std::vector<std::vector<T>> branch_cuts(ss.curves.size());
            {
                struct Span
                {
                    std::size_t curve;
                    T a, b;
                };
                std::vector<Span> spans;
                for (std::size_t i = 0; i < ss.curves.size(); ++i)
                    for (const auto &sp : ss.curves[i].branch_spans)
                        spans.push_back({i, sp[0], sp[1]});
                IntersectOptions<T> o;
                o.tol = T(100) * tol; // two interpolations of the curves, each within its deviation
                for (std::size_t i = 0; i < spans.size(); ++i)
                    for (std::size_t j = i + 1; j < spans.size(); ++j)
                    {
                        const auto &x = spans[i], &y = spans[j];
                        const auto &cx = *ss.curves[x.curve].curve, &cy = *ss.curves[y.curve].curve;
                        const auto r = intersect(cx, x.a, x.b, cy, y.a, y.b, o);
                        res.complete = res.complete && r.complete;
                        for (const auto &h : r.hits)
                        {
                            branch_cuts[x.curve].push_back(h.t1);
                            branch_cuts[y.curve].push_back(h.t2);
                        }
                    }
            }
            for (std::size_t ci = 0; ci < ss.curves.size(); ++ci)
            {
                const auto &ic = ss.curves[ci];
                std::vector<T> ts = branch_cuts[ci];
                std::vector<std::array<T, 2>> along;
                const bool done_a = boundary_crossings(m, fa, *ic.pc1, ic.a, ic.b, ts, along);
                const bool done_b = tool || boundary_crossings(m, fb, *ic.pc2, ic.a, ic.b, ts, along);
                res.complete = res.complete && done_a && done_b;
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
                    // (a point within the tolerance of two surfaces crossing at right angles is up
                    // to sqrt(2) tol from their common border: twice the tolerance to find it On)
                    const auto la = classify(m, fa, ic.pc1->value(tm), T(2) * tol), lb = inside_b(ic.pc2->value(tm));
                    if (la == Location::Out || lb == Location::Out)
                        continue;
                    // the section runs along edges of the shapes: those of one face are reported once,
                    // and the edges of both with the face of the other shape they lie on
                    auto on_boundary = [&](std::span<const EdgeId> own, std::span<const EdgeId> others) {
                        for (auto e : own.empty() ? others : own)
                            if (!std::ranges::contains(res.boundary_edges, e))
                                res.boundary_edges.push_back(e);
                        auto lies_on = [&](std::span<const EdgeId> es, FaceId other) {
                            for (auto e : es)
                                if (other.valid() && !std::ranges::contains(res.edges_on_faces, std::pair{e, other}))
                                    res.edges_on_faces.emplace_back(e, other);
                        };
                        lies_on(own, fb);
                        lies_on(others, fa);
                    };
                    if ((la == Location::On || lb == Location::On) &&
                        std::ranges::any_of(along, [&](const auto &iv) { return iv[0] <= tm && tm <= iv[1]; }))
                    {
                        // along the boundary of a face (of both when their boundaries run along each other)
                        const auto pm = ic.curve->value(tm);
                        std::vector<EdgeId> ea, eb;
                        if (auto e = la == Location::On ? boundary_edge(m, fa, pm, tol) : std::nullopt)
                            ea.push_back(*e);
                        if (auto e = lb == Location::On && !tool ? boundary_edge(m, fb, pm, tol) : std::nullopt)
                            eb.push_back(*e);
                        on_boundary(ea, eb);
                        continue;
                    }
                    // from a vertex of the boundary to another within the confusion radius of the edges
                    // between them (faces a few tolerances apart), the piece is these edges: a new edge
                    // would share their vertices
                    const auto ea = edges_along(m, fa, *ic.curve, t0, t1, tol);
                    const auto eb = tool ? std::vector<EdgeId>{} : edges_along(m, fb, *ic.curve, t0, t1, tol);
                    on_boundary(ea, eb);
                    // a piece near the boundary but not along it is a contact at a point
                    if (!ea.empty() || !eb.empty() || la == Location::On || lb == Location::On)
                        continue;
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
                    if (!std::ranges::contains(res.points, v))
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
            auto root = [&](this auto &self, std::size_t x) -> std::size_t { return parent[x] == x ? x : parent[x] = self(parent[x]); };
            auto all = res.edges | std::views::transform(&SectionEdge<T>::edge) | std::ranges::to<std::vector>();
            all.append_range(res.boundary_edges);
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
