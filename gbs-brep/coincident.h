#pragma once
/**
 * @file coincident.h
 * @brief Imprint of two coincident faces (on the same surface) on each other, before a boolean
 * operation splits them; imprint of an edge lying on a face.
 *
 * 1. The boundary edges of both faces are crossed (curve-curve in 3D): an edge is split where an
 *    edge of the other face crosses it, starts or ends on it; vertices closer than the tolerance
 *    are merged.
 * 2. Overlapping edges, which have the same vertices after step 1, are merged into one: the
 *    co-edges of the absorbed edge use the kept one, with a pcurve computed on their face.
 * 3. The boundary edges of each face lying inside the other face are cuts of the other face, with
 *    their pcurve on it.
 *
 * The pcurves of the edges on the faces come from the closed forms on primitive surfaces and from
 * the projection of the edges on the other surfaces (`pcurve.h`).
 */
#include <gbs-brep/classify.h>
#include <gbs-brep/closure.h>
#include <gbs-brep/edit.h>
#include <gbs-brep/pcurve.h>
#include <gbs-brep/explore.h>
#include <gbs-brep/split.h>
#include <gbs/intersect/curves.h>
#include <gbs/intersect/surface_surface.h>

#include <algorithm>
#include <cmath>
#include <expected>
#include <memory>
#include <numbers>
#include <optional>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs::brep
{
    namespace coincident_detail
    {
        /**
         * Pcurve of the edge e on the face f, whose surface closes as `cl`, and its deviation: on
         * a primitive, exact when linear and interpolated otherwise; on another surface, projected
         * and interpolated within the approximation tolerance of the pcurves (the parametrization
         * of a rational surface may have kinks that a smooth pcurve only approaches). The error
         * says why there is none: the edge leaves the surface (`OffSurface`), crosses its seam
         * (`CrossesSeam`, split it first) or the projection did not converge.
         */
        template <std::floating_point T>
        auto edge_pcurve(const Model<T> &m, EdgeId e, FaceId f, T tol, const SurfaceClosure &cl) -> std::expected<PCurveFit<T>, PCurveErrc>
        {
            const auto &srf = *m.face(f).surface;
            const auto &ed = m.edge(e);
            if (auto exact = plane_pcurve(ed, srf, std::max(tol, ed.tol)))
                return *exact;
            const auto *es = dynamic_cast<const ElementarySurface<T> *>(&srf);
            if (!es)
            {
                MakeFaceOptions<T> opts;
                opts.tol = tol;
                opts.pcurve_tol = std::max({tol, ed.tol, brep_pcurve_approx_tol<T>});
                return project_pcurve(ed, srf, cl, std::optional<point<T, 2>>{}, opts);
            }
            auto box = classify_detail::face_uv_box(m, f);
            // the box of the face, grown a little: the edge may run along its border
            const T du = T(1e-6) * (T(1) + box[1] - box[0]), dv = T(1e-6) * (T(1) + box[3] - box[2]);
            box = {box[0] - du, box[1] + du, box[2] - dv, box[3] + dv};
            const auto [pc, dev, exact] = intersect_detail::primitive_pcurve(*ed.curve, ed.u1, ed.u2, *es, box, tol);
            if (!pc || dev > std::max(tol, ed.tol))
                return std::unexpected(PCurveErrc::OffSurface);
            return PCurveFit<T>{pc, dev, exact};
        }

        /// Raises the tolerance of the edge e (and of its vertices) to the deviation of a pcurve.
        template <std::floating_point T>
        void cover_deviation(Model<T> &m, EdgeId e, T dev)
        {
            auto &ed = m.edge(e);
            if (dev <= ed.tol)
                return;
            ed.tol = dev;
            for (auto v : {ed.v1, ed.v2})
                m.vertex(v).tol = std::max(m.vertex(v).tol, dev);
        }

        /// Periods of the parameters of a surface closing as `cl` (0 for a parameter that does not).
        template <std::floating_point T>
        auto periods(const Surface<T, 3> &s, const SurfaceClosure &cl) -> std::array<T, 2>
        {
            const T pi2 = T(2) * std::numbers::pi_v<T>;
            if (const auto *es = dynamic_cast<const ElementarySurface<T> *>(&s))
            {
                const auto per = intersect_detail::periodic(*es);
                return {per[0] ? pi2 : T(0), per[1] ? pi2 : T(0)};
            }
            const auto b = s.bounds();
            return {cl.closed_u ? b[1] - b[0] : T(0), cl.closed_v ? b[3] - b[2] : T(0)};
        }

        /// Parameter of the edge e where it crosses the seam of the surface s (closing as `cl`).
        template <std::floating_point T>
        auto seam_crossing(const Model<T> &m, EdgeId e, const Surface<T, 3> &s, const SurfaceClosure &cl, T tol) -> std::optional<T>
        {
            const auto &ed = m.edge(e);
            const auto b = s.bounds();
            for (int i : {0, 1})
            {
                if (!(i == 0 ? cl.closed_u : cl.closed_v))
                    continue;
                const T period = b[2 * i + 1] - b[2 * i];
                // the parameter jumps by about a period where the edge goes across the seam
                auto par = [&](T t) { return detail::project_point_global(s, ed.curve->value(t), tol * T(1e-3))[i]; };
                constexpr int n = 64;
                T t0 = ed.u1, x0 = par(t0);
                for (int k = 1; k <= n; ++k)
                {
                    const T t1 = ed.u1 + (ed.u2 - ed.u1) * T(k) / T(n), x1 = par(t1);
                    if (std::abs(x1 - x0) > period / T(2))
                    {
                        T a = t0, xa = x0, c = t1;
                        for (int it = 0; it < 60; ++it)
                        {
                            const T mid = (a + c) / T(2), xm = par(mid);
                            if (std::abs(xm - xa) > period / T(2))
                                c = mid;
                            else
                            {
                                a = mid;
                                xa = xm;
                            }
                        }
                        return (a + c) / T(2);
                    }
                    t0 = t1;
                    x0 = x1;
                }
            }
            return std::nullopt;
        }

        /// Splits the edge e where it crosses the seam of the face f, at the vertex `vertex_at` gives; the new half.
        template <std::floating_point T, typename VertexAt>
        auto split_at_seam(Model<T> &m, EdgeId e, FaceId f, T tol, const SurfaceClosure &cl, VertexAt &vertex_at) -> std::optional<EdgeId>
        {
            const auto t = seam_crossing(m, e, *m.face(f).surface, cl, tol);
            if (!t)
                return std::nullopt;
            const auto v = vertex_at(m.edge(e).curve->value(*t));
            return split_detail::split_at_vertex(m, e, *t, v);
        }

        /**
         * Cuts of the face `into` made by the edges `edges` (of another face lying on it): each
         * edge gets its pcurve on `into` and is a cut when its middle is inside. An edge crossing
         * the seam of `into` is split there first; an edge leaving the surface (a B-spline ends,
         * or the edge only touches it) is not a cut. False when the projection does not converge.
         */
        template <std::floating_point T, typename VertexAt>
        bool cuts_on(Model<T> &m, std::vector<EdgeId> edges, FaceId into, T tol, VertexAt &vertex_at, std::vector<FaceCut<T>> &out)
        {
            const auto cl = surface_closure(*m.face(into).surface, tol);
            for (std::size_t k = 0; k < edges.size(); ++k)
            {
                const auto e = edges[k];
                if (!m.alive(e))
                    continue;
                const auto pc = edge_pcurve(m, e, into, tol, cl);
                if (!pc)
                {
                    if (pc.error() == PCurveErrc::CrossesSeam)
                    {
                        const auto second = split_at_seam(m, e, into, tol, cl, vertex_at);
                        if (!second)
                            return false;
                        edges.push_back(e);
                        edges.push_back(*second);
                        continue;
                    }
                    // off the surface somewhere (a B-spline ends, or the edge only touches it): not a cut
                    if (pc.error() == PCurveErrc::OffSurface)
                        continue;
                    return false;
                }
                cover_deviation(m, e, pc->deviation);
                const auto &ed = m.edge(e);
                if (classify(m, into, pc->pcurve->value((ed.u1 + ed.u2) / T(2)), tol) == Location::In)
                    out.push_back({e, pc->pcurve});
            }
            return true;
        }

        /// Replaces the vertex `from` by `to` in every edge of the model.
        template <std::floating_point T>
        void replace_vertex(Model<T> &m, VertexId from, VertexId to)
        {
            if (from == to)
                return;
            for (auto e : m.template ids<EdgeId>())
            {
                auto &ed = m.edge(e);
                if (ed.v1 == from)
                    ed.v1 = to;
                if (ed.v2 == from)
                    ed.v2 = to;
            }
            auto &vt = m.vertex(to);
            vt.tol = std::max(vt.tol, norm(vt.pnt - m.vertex(from).pnt) + m.vertex(from).tol);
            m.erase(from);
        }

        /// Replaces the edge `from` by `to` (same ends) in every wire of the model, with pcurves on their faces.
        template <std::floating_point T>
        bool replace_edge(Model<T> &m, EdgeId from, EdgeId to, T tol)
        {
            const auto &ef = m.edge(from), &et = m.edge(to);
            const bool same = ef.v1 == et.v1;
            const auto mid = et.curve->value((et.u1 + et.u2) / T(2));
            const T t_old = intersect_detail::project(*ef.curve, ef.u1, ef.u2, mid).first;
            for (auto f : m.template ids<FaceId>())
            {
                std::optional<SurfaceClosure> cl; // computed once per face using `from`
                for (auto w : m.face(f).wires)
                    for (auto &ce : m.wire(w).coedges)
                        if (ce.edge == from)
                        {
                            if (!cl)
                                cl = surface_closure(*m.face(f).surface, tol);
                            const auto r = edge_pcurve(m, to, f, tol, *cl);
                            if (!r || !ce.pcurve)
                                return false;
                            cover_deviation(m, to, r->deviation);
                            // keep the side of the seam of the former pcurve: shift by whole periods
                            const auto per = periods(*m.face(f).surface, *cl);
                            const auto a = ce.pcurve->value(t_old), b = r->pcurve->value((et.u1 + et.u2) / T(2));
                            point<T, 2> shift{};
                            for (int i = 0; i < 2; ++i)
                                if (per[i] > T(0))
                                    shift[i] = per[i] * std::round((a[i] - b[i]) / per[i]);
                            auto pc = shifted_pcurve(r->pcurve, shift);
                            if (!pc)
                                return false;
                            ce = CoEdge<T>{to, same ? ce.orient : reverse(ce.orient), pc};
                        }
            }
            auto &t = m.edge(to);
            t.tol = std::max(t.tol, ef.tol);
            m.erase(from);
            return true;
        }

        template <std::floating_point T>
        auto boundary_edges(const Model<T> &m, FaceId f) -> std::vector<EdgeId>
        {
            std::vector<EdgeId> r;
            for (auto w : m.face(f).wires)
                for (const auto &ce : m.wire(w).coedges)
                    if (!m.edge(ce.edge).degenerate && !std::ranges::contains(r, ce.edge))
                        r.push_back(ce.edge);
            return r;
        }

        /// The vertex at p: an end of the edge, or nullopt.
        template <std::floating_point T>
        auto end_at(const Model<T> &m, EdgeId e, const point<T, 3> &p, T tol) -> std::optional<VertexId>
        {
            const auto &ed = m.edge(e);
            for (auto v : {ed.v1, ed.v2})
                if (norm(m.vertex(v).pnt - p) <= std::max(tol, m.vertex(v).tol))
                    return v;
            return std::nullopt;
        }

        /**
         * Crossing of the edges ea and eb: the first point where they cross, or where one starts
         * or ends on the other, that needs a split or a merge of vertices is applied. Returns true
         * when something changed (`ok` false on failure); `second_a` gets the new half of ea when
         * ea is split.
         */
        template <std::floating_point T, typename VertexAt>
        bool imprint_pair(Model<T> &m, EdgeId ea, EdgeId eb, T tol, VertexAt &vertex_at, bool &ok, std::optional<EdgeId> &second_a)
        {
            IntersectOptions<T> o;
            o.tol = tol;
            const auto &a = m.edge(ea), &b = m.edge(eb);
            const auto r = intersect(*a.curve, a.u1, a.u2, *b.curve, b.u1, b.u2, o);
            if (!r.complete)
            {
                ok = false; // capped: a crossing may be missing
                return false;
            }
            std::vector<std::pair<point<T, 3>, std::array<T, 2>>> pts;
            for (const auto &h : r.hits)
                pts.push_back({h.p, {h.t1, h.t2}});
            for (const auto &ov : r.overlaps)
            {
                pts.push_back({a.curve->value(ov.t1a), {ov.t1a, ov.t2a}});
                pts.push_back({a.curve->value(ov.t1b), {ov.t1b, ov.t2b}});
            }
            for (const auto &[p, t] : pts)
            {
                const auto va = end_at(m, ea, p, tol), vb = end_at(m, eb, p, tol);
                if (va && vb)
                {
                    if (*va != *vb)
                    {
                        replace_vertex(m, *vb, *va);
                        return true;
                    }
                    continue;
                }
                // the vertex of the crossing: an existing end, or a new one
                const auto v = va ? *va : (vb ? *vb : vertex_at(p));
                if (!va)
                {
                    second_a = split_detail::split_at_vertex(m, ea, t[0], v);
                    if (!second_a)
                    {
                        ok = false;
                        return false;
                    }
                }
                if (!vb && !split_detail::split_at_vertex(m, eb, t[1], v))
                {
                    ok = false;
                    return false;
                }
                return true;
            }
            return false;
        }

        /**
         * One step of step 1 on the faces fa, fb: the first crossing of two of their boundary
         * edges that needs a split or a merge of vertices is applied. Returns false when there is
         * none left (or on failure, with `ok` set to false).
         */
        template <std::floating_point T, typename VertexAt>
        bool imprint_step(Model<T> &m, FaceId fa, FaceId fb, T tol, VertexAt &vertex_at, bool &ok)
        {
            const auto ea_list = boundary_edges(m, fa), eb_list = boundary_edges(m, fb);
            for (auto ea : ea_list)
                for (auto eb : eb_list)
                {
                    if (ea == eb)
                        continue;
                    std::optional<EdgeId> second;
                    if (imprint_pair(m, ea, eb, tol, vertex_at, ok, second))
                        return true;
                    if (!ok)
                        return false;
                }
            return false;
        }

        /// Merges the edge `e` with an edge of the boundary of `f` having the same ends and middle; true if merged.
        template <std::floating_point T>
        auto merge_overlapping(Model<T> &m, EdgeId e, FaceId f, T tol) -> std::optional<bool>
        {
            for (auto b : boundary_edges(m, f))
            {
                if (b == e)
                    continue;
                const auto &a = m.edge(e), &eb = m.edge(b);
                const bool ends = (a.v1 == eb.v1 && a.v2 == eb.v2) || (a.v1 == eb.v2 && a.v2 == eb.v1);
                if (!ends)
                    continue;
                const auto mid = eb.curve->value((eb.u1 + eb.u2) / T(2));
                if (intersect_detail::project(*a.curve, a.u1, a.u2, mid).second > std::max(tol, a.tol))
                    continue; // two different arcs between the same vertices
                if (!replace_edge(m, b, e, tol))
                    return std::nullopt;
                return true;
            }
            return false;
        }
    }

    /**
     * @brief Imprints the coincident faces fa and fb on each other (see the file comment).
     * `vertex_at(p)` gives the vertex at a crossing (an existing one within the tolerance, or a
     * new one). Returns the cuts of each face, or nullopt when a pcurve cannot be computed.
     */
    template <std::floating_point T, typename VertexAt>
    [[nodiscard]] auto imprint_coincident(Model<T> &m, FaceId fa, FaceId fb, T tol, VertexAt vertex_at)
        -> std::optional<std::pair<std::vector<FaceCut<T>>, std::vector<FaceCut<T>>>>
    {
        using namespace coincident_detail;
        // 1. crossings: split, merge vertices; one change per step, a bound on the steps that the
        //    crossings of the edges do not reach (reached: nullopt, rather than an imprint left
        //    unfinished)
        bool ok = true;
        const std::size_t max_steps = 1000 + 8 * boundary_edges(m, fa).size() * boundary_edges(m, fb).size();
        for (std::size_t steps = 0; imprint_step(m, fa, fb, tol, vertex_at, ok);)
            if (++steps == max_steps)
                return std::nullopt;
        if (!ok)
            return std::nullopt;
        // 2. overlapping edges: same ends and same middle
        for (bool merged = true; merged;)
        {
            merged = false;
            for (auto ea : boundary_edges(m, fa))
            {
                const auto r = merge_overlapping(m, ea, fb, tol);
                if (!r)
                    return std::nullopt;
                if ((merged = *r))
                    break;
            }
        }
        // 3. edges of one face inside the other: cuts
        std::pair<std::vector<FaceCut<T>>, std::vector<FaceCut<T>>> cuts;
        auto inner = [&](FaceId from, FaceId into, std::vector<FaceCut<T>> &out) {
            const auto own = boundary_edges(m, into);
            std::vector<EdgeId> edges;
            for (auto e : boundary_edges(m, from))
                if (!std::ranges::contains(own, e))
                    edges.push_back(e);
            return cuts_on(m, std::move(edges), into, tol, vertex_at, out);
        };
        if (!inner(fb, fa, cuts.first) || !inner(fa, fb, cuts.second))
            return std::nullopt;
        return cuts;
    }

    /**
     * @brief Imprints the edge `e`, lying on the face `f` of another shape, on that
     * face: the edge and the boundary of the face are split where they cross or where one ends on
     * the other, coincident vertices are merged, the pieces of the edge running along the
     * boundary replace the boundary edges they overlap, and the pieces inside the face are
     * returned as its cuts. `vertex_at(p)` gives the vertex at a crossing. nullopt when a pcurve
     * cannot be computed on the face.
     */
    template <std::floating_point T, typename VertexAt>
    [[nodiscard]] auto imprint_edge(Model<T> &m, EdgeId e, FaceId f, T tol, VertexAt vertex_at) -> std::optional<std::vector<FaceCut<T>>>
    {
        using namespace coincident_detail;
        std::vector<EdgeId> pieces{e};
        bool ok = true;
        const std::size_t max_steps = 1000 + 8 * boundary_edges(m, f).size();
        for (std::size_t steps = 0;; ++steps)
        {
            if (steps == max_steps)
                return std::nullopt; // not finished
            bool changed = false;
            for (std::size_t k = 0; k < pieces.size() && !changed; ++k)
            {
                if (!m.alive(pieces[k]))
                    continue;
                for (auto b : boundary_edges(m, f))
                {
                    if (b == pieces[k] || !m.alive(pieces[k]))
                        continue;
                    std::optional<EdgeId> second;
                    changed = imprint_pair(m, pieces[k], b, tol, vertex_at, ok, second);
                    if (!ok)
                        return std::nullopt;
                    if (second)
                        pieces.push_back(*second);
                    if (changed)
                        break;
                }
            }
            if (!changed)
                break;
        }
        std::vector<EdgeId> inside;
        for (auto p : pieces)
        {
            if (!m.alive(p))
                continue;
            const auto merged = merge_overlapping(m, p, f, tol);
            if (!merged)
                return std::nullopt;
            if (*merged)
                continue;
            if (!std::ranges::contains(boundary_edges(m, f), p))
                inside.push_back(p);
        }
        std::vector<FaceCut<T>> cuts;
        if (!cuts_on(m, std::move(inside), f, tol, vertex_at, cuts))
            return std::nullopt;
        return cuts;
    }
}
