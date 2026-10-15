#pragma once
/**
 * @file split.h
 * @brief Split of faces by edges crossing them, and of shapes by a surface or another shape.
 *
 * A face is split in its parameter space:
 * 1. The cuts are edges with a pcurve on the face (section edges). Where a cut ends inside a
 *    boundary edge, that edge is split at the vertex of the cut (`split_edge`): the neighbour face
 *    gets the same pieces and the shell stays closed.
 * 2. Half-edges: each boundary co-edge once, in the sense that keeps the face on its left; each
 *    cut twice. Nodes are the vertices at a position in (u, v): the two sides of a seam are two
 *    nodes of the same vertex.
 * 3. Cuts with a free end (dangling) are removed and reported.
 * 4. Loops: from each half-edge, the next one is the first outgoing half-edge clockwise from the
 *    way back, which keeps the smallest region on the left. Loops of positive area are the outer
 *    boundaries of the new faces, those of negative area their holes, each given to the smallest
 *    outer loop containing it.
 * 5. The new faces share the surface of the face; they replace it in its shells, with its sense.
 */
#include <gbs-brep/builders.h>
#include <gbs-brep/classify.h>
#include <gbs-brep/closure.h>
#include <gbs-brep/edit.h>
#include <gbs-brep/section.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <numbers>
#include <optional>
#include <span>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs::brep
{
    /// An edge cutting a face, with its pcurve on the face (parametrized like the edge).
    template <std::floating_point T>
    struct FaceCut
    {
        EdgeId edge;
        std::shared_ptr<Curve<T, 2>> pcurve;
    };

    struct SplitFaceResult
    {
        std::vector<FaceId> faces;     ///< the pieces (the face itself when nothing splits it)
        std::vector<EdgeId> dangling;  ///< cuts with a free end, not used
    };

    namespace split_detail
    {
        template <std::floating_point T>
        struct HalfEdge
        {
            EdgeId edge;
            std::shared_ptr<Curve<T, 2>> pc;
            T t0{}, t1{};         ///< from t0 to t1 on the edge (and its pcurve)
            VertexId from, to;
            std::size_t n0{}, n1{}; ///< nodes
            bool cut{};
            bool alive{true};
        };

        /// Direction leaving the start of a pcurve between t0 and t1, at two scales.
        template <std::floating_point T>
        auto departure(const Curve<T, 2> &pc, T t0, T t1, T frac) -> T
        {
            const auto a = pc.value(t0), b = pc.value(t0 + frac * (t1 - t0));
            return std::atan2(b[1] - a[1], b[0] - a[0]);
        }

        /// Splits the edge e at u and puts the vertex v at the split (instead of a new one); the new second half.
        template <std::floating_point T>
        auto split_at_vertex(Model<T> &m, EdgeId e, T u, VertexId v) -> std::optional<EdgeId>
        {
            auto r = split_edge(m, e, u);
            if (!r)
                return std::nullopt;
            const auto created = r->vertex;
            m.edge(r->first).v2 = v;
            m.edge(r->second).v1 = v;
            auto &vx = m.vertex(v);
            vx.tol = std::max(vx.tol, norm(vx.pnt - m.vertex(created).pnt) + m.vertex(created).tol);
            m.erase(created);
            return r->second;
        }
    }

    /**
     * @brief Splits the face `f` by the edges `cuts` (each with its pcurve on the face).
     *
     * Boundary edges are split where cuts end on them (in every wire using them). The face is
     * replaced in its shells by its pieces, with the same sense, and erased; a face that the cuts
     * do not split is kept.
     */
    template <std::floating_point T>
    [[nodiscard]] auto split_face(Model<T> &m, FaceId f, std::span<const FaceCut<T>> cuts, std::type_identity_t<T> tol = brep_default_tolerance<T>)
        -> BuildResult<SplitFaceResult>
    {
        using namespace split_detail;
        if (!m.alive(f))
            return detail::fail(BuildErrc::InvalidId, "split_face: dead or invalid face", {f});
        for (const auto &c : cuts)
            if (!m.alive(c.edge) || !c.pcurve)
                return detail::fail(BuildErrc::InvalidId, "split_face: dead cut edge or cut without pcurve", {c.edge});
        SplitFaceResult r;
        if (cuts.empty())
        {
            r.faces.push_back(f);
            return r;
        }

        // ---- 0. tolerances of the cuts: vertex >= edge >= face
        for (const auto &c : cuts)
        {
            auto &e = m.edge(c.edge);
            e.tol = std::max(e.tol, m.face(f).tol);
            for (auto v : {e.v1, e.v2})
                m.vertex(v).tol = std::max(m.vertex(v).tol, e.tol);
        }

        // ---- 1. boundary edges split where cuts end on them
        auto boundary_vertices = [&] {
            std::vector<VertexId> vs;
            for (auto w : m.face(f).wires)
                for (const auto &ce : m.wire(w).coedges)
                {
                    vs.push_back(m.edge(ce.edge).v1);
                    vs.push_back(m.edge(ce.edge).v2);
                }
            return vs;
        };
        for (const auto &c : cuts)
            for (auto v : {m.edge(c.edge).v1, m.edge(c.edge).v2})
            {
                const auto bv = boundary_vertices();
                if (std::ranges::find(bv, v) != bv.end())
                    continue;
                const auto p = m.vertex(v).pnt;
                std::optional<std::pair<EdgeId, T>> on;
                for (auto w : m.face(f).wires)
                    for (const auto &ce : m.wire(w).coedges)
                    {
                        const auto &e = m.edge(ce.edge);
                        if (e.degenerate)
                            continue;
                        const auto [u, d] = intersect_detail::project(*e.curve, e.u1, e.u2, p);
                        if (d <= std::max({tol, e.tol, m.vertex(v).tol}) && !on)
                            on = std::pair{ce.edge, u};
                    }
                if (on)
                    split_at_vertex(m, on->first, on->second, v);
            }

        // ---- 2. half-edges and nodes
        const auto &face = m.face(f);
        const auto box = classify_detail::face_uv_box(m, f);
        const T eps_uv = T(1e-3) * std::hypot(box[1] - box[0], box[3] - box[2]);
        std::vector<HalfEdge<T>> hes;
        // the sense of the boundary: the outer wire must turn counter-clockwise in (u, v)
        auto loop_polygon = [&](const std::vector<HalfEdge<T>> &loop) {
            std::vector<point<T, 2>> poly;
            for (const auto &h : loop)
                for (int k = 0; k < 8; ++k)
                    poly.push_back(h.pc->value(h.t0 + (h.t1 - h.t0) * T(k) / T(8)));
            return poly;
        };
        bool flip = false;
        {
            std::vector<HalfEdge<T>> outer;
            for (const auto &ce : m.wire(face.wires.front()).coedges)
            {
                const auto &e = m.edge(ce.edge);
                const bool fwd = ce.orient == Orientation::Forward;
                outer.push_back({ce.edge, ce.pcurve, fwd ? e.u1 : e.u2, fwd ? e.u2 : e.u1});
            }
            flip = signed_area(loop_polygon(outer)) < T(0);
        }
        for (auto w : face.wires)
            for (const auto &ce : m.wire(w).coedges)
            {
                const auto &e = m.edge(ce.edge);
                if (!ce.pcurve)
                    return detail::fail(BuildErrc::InvalidId, "split_face: boundary co-edge without pcurve", {ce.edge});
                bool fwd = ce.orient == Orientation::Forward;
                if (flip)
                    fwd = !fwd;
                hes.push_back({ce.edge, ce.pcurve, fwd ? e.u1 : e.u2, fwd ? e.u2 : e.u1, fwd ? e.v1 : e.v2, fwd ? e.v2 : e.v1});
            }
        const auto n_boundary = hes.size();
        for (const auto &c : cuts)
        {
            const auto &e = m.edge(c.edge);
            hes.push_back({c.edge, c.pcurve, e.u1, e.u2, e.v1, e.v2, 0, 0, true});
            hes.push_back({c.edge, c.pcurve, e.u2, e.u1, e.v2, e.v1, 0, 0, true});
        }
        struct Node
        {
            VertexId v;
            point<T, 2> uv;
        };
        std::vector<Node> nodes;
        auto node_at = [&](VertexId v, const point<T, 2> &uv) {
            for (std::size_t k = 0; k < nodes.size(); ++k)
                if (nodes[k].v == v && norm(nodes[k].uv - uv) <= eps_uv)
                    return k;
            nodes.push_back({v, uv});
            return nodes.size() - 1;
        };
        for (auto &h : hes)
        {
            h.n0 = node_at(h.from, h.pc->value(h.t0));
            h.n1 = node_at(h.to, h.pc->value(h.t1));
        }

        // ---- 3. dangling cuts: a node reached by a single cut
        for (bool changed = true; changed;)
        {
            changed = false;
            std::vector<std::size_t> degree(nodes.size(), 0);
            for (const auto &h : hes)
                if (h.alive)
                    ++degree[h.n0];
            for (std::size_t k = n_boundary; k < hes.size(); k += 2)
                if (hes[k].alive && (degree[hes[k].n0] <= 1 || degree[hes[k].n1] <= 1))
                {
                    hes[k].alive = hes[k + 1].alive = false;
                    r.dangling.push_back(hes[k].edge);
                    changed = true;
                }
        }

        // ---- 4. loops
        std::vector<std::vector<std::size_t>> out(nodes.size());
        for (std::size_t k = 0; k < hes.size(); ++k)
            if (hes[k].alive)
                out[hes[k].n0].push_back(k);
        auto leave = [&](std::size_t k, T frac) { return departure(*hes[k].pc, hes[k].t0, hes[k].t1, frac); };
        auto back = [&](std::size_t k, T frac) { return departure(*hes[k].pc, hes[k].t1, hes[k].t0, frac); };
        const T two_pi = T(2) * std::numbers::pi_v<T>;
        auto next_of = [&](std::size_t k) -> std::optional<std::size_t> {
            const auto &cands = out[hes[k].n1];
            std::optional<std::size_t> best;
            T best_turn{};
            for (std::size_t c : cands)
            {
                // clockwise turn from the way back to the candidate, in (0, 2 pi]; ties broken farther along
                T turn{};
                for (T frac : {T(1e-6), T(1e-2), T(0.2)})
                {
                    turn = std::fmod(back(k, frac) - leave(c, frac) + T(2) * two_pi, two_pi);
                    if (turn > T(1e-7) && turn < two_pi - T(1e-7))
                        break;
                }
                if (turn <= T(1e-7))
                    turn = two_pi; // the way back itself: last resort
                if (!best || turn < best_turn)
                {
                    best = c;
                    best_turn = turn;
                }
            }
            return best;
        };
        std::vector<bool> used(hes.size(), false);
        std::vector<std::vector<std::size_t>> loops;
        for (std::size_t s = 0; s < hes.size(); ++s)
        {
            if (!hes[s].alive || used[s])
                continue;
            std::vector<std::size_t> loop;
            std::size_t k = s;
            bool ok = true;
            for (std::size_t guard = 0; guard <= hes.size(); ++guard)
            {
                if (used[k])
                {
                    ok = k == s;
                    break;
                }
                used[k] = true;
                loop.push_back(k);
                const auto nx = next_of(k);
                if (!nx)
                {
                    ok = false;
                    break;
                }
                k = *nx;
            }
            if (!ok)
                return detail::fail(BuildErrc::Disconnected, "split_face: the cuts and the boundary do not form closed loops", {f});
            loops.push_back(std::move(loop));
        }

        // ---- 5. outer loops and holes
        std::vector<std::vector<point<T, 2>>> polys;
        std::vector<T> areas;
        for (const auto &lp : loops)
        {
            std::vector<HalfEdge<T>> l;
            for (auto k : lp)
                l.push_back(hes[k]);
            polys.push_back(loop_polygon(l));
            areas.push_back(signed_area(polys.back()));
        }
        std::vector<std::size_t> outers;
        for (std::size_t k = 0; k < loops.size(); ++k)
            if (areas[k] > T(0))
                outers.push_back(k);
        if (outers.size() <= 1 && r.dangling.size() * 2 == hes.size() - n_boundary)
        {
            r.faces.push_back(f); // nothing splits the face
            return r;
        }
        std::vector<std::vector<std::size_t>> holes(loops.size());
        for (std::size_t k = 0; k < loops.size(); ++k)
        {
            if (areas[k] > T(0))
                continue;
            std::optional<std::size_t> owner;
            for (auto o : outers)
                if (inside(polys[o], polys[k].front()) && (!owner || areas[o] < areas[*owner]))
                    owner = o;
            if (!owner)
                return detail::fail(BuildErrc::HoleOutsideOuter, "split_face: a loop of the split is outside every outer loop", {f});
            holes[*owner].push_back(k);
        }

        // ---- 6. faces
        auto make_wire_of = [&](std::size_t l) {
            Wire<T> w;
            w.closed = true;
            for (auto k : loops[l])
            {
                const auto &h = hes[k];
                w.coedges.push_back({h.edge, h.t0 < h.t1 ? Orientation::Forward : Orientation::Reversed, h.pc});
            }
            return m.add(std::move(w));
        };
        const auto old = m.face(f);
        const auto name = std::string(m.name(f));
        for (auto o : outers)
        {
            Face<T> nf{old.surface, {make_wire_of(o)}, old.tol, false};
            for (auto h : holes[o])
                nf.wires.push_back(make_wire_of(h));
            const auto id = m.add(std::move(nf));
            if (!name.empty())
                m.setName(id, name);
            r.faces.push_back(id);
        }
        for (auto sid : m.template ids<ShellId>())
        {
            auto &uses = m.shell(sid).faces;
            for (std::size_t i = 0; i < uses.size(); ++i)
                if (uses[i].face == f)
                {
                    const auto orient = uses[i].orient;
                    uses.erase(uses.begin() + static_cast<std::ptrdiff_t>(i));
                    for (std::size_t k = 0; k < r.faces.size(); ++k)
                        uses.insert(uses.begin() + static_cast<std::ptrdiff_t>(i + k), FaceUse{r.faces[k], orient});
                    i += r.faces.size() - 1;
                }
        }
        for (auto w : old.wires)
            m.erase(w);
        m.erase(f);
        return r;
    }

    template <std::floating_point T>
    [[nodiscard]] auto split_face(Model<T> &m, FaceId f, const std::vector<FaceCut<T>> &cuts, std::type_identity_t<T> tol = brep_default_tolerance<T>)
        -> BuildResult<SplitFaceResult>
    {
        return split_face(m, f, std::span<const FaceCut<T>>{cuts}, tol);
    }

    template <std::floating_point T>
    struct SplitResult
    {
        std::vector<std::pair<FaceId, std::vector<FaceId>>> faces; ///< each split face and its pieces
        std::vector<EdgeId> dangling;
        SectionResult<T> section;
    };

    namespace split_detail
    {
        template <std::floating_point T>
        auto split_by_section(Model<T> &m, SectionResult<T> sec, T tol) -> BuildResult<SplitResult<T>>
        {
            std::map<std::uint32_t, std::vector<FaceCut<T>>> by_face;
            for (const auto &se : sec.edges)
                by_face[se.face_a.index].push_back({se.edge, se.pcurve_a});
            SplitResult<T> r;
            for (auto &[index, cuts] : by_face)
            {
                const FaceId f{index};
                auto s = split_face(m, f, cuts, tol);
                if (!s)
                    return std::unexpected(s.error());
                r.dangling.insert(r.dangling.end(), s->dangling.begin(), s->dangling.end());
                if (s->faces.size() > 1 || s->faces.front() != f)
                    r.faces.emplace_back(f, s->faces);
            }
            r.section = std::move(sec);
            return r;
        }
    }

    /**
     * @brief Splits the faces of `shape` along its section by `tool` (a shape): the shells and
     * solids keep their uses, a split face being replaced by its pieces. `tool` is not modified,
     * except that the section reuses its vertices.
     */
    template <std::floating_point T>
    [[nodiscard]] auto split(Model<T> &m, const ShapeId &shape, const ShapeId &tool, std::type_identity_t<SectionOptions<T>> opts = {})
        -> BuildResult<SplitResult<T>>
    {
        auto sec = section(m, shape, tool, opts);
        if (!sec)
            return std::unexpected(sec.error());
        return split_detail::split_by_section(m, std::move(*sec), opts.tol);
    }

    /// Splits the faces of `shape` along its section by the surface `tool` (on `tool_box`, its bounds by default).
    template <std::floating_point T>
    [[nodiscard]] auto split(Model<T> &m, const ShapeId &shape, std::type_identity_t<std::shared_ptr<Surface<T, 3>>> tool,
                             std::optional<std::array<T, 4>> tool_box = std::nullopt, std::type_identity_t<SectionOptions<T>> opts = {})
        -> BuildResult<SplitResult<T>>
    {
        auto sec = section(m, shape, std::move(tool), tool_box, opts);
        if (!sec)
            return std::unexpected(sec.error());
        return split_detail::split_by_section(m, std::move(*sec), opts.tol);
    }
}
