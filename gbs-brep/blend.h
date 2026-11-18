#pragma once
/**
 * @file blend.h
 * @brief Fillets and chamfers of edges of a solid (`fillet`, `chamfer`): the faces of each edge
 * cut along the contacts of its blend, the faces at its ends cut by the blend, the face of the
 * blend between them; a new solid (stage 4, `docs/sources/design/brep_palier4.md`, § 4.2 and
 * § 4.5; note `brep4_pr06_fillet_topology.md`).
 *
 * The solid is copied, and the copy edited:
 * 1. **Geometry.** The blend of each edge (`fillet_geometry`, `chamfer_geometry`): its surface
 *    and its contacts, continued past the vertices of the edge.
 * 2. **Ends.** At each vertex of an open edge, the blend ends:
 *    - **free**, at a vertex of three edges where no other edge is blended: the third face F3
 *      cuts it. The contacts meet F3 on the other two edges (Gauss-Newton), which are split
 *      there; the **end curve** is the intersection of the blend and of F3 between them, marched
 *      across the blend (Gauss-Newton at each point, the points doubled until the interpolated
 *      curve and its pcurves are within tol / 10 of the exact points between them);
 *    - **joined** to the blend of the next edge of a **chain**, a blended edge going on
 *      tangentially at the vertex: both blends end on the same arc across them at the vertex
 *      (their contacts meet there, checked); a contact passing onto the next face splits the edge
 *      between both faces (that edge must run through it);
 *    - **vanishing**, where the faces become tangent at the vertex (`fillet_geometry` ends the
 *      fillet there): its contacts meet at the vertex, the face of the fillet ends in a point;
 *    - a closed edge has no end: the face of the blend runs its whole loop, across a seam.
 * 3. **Faces.** In the wire of each face of an edge, the run of co-edges between the ends of the
 *    contacts that holds the blended edge is replaced by the contacts (chained along a chain of
 *    edges, from one end to the other; a closed chain clear of the boundary replaces a whole
 *    wire made of blended edges); in the wire of F3, the run through the vertex by the end curve.
 *    Convex and concave edges alike: the face loses the band along the edge (convex) or the
 *    corner it gets back (concave, F3).
 * 4. The face of each blend, on its surface, bounded by its contacts and its ends, oriented as the
 *    faces of its edge; the edges and vertices left out erased; the solid checked (`check`).
 *
 * Refused (`NonGeneric`): two blended edges meeting at a vertex but tangentially (the vertex
 * blends of PR 7), a free end at a vertex of more than three edges, a contact that does not meet
 * the face at an end within the blend continued (`FilletOptions::extension`) or meets it beyond
 * the edge between them (a blend wider than a face next to it), a contact ending inside a face,
 * a result that is not a valid solid; the refusals of the geometry. On failure the model is left
 * as it was.
 */
#include <gbs-brep/check.h>
#include <gbs-brep/closure.h>
#include <gbs-brep/edit.h>
#include <gbs-brep/fillet.h>
#include <gbs-brep/pcurve.h>
#include <gbs-brep/split.h>
#include <gbs/bscinterp.h>
#include <gbs/curveonsurface.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs::brep
{
    namespace blend_detail
    {
        template <std::floating_point T>
        using Vec = thick_detail::Vec<T>;
        template <std::floating_point T>
        using Mat = thick_detail::Mat<T>;

        template <std::floating_point T>
        [[noreturn]] void refuse(const std::string &op, const std::string &why, std::vector<ShapeId> shapes = {})
        {
            throw BuildFailure(BuildError{BuildErrc::NonGeneric, op + ": " + why, std::move(shapes)});
        }

        /// How a blend ends at a vertex of its edge.
        enum class EndKind
        {
            Free,   ///< cut by the third face at the vertex
            Joint,  ///< joined to the blend of the next edge of a chain
            Vanish, ///< the contacts meet at the vertex
            Loop    ///< a closed edge: the seam of the face of the blend
        };

        /// An end of a blend: its contacts there and the edge across it.
        template <std::floating_point T>
        struct End
        {
            EndKind kind{EndKind::Free};
            VertexId vertex;                       ///< the vertex of the edge
            T t{};                                 ///< its parameter on the blend
            std::array<VertexId, 2> v;             ///< the ends of the contacts
            std::array<T, 2> tc{};                 ///< the parameters of the contacts there
            EdgeId edge;                           ///< the edge across (none: a vanishing end whose contacts meet in (u, v))
            std::shared_ptr<Curve<T, 2>> on_blend; ///< its pcurve on the blend, parametrized like it
            bool forward{true};                    ///< it runs from the contact 0 of the blend to its contact 1
            // a free end: the third face, the edges cut and the end curve on it
            FaceId f3;
            std::shared_ptr<Curve<T, 2>> on_f3;
        };

        template <std::floating_point T>
        struct Blend
        {
            EdgeBlend<T> g;
            std::array<End<T>, 2> ends; ///< at a (lower parameter) and at b
            std::array<EdgeId, 2> contacts;
        };

        /// The solid being edited: its shells and the orientation of each face in them.
        template <std::floating_point T>
        struct SolidView
        {
            const Model<T> &m;
            std::vector<ShellId> shells;

            SolidView(const Model<T> &m_, SolidId s) : m{m_}
            {
                shells.push_back(m.solid(s).outer);
                shells.insert(shells.end(), m.solid(s).voids.begin(), m.solid(s).voids.end());
            }

            auto faces() const -> std::vector<FaceUse>
            {
                std::vector<FaceUse> r;
                for (auto sh : shells)
                    r.insert(r.end(), m.shell(sh).faces.begin(), m.shell(sh).faces.end());
                return r;
            }
            auto orient(FaceId f) const -> Orientation
            {
                for (const auto &fu : faces())
                    if (fu.face == f)
                        return fu.orient;
                return Orientation::Forward;
            }
            auto shell_of(FaceId f) const -> ShellId
            {
                for (auto sh : shells)
                    for (const auto &fu : m.shell(sh).faces)
                        if (fu.face == f)
                            return sh;
                return {};
            }
            /// The faces using the edge e (a seam: twice).
            auto users(EdgeId e) const -> std::vector<FaceId>
            {
                std::vector<FaceId> r;
                for (const auto &fu : faces())
                    for (auto w : m.face(fu.face).wires)
                        for (const auto &ce : m.wire(w).coedges)
                            if (ce.edge == e)
                                r.push_back(fu.face);
                return r;
            }
            /// The edges at the vertex v, but degenerate ones, each once.
            auto edges_at(VertexId v) const -> std::vector<EdgeId>
            {
                std::set<EdgeId> r;
                for (const auto &fu : faces())
                    for (auto w : m.face(fu.face).wires)
                        for (const auto &ce : m.wire(w).coedges)
                        {
                            const auto &e = m.edge(ce.edge);
                            if (!e.degenerate && (e.v1 == v || e.v2 == v))
                                r.insert(ce.edge);
                        }
                return {r.begin(), r.end()};
            }
            /// The pcurve of the edge e on the face f (its first use).
            auto pcurve(FaceId f, EdgeId e) const -> std::shared_ptr<Curve<T, 2>>
            {
                for (auto w : m.face(f).wires)
                    for (const auto &ce : m.wire(w).coedges)
                        if (ce.edge == e)
                            return ce.pcurve;
                return nullptr;
            }
            /// The outward normal of the face f at uv of its surface s.
            auto outward(FaceId f, const Surface<T, 3> &s, const point<T, 2> &uv) const -> point<T, 3>
            {
                const auto n = normal_direction(s, uv[0], uv[1]);
                return orient(f) == Orientation::Forward ? n : T(-1) * n;
            }
        };

        /// The point of the blend across it at t, s from the contact 0 (s = 0) to the contact 1, and its (u, v).
        template <std::floating_point T>
        auto across(const EdgeBlend<T> &g, T t, T s) -> point<T, 2>
        {
            const auto p = g.contacts[0].on_blend->value(t), q = g.contacts[1].on_blend->value(t);
            return p + s * (q - p);
        }

        template <std::floating_point T>
        auto across_dt(const EdgeBlend<T> &g, T t, T s) -> point<T, 2>
        {
            const auto p = g.contacts[0].on_blend->value(t, 1), q = g.contacts[1].on_blend->value(t, 1);
            return p + s * (q - p);
        }

        /// The periods of a surface (0: not periodic).
        template <std::floating_point T>
        auto periods(const Surface<T, 3> &s) -> std::array<T, 2>
        {
            if (const auto *e = dynamic_cast<const ElementarySurface<T> *>(&s))
                return e->periods();
            return {T(0), T(0)};
        }

        /// The vertex tolerance raised to cover the point p within tol.
        template <std::floating_point T>
        void cover(Model<T> &m, VertexId v, const point<T, 3> &p, T tol)
        {
            auto &vx = m.vertex(v);
            vx.tol = std::max({vx.tol, tol, norm(vx.pnt - p) * (T(1) + T(1e-9))});
        }

        /// The tolerance of the edge raised to its gaps to the surfaces of its co-edges, its vertices to cover it.
        template <std::floating_point T>
        void settle(Model<T> &m, EdgeId eid, std::type_identity_t<std::span<const std::pair<const Surface<T, 3> *, std::shared_ptr<Curve<T, 2>>>>> pcurves)
        {
            auto &e = m.edge(eid);
            if (e.degenerate)
            {
                for (const auto &[s, pc] : pcurves)
                    for (int i = 0; i <= 8; ++i)
                    {
                        const auto uv = pc->value(e.u1 + (e.u2 - e.u1) * T(i) / T(8));
                        cover(m, e.v1, s->value(uv[0], uv[1]), e.tol);
                    }
                return;
            }
            for (const auto &[s, pc] : pcurves)
                for (int i = 0; i <= 16; ++i)
                {
                    const T u = e.u1 + (e.u2 - e.u1) * T(i) / T(16);
                    const auto uv = pc->value(u);
                    e.tol = std::max(e.tol, norm(s->value(uv[0], uv[1]) - e.curve->value(u)) * (T(1) + T(1e-6)));
                }
            for (const auto &[v, u] : {std::pair{e.v1, e.u1}, std::pair{e.v2, e.u2}})
            {
                cover(m, v, e.curve->value(u), e.tol);
                for (const auto &[s, pc] : pcurves)
                {
                    const auto uv = pc->value(u);
                    cover(m, v, s->value(uv[0], uv[1]), e.tol);
                }
            }
        }

        /// The co-edge's start and end vertices, its sense applied.
        template <std::floating_point T>
        auto ends_of(const Model<T> &m, const CoEdge<T> &ce) -> std::array<VertexId, 2>
        {
            const auto &e = m.edge(ce.edge);
            return ce.orient == Orientation::Forward ? std::array{e.v1, e.v2} : std::array{e.v2, e.v1};
        }

        /// The co-edges run the other way.
        template <std::floating_point T>
        auto reversed(std::vector<CoEdge<T>> c) -> std::vector<CoEdge<T>>
        {
            std::ranges::reverse(c);
            for (auto &ce : c)
                ce.orient = reverse(ce.orient);
            return c;
        }

        /**
         * In a wire of the face f, the run of co-edges from the vertex a to the vertex b (a == b:
         * a loop through a) holding a doomed edge, the rest of the wire holding none, replaced by
         * `chain` (co-edges from a to b; tried the other way round too, and a closed chain both
         * ways), its pcurves moved by periods to follow on. The edges of the run are added to
         * `removed`. False when there is no such run.
         */
        template <std::floating_point T>
        bool replace_run(Model<T> &m, FaceId f, const std::vector<CoEdge<T>> &chain, VertexId a, VertexId b, const std::set<EdgeId> &doomed,
                         std::set<EdgeId> &removed)
        {
            const auto &srf = *m.face(f).surface;
            const auto per = periods(srf);
            const auto box = classify_detail::face_uv_box(m, f);
            const std::array<T, 2> size{std::max(box[1] - box[0], T(1e-300)), std::max(box[3] - box[2], T(1e-300))};
            for (auto wid : m.face(f).wires)
            {
                const auto cs = m.wire(wid).coedges;
                const std::size_t n = cs.size();
                auto has_doomed = [&](std::size_t from, std::size_t count) {
                    for (std::size_t i = 0; i < count; ++i)
                        if (doomed.contains(cs[(from + i) % n].edge))
                            return true;
                    return false;
                };
                for (const bool flip : {false, true})
                {
                    const auto run_chain = flip ? reversed(chain) : chain;
                    const VertexId from = flip ? b : a, to = flip ? a : b;
                    for (std::size_t ia = 0; ia < n; ++ia)
                    {
                        if (ends_of(m, cs[ia])[0] != from)
                            continue;
                        for (std::size_t ib = 0; ib < n; ++ib)
                        {
                            if (ends_of(m, cs[ib])[0] != to || (ia == ib && from != to))
                                continue;
                            const std::size_t count = ia == ib ? n : (ib + n - ia) % n;
                            if (!has_doomed(ia, count) || has_doomed(ib, n - count))
                                continue;
                            // the chain from the (u, v) where the run starts, moved by periods
                            const auto uv_in = coedge_uv(m, cs[ia], m.edge(cs[ia].edge).u1);
                            const auto uv_out = coedge_uv(m, cs[ib], m.edge(cs[ib].edge).u1);
                            auto cur = uv_in;
                            std::vector<CoEdge<T>> placed;
                            bool ok = true;
                            for (auto ce : run_chain)
                            {
                                const auto &e = m.edge(ce.edge);
                                const auto st = coedge_uv(m, ce, e.u1);
                                point<T, 2> shift{};
                                for (std::size_t k = 0; k < 2; ++k)
                                    if (per[k] > T(0))
                                        shift[k] = per[k] * std::round((cur[k] - st[k]) / per[k]);
                                ce.pcurve = shifted_pcurve(ce.pcurve, shift);
                                if (!ce.pcurve)
                                {
                                    ok = false;
                                    break;
                                }
                                cur = coedge_uv(m, ce, e.u2);
                                placed.push_back(ce);
                            }
                            if (!ok)
                                continue;
                            // the chain ends where the run did (not a period away)
                            const auto in0 = coedge_uv(m, placed.front(), m.edge(placed.front().edge).u1);
                            for (const auto &[p, q] : {std::pair{in0, uv_in}, std::pair{cur, uv_out}})
                                for (std::size_t k = 0; k < 2; ++k)
                                    if (std::abs(p[k] - q[k]) > T(1e-4) * size[k])
                                        ok = false;
                            if (!ok)
                                continue;
                            std::vector<CoEdge<T>> nc;
                            for (std::size_t i = 0; i < n - (ia == ib ? n : count); ++i)
                                nc.push_back(cs[(ib + i) % n]);
                            nc.insert(nc.end(), placed.begin(), placed.end());
                            for (std::size_t i = 0; i < count; ++i)
                                removed.insert(cs[(ia + i) % n].edge);
                            m.wire(wid).coedges = std::move(nc);
                            m.face(f).natural_bounds = false;
                            return true;
                        }
                    }
                }
            }
            return false;
        }

        /// A contact of a blend on a face: its edge, its pcurve on the face, and the blended edge.
        template <std::floating_point T>
        struct Cut
        {
            EdgeId edge;
            std::shared_ptr<Curve<T, 2>> pcurve;
            EdgeId blended;
        };

        /// The point p on the boundary of the face f: a vertex of it, or a point of one of its edges (split there, at v).
        template <std::floating_point T>
        void put_on_boundary(Model<T> &m, FaceId f, VertexId v, T tol, const std::string &op)
        {
            for (auto w : m.face(f).wires)
                for (const auto &ce : m.wire(w).coedges)
                    if (ends_of(m, ce)[0] == v)
                        return;
            const auto p = m.vertex(v).pnt;
            std::optional<EdgeId> best;
            T gap = std::numeric_limits<T>::infinity(), at{};
            for (auto w : m.face(f).wires)
                for (const auto &ce : m.wire(w).coedges)
                {
                    const auto &e = m.edge(ce.edge);
                    if (e.degenerate || !e.curve)
                        continue;
                    const std::array<T, 2> range{e.u1, e.u2};
                    T u0 = e.u1, d0 = std::numeric_limits<T>::infinity();
                    for (int i = 0; i <= 64; ++i)
                    {
                        const T u = e.u1 + (e.u2 - e.u1) * T(i) / T(64);
                        if (const T d = norm(e.curve->value(u) - p); d < d0)
                            d0 = d, u0 = u;
                    }
                    const auto [u, d] = newton_point_curve(*e.curve, p, u0, range, T(1e-14) * (T(1) + std::abs(e.u2 - e.u1)));
                    if (std::min(d, d0) < gap)
                        gap = std::min(d, d0), at = d < d0 ? u : u0, best = ce.edge;
                }
            const T reach = T(100) * tol;
            if (!best || gap > reach)
                refuse<T>(op, "a contact ends inside a face (a chain of edges whose contacts leave their faces apart from an edge)", {f});
            const auto &e = m.edge(*best);
            const T eps = T(1e-9) * (e.u2 - e.u1);
            if (!(at > e.u1 + eps && at < e.u2 - eps))
                refuse<T>(op, "a contact ends at a vertex of a face", {f, *best});
            cover(m, v, e.curve->value(at), e.tol);
            if (!split_detail::split_at_vertex(m, *best, at, v))
                refuse<T>(op, "an edge that cannot be split where a contact ends", {*best});
        }

        /// The contacts on the face f chained, from end to end (or round a loop), and the blended edges along them.
        template <std::floating_point T>
        void cut_face(Model<T> &m, FaceId f, const std::vector<Cut<T>> &cuts, T tol, std::set<EdgeId> &removed, const std::string &op)
        {
            std::map<VertexId, std::vector<std::size_t>> at;
            for (std::size_t i = 0; i < cuts.size(); ++i)
                for (auto v : {m.edge(cuts[i].edge).v1, m.edge(cuts[i].edge).v2})
                    at[v].push_back(i);
            for (const auto &[v, is] : at)
                if (is.size() > 2)
                    refuse<T>(op, "more than two contacts meet at a point of a face", {f});
            std::vector<bool> used(cuts.size(), false);
            auto chain_from = [&](std::size_t i0, VertexId start) {
                std::vector<CoEdge<T>> chain;
                std::set<EdgeId> blended;
                VertexId v = start;
                std::size_t i = i0;
                for (;;)
                {
                    used[i] = true;
                    const auto &e = m.edge(cuts[i].edge);
                    const bool fwd = e.v1 == v;
                    chain.push_back(CoEdge<T>{cuts[i].edge, fwd ? Orientation::Forward : Orientation::Reversed, cuts[i].pcurve});
                    blended.insert(cuts[i].blended);
                    v = fwd ? e.v2 : e.v1;
                    const auto &next = at.at(v);
                    const auto it = std::ranges::find_if(next, [&](std::size_t j) { return !used[j]; });
                    if (it == next.end())
                        break;
                    i = *it;
                }
                return std::tuple{chain, blended, v};
            };
            std::map<EdgeId, EdgeId> blended_of;
            for (const auto &c : cuts)
                blended_of[c.edge] = c.blended;
            // a vertex on the boundary of the face: one of its vertices, or on one of its edges but the blended ones
            auto on_boundary = [&](VertexId v, const std::set<EdgeId> &blended) {
                const auto p = m.vertex(v).pnt;
                for (auto w : m.face(f).wires)
                    for (const auto &ce : m.wire(w).coedges)
                    {
                        if (ends_of(m, ce)[0] == v)
                            return true;
                        const auto &e = m.edge(ce.edge);
                        if (e.degenerate || !e.curve || blended.contains(ce.edge))
                            continue;
                        T d0 = std::numeric_limits<T>::infinity(), u0 = e.u1;
                        for (int j = 0; j <= 64; ++j)
                        {
                            const T u = e.u1 + (e.u2 - e.u1) * T(j) / T(64);
                            if (const T d = norm(e.curve->value(u) - p); d < d0)
                                d0 = d, u0 = u;
                        }
                        const T d = newton_point_curve(*e.curve, p, u0, std::array{e.u1, e.u2}, T(1e-14) * (T(1) + std::abs(e.u2 - e.u1)))[1];
                        if (std::min(d, d0) <= T(100) * tol)
                            return true;
                    }
                return false;
            };
            // a chain cut where it touches the boundary (a seam it crosses, an edge between two
            // faces of a chain), each run replacing the band along its blended edges
            auto process = [&](std::vector<CoEdge<T>> chain, bool closed) {
                std::set<EdgeId> blended;
                for (const auto &ce : chain)
                    blended.insert(blended_of.at(ce.edge));
                const std::size_t n = chain.size();
                std::vector<std::size_t> on;
                for (std::size_t k = closed ? 0 : 1; k < n; ++k)
                    if (on_boundary(ends_of(m, chain[k])[0], blended))
                        on.push_back(k);
                if (closed && on.empty())
                {
                    // a wire made of blended edges only, replaced as a whole
                    bool done = false;
                    for (auto w : m.face(f).wires)
                    {
                        auto &cs = m.wire(w).coedges;
                        if (!std::ranges::all_of(cs, [&](const CoEdge<T> &ce) { return blended.contains(ce.edge); }))
                            continue;
                        // the loop turning the same way as the wire in (u, v), on its periods
                        const T area = uv_signed_area(m, std::span<const CoEdge<T>>{cs});
                        auto loop = uv_signed_area(m, std::span<const CoEdge<T>>{chain}) * area < T(0) ? reversed(chain) : chain;
                        for (auto ce : cs)
                            removed.insert(ce.edge);
                        const auto per = periods(*m.face(f).surface);
                        const auto ref = coedge_uv(m, cs.front(), m.edge(cs.front().edge).u1);
                        auto cur = ref;
                        for (std::size_t k = 0; k < loop.size(); ++k)
                        {
                            auto &ce = loop[k];
                            const auto st = coedge_uv(m, ce, m.edge(ce.edge).u1);
                            point<T, 2> shift{};
                            for (std::size_t d = 0; d < 2; ++d)
                                if (per[d] > T(0))
                                    shift[d] = per[d] * std::round((cur[d] - st[d]) / per[d]);
                            if (k == 0)
                                for (std::size_t d = 0; d < 2; ++d)
                                    if (per[d] > T(0)) // the loop in the window of the face, near where the wire was
                                        shift[d] = per[d] * std::round((ref[d] - st[d]) / per[d]);
                            ce.pcurve = shifted_pcurve(ce.pcurve, shift);
                            if (!ce.pcurve)
                                refuse<T>(op, "a contact that cannot be moved by a period", {f});
                            cur = coedge_uv(m, ce, m.edge(ce.edge).u2);
                        }
                        cs = std::move(loop);
                        m.face(f).natural_bounds = false;
                        done = true;
                        break;
                    }
                    if (!done)
                        refuse<T>(op, "a loop of contacts inside a face that bounds no wire of blended edges", {f});
                    return;
                }
                if (closed)
                {
                    const auto first = static_cast<std::ptrdiff_t>(on.front());
                    std::ranges::rotate(chain, chain.begin() + first);
                    for (auto &k : on)
                        k -= static_cast<std::size_t>(first);
                }
                else
                    on.insert(on.begin(), 0);
                on.push_back(n);
                for (std::size_t i = 0; i + 1 < on.size(); ++i)
                {
                    const std::vector<CoEdge<T>> sub(chain.begin() + static_cast<std::ptrdiff_t>(on[i]), chain.begin() + static_cast<std::ptrdiff_t>(on[i + 1]));
                    std::set<EdgeId> along;
                    for (const auto &ce : sub)
                        along.insert(blended_of.at(ce.edge));
                    const auto va = ends_of(m, sub.front())[0], vb = ends_of(m, sub.back())[1];
                    put_on_boundary(m, f, va, tol, op);
                    put_on_boundary(m, f, vb, tol, op);
                    if (!replace_run(m, f, sub, va, vb, along, removed))
                        refuse<T>(op, "the contacts do not bound a band of a face along its blended edges", {f});
                }
            };
            for (const auto &[v, is] : at)
                if (is.size() == 1 && !used[is.front()])
                    process(std::get<0>(chain_from(is.front(), v)), false);
            for (std::size_t i = 0; i < cuts.size(); ++i)
                if (!used[i])
                    process(std::get<0>(chain_from(i, m.edge(cuts[i].edge).v1)), true);
        }

        /// The parameter t of the curve c where it meets the surface s (Gauss-Newton from t0, uv0), with its (u, v).
        template <std::floating_point T>
        auto meet(const Curve<T, 3> &c, const Surface<T, 3> &s, T t0, const point<T, 2> &uv0, T tol) -> std::optional<std::pair<T, point<T, 2>>>
        {
            Vec<T> x(3);
            x << t0, uv0[0], uv0[1];
            const T res = thick_detail::gauss_newton<T>(x, [&](const Vec<T> &y, Vec<T> &f, Mat<T> &J) {
                const auto p = c.value(y[0]), dp = c.value(y[0], 1);
                const auto q = s.value(y[1], y[2]), su = s.value(y[1], y[2], 1, 0), sv = s.value(y[1], y[2], 0, 1);
                f.resize(3);
                J.resize(3, 3);
                for (int i = 0; i < 3; ++i)
                {
                    f[i] = p[std::size_t(i)] - q[std::size_t(i)];
                    J.row(i) << dp[std::size_t(i)], -su[std::size_t(i)], -sv[std::size_t(i)];
                }
            });
            if (!(res <= tol / T(10)))
                return std::nullopt;
            return std::pair{x[0], point<T, 2>{x[1], x[2]}};
        }

        /// The point of the blend across it at s that lies on the surface s3: its t on the blend and its (u, v) on s3.
        template <std::floating_point T>
        auto meet_across(const EdgeBlend<T> &g, const Surface<T, 3> &s3, T s, T t0, const point<T, 2> &uv0, T tol)
            -> std::optional<std::pair<T, point<T, 2>>>
        {
            Vec<T> x(3);
            x << t0, uv0[0], uv0[1];
            const auto &b = *g.surface;
            const T res = thick_detail::gauss_newton<T>(x, [&](const Vec<T> &y, Vec<T> &f, Mat<T> &J) {
                const auto w = across(g, y[0], s), dw = across_dt(g, y[0], s);
                const auto p = b.value(w[0], w[1]);
                const auto dp = dw[0] * b.value(w[0], w[1], 1, 0) + dw[1] * b.value(w[0], w[1], 0, 1);
                const auto q = s3.value(y[1], y[2]), su = s3.value(y[1], y[2], 1, 0), sv = s3.value(y[1], y[2], 0, 1);
                f.resize(3);
                J.resize(3, 3);
                for (int i = 0; i < 3; ++i)
                {
                    f[i] = p[std::size_t(i)] - q[std::size_t(i)];
                    J.row(i) << dp[std::size_t(i)], -su[std::size_t(i)], -sv[std::size_t(i)];
                }
            });
            if (!(res <= tol / T(10)) || !(x[0] >= g.a - T(1e-9) * (g.b - g.a)) || !(x[0] <= g.b + T(1e-9) * (g.b - g.a)))
                return std::nullopt;
            return std::pair{x[0], point<T, 2>{x[1], x[2]}};
        }

        /// The face at a free end of a blend: the edges between it and the faces of the edge, its surface continued.
        template <std::floating_point T>
        struct ThirdFace
        {
            std::array<EdgeId, 2> side;
            FaceId face;
            std::shared_ptr<Surface<T, 3>> surface;
            point<T, 2> uv; ///< of the vertex
        };

        /// The third face at the vertex v, the end of the blend g at a vertex of three edges.
        template <std::floating_point T>
        auto third_face(const Model<T> &m, const SolidView<T> &sv, const EdgeBlend<T> &g, VertexId v, T tol, const std::string &op) -> ThirdFace<T>
        {
            const auto at = sv.edges_at(v);
            if (at.size() != 3)
                refuse<T>(op, "a blend ending at a vertex of more than three edges", {g.edge, v});
            ThirdFace<T> tf;
            std::array<FaceId, 2> third;
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto f = g.contacts[k].face;
                std::optional<EdgeId> e;
                for (auto x : at)
                    if (x != g.edge && std::ranges::contains(sv.users(x), f))
                        e = x;
                if (!e)
                    refuse<T>(op, "a blend ending at a vertex where a face of the edge has no other edge", {g.edge, v});
                tf.side[k] = *e;
                const auto us = sv.users(*e);
                if (us.size() != 2)
                    refuse<T>(op, "an edge at the end of a blend not shared by two faces", {*e});
                third[k] = us[0] == f ? us[1] : us[0];
            }
            if (tf.side[0] == tf.side[1] || third[0] != third[1])
                refuse<T>(op, "a blend ending at a vertex without a third face", {g.edge, v});
            // an edge going on tangentially: a chain blended in part
            auto tangent_at = [&](EdgeId eid) {
                const auto &e = m.edge(eid);
                return fillet_detail::unit<T>(e.curve->value(e.v1 == v ? e.u1 : e.u2, 1));
            };
            for (std::size_t k = 0; k < 2; ++k)
                if (std::abs(tangent_at(g.edge) * tangent_at(tf.side[k])) > T(1) - T(1e-6))
                    refuse<T>(op, "a chain of tangent edges blended in part (blend the whole chain)", {g.edge, tf.side[k]});
            tf.face = third[0];
            // the surface of the third face continued (the corner of a concave edge is outside it)
            const T reach = T(4) * (norm(g.contacts[0].curve->value(g.a) - g.contacts[1].curve->value(g.a)) +
                                    norm(g.contacts[0].curve->value(g.b) - g.contacts[1].curve->value(g.b))) +
                            T(10) * tol;
            tf.surface = thick_detail::extended_surface(m.face(tf.face).surface, classify_detail::face_uv_box(m, tf.face), reach, tol).surface;
            const auto pc = sv.pcurve(tf.face, tf.side[0]);
            const auto &es = m.edge(tf.side[0]);
            tf.uv = pc->value(es.v1 == v ? es.u1 : es.u2);
            return tf;
        }

        /// Where the contact k of the blend meets the third face, within the blend.
        template <std::floating_point T>
        auto meet_contact(const EdgeBlend<T> &g, std::size_t k, const ThirdFace<T> &tf, T t, T tol) -> std::optional<std::pair<T, point<T, 2>>>
        {
            const auto x = meet(*g.contacts[k].curve, *tf.surface, t, tf.uv, tol);
            if (!x || !(x->first >= g.a - T(1e-9) * (g.b - g.a) && x->first <= g.b + T(1e-9) * (g.b - g.a)))
                return std::nullopt;
            return x;
        }

        /**
         * The free end of the blend at its end j (at the vertex v): the third face, the contacts
         * met on the edges between it and the faces of the edge (split there), the end curve from
         * one to the other across the blend.
         */
        template <std::floating_point T>
        void free_end(Model<T> &m, const SolidView<T> &sv, Blend<T> &bl, std::size_t j, T tol, const std::string &op)
        {
            auto &end = bl.ends[j];
            const auto &g = bl.g;
            const VertexId v = end.vertex;
            const auto tf = third_face(m, sv, g, v, tol, op);
            const auto side = tf.side;
            const FaceId f3 = tf.face;
            end.f3 = f3;
            const auto &s3 = tf.surface;
            // the contacts on the third face, on the edges next to it
            std::array<point<T, 2>, 2> uv3;
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto x = meet_contact(g, k, tf, end.t, tol);
                if (!x)
                    refuse<T>(op, "a contact does not meet the face at an end of its edge within the blend continued", {g.edge, f3});
                end.tc[k] = std::clamp(x->first, g.a, g.b);
                uv3[k] = x->second;
                const auto p = g.contacts[k].curve->value(end.tc[k]);
                const auto &e = m.edge(side[k]);
                T u0 = e.u1, d0 = std::numeric_limits<T>::infinity();
                for (int i = 0; i <= 64; ++i)
                {
                    const T u = e.u1 + (e.u2 - e.u1) * T(i) / T(64);
                    if (const T d = norm(e.curve->value(u) - p); d < d0)
                        d0 = d, u0 = u;
                }
                const auto [u, d] = newton_point_curve(*e.curve, p, u0, std::array{e.u1, e.u2}, T(1e-14) * (T(1) + std::abs(e.u2 - e.u1)));
                const T eps = T(1e-9) * (e.u2 - e.u1);
                if (!(d <= T(100) * tol) || !(u > e.u1 + eps && u < e.u2 - eps))
                    refuse<T>(op, "a blend wider than a face at an end of its edge", {g.edge, side[k]});
                end.v[k] = m.add(Vertex<T>{p, tol});
                cover(m, end.v[k], e.curve->value(u), std::max(tol, e.tol));
                if (!split_detail::split_at_vertex(m, side[k], u, end.v[k]))
                    refuse<T>(op, "an edge that cannot be split where a contact ends", {side[k]});
            }
            // the end curve across the blend, from the contact 0 to the contact 1
            struct Pt
            {
                T t;
                point<T, 2> uv3;
                point<T, 3> p;
            };
            auto solve = [&](T s, const Pt &from) -> Pt {
                const auto x = meet_across(g, *s3, s, from.t, from.uv3, tol);
                if (!x)
                    refuse<T>(op, "the blend does not meet the face at an end of its edge", {g.edge, f3});
                const auto w = across(g, x->first, s);
                return {x->first, x->second, g.surface->value(w[0], w[1])};
            };
            Pt first{end.tc[0], uv3[0], {}};
            first = solve(T(0), first);
            std::shared_ptr<BSCurve<T, 3>> c3;
            std::shared_ptr<BSCurve<T, 2>> cb, cf;
            T dev = std::numeric_limits<T>::infinity();
            for (std::size_t n = 8; n <= 512; n *= 2)
            {
                std::vector<Pt> pts{first};
                for (std::size_t i = 1; i <= n; ++i)
                    pts.push_back(solve(T(i) / T(n), pts.back()));
                std::vector<T> s(n + 1);
                points_vector<T, 3> q3;
                points_vector<T, 2> qb, qf;
                for (std::size_t i = 0; i <= n; ++i)
                {
                    s[i] = T(i) / T(n);
                    q3.push_back(pts[i].p);
                    qb.push_back(across(g, pts[i].t, s[i]));
                    qf.push_back(pts[i].uv3);
                }
                auto k3 = std::make_shared<BSCurve<T, 3>>(interpolate(q3, s, 3));
                auto kb = std::make_shared<BSCurve<T, 2>>(interpolate(qb, s, 3));
                auto kf = std::make_shared<BSCurve<T, 2>>(interpolate(qf, s, 3));
                T d{};
                for (std::size_t i = 0; i < n; ++i)
                {
                    const T sm = (s[i] + s[i + 1]) / T(2);
                    const auto x = solve(sm, pts[i]);
                    const auto wb = kb->value(sm), wf = kf->value(sm);
                    d = std::max({d, norm(k3->value(sm) - x.p), norm(g.surface->value(wb[0], wb[1]) - x.p), norm(s3->value(wf[0], wf[1]) - x.p)});
                }
                if (d < dev)
                    dev = d, c3 = k3, cb = kb, cf = kf;
                if (dev <= tol / T(10))
                    break;
            }
            if (!(dev <= tol))
                refuse<T>(op, "the end curve of a blend cannot be approximated within the tolerance", {g.edge, f3});
            Edge<T> ne;
            ne.curve = c3;
            ne.u1 = T(0);
            ne.u2 = T(1);
            ne.v1 = end.v[0];
            ne.v2 = end.v[1];
            ne.tol = std::max(tol, dev);
            end.edge = m.add(std::move(ne));
            end.on_blend = cb;
            end.on_f3 = cf;
            end.forward = true;
            const std::array<std::pair<const Surface<T, 3> *, std::shared_ptr<Curve<T, 2>>>, 2> pcs{std::pair{g.surface.get(), cb},
                                                                                                      std::pair{m.face(f3).surface.get(), cf}};
            settle(m, end.edge, std::span{pcs});
        }

        /// Two blends joined at the vertex v: the arc across both there, their contacts meeting at its ends.
        template <std::floating_point T>
        void joint(Model<T> &m, Blend<T> &b1, std::size_t j1, Blend<T> &b2, std::size_t j2, T tol, T approx_tol, const std::string &op)
        {
            auto &e1 = b1.ends[j1], &e2 = b2.ends[j2];
            const auto &g1 = b1.g, &g2 = b2.g;
            const T limit = T(10) * tol + (g1.exact ? T(0) : approx_tol) + (g2.exact ? T(0) : approx_tol);
            std::array<std::size_t, 2> k2{};
            std::array<T, 2> gap{};
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto p = g1.contacts[k].curve->value(e1.t);
                const T d0 = norm(p - g2.contacts[0].curve->value(e2.t)), d1 = norm(p - g2.contacts[1].curve->value(e2.t));
                k2[k] = d0 <= d1 ? 0 : 1;
                gap[k] = std::min(d0, d1);
            }
            const auto pc1 = fillet_detail::segment(across(g1, e1.t, T(0)), across(g1, e1.t, T(1)), T(0), T(1));
            const auto pc2 = fillet_detail::segment(across(g2, e2.t, T(k2[0])), across(g2, e2.t, T(k2[1])), T(0), T(1));
            T arc_gap{};
            for (int i = 1; i < 8; ++i)
            {
                const T s = T(i) / T(8);
                const auto w1 = pc1->value(s), w2 = pc2->value(s);
                arc_gap = std::max(arc_gap, norm(g1.surface->value(w1[0], w1[1]) - g2.surface->value(w2[0], w2[1])));
            }
            if (k2[0] == k2[1] || gap[0] > limit || gap[1] > limit || arc_gap > limit)
                refuse<T>(op, "blended edges meeting at a vertex but tangentially, or blends that do not join there (vertex blends: PR 7)",
                          {g1.edge, g2.edge, e1.vertex});
            for (std::size_t k = 0; k < 2; ++k)
            {
                e1.v[k] = m.add(Vertex<T>{g1.contacts[k].curve->value(e1.t), std::max(tol, gap[k])});
                cover(m, e1.v[k], g2.contacts[k2[k]].curve->value(e2.t), tol);
                e2.v[k2[k]] = e1.v[k];
                e1.tc[k] = e1.t;
                e2.tc[k] = e2.t;
            }
            Edge<T> ne;
            ne.curve = std::make_shared<CurveOnSurface<T, 3>>(pc1, g1.surface);
            ne.u1 = T(0);
            ne.u2 = T(1);
            ne.v1 = e1.v[0];
            ne.v2 = e1.v[1];
            ne.tol = tol;
            const auto eid = m.add(std::move(ne));
            e1.kind = e2.kind = EndKind::Joint;
            e1.edge = e2.edge = eid;
            e1.on_blend = pc1;
            e2.on_blend = pc2;
            e1.forward = true;
            e2.forward = k2[0] == 0;
            const std::array<std::pair<const Surface<T, 3> *, std::shared_ptr<Curve<T, 2>>>, 2> pcs{std::pair{g1.surface.get(), pc1},
                                                                                                      std::pair{g2.surface.get(), pc2}};
            settle(m, eid, std::span{pcs});
        }

        /// The contact k of the blend as an edge between the ends of the blend.
        template <std::floating_point T>
        auto contact_edge(Model<T> &m, const Blend<T> &bl, std::size_t k, T tol) -> EdgeId
        {
            const auto &ct = bl.g.contacts[k];
            Edge<T> e;
            e.curve = ct.curve;
            e.u1 = bl.ends[0].tc[k];
            e.u2 = bl.ends[1].tc[k];
            e.v1 = bl.ends[0].v[k];
            e.v2 = bl.ends[1].v[k];
            e.tol = std::max(tol, ct.deviation * (T(1) + T(1e-6)));
            return m.add(std::move(e));
        }

        /// The face of the blend, bounded by its contacts and its ends, its use oriented as the faces of the edge.
        template <std::floating_point T>
        auto blend_face(Model<T> &m, const SolidView<T> &sv, const Blend<T> &bl, T tol) -> FaceUse
        {
            const auto &g = bl.g;
            std::vector<CoEdge<T>> cs;
            cs.push_back(CoEdge<T>{bl.contacts[0], Orientation::Forward, g.contacts[0].on_blend});
            const auto &hi = bl.ends[1], &lo = bl.ends[0];
            if (hi.edge.valid())
                cs.push_back(CoEdge<T>{hi.edge, hi.forward ? Orientation::Forward : Orientation::Reversed, hi.on_blend});
            cs.push_back(CoEdge<T>{bl.contacts[1], Orientation::Reversed, g.contacts[1].on_blend});
            if (lo.edge.valid())
                cs.push_back(CoEdge<T>{lo.edge, lo.forward ? Orientation::Reversed : Orientation::Forward, lo.on_blend});
            if (uv_signed_area(m, std::span<const CoEdge<T>>{cs}) < T(0))
                cs = reversed(cs);
            Wire<T> w;
            w.coedges = std::move(cs);
            w.closed = true;
            const auto wid = m.add(std::move(w));
            const auto fid = m.add(Face<T>{g.surface, {wid}, tol, false});
            // oriented as the faces of the edge: its normal along the sum of theirs, at the middle
            const T tm = (bl.ends[0].tc[0] + bl.ends[1].tc[0]) / T(2);
            point<T, 3> out{};
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto &ct = g.contacts[k];
                const T t = (bl.ends[0].tc[k] + bl.ends[1].tc[k]) / T(2);
                out = out + sv.outward(ct.face, *ct.support, ct.pcurve->value(t));
            }
            const auto w0 = across(g, tm, T(0.5));
            const auto n = normal_direction(*g.surface, w0[0], w0[1]);
            return FaceUse{fid, n * out >= T(0) ? Orientation::Forward : Orientation::Reversed};
        }

        /**
         * The blends of the edges of the solid s (their geometry from `geometry`, on a copy):
         * the new solid (file header).
         */
        template <std::floating_point T, typename Geometry>
        auto blend(Model<T> &m, SolidId s, std::span<const EdgeId> edges, Geometry &&geometry, const FilletOptions<T> &o, const std::string &op)
            -> SolidId
        {
            if (!m.alive(s))
                fillet_detail::refuse<T>(BuildErrc::InvalidId, op, "dead or invalid solid", {s});
            if (edges.empty())
                fillet_detail::refuse<T>(BuildErrc::EmptyInput, op, "no edge given", {s});
            for (std::size_t i = 0; i < edges.size(); ++i)
                if (std::ranges::contains(edges.first(i), edges[i]))
                    fillet_detail::refuse<T>(BuildErrc::InvalidId, op, "an edge given twice", {edges[i]});
            const T tol = o.tol;
            // ---- the copy, and the geometry of each blend on it
            detail::ShapeCopier<T> copier{m, Placement<T>{}};
            const auto ns = std::get<SolidId>(copier.copy(ShapeId{s}));
            std::vector<Blend<T>> bls;
            for (auto e : edges)
            {
                if (!m.alive(e))
                    fillet_detail::refuse<T>(BuildErrc::InvalidId, op, "dead or invalid edge", {e});
                const auto ce = copier.copy(e);
                Blend<T> bl;
                bl.g = geometry(ns, ce, copier, o);
                bls.push_back(std::move(bl));
            }
            const SolidView<T> sv{m, ns};
            // a free end met past the blend continued (an edge oblique to the face at its end): the
            // blend continued further, by 2, 4 and 8 times as much, while its geometry can be
            for (auto &bl : bls)
            {
                if (bl.g.closed)
                    continue;
                for (int attempt = 0; attempt < 3; ++attempt)
                {
                    const auto &g = bl.g;
                    const auto &ed = m.edge(g.edge);
                    bool short_end = false;
                    for (const auto &[v, t] : {std::pair{ed.v1, g.t1}, std::pair{ed.v2, g.t2}})
                    {
                        const bool joined = std::ranges::any_of(bls, [&](const Blend<T> &b2) {
                            const auto &e2 = m.edge(b2.g.edge);
                            return b2.g.edge != g.edge && (e2.v1 == v || e2.v2 == v);
                        });
                        const bool vanishing = (t == g.a || t == g.b) && norm(g.contacts[0].curve->value(t) - g.contacts[1].curve->value(t)) <= T(10) * tol;
                        if (joined || vanishing)
                            continue;
                        const auto tf = third_face(m, sv, g, v, tol, op);
                        for (std::size_t k = 0; k < 2; ++k)
                            short_end = short_end || !meet_contact(g, k, tf, t, tol);
                    }
                    if (!short_end)
                        break;
                    auto o2 = o;
                    o2.extension = o.extension * std::pow(T(2), T(attempt + 1));
                    try
                    {
                        bl.g = geometry(ns, bl.g.edge, copier, o2);
                    }
                    catch (const BuildFailure &)
                    {
                        break;
                    }
                }
            }
            std::set<EdgeId> blended;
            for (const auto &bl : bls)
                blended.insert(bl.g.edge);

            // ---- the ends: lo at a, hi at b
            for (auto &bl : bls)
            {
                const auto &g = bl.g;
                const auto &ed = m.edge(g.edge);
                if (g.closed)
                {
                    for (std::size_t j = 0; j < 2; ++j)
                    {
                        auto &end = bl.ends[j];
                        end.kind = EndKind::Loop;
                        end.vertex = ed.v1;
                        end.t = j == 0 ? g.a : g.b;
                        end.tc = {end.t, end.t};
                    }
                    for (std::size_t k = 0; k < 2; ++k)
                        bl.ends[0].v[k] = bl.ends[1].v[k] = m.add(Vertex<T>{g.contacts[k].curve->value(g.a), tol});
                    const auto pa = fillet_detail::segment(across(g, g.a, T(0)), across(g, g.a, T(1)), T(0), T(1));
                    const auto pb = fillet_detail::segment(across(g, g.b, T(0)), across(g, g.b, T(1)), T(0), T(1));
                    Edge<T> seam;
                    seam.curve = std::make_shared<CurveOnSurface<T, 3>>(pa, g.surface);
                    seam.u1 = T(0);
                    seam.u2 = T(1);
                    seam.v1 = bl.ends[0].v[0];
                    seam.v2 = bl.ends[0].v[1];
                    seam.tol = tol;
                    const auto sid = m.add(std::move(seam));
                    bl.ends[0].edge = bl.ends[1].edge = sid;
                    bl.ends[0].on_blend = pa;
                    bl.ends[1].on_blend = pb;
                    const std::array<std::pair<const Surface<T, 3> *, std::shared_ptr<Curve<T, 2>>>, 2> pcs{std::pair{g.surface.get(), pa},
                                                                                                              std::pair{g.surface.get(), pb}};
                    settle(m, sid, std::span{pcs});
                    continue;
                }
                const bool along = g.t1 < g.t2;
                bl.ends[0].vertex = along ? ed.v1 : ed.v2;
                bl.ends[1].vertex = along ? ed.v2 : ed.v1;
                bl.ends[0].t = std::min(g.t1, g.t2);
                bl.ends[1].t = std::max(g.t1, g.t2);
            }
            // joints first (both blends at once), then the vanishing and free ends
            for (std::size_t i = 0; i < bls.size(); ++i)
            {
                if (bls[i].g.closed)
                    continue;
                for (std::size_t j = 0; j < 2; ++j)
                {
                    auto &end = bls[i].ends[j];
                    if (end.edge.valid() || end.kind == EndKind::Vanish)
                        continue;
                    std::vector<std::pair<std::size_t, std::size_t>> others;
                    for (std::size_t i2 = 0; i2 < bls.size(); ++i2)
                    {
                        const auto &e2 = m.edge(bls[i2].g.edge);
                        if (i2 == i || !(e2.v1 == end.vertex || e2.v2 == end.vertex))
                            continue;
                        if (bls[i2].g.closed)
                            refuse<T>(op, "a blended closed edge at a vertex of another blended edge (vertex blends: PR 7)", {bls[i2].g.edge, end.vertex});
                        for (std::size_t j2 = 0; j2 < 2; ++j2)
                            if (bls[i2].ends[j2].vertex == end.vertex)
                                others.emplace_back(i2, j2);
                    }
                    if (others.size() > 1)
                        refuse<T>(op, "more than two blended edges at a vertex (vertex blends: PR 7)", {end.vertex});
                    if (others.size() == 1)
                    {
                        const auto [i2, j2] = others.front();
                        joint(m, bls[i], j, bls[i2], j2, tol, o.approx_tol, op);
                    }
                }
            }
            for (auto &bl : bls)
            {
                if (bl.g.closed)
                    continue;
                const auto &g = bl.g;
                for (std::size_t j = 0; j < 2; ++j)
                {
                    auto &end = bl.ends[j];
                    if (end.kind == EndKind::Joint)
                        continue;
                    const T limit = j == 0 ? g.a : g.b;
                    if (end.t == limit && norm(g.contacts[0].curve->value(end.t) - g.contacts[1].curve->value(end.t)) <= T(10) * tol)
                    {
                        // vanishing: the contacts meet at the vertex
                        end.kind = EndKind::Vanish;
                        end.v = {end.vertex, end.vertex};
                        end.tc = {end.t, end.t};
                        for (std::size_t k = 0; k < 2; ++k)
                            cover(m, end.vertex, g.contacts[k].curve->value(end.t), tol);
                        const auto p = across(g, end.t, T(0)), q = across(g, end.t, T(1));
                        if (norm(q - p) > T(1e-12) * (T(1) + norm(p)))
                        {
                            end.on_blend = fillet_detail::segment(p, q, T(0), T(1));
                            end.edge = m.add(Edge<T>{nullptr, T(0), T(1), end.vertex, end.vertex, m.vertex(end.vertex).tol, true, true});
                            const std::array<std::pair<const Surface<T, 3> *, std::shared_ptr<Curve<T, 2>>>, 1> pcs{std::pair{g.surface.get(), end.on_blend}};
                            settle(m, end.edge, std::span{pcs});
                        }
                        continue;
                    }
                    end.kind = EndKind::Free;
                    free_end(m, sv, bl, j, tol, op);
                }
            }

            // ---- the contacts, and the faces of the edges cut along them
            std::map<FaceId, std::vector<Cut<T>>> cuts;
            for (auto &bl : bls)
                for (std::size_t k = 0; k < 2; ++k)
                {
                    if (!(bl.ends[1].tc[k] > bl.ends[0].tc[k]))
                        refuse<T>(op, "a contact reduced to a point or turned back", {bl.g.edge});
                    bl.contacts[k] = contact_edge(m, bl, k, tol);
                    const auto &ct = bl.g.contacts[k];
                    const std::array<std::pair<const Surface<T, 3> *, std::shared_ptr<Curve<T, 2>>>, 2> pcs{
                        std::pair{ct.support.get(), ct.pcurve}, std::pair{bl.g.surface.get(), ct.on_blend}};
                    settle(m, bl.contacts[k], std::span{pcs});
                    cuts[ct.face].push_back(Cut<T>{bl.contacts[k], ct.pcurve, bl.g.edge});
                }
            std::set<EdgeId> removed;
            for (const auto &[f, cs] : cuts)
                cut_face(m, f, cs, tol, removed, op);
            // the third faces at the free ends: the run through the vertex replaced by the end curve
            for (const auto &bl : bls)
                for (const auto &end : bl.ends)
                    if (end.kind == EndKind::Free)
                    {
                        std::set<EdgeId> doomed;
                        for (auto w : m.face(end.f3).wires)
                            for (const auto &ce : m.wire(w).coedges)
                                if (const auto &e = m.edge(ce.edge); !e.degenerate && (e.v1 == end.vertex || e.v2 == end.vertex))
                                    doomed.insert(ce.edge);
                        const std::vector<CoEdge<T>> chain{CoEdge<T>{end.edge, Orientation::Forward, end.on_f3}};
                        if (!replace_run(m, end.f3, chain, end.v[0], end.v[1], doomed, removed))
                            refuse<T>(op, "the end curve of a blend does not cut the corner of the face at the end", {bl.g.edge, end.f3});
                    }

            // ---- the faces of the blends, in the shell of their edge
            for (const auto &bl : bls)
            {
                const auto sh = sv.shell_of(bl.g.contacts[0].face);
                const auto fu = blend_face(m, sv, bl, tol);
                m.shell(sh).faces.push_back(fu);
            }

            // ---- the edges and vertices left out erased
            std::set<EdgeId> used;
            std::set<VertexId> used_v;
            for (const auto &fu : sv.faces())
                for (auto w : m.face(fu.face).wires)
                    for (const auto &ce : m.wire(w).coedges)
                    {
                        used.insert(ce.edge);
                        used_v.insert(m.edge(ce.edge).v1);
                        used_v.insert(m.edge(ce.edge).v2);
                    }
            for (auto e : removed)
            {
                if (used.contains(e) || !m.alive(e))
                    continue;
                for (auto v : {m.edge(e).v1, m.edge(e).v2})
                    if (!used_v.contains(v) && m.alive(v))
                        m.erase(v);
                m.erase(e);
            }
            const auto report = check(m, ShapeId{ns});
            if (!report.ok())
            {
                const auto &first = report.entries.front();
                refuse<T>(op, std::string("the result is not a valid solid: ") + to_string(first.issue) + " " + first.detail, {first.shape});
            }
            return ns;
        }

        template <std::floating_point T, typename F>
        auto guarded(Model<T> &m, const std::string &op, F &&f) -> BuildResult<SolidId>
        {
            const auto size = detail::arena_size(m);
            try
            {
                return f();
            }
            catch (const BuildFailure &e)
            {
                detail::rollback(m, size);
                return std::unexpected(e.error());
            }
            catch (const std::exception &e)
            {
                detail::rollback(m, size);
                return detail::fail(BuildErrc::NonGeneric, op + ": " + e.what());
            }
        }
    }

    /**
     * @brief The solid `s` with its edges `edges` filleted at the radius r (file header): a new
     * solid, `s` unchanged. The edges blended together meet only along chains (tangent at their
     * common vertices); on failure the model is left as it was.
     */
    template <std::floating_point T>
    [[nodiscard]] auto fillet(Model<T> &m, SolidId s, std::span<const EdgeId> edges, std::type_identity_t<T> r, const FilletOptions<T> &o = {})
        -> BuildResult<SolidId>
    {
        const std::string op = "fillet";
        return blend_detail::guarded(m, op, [&] {
            return blend_detail::blend(
                m, s, edges, [&](SolidId ns, EdgeId e, detail::ShapeCopier<T> &, const FilletOptions<T> &oe) { return fillet_detail::fillet(m, ns, e, r, oe); }, o, op);
        });
    }

    template <std::floating_point T>
    [[nodiscard]] auto fillet(Model<T> &m, SolidId s, const std::vector<EdgeId> &edges, std::type_identity_t<T> r, const FilletOptions<T> &o = {})
        -> BuildResult<SolidId>
    {
        return fillet(m, s, std::span<const EdgeId>{edges}, r, o);
    }

    /**
     * @brief The solid `s` with its edges `edges` chamfered at the distance d1 on the face `f`
     * and d2 on the other face of each edge (every edge a boundary of `f`), as `fillet`.
     */
    template <std::floating_point T>
    [[nodiscard]] auto chamfer(Model<T> &m, SolidId s, std::span<const EdgeId> edges, FaceId f, std::type_identity_t<T> d1, std::type_identity_t<T> d2,
                               const FilletOptions<T> &o = {}) -> BuildResult<SolidId>
    {
        const std::string op = "chamfer";
        return blend_detail::guarded(m, op, [&] {
            if (!m.alive(f))
                fillet_detail::refuse<T>(BuildErrc::InvalidId, op, "dead or invalid face", {f});
            return blend_detail::blend(
                m, s, edges, [&](SolidId ns, EdgeId e, detail::ShapeCopier<T> &c, const FilletOptions<T> &oe) {
                    return fillet_detail::chamfer(m, ns, e, c.copy(f), d1, d2, oe);
                },
                o, op);
        });
    }

    template <std::floating_point T>
    [[nodiscard]] auto chamfer(Model<T> &m, SolidId s, const std::vector<EdgeId> &edges, FaceId f, std::type_identity_t<T> d1, std::type_identity_t<T> d2,
                               const FilletOptions<T> &o = {}) -> BuildResult<SolidId>
    {
        return chamfer(m, s, std::span<const EdgeId>{edges}, f, d1, d2, o);
    }

    /**
     * @brief The solid `s` with its edges `edges` chamfered at the distance d on both faces of
     * each edge, as `fillet`.
     */
    template <std::floating_point T>
    [[nodiscard]] auto chamfer(Model<T> &m, SolidId s, std::span<const EdgeId> edges, std::type_identity_t<T> d, const FilletOptions<T> &o = {})
        -> BuildResult<SolidId>
    {
        const std::string op = "chamfer";
        return blend_detail::guarded(m, op, [&] {
            return blend_detail::blend(
                m, s, edges,
                [&](SolidId ns, EdgeId e, detail::ShapeCopier<T> &, const FilletOptions<T> &oe) {
                    // the face of d1 any face of the edge
                    const auto uses = shell_edge_uses(m, m.solid(ns).outer);
                    FaceId f{};
                    if (const auto it = uses.find(e); it != uses.end() && !it->second.empty())
                        f = it->second.front().face;
                    else
                        for (auto v : m.solid(ns).voids)
                            if (const auto u = shell_edge_uses(m, v); u.contains(e) && !u.at(e).empty())
                                f = u.at(e).front().face;
                    return fillet_detail::chamfer(m, ns, e, f, d, d, oe);
                },
                o, op);
        });
    }

    template <std::floating_point T>
    [[nodiscard]] auto chamfer(Model<T> &m, SolidId s, const std::vector<EdgeId> &edges, std::type_identity_t<T> d, const FilletOptions<T> &o = {})
        -> BuildResult<SolidId>
    {
        return chamfer(m, s, std::span<const EdgeId>{edges}, d, o);
    }
}
