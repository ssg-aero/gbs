#pragma once
/**
 * @file blend.h
 * @brief Fillets and chamfers of edges of a solid (`fillet`, `chamfer`): the faces of each edge
 * cut along the contacts of its blend, the faces at its ends cut by the blend, the face of the
 * blend between them; a new solid (stage 4, `docs/sources/design/brep_palier4.md`, § 4.2,
 * § 4.5 and § 4.6; notes `brep4_pr06_fillet_topology.md`, `brep4_pr07_fillet_corners.md`).
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
 *    - at a **corner**, a vertex of three edges all filleted (`fillet`): the three centre lines
 *      meet at the centre of a sphere of the radius, at r from the three faces; each fillet ends
 *      on its arc through that centre, a great circle of the sphere, and the face of the sphere
 *      bounded by the three arcs fills the corner (convex or concave alike), its poles and its
 *      seam clear of the arcs (`corner_frame`); the contacts of two fillets on their common face
 *      meet at the foot of the centre on it;
 *    - a closed edge has no end: the face of the blend runs its whole loop, across a seam.
 * 3. **Faces.** In the wire of each face of an edge, the run of co-edges between the ends of the
 *    contacts that holds the blended edge is replaced by the contacts (chained along a chain of
 *    edges, from one end to the other; a closed chain clear of the boundary replaces a whole
 *    wire made of blended edges); in the wire of F3, the run through the vertex by the end curve.
 *    Convex and concave edges alike: the face loses the band along the edge (convex) or the
 *    corner it gets back (concave, F3).
 * 4. The face of each blend, on its surface, bounded by its contacts and its ends, oriented as the
 *    faces of its edge, and the faces of the corners; the edges and vertices left out erased, the
 *    edges at least as loose as the faces they bound; the solid checked (`check`).
 *
 * A contact that does not meet the face at a free end within the blend continued is tried again
 * with the blend continued 2, 4 and 8 times further (from one radius or distance at least,
 * `FilletOptions::extension` being 0), keeping the last geometry that could be computed.
 *
 * Refused (`NonGeneric`): two blended edges meeting at a vertex but not tangentially, three
 * blended edges at a vertex but fillets at a vertex of three edges (chamfers, a fourth edge),
 * a corner of convex and concave edges, of fillets whose centre lines do not meet within the
 * blends continued or whose contacts do not meet on a face, more than three blended edges at a
 * vertex, a free end at a vertex of more than three edges, a contact that does not meet
 * the face at an end within the blend continued or meets it beyond the edge between them (a blend
 * wider than a face next to it) or at a vertex of it (blends meeting there, as a full round), a
 * contact ending inside a face or at a vertex of it, a result that is not a valid solid; the
 * refusals of the geometry. On failure the model is left as it was.
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
#include <initializer_list>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <variant>
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
            Loop,   ///< a closed edge: the seam of the face of the blend
            Corner, ///< three fillets at a vertex: the sphere of the corner
            Lune    ///< a kink of the radius law between two pieces of a fillet: the lune of the ball
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
            // a free end: the third face, the edges cut and the end curve on it (a corner: the
            // end curve on the sphere)
            FaceId f3;
            std::shared_ptr<Curve<T, 2>> on_f3;
        };

        template <std::floating_point T>
        struct Blend
        {
            EdgeBlend<T> g;
            std::array<End<T>, 2> ends; ///< at a (lower parameter) and at b
            std::array<EdgeId, 2> contacts;
            /// a piece of a fillet in pieces (`EdgeBlend::pieces`): `g` the piece alone; its ends
            /// at a vertex of its edge, or at a kink of the law (`kink`); the lune of the law at an
            /// end at a vertex (a kink there: between the blends of two edges of a chain)
            bool part{};
            std::array<bool, 2> kink{};
            std::array<std::optional<BlendLune<T>>, 2> lune_at;
        };

        /// The face of a lune, oriented later from the face of the blend `owner` along their common arc.
        struct LuneFace
        {
            ShellId shell;
            FaceId face;
            EdgeId arc;             ///< its arc on the blend `owner`
            Orientation in_wire{};  ///< the use of `arc` in its wire
            std::size_t owner{};    ///< the index of that blend
        };

        /// The blend g as one surface: its piece i alone (its surface, its contacts on it, on [a, b] of the piece).
        template <std::floating_point T>
        auto piece_view(const EdgeBlend<T> &g, std::size_t i) -> EdgeBlend<T>
        {
            EdgeBlend<T> v = g;
            const auto &p = g.pieces[i];
            v.pieces.clear();
            v.lunes.clear();
            v.surface = p.surface;
            v.a = p.a, v.b = p.b;
            v.deviation = p.deviation;
            for (std::size_t k = 0; k < 2; ++k)
                v.contacts[k].on_blend = p.on_blend[k];
            return v;
        }

        /// The blend g as one surface around t: itself, or its piece holding t.
        template <std::floating_point T>
        auto view_at(const EdgeBlend<T> &g, T t) -> EdgeBlend<T>
        {
            if (g.pieces.empty())
                return g;
            std::size_t best = 0;
            T gap = std::numeric_limits<T>::infinity();
            for (std::size_t i = 0; i < g.pieces.size(); ++i)
            {
                const T d = t < g.pieces[i].a ? g.pieces[i].a - t : t > g.pieces[i].b ? t - g.pieces[i].b : T(0);
                if (d < gap)
                    gap = d, best = i;
            }
            return piece_view(g, best);
        }

        /**
         * The solid being edited: its shells, the orientation and the shell of each of its faces
         * (fixed), and an index of its topology (the faces of each edge, the edges at each vertex,
         * the first pcurve of each edge on each face), rebuilt on the first query after an edit
         * (`invalidate`) rather than scanning every wire at each query.
         */
        template <std::floating_point T>
        struct SolidView
        {
            const Model<T> &m;
            std::vector<ShellId> shells;
            std::map<FaceId, std::pair<ShellId, Orientation>> uses; ///< the faces of the solid as it was given

            SolidView(const Model<T> &m_, SolidId s) : m{m_}
            {
                shells.push_back(m.solid(s).outer);
                shells.insert(shells.end(), m.solid(s).voids.begin(), m.solid(s).voids.end());
                for (auto sh : shells)
                    for (const auto &fu : m.shell(sh).faces)
                        uses.emplace(fu.face, std::pair{sh, fu.orient});
            }

            /// The faces of the shells as they are now (the faces of the blends included).
            auto faces() const -> std::vector<FaceUse>
            {
                std::vector<FaceUse> r;
                for (auto sh : shells)
                    r.insert(r.end(), m.shell(sh).faces.begin(), m.shell(sh).faces.end());
                return r;
            }
            auto orient(FaceId f) const -> Orientation
            {
                const auto it = uses.find(f);
                return it == uses.end() ? Orientation::Forward : it->second.second;
            }
            auto shell_of(FaceId f) const -> ShellId
            {
                const auto it = uses.find(f);
                return it == uses.end() ? ShellId{} : it->second.first;
            }
            /// The index out of date: the model was edited (an edge split, a wire changed).
            void invalidate() const { stale_ = true; }
            /// The faces using the edge e (a seam: twice).
            auto users(EdgeId e) const -> std::vector<FaceId>
            {
                index();
                const auto it = users_.find(e);
                return it == users_.end() ? std::vector<FaceId>{} : it->second;
            }
            /// The edges at the vertex v, but degenerate ones, each once.
            auto edges_at(VertexId v) const -> std::vector<EdgeId>
            {
                index();
                const auto it = at_.find(v);
                return it == at_.end() ? std::vector<EdgeId>{} : std::vector<EdgeId>{it->second.begin(), it->second.end()};
            }
            /// The pcurve of the edge e on the face f (its first use).
            auto pcurve(FaceId f, EdgeId e) const -> std::shared_ptr<Curve<T, 2>>
            {
                index();
                const auto it = pc_.find({f, e});
                return it == pc_.end() ? nullptr : it->second;
            }
            /// The outward normal of the face f at uv of its surface s.
            auto outward(FaceId f, const Surface<T, 3> &s, const point<T, 2> &uv) const -> point<T, 3>
            {
                const auto n = normal_direction(s, uv[0], uv[1]);
                return orient(f) == Orientation::Forward ? n : T(-1) * n;
            }

        private:
            mutable bool stale_{true};
            mutable std::map<EdgeId, std::vector<FaceId>> users_;
            mutable std::map<VertexId, std::set<EdgeId>> at_;
            mutable std::map<std::pair<FaceId, EdgeId>, std::shared_ptr<Curve<T, 2>>> pc_;

            void index() const
            {
                if (!stale_)
                    return;
                users_.clear(), at_.clear(), pc_.clear();
                for (const auto &fu : faces())
                    for (auto w : m.face(fu.face).wires)
                        for (const auto &ce : m.wire(w).coedges)
                        {
                            users_[ce.edge].push_back(fu.face);
                            pc_.try_emplace({fu.face, ce.edge}, ce.pcurve);
                            const auto &e = m.edge(ce.edge);
                            if (!e.degenerate)
                                at_[e.v1].insert(ce.edge), at_[e.v2].insert(ce.edge);
                        }
                stale_ = false;
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

        /// A pcurve of an edge on a surface, for `settle`.
        template <std::floating_point T>
        using OnSurface = std::pair<const Surface<T, 3> *, std::shared_ptr<Curve<T, 2>>>;

        /**
         * A new edge on the curve c (none: a degenerate edge) over [u1, u2] from v1 to v2, its
         * tolerance tol raised to its gaps to its pcurves and its vertices to cover it (`settle`).
         */
        template <std::floating_point T>
        auto add_edge(Model<T> &m, std::shared_ptr<Curve<T, 3>> c, T u1, T u2, VertexId v1, VertexId v2, T tol,
                      std::type_identity_t<std::initializer_list<OnSurface<T>>> pcs) -> EdgeId
        {
            Edge<T> e;
            e.degenerate = !c;
            e.curve = std::move(c);
            e.u1 = u1, e.u2 = u2, e.v1 = v1, e.v2 = v2, e.tol = tol;
            const auto id = m.add(std::move(e));
            const std::vector<OnSurface<T>> list(pcs);
            settle(m, id, std::span<const OnSurface<T>>{list});
            return id;
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
                    if (const auto [u, d] = fillet_detail::edge_point(e, p); d < gap)
                        gap = d, at = u, best = ce.edge;
                }
            const T reach = T(100) * tol;
            if (!best || gap > reach)
                refuse<T>(op, "a contact ends inside a face (a chain of edges whose contacts leave their faces apart from an edge)", {f});
            const auto &e = m.edge(*best);
            const T eps = T(1e-9) * (e.u2 - e.u1);
            if (!(at > e.u1 + eps && at < e.u2 - eps))
                refuse<T>(op, "a contact ends at a vertex of a face (blends meeting there: not handled)", {f, *best});
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
                        if (fillet_detail::edge_point(e, p)[1] <= T(100) * tol)
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
                const auto [u, d] = fillet_detail::edge_point(e, p);
                const T eps = T(1e-9) * (e.u2 - e.u1);
                if (!(d <= T(100) * tol))
                    refuse<T>(op, "a blend wider than a face at an end of its edge", {g.edge, side[k]});
                if (!(u > e.u1 + eps && u < e.u2 - eps))
                    refuse<T>(op, "a contact meets the face at an end of its edge at a vertex of it (blends meeting there, as a full round: not handled)",
                              {g.edge, side[k]});
                end.v[k] = m.add(Vertex<T>{p, tol});
                cover(m, end.v[k], e.curve->value(u), std::max(tol, e.tol));
                if (!split_detail::split_at_vertex(m, side[k], u, end.v[k]))
                    refuse<T>(op, "an edge that cannot be split where a contact ends", {side[k]});
                sv.invalidate();
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
            // n + 1 points across, n doubled from 8: the points at the middles, solved for the check,
            // are the new points
            std::size_t n = 8;
            std::vector<Pt> pts{first};
            std::vector<T> s{T(0)};
            for (std::size_t i = 1; i <= n; ++i)
            {
                s.push_back(T(i) / T(n));
                pts.push_back(solve(s.back(), pts.back()));
            }
            for (;;)
            {
                points_vector<T, 3> q3;
                points_vector<T, 2> qb, qf;
                for (std::size_t i = 0; i <= n; ++i)
                {
                    q3.push_back(pts[i].p);
                    qb.push_back(across(g, pts[i].t, s[i]));
                    qf.push_back(pts[i].uv3);
                }
                auto k3 = std::make_shared<BSCurve<T, 3>>(interpolate(q3, s, 3));
                auto kb = std::make_shared<BSCurve<T, 2>>(interpolate(qb, s, 3));
                auto kf = std::make_shared<BSCurve<T, 2>>(interpolate(qf, s, 3));
                std::vector<Pt> mid;
                T d{};
                for (std::size_t i = 0; i < n; ++i)
                {
                    const T sm = (s[i] + s[i + 1]) / T(2);
                    mid.push_back(solve(sm, pts[i]));
                    const auto &x = mid.back();
                    const auto wb = kb->value(sm), wf = kf->value(sm);
                    d = std::max({d, norm(k3->value(sm) - x.p), norm(g.surface->value(wb[0], wb[1]) - x.p), norm(s3->value(wf[0], wf[1]) - x.p)});
                }
                if (d < dev)
                    dev = d, c3 = k3, cb = kb, cf = kf;
                if (dev <= tol / T(10) || 2 * n > 512)
                    break;
                std::vector<Pt> p2;
                std::vector<T> s2;
                for (std::size_t i = 0; i <= n; ++i)
                {
                    p2.push_back(pts[i]), s2.push_back(s[i]);
                    if (i < n)
                        p2.push_back(mid[i]), s2.push_back((s[i] + s[i + 1]) / T(2));
                }
                pts = std::move(p2), s = std::move(s2), n *= 2;
            }
            if (!(dev <= tol))
                refuse<T>(op, "the end curve of a blend cannot be approximated within the tolerance", {g.edge, f3});
            end.edge = add_edge(m, std::shared_ptr<Curve<T, 3>>{c3}, T(0), T(1), end.v[0], end.v[1], std::max(tol, dev),
                                {{g.surface.get(), cb}, {m.face(f3).surface.get(), cf}});
            end.on_blend = cb;
            end.on_f3 = cf;
            end.forward = true;
        }

        /// Two blends joined at the vertex v: the arc across both there, their contacts meeting at its ends.
        template <std::floating_point T>
        auto lune_face(Model<T> &m, const SolidView<T> &sv, Blend<T> &b1, std::size_t j1, Blend<T> &b2, std::size_t j2, const SphericalSurface<T> &sphere,
                       const std::shared_ptr<SphericalSurface<T>> &owner, T tol, const std::string &op) -> LuneFace;

        template <std::floating_point T>
        void joint(Model<T> &m, const SolidView<T> &sv, Blend<T> &b1, std::size_t j1, Blend<T> &b2, std::size_t j2, T tol, T approx_tol, const std::string &op,
                   std::vector<LuneFace> &lunes, std::size_t i1)
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
            if (k2[0] == k2[1] || gap[0] > limit || gap[1] > limit)
                refuse<T>(op, "blended edges meeting at a vertex but not tangentially, or blends that do not join there (not handled)", {g1.edge, g2.edge, e1.vertex});
            std::optional<BlendCrease<T>> crease;
            if (arc_gap > limit)
            {
                // the arcs of the blends differ at the vertex (a kink of the radius law there, or
                // a radius varying where the curvature of the edges jumps: the centre line has a
                // corner): the balls next to the vertex cover the middle of the lune of its ball
                // between the arcs (the envelopes on either side meet along a crease), or not (the
                // lune is part of the fillet)
                const auto c = g1.spine->value(e1.t);
                const auto p1 = g1.contacts[0].curve->value(e1.t), p2 = g1.contacts[1].curve->value(e1.t);
                const T r = norm(p1 - c);
                const auto a1 = CurveOnSurface<T, 3>(pc1, g1.surface).value(T(0.5)), a2 = CurveOnSurface<T, 3>(pc2, g2.surface).value(T(0.5));
                const auto q = c + r * fillet_detail::unit((a1 + a2) / T(2) - c);
                auto covers = [&](const EdgeBlend<T> &g, T t) {
                    const auto cc = g.spine->value(t);
                    return norm(q - cc) < norm(g.contacts[0].curve->value(t) - cc) - T(1e-12) * (T(1) + r);
                };
                const T d1 = T(1e-4) * (g1.b - g1.a), d2 = T(1e-4) * (g2.b - g2.a);
                const bool covered = covers(g1, j1 == 1 ? e1.t - d1 : e1.t + d1) || covers(g2, j2 == 1 ? e2.t - d2 : e2.t + d2);
                if (!covered)
                {
                    // the lune of the ball: its sphere, its seam and poles away from the lune
                    std::shared_ptr<SphericalSurface<T>> sphere;
                    if (const auto &l = b1.lune_at[j1] ? b1.lune_at[j1] : b2.lune_at[j2])
                        sphere = l->sphere;
                    else
                    {
                        constexpr T pi = std::numbers::pi_v<T>;
                        const auto x = fillet_detail::unit((a1 + a2) / T(2) - c), z = fillet_detail::unit(cross(x, p2 - p1));
                        sphere = std::make_shared<SphericalSurface<T>>(r, ax2<T, 3>{c, z, x}, std::array<T, 4>{-pi, pi, -pi / T(2), pi / T(2)});
                    }
                    for (std::size_t k = 0; k < 2; ++k)
                    {
                        e1.v[k] = m.add(Vertex<T>{g1.contacts[k].curve->value(e1.t), std::max(tol, gap[k])});
                        cover(m, e1.v[k], g2.contacts[k2[k]].curve->value(e2.t), tol);
                        e2.v[k2[k]] = e1.v[k];
                    }
                    lunes.push_back(lune_face(m, sv, b1, j1, b2, j2, *sphere, sphere, tol, op));
                    lunes.back().owner = i1;
                    return;
                }
                crease = fillet_detail::crease_curve(g1.surface, e1.t, g2.surface, e2.t, k2[0] == 0, tol);
                if (!crease)
                    refuse<T>(op, "the fillets of two edges of a chain do not meet along a crease at their vertex", {g1.edge, g2.edge, e1.vertex});
            }
            for (std::size_t k = 0; k < 2; ++k)
            {
                e1.v[k] = m.add(Vertex<T>{g1.contacts[k].curve->value(e1.t), std::max(tol, gap[k])});
                cover(m, e1.v[k], g2.contacts[k2[k]].curve->value(e2.t), tol);
                e2.v[k2[k]] = e1.v[k];
                e1.tc[k] = e1.t;
                e2.tc[k] = e2.t;
            }
            // (along the crease from contact 0 of b1: on b2 its pcurve runs from the contact of b2
            // that meets contact 0 of b1, as the arc pc2)
            const auto q1 = crease ? crease->on_pieces[0] : pc1, q2 = crease ? crease->on_pieces[1] : pc2;
            const auto eid = add_edge(m, std::shared_ptr<Curve<T, 3>>{std::make_shared<CurveOnSurface<T, 3>>(q1, g1.surface)}, T(0), T(1), e1.v[0], e1.v[1],
                                      std::max(tol, crease ? crease->deviation * (T(1) + T(1e-6)) : T(0)), {{g1.surface.get(), q1}, {g2.surface.get(), q2}});
            e1.kind = e2.kind = EndKind::Joint;
            e1.edge = e2.edge = eid;
            e1.on_blend = q1;
            e2.on_blend = q2;
            e1.forward = true;
            e2.forward = k2[0] == 0;
        }

        /**
         * The frame at o of the sphere of a corner, its poles and its seam clear of the patch of
         * the corner: the spherical triangle of the feet f (directions from o; its sides the arcs,
         * great circles), `dirs` the directions of points along the arcs. The pole z: of 2048
         * directions spread over the sphere (Fibonacci), the farthest from the arcs with neither
         * it nor its opposite in the patch; the seam (u = +-pi) in the middle of the widest gap
         * the arcs leave round z, x opposite. (x toward the sum of the feet, z across, held for a
         * box but nothing kept the pole out of the patch nor the seam off the arcs.) At a sharp
         * corner the patch nearly fills a hemisphere and every pole is a few degrees from the
         * arcs (4 at a spike of faces 160 degrees apart): their pcurves are interpolated finer
         * there (`corner`). None when no pole is more than 3 degrees from the arcs.
         */
        template <std::floating_point T>
        auto corner_frame(const point<T, 3> &o, const std::array<point<T, 3>, 3> &f, const std::vector<point<T, 3>> &dirs) -> std::optional<ax2<T, 3>>
        {
            constexpr T pi = std::numbers::pi_v<T>;
            auto inside = [&](const point<T, 3> &q) {
                for (std::size_t i = 0; i < 3; ++i)
                {
                    const auto n = cross(f[i], f[(i + 1) % 3]);
                    if (!((q * n) * (f[(i + 2) % 3] * n) > T(0)))
                        return false;
                }
                return true;
            };
            constexpr int count = 2048;
            const T turn = pi * (T(3) - std::sqrt(T(5)));
            point<T, 3> z{};
            T best = T(1);
            for (int i = 0; i < count; ++i)
            {
                const T h = T(1) - (T(2) * T(i) + T(1)) / T(count), s = std::sqrt(T(1) - h * h);
                const point<T, 3> q{s * std::cos(turn * T(i)), s * std::sin(turn * T(i)), h};
                if (inside(q) || inside(T(-1) * q))
                    continue;
                T near{};
                for (const auto &p : dirs)
                    near = std::max(near, std::abs(p * q));
                if (near < best)
                    best = near, z = q;
            }
            if (!(best < std::cos(T(3) * pi / T(180))))
                return std::nullopt;
            auto e1 = cross(z, std::abs(z[0]) < T(0.9) ? point<T, 3>{T(1), T(0), T(0)} : point<T, 3>{T(0), T(1), T(0)});
            e1 = e1 / norm(e1);
            const auto e2 = cross(z, e1);
            std::vector<T> lon;
            for (const auto &p : dirs)
                lon.push_back(std::atan2(p * e2, p * e1));
            std::ranges::sort(lon);
            T gap = lon.front() + T(2) * pi - lon.back(), seam = lon.back() + gap / T(2);
            for (std::size_t i = 1; i < lon.size(); ++i)
                if (lon[i] - lon[i - 1] > gap)
                    gap = lon[i] - lon[i - 1], seam = (lon[i] + lon[i - 1]) / T(2);
            const auto x = T(-1) * (std::cos(seam) * e1 + std::sin(seam) * e2);
            return ax2<T, 3>{o, z, x};
        }

        /**
         * Three fillets at a vertex of three edges (the corner of a box), at the ends j of the
         * blends i of `at`: their centre lines meet at the centre c of a sphere of their radius,
         * at r from the three faces (Gauss-Newton on the parameters of the three lines). Each
         * fillet ends on its arc through c (a great circle of the sphere); the contacts of two
         * fillets on their common face meet at the foot of c on it. The face of the sphere,
         * bounded by the three arcs (`corner_frame`: its poles and its seam clear of them), the
         * solid inside it at a convex corner and outside at a concave one, is returned with the
         * shell it goes in. Refused: centre lines that do not meet within the blends continued
         * (the centre past a vertex by more than `FilletOptions::extension`), contacts that do not
         * meet on a face, a face not shared by two of the fillets, edges convex and concave at the
         * corner.
         */
        template <std::floating_point T>
        auto corner(Model<T> &m, const SolidView<T> &sv, std::vector<Blend<T>> &bls, const std::array<std::pair<std::size_t, std::size_t>, 3> &at, T tol,
                    T approx_tol, const std::string &op) -> std::pair<ShellId, FaceUse>
        {
            const VertexId vx = bls[at[0].first].ends[at[0].second].vertex;
            std::vector<ShapeId> involved{vx};
            T limit = T(10) * tol;
            for (const auto &[i, j] : at)
            {
                involved.push_back(bls[i].g.edge);
                limit += bls[i].g.exact ? T(0) : approx_tol;
            }
            for (const auto &[i, j] : at)
                if (bls[i].g.convex != bls[at[0].first].g.convex)
                    refuse<T>(op, "a corner of convex and concave blended edges (not handled)", involved);
            // the centre: the three centre lines meeting, from the vertex
            Vec<T> x(3);
            for (std::size_t q = 0; q < 3; ++q)
                x[Eigen::Index(q)] = bls[at[q].first].ends[at[q].second].t;
            auto spine = [&](std::size_t q, T t) {
                const auto &g = bls[at[q].first].g;
                return g.spine->value(std::clamp(t, g.a, g.b));
            };
            const T res = thick_detail::gauss_newton<T>(x, [&](const Vec<T> &y, Vec<T> &f, Mat<T> &J) {
                f.resize(6);
                J = Mat<T>::Zero(6, 3);
                for (std::size_t q = 0; q < 2; ++q)
                {
                    const auto &ga = bls[at[q].first].g, &gb = bls[at[q + 1].first].g;
                    const auto d = ga.spine->value(y[Eigen::Index(q)]) - gb.spine->value(y[Eigen::Index(q + 1)]);
                    const auto da = ga.spine->value(y[Eigen::Index(q)], 1), db = gb.spine->value(y[Eigen::Index(q + 1)], 1);
                    for (int c = 0; c < 3; ++c)
                    {
                        f[Eigen::Index(3 * q) + c] = d[std::size_t(c)];
                        J(Eigen::Index(3 * q) + c, Eigen::Index(q)) = da[std::size_t(c)];
                        J(Eigen::Index(3 * q) + c, Eigen::Index(q + 1)) = -db[std::size_t(c)];
                    }
                }
            });
            point<T, 3> c{};
            for (std::size_t q = 0; q < 3; ++q)
                c = c + spine(q, x[Eigen::Index(q)]) / T(3);
            for (std::size_t q = 0; q < 3; ++q)
            {
                const auto &g = bls[at[q].first].g;
                const T t = x[Eigen::Index(q)];
                if (!(res <= limit) || !(t >= g.a && t <= g.b) || norm(spine(q, t) - c) > limit)
                    refuse<T>(op, "fillets at a corner whose centre lines do not meet within the blends continued", involved);
            }
            // the feet of the centre on the faces: where the contacts of two fillets meet
            std::map<FaceId, std::vector<point<T, 3>>> feet;
            T r{};
            for (std::size_t q = 0; q < 3; ++q)
            {
                const auto &g = bls[at[q].first].g;
                for (std::size_t k = 0; k < 2; ++k)
                {
                    const auto p = g.contacts[k].curve->value(x[Eigen::Index(q)]);
                    feet[g.contacts[k].face].push_back(p);
                    r += norm(p - c) / T(6);
                }
            }
            if (feet.size() != 3 || std::ranges::any_of(feet, [](const auto &f) { return f.second.size() != 2; }))
                refuse<T>(op, "fillets at a corner whose faces are not shared two by two", involved);
            std::map<FaceId, VertexId> foot;
            std::array<point<T, 3>, 3> fdir; // the directions of the feet from c
            std::size_t nf = 0;
            for (const auto &[f, ps] : feet)
            {
                const T gap = norm(ps[0] - ps[1]);
                if (gap > limit || std::abs(norm(ps[0] - c) - r) > limit || std::abs(norm(ps[1] - c) - r) > limit)
                    refuse<T>(op, "fillets at a corner whose contacts do not meet on a face", involved);
                foot[f] = m.add(Vertex<T>{(ps[0] + ps[1]) / T(2), std::max(tol, gap)});
                const auto e = (ps[0] + ps[1]) / T(2) - c;
                fdir[nf++] = e / norm(e);
            }
            // the arcs: the ends of the fillets on their surfaces, and the directions of their points
            std::array<std::shared_ptr<CurveOnSurface<T, 3>>, 3> curves;
            std::vector<point<T, 3>> dirs;
            for (std::size_t q = 0; q < 3; ++q)
            {
                auto &end = bls[at[q].first].ends[at[q].second];
                const auto &g = bls[at[q].first].g;
                const T t = x[Eigen::Index(q)];
                end.on_blend = fillet_detail::segment(across(g, t, T(0)), across(g, t, T(1)), T(0), T(1));
                curves[q] = std::make_shared<CurveOnSurface<T, 3>>(end.on_blend, g.surface);
                for (int i = 0; i <= 32; ++i)
                {
                    const auto e = curves[q]->value(T(i) / T(32)) - c;
                    dirs.push_back(e / norm(e));
                }
            }
            // the sphere, its poles and its seam clear of the corner
            const auto frame = corner_frame(c, fdir, dirs);
            if (!frame)
                refuse<T>(op, "a corner whose sphere has no pole clear of its arcs", involved);
            constexpr T pi = std::numbers::pi_v<T>;
            const auto sphere = std::make_shared<SphericalSurface<T>>(r, *frame, std::array<T, 4>{-pi, pi, -pi / T(2), pi / T(2)});
            // the arcs on the sphere too (`sphere_pcurve`: interpolated until within tol / 10
            // between their points), the edge as loose as the gap found at the points and between
            // them
            std::vector<CoEdge<T>> arcs;
            for (std::size_t q = 0; q < 3; ++q)
            {
                auto &bl = bls[at[q].first];
                auto &end = bl.ends[at[q].second];
                const auto &g = bl.g;
                const T t = x[Eigen::Index(q)];
                const auto &curve = *curves[q];
                end.kind = EndKind::Corner;
                end.t = t;
                end.tc = {t, t};
                for (std::size_t k = 0; k < 2; ++k)
                {
                    end.v[k] = foot.at(g.contacts[k].face);
                    cover(m, end.v[k], g.contacts[k].curve->value(t), tol);
                }
                const auto [on, gap] = fillet_detail::sphere_pcurve(*sphere, curve, tol);
                if (gap > limit)
                    refuse<T>(op, "fillets at a corner whose radius varies there (their ends are not on the sphere)", involved);
                // a radius dropping toward the corner: the end of the fillet a circle of its ball
                // shifted toward the corner, the balls before it covering the sphere there
                {
                    const T dt = T(1e-6) * (g.b - g.a), toward = at[q].second == 1 ? T(1) : T(-1);
                    const auto tc = fillet_detail::unit(g.spine->value(t + toward * dt) - g.spine->value(t - toward * dt));
                    const T l = ((g.contacts[0].curve->value(t) - c) * tc + (g.contacts[1].curve->value(t) - c) * tc) / T(2);
                    if (l > limit)
                        refuse<T>(op, "fillets at a corner whose radius drops toward it (the balls next to the corner cover its sphere)", involved);
                }
                end.on_f3 = on;
                end.edge = add_edge(m, std::shared_ptr<Curve<T, 3>>{curves[q]}, T(0), T(1), end.v[0], end.v[1], tol,
                                    {{g.surface.get(), end.on_blend}, {sphere.get(), end.on_f3}});
                auto &e = m.edge(end.edge);
                e.tol = std::max(e.tol, gap * (T(1) + T(1e-6)));
                for (auto v : {e.v1, e.v2})
                    cover(m, v, m.vertex(v).pnt, e.tol);
                end.forward = true;
                arcs.push_back(CoEdge<T>{end.edge, Orientation::Forward, end.on_f3});
            }
            // the face of the sphere: the three arcs chained
            std::vector<CoEdge<T>> cs{arcs[0]};
            std::vector<bool> used{true, false, false};
            for (int n = 0; n < 2; ++n)
            {
                const auto last = ends_of(m, cs.back())[1];
                bool found = false;
                for (std::size_t q = 0; q < 3 && !found; ++q)
                {
                    if (used[q])
                        continue;
                    const auto &e = m.edge(arcs[q].edge);
                    if (e.v1 == last || e.v2 == last)
                    {
                        used[q] = true;
                        found = true;
                        cs.push_back(CoEdge<T>{arcs[q].edge, e.v1 == last ? Orientation::Forward : Orientation::Reversed, arcs[q].pcurve});
                    }
                }
                if (!found)
                    refuse<T>(op, "the arcs at a corner do not close", involved);
            }
            if (ends_of(m, cs.back())[1] != ends_of(m, cs.front())[0])
                refuse<T>(op, "the arcs at a corner do not close", involved);
            if (uv_signed_area(m, std::span<const CoEdge<T>>{cs}) < T(0))
                cs = reversed(cs);
            Wire<T> wire;
            wire.coedges = std::move(cs);
            wire.closed = true;
            const auto wid = m.add(std::move(wire));
            const auto fid = m.add(Face<T>{sphere, {wid}, tol, false});
            // the solid inside the sphere at a convex corner, outside at a concave one: the face
            // forward where the normal of the sphere (at the middle of an arc, clear of the poles)
            // points out of it at a convex corner, in at a concave one
            const auto pm = curves[0]->value(T(0.5));
            const auto um = sphere->parameters(pm);
            const bool out = normal_direction(*sphere, um[0], um[1]) * (pm - c) > T(0);
            const auto &g0 = bls[at[0].first].g;
            return {sv.shell_of(g0.contacts[0].face), FaceUse{fid, out == g0.convex ? Orientation::Forward : Orientation::Reversed}};
        }

        /**
         * The lune of the ball `sphere` between the end j1 of b1 and the end j2 of b2, their
         * contacts already on the same vertices (`End::v`): on each, the arc across it at its end
         * (two characteristic circles of the ball, both through the contacts), an edge with its
         * pcurves on its blend and on the sphere (`sphere_pcurve`, the edge as loose as the gap
         * found); the face of the sphere between both arcs, oriented as the blends (the solid
         * inside the ball along a convex edge, outside along a concave one).
         */
        template <std::floating_point T>
        auto lune_face(Model<T> &m, const SolidView<T> &sv, Blend<T> &b1, std::size_t j1, Blend<T> &b2, std::size_t j2, const SphericalSurface<T> &sphere,
                       const std::shared_ptr<SphericalSurface<T>> &owner, T tol, const std::string &op) -> LuneFace
        {
            std::array<CoEdge<T>, 2> arcs;
            std::shared_ptr<CurveOnSurface<T, 3>> first;
            for (std::size_t q = 0; q < 2; ++q)
            {
                auto &bl = q == 0 ? b1 : b2;
                auto &end = bl.ends[q == 0 ? j1 : j2];
                const auto &g = bl.g;
                end.kind = EndKind::Lune;
                end.tc = {end.t, end.t};
                end.on_blend = fillet_detail::segment(across(g, end.t, T(0)), across(g, end.t, T(1)), T(0), T(1));
                const auto curve = std::make_shared<CurveOnSurface<T, 3>>(end.on_blend, g.surface);
                if (q == 0)
                    first = curve;
                const auto [on, gap] = fillet_detail::sphere_pcurve(sphere, *curve, tol);
                end.on_f3 = on;
                end.edge = add_edge(m, std::shared_ptr<Curve<T, 3>>{curve}, T(0), T(1), end.v[0], end.v[1], tol, {{g.surface.get(), end.on_blend}, {&sphere, on}});
                auto &e = m.edge(end.edge);
                e.tol = std::max(e.tol, gap * (T(1) + T(1e-6)));
                for (auto v : {e.v1, e.v2})
                    cover(m, v, m.vertex(v).pnt, e.tol);
                end.forward = true;
                arcs[q] = CoEdge<T>{end.edge, Orientation::Forward, on};
            }
            // the two arcs end to end round the lune
            std::vector<CoEdge<T>> cs{arcs[0]};
            const auto last = ends_of(m, arcs[0])[1];
            const auto &e2 = m.edge(arcs[1].edge);
            if (!(e2.v1 == last || e2.v2 == last) || (e2.v1 == last ? e2.v2 : e2.v1) != ends_of(m, arcs[0])[0])
                refuse<T>(op, "the arcs of a lune do not close", {b1.g.edge, b2.g.edge});
            cs.push_back(CoEdge<T>{arcs[1].edge, e2.v1 == last ? Orientation::Forward : Orientation::Reversed, arcs[1].pcurve});
            if (uv_signed_area(m, std::span<const CoEdge<T>>{cs}) < T(0))
                cs = reversed(cs);
            const auto in_wire = std::ranges::find(cs, arcs[0].edge, &CoEdge<T>::edge)->orient;
            Wire<T> wire;
            wire.coedges = std::move(cs);
            wire.closed = true;
            const auto wid = m.add(std::move(wire));
            const auto fid = m.add(Face<T>{owner, {wid}, tol, false});
            (void)first;
            return LuneFace{sv.shell_of(b1.g.contacts[0].face), fid, arcs[0].edge, in_wire, 0};
        }

        /// The contact k of the blend as an edge between the ends of the blend.
        template <std::floating_point T>
        auto contact_edge(Model<T> &m, const Blend<T> &bl, std::size_t k, T tol) -> EdgeId
        {
            const auto &ct = bl.g.contacts[k];
            return add_edge(m, ct.curve, bl.ends[0].tc[k], bl.ends[1].tc[k], bl.ends[0].v[k], bl.ends[1].v[k], std::max(tol, ct.deviation * (T(1) + T(1e-6))),
                            {{ct.support.get(), ct.pcurve}, {bl.g.surface.get(), ct.on_blend}});
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
         * the new solid (file header). `fillets`: the blends are fillets, three of them at a
         * vertex of three edges meeting on a sphere (`corner`).
         */
        template <std::floating_point T, typename Geometry>
        auto blend(Model<T> &m, SolidId s, std::span<const EdgeId> edges, Geometry &&geometry, const FilletOptions<T> &o, const std::string &op,
                   bool fillets) -> SolidId
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
                        const auto gv = view_at(g, t);
                        const auto tf = third_face(m, sv, gv, v, tol, op);
                        for (std::size_t k = 0; k < 2; ++k)
                            short_end = short_end || !meet_contact(gv, k, tf, t, tol);
                    }
                    if (!short_end)
                        break;
                    // (from an extension of one radius at least: none asked would stay none)
                    auto o2 = o;
                    o2.extension = std::max(o.extension, T(1)) * std::pow(T(2), T(attempt + 1));
                    try
                    {
                        bl.g = geometry(ns, bl.g.edge, copier, o2);
                    }
                    catch (const std::exception &)
                    {
                        // the blend as it was (a support that cannot be continued that far)
                        break;
                    }
                }
            }
            // ---- fillets in pieces (kinks of a radius law): a blend per piece within its edge,
            // a lune of the ball between two pieces (none between the last piece of a loop and its
            // first where the law does not kink there: they share their arc)
            struct Between
            {
                std::size_t before, after;
                std::optional<BlendLune<T>> lune;
                std::optional<BlendCrease<T>> crease;
            };
            std::vector<Between> between;
            {
                std::vector<Blend<T>> parts;
                for (auto &bl : bls)
                {
                    const auto &g = bl.g;
                    if (g.pieces.empty())
                    {
                        parts.push_back(std::move(bl));
                        continue;
                    }
                    const T lo = g.closed ? g.a : std::min(g.t1, g.t2), hi = g.closed ? g.b : std::max(g.t1, g.t2);
                    const T eps = T(1e-9) * (g.b - g.a);
                    auto lune_near = [&](T t) -> std::optional<BlendLune<T>> {
                        for (const auto &l : g.lunes)
                            if (std::abs(l.t - t) <= eps)
                                return l;
                        return std::nullopt;
                    };
                    auto crease_near = [&](T t) -> std::optional<BlendCrease<T>> {
                        for (const auto &c : g.creases)
                            if (std::abs(c.t - t) <= eps)
                                return c;
                        return std::nullopt;
                    };
                    std::optional<std::size_t> prev;
                    const std::size_t first = parts.size();
                    for (std::size_t i = 0; i < g.pieces.size(); ++i)
                    {
                        const auto &pc = g.pieces[i];
                        if (pc.b <= lo + eps || pc.a >= hi - eps)
                            continue; // (a piece past a vertex: the blend continued only)
                        Blend<T> b;
                        b.g = piece_view(g, i);
                        b.part = true;
                        b.kink = g.closed ? std::array<bool, 2>{true, true} : std::array<bool, 2>{pc.a > lo + eps, pc.b < hi - eps};
                        if (!b.kink[0])
                            b.lune_at[0] = lune_near(lo);
                        if (!b.kink[1])
                            b.lune_at[1] = lune_near(hi);
                        if (prev)
                            between.push_back({*prev, parts.size(), lune_near(pc.a), crease_near(pc.a)});
                        prev = parts.size();
                        parts.push_back(std::move(b));
                    }
                    if (g.closed && prev)
                    {
                        auto l = lune_near(g.b);
                        auto c = crease_near(g.b);
                        between.push_back({*prev, first, l ? l : lune_near(g.a), c ? c : crease_near(g.a)});
                    }
                }
                bls = std::move(parts);
            }
            std::set<EdgeId> blended;
            for (const auto &bl : bls)
                blended.insert(bl.g.edge);

            // ---- the ends: lo at a, hi at b
            for (auto &bl : bls)
            {
                const auto &g = bl.g;
                const auto &ed = m.edge(g.edge);
                if (bl.part)
                {
                    const bool along = g.t1 < g.t2;
                    for (std::size_t j = 0; j < 2; ++j)
                    {
                        auto &end = bl.ends[j];
                        if (bl.kink[j])
                        {
                            end.kind = EndKind::Lune;
                            end.t = j == 0 ? g.a : g.b;
                            continue;
                        }
                        end.vertex = (j == 0) == along ? ed.v1 : ed.v2;
                        end.t = j == 0 ? std::min(g.t1, g.t2) : std::max(g.t1, g.t2);
                    }
                    continue;
                }
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
                    const auto sid = add_edge(m, std::shared_ptr<Curve<T, 3>>{std::make_shared<CurveOnSurface<T, 3>>(pa, g.surface)}, T(0), T(1),
                                              bl.ends[0].v[0], bl.ends[0].v[1], tol, {{g.surface.get(), pa}, {g.surface.get(), pb}});
                    bl.ends[0].edge = bl.ends[1].edge = sid;
                    bl.ends[0].on_blend = pa;
                    bl.ends[1].on_blend = pb;
                    continue;
                }
                const bool along = g.t1 < g.t2;
                bl.ends[0].vertex = along ? ed.v1 : ed.v2;
                bl.ends[1].vertex = along ? ed.v2 : ed.v1;
                bl.ends[0].t = std::min(g.t1, g.t2);
                bl.ends[1].t = std::max(g.t1, g.t2);
            }
            // the pieces of a fillet: their lunes, or their common arc
            std::vector<std::pair<ShellId, FaceUse>> corners;
            std::vector<LuneFace> lunes;
            for (const auto &bw : between)
            {
                auto &b1 = bls[bw.before], &b2 = bls[bw.after];
                auto &e1 = b1.ends[1], &e2 = b2.ends[0];
                for (std::size_t k = 0; k < 2; ++k)
                {
                    e1.v[k] = m.add(Vertex<T>{b1.g.contacts[k].curve->value(e1.t), tol});
                    cover(m, e1.v[k], b2.g.contacts[k].curve->value(e2.t), tol);
                    e2.v[k] = e1.v[k];
                }
                if (bw.lune)
                {
                    lunes.push_back(lune_face(m, sv, b1, 1, b2, 0, *bw.lune->sphere, bw.lune->sphere, tol, op));
                    lunes.back().owner = bw.before;
                    continue;
                }
                // a crease (a peak of the law), or the arc common to the last piece of a loop and its first
                const auto pa = bw.crease ? bw.crease->on_pieces[1] : fillet_detail::segment(across(b2.g, e2.t, T(0)), across(b2.g, e2.t, T(1)), T(0), T(1));
                const auto pb = bw.crease ? bw.crease->on_pieces[0] : fillet_detail::segment(across(b1.g, e1.t, T(0)), across(b1.g, e1.t, T(1)), T(0), T(1));
                const auto sid = add_edge(m, std::shared_ptr<Curve<T, 3>>{std::make_shared<CurveOnSurface<T, 3>>(pa, b2.g.surface)}, T(0), T(1), e1.v[0], e1.v[1],
                                          std::max(tol, bw.crease ? bw.crease->deviation * (T(1) + T(1e-6)) : T(0)), {{b2.g.surface.get(), pa}, {b1.g.surface.get(), pb}});
                e1.kind = e2.kind = EndKind::Loop;
                e1.edge = e2.edge = sid;
                e1.on_blend = pb;
                e2.on_blend = pa;
                e1.tc = {e1.t, e1.t};
                e2.tc = {e2.t, e2.t};
            }
            // joints and corners (all the blends at once), then the vanishing and free ends
            for (std::size_t i = 0; i < bls.size(); ++i)
            {
                if (bls[i].g.closed)
                    continue;
                for (std::size_t j = 0; j < 2; ++j)
                {
                    auto &end = bls[i].ends[j];
                    if (end.edge.valid() || end.kind == EndKind::Vanish || end.kind == EndKind::Lune)
                        continue;
                    std::vector<std::pair<std::size_t, std::size_t>> others;
                    for (std::size_t i2 = 0; i2 < bls.size(); ++i2)
                    {
                        const auto &e2 = m.edge(bls[i2].g.edge);
                        if (i2 == i || !(e2.v1 == end.vertex || e2.v2 == end.vertex))
                            continue;
                        if (bls[i2].g.closed)
                            refuse<T>(op, "a blended closed edge at a vertex of another blended edge (not handled)", {bls[i2].g.edge, end.vertex});
                        for (std::size_t j2 = 0; j2 < 2; ++j2)
                            if (bls[i2].ends[j2].vertex == end.vertex)
                                others.emplace_back(i2, j2);
                    }
                    if (others.size() == 2 && fillets && sv.edges_at(end.vertex).size() == 3)
                    {
                        corners.push_back(corner(m, sv, bls, {std::pair{i, j}, others[0], others[1]}, tol, o.approx_tol, op));
                        continue;
                    }
                    if (others.size() > 1)
                        refuse<T>(op, "more than two blended edges at a vertex, or three that are not fillets at a vertex of three edges", {end.vertex});
                    if (others.size() == 1)
                    {
                        const auto [i2, j2] = others.front();
                        joint(m, sv, bls[i], j, bls[i2], j2, tol, o.approx_tol, op, lunes, i);
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
                    if (end.kind == EndKind::Joint || end.kind == EndKind::Corner || end.kind == EndKind::Lune || end.kind == EndKind::Loop)
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
                            end.edge = add_edge(m, std::shared_ptr<Curve<T, 3>>{}, T(0), T(1), end.vertex, end.vertex, m.vertex(end.vertex).tol,
                                                {{g.surface.get(), end.on_blend}});
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
            std::vector<FaceUse> faces_of;
            for (const auto &bl : bls)
            {
                const auto sh = sv.shell_of(bl.g.contacts[0].face);
                const auto fu = blend_face(m, sv, bl, tol);
                m.shell(sh).faces.push_back(fu);
                faces_of.push_back(fu);
            }
            // a lune oriented against the face of its blend along their common arc
            for (const auto &l : lunes)
            {
                const auto &fu = faces_of[l.owner];
                std::optional<Orientation> used;
                for (auto w : m.face(fu.face).wires)
                    for (const auto &ce : m.wire(w).coedges)
                        if (ce.edge == l.arc)
                            used = ce.orient;
                if (!used)
                    refuse<T>(op, "a lune whose arc does not bound its blend", {l.face});
                // along the arc as seen from outside: forward in the face of the blend, backward in the lune
                const bool blend_along = (*used == Orientation::Forward) == (fu.orient == Orientation::Forward);
                const bool lune_forward = (l.in_wire == Orientation::Forward) != blend_along;
                m.shell(l.shell).faces.push_back(FaceUse{l.face, lune_forward ? Orientation::Forward : Orientation::Reversed});
            }
            for (const auto &[sh, fu] : corners)
                m.shell(sh).faces.push_back(fu);

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
            // the new edges at least as loose as the faces they bound (a face sewn or cut looser
            // than tol: a skinned blade after a boolean), their vertices as their edges
            for (const auto &fu : sv.faces())
                for (auto w : m.face(fu.face).wires)
                    for (const auto &ce : m.wire(w).coedges)
                    {
                        auto &e = m.edge(ce.edge);
                        e.tol = std::max(e.tol, m.face(fu.face).tol);
                        for (auto v : {e.v1, e.v2})
                            m.vertex(v).tol = std::max(m.vertex(v).tol, e.tol);
                    }
            const auto report = check(m, ShapeId{ns});
            if (!report.ok())
            {
                const auto &first = report.entries.front();
                refuse<T>(op, std::string("the result is not a valid solid: ") + to_string(first.issue) + " " + first.detail, {first.shape});
            }
            return ns;
        }

        /**
         * The edges given grouped into chains along which a radius law runs: edges meeting two by
         * two at a vertex where no other edge given meets them, tangent there; each chain in order
         * from one end to the other (`ChainEdge::forward`: the edge runs along it), or round a loop
         * (`closed`; a closed edge alone is a loop).
         */
        template <std::floating_point T>
        auto law_chains(const Model<T> &m, std::span<const EdgeId> edges) -> std::vector<std::pair<std::vector<fillet_detail::ChainEdge>, bool>>
        {
            std::map<VertexId, std::vector<EdgeId>> at;
            for (auto e : edges)
            {
                if (!m.alive(e) || m.edge(e).degenerate || !m.edge(e).curve)
                    continue;
                const auto &ed = m.edge(e);
                at[ed.v1].push_back(e);
                if (ed.v2 != ed.v1)
                    at[ed.v2].push_back(e);
            }
            auto dir = [&](EdgeId e, VertexId v) {
                const auto &ed = m.edge(e);
                const auto d = ed.curve->value(ed.v1 == v ? ed.u1 : ed.u2, 1);
                return d / norm(d);
            };
            auto next = [&](EdgeId e, VertexId v) -> std::optional<EdgeId> {
                const auto &l = at[v];
                if (l.size() != 2 || l[0] == l[1])
                    return std::nullopt;
                const EdgeId o = l[0] == e ? l[1] : l[0];
                if (!(std::abs(dir(e, v) * dir(o, v)) > T(1) - T(1e-6)))
                    return std::nullopt;
                return o;
            };
            auto other = [&](EdgeId e, VertexId v) { return m.edge(e).v1 == v ? m.edge(e).v2 : m.edge(e).v1; };
            std::set<EdgeId> done;
            std::vector<std::pair<std::vector<fillet_detail::ChainEdge>, bool>> r;
            for (auto e : edges)
            {
                if (done.contains(e) || !m.alive(e) || m.edge(e).degenerate || !m.edge(e).curve)
                    continue;
                const auto &ed = m.edge(e);
                std::vector<fillet_detail::ChainEdge> chain{{e, true}};
                done.insert(e);
                bool closed = ed.v1 == ed.v2;
                if (!closed)
                {
                    // forward from v2, then backward from v1
                    VertexId v = ed.v2;
                    for (auto cur = e;;)
                    {
                        const auto n = next(cur, v);
                        if (!n)
                            break;
                        if (*n == e)
                        {
                            closed = true;
                            break;
                        }
                        chain.push_back({*n, m.edge(*n).v1 == v});
                        done.insert(*n);
                        v = other(*n, v), cur = *n;
                    }
                    if (!closed)
                    {
                        v = ed.v1;
                        for (auto cur = e;;)
                        {
                            const auto n = next(cur, v);
                            if (!n || done.contains(*n))
                                break;
                            chain.insert(chain.begin(), {*n, m.edge(*n).v2 == v});
                            done.insert(*n);
                            v = other(*n, v), cur = *n;
                        }
                    }
                }
                r.emplace_back(std::move(chain), closed);
            }
            return r;
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
     * solid, `s` unchanged. The edges filleted together meet along chains (tangent at their
     * common vertices) or three at a vertex of three edges (a corner: a sphere); on failure the
     * model is left as it was.
     */
    template <std::floating_point T>
    [[nodiscard]] auto fillet(Model<T> &m, SolidId s, std::span<const EdgeId> edges, std::type_identity_t<T> r, const FilletOptions<T> &o = {})
        -> BuildResult<SolidId>
    {
        const std::string op = "fillet";
        return blend_detail::guarded(m, op, [&] {
            return blend_detail::blend(
                m, s, edges, [&](SolidId ns, EdgeId e, detail::ShapeCopier<T> &, const FilletOptions<T> &oe) { return fillet_detail::fillet(m, ns, e, r, oe); }, o, op, true);
        });
    }

    template <std::floating_point T>
    [[nodiscard]] auto fillet(Model<T> &m, SolidId s, const std::vector<EdgeId> &edges, std::type_identity_t<T> r, const FilletOptions<T> &o = {})
        -> BuildResult<SolidId>
    {
        return fillet(m, s, std::span<const EdgeId>{edges}, r, o);
    }

    /**
     * @brief The solid `s` with its edges `edges` filleted at a radius following `law` (stage 5,
     * `docs/sources/design/brep_palier5_variable_fillets.md`): a new solid, `s` unchanged. The law
     * runs along each chain of tangent edges given (its points on the edges of the chain, by the
     * arc length along it); a chain of the same radius everywhere is a fillet of constant radius.
     * Where the slope of the law jumps (a point of a linear law), the fillet comes in pieces with
     * the lune of the ball between them (also between the fillets of two edges of a chain); the
     * ends, chains, loops and corners of `fillet` (a corner only where the three radii and their
     * arcs agree on the sphere). On failure the model is left as it was.
     */
    template <std::floating_point T>
    [[nodiscard]] auto fillet(Model<T> &m, SolidId s, std::span<const EdgeId> edges, const RadiusLaw<T> &law, const FilletOptions<T> &o = {})
        -> BuildResult<SolidId>
    {
        const std::string op = "fillet";
        return blend_detail::guarded(m, op, [&] {
            for (const auto &p : law.points)
                if (!std::ranges::contains(edges, p.edge))
                    fillet_detail::refuse<T>(BuildErrc::InvalidBounds, op, "a point of the law on an edge not filleted", {p.edge});
            // the law on each edge, from the chains of the edges given (constant: a radius)
            std::map<EdgeId, std::variant<T, fillet_detail::EdgeRadius<T>>> radius;
            for (const auto &[chain, closed] : blend_detail::law_chains(m, edges))
            {
                const auto r = fillet_detail::chain_radius(m, chain, closed, law, o.tol, op);
                for (std::size_t i = 0; i < chain.size(); ++i)
                {
                    if (r)
                        radius.emplace(chain[i].edge, (*r)[i]);
                    else
                        radius.emplace(chain[i].edge, std::ranges::find_if(law.points, [&](const RadiusPoint<T> &p) {
                                                          return std::ranges::any_of(chain, [&](const auto &ce) { return ce.edge == p.edge; });
                                                      })->radius);
                }
            }
            std::map<EdgeId, EdgeId> original; // the copies of the edges
            return blend_detail::blend(
                m, s, edges,
                [&](SolidId ns, EdgeId e, detail::ShapeCopier<T> &c, const FilletOptions<T> &oe) {
                    if (original.empty())
                        for (auto x : edges)
                            original[c.copy(x)] = x;
                    const auto &r = radius.at(original.at(e));
                    if (const auto *v = std::get_if<T>(&r))
                        return fillet_detail::fillet(m, ns, e, *v, oe);
                    return fillet_detail::fillet_law(m, ns, e, std::get<fillet_detail::EdgeRadius<T>>(r), oe);
                },
                o, op, true);
        });
    }

    template <std::floating_point T>
    [[nodiscard]] auto fillet(Model<T> &m, SolidId s, const std::vector<EdgeId> &edges, const RadiusLaw<T> &law, const FilletOptions<T> &o = {})
        -> BuildResult<SolidId>
    {
        return fillet(m, s, std::span<const EdgeId>{edges}, law, o);
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
                o, op, false);
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
                o, op, false);
        });
    }

    template <std::floating_point T>
    [[nodiscard]] auto chamfer(Model<T> &m, SolidId s, const std::vector<EdgeId> &edges, std::type_identity_t<T> d, const FilletOptions<T> &o = {})
        -> BuildResult<SolidId>
    {
        return chamfer(m, s, std::span<const EdgeId>{edges}, d, o);
    }
}
