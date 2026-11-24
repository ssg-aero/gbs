#pragma once
/**
 * @file blend.h
 * @brief Fillets and chamfers of edges of a solid (`fillet`, `chamfer`): the faces of each edge
 * cut along the contacts of its blend, the faces at its ends cut by the blend, the face of the
 * blend between them; a new solid (stage 4, `docs/sources/design/brep_palier4.md`, § 4.2,
 * § 4.5 and § 4.6; notes `brep4_pr06_fillet_topology.md`, `brep4_pr07_fillet_corners.md`;
 * stage 6, `brep_palier6_vertex_blends.md`, notes `brep6_pr01_mitres.md`, `brep6_pr03_vertex_blends.md`;
 * pockets, `brep_swallowed_faces.md`, note `brep_swallowed_faces_pr01.md`).
 *
 * **Pockets** (fillets: by radius, law or chord). Where the faces of two edges filleted meet
 * again past their ends round a few small faces (a pocket: the sliver left on the top of a
 * cylinder cut by a tilted plane, a flat across an edge), the fillet may pass beside them without
 * touching them. The solid is then filleted with each pocket filled first, on a copy: the
 * surfaces of both faces continued, their intersection from one end to the other (the bridge) as
 * an edge in place of the pocket, filleted with the others; the result is kept when every point
 * of the pocket is on the side of it where the bridge is (outside it; inside for a concave
 * fillet) and no face of the pocket meets it: the pocket and its edges, given or not, are
 * swallowed (with the points of a law and the setbacks on them: without effect). Otherwise the
 * solid is filleted as it is; a fillet found to cut through a pocket then, if refused, is
 * refused as an overflow (`through_pockets`).
 *
 * The solid is copied, and the copy edited:
 * 1. **Geometry.** The blend of each edge (`fillet_geometry`, `chamfer_geometry`): its surface
 *    and its contacts, continued past the vertices of the edge.
 * 2. **Ends.** At each vertex of an open edge, the blend ends:
 *    - **free**, at a vertex of three edges where no other edge is blended: the third face F3
 *      cuts it (at a vertex of n > 3 edges, the n - 2 faces round it that do not bound the edge,
 *      in turn: across each, from an edge round the vertex to the next, split where the blend
 *      meets it). The contacts meet F3 on the other two edges (Gauss-Newton), which are split
 *      there; the **end curve** is the intersection of the blend and of F3 between them, marched
 *      across the blend (Gauss-Newton at each point, the points doubled until the interpolated
 *      curve and its pcurves are within tol / 10 of the exact points between them);
 *    - **joined** to the blend of the next edge of a **chain**, a blended edge going on
 *      tangentially at the vertex: both blends end on the same arc across them at the vertex
 *      (their contacts meet there, checked; the arc of the second projected from that of the
 *      first where both surfaces put it on other parameters); a contact passing onto the next
 *      face splits the edge between both faces (that edge must run through it);
 *    - **vanishing**, where the faces become tangent at the vertex (`fillet_geometry` ends the
 *      fillet there): its contacts meet at the vertex, the face of the fillet ends in a point;
 *    - at a **corner**, a vertex of n edges all filleted (`fillet`): the n centre lines meet at
 *      the centre of a sphere of the radius, at r from the n faces; each fillet ends on its arc
 *      through that centre, a great circle of the sphere, and the face of the sphere bounded by
 *      the n arcs fills the corner (convex or concave alike), its poles and its seam clear of the
 *      arcs (`corner_frame`); the contacts of two fillets on their common face meet at the foot of
 *      the centre on it;
 *    - in a **mitre**, two blends at a vertex of three edges, not tangent there (the third edge
 *      sharp): they meet along a curve from where their contacts cross on the face between them;
 *      the narrower ends on it, the wider goes on past it, cut by the face beyond the sharp edge
 *      (`mitre`); three chamfers at a vertex of three edges meet two by two along lines through
 *      their common point (`mitre3`);
 *    - at a **vertex blend**, n fillets at a vertex of n edges whose sphere is not exact (radii
 *      that differ, a radius dropping toward the corner, convex and concave edges) or with
 *      setbacks asked (`FilletOptions::setbacks`): each fillet set back from where the contacts
 *      of its neighbours cross its own, the hole between them filled by one B-spline face tangent
 *      to the fillets and the faces round it (`vertex_blend`, `fill_n_sided`);
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
 * A contact that does not meet the face at a free end within the blend continued (at a free end
 * at a vertex of more than three edges, the faces next to its own round it; in a mitre, the face
 * across the sharp edge: `reaches`) is tried again with the blend continued 2, 4 and 8 times
 * further (from one radius or distance at least, `FilletOptions::extension` being 0), keeping the
 * last geometry that could be computed.
 *
 * Refused (`NonGeneric`): a fillet cutting through a pocket (overflow), two blended edges meeting not tangentially at a vertex of more than
 * three edges, blended edges at a vertex where some are left sharp but in a mitre of two, a mitre
 * of convex and concave edges, a vertex blend whose contacts do not cross or whose hole cannot be
 * filled within the tolerances, a setback shorter than where the contacts cross, so little
 * longer that its joining curves would be too short, or not at a vertex blend, a mitre whose contacts do not cross or whose blends do not meet within the blends
 * continued, a free end at a vertex of more than three edges meeting the edges round it out of
 * turn, a contact that does not meet
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
#include <gbs/bscanalysis.h>
#include <gbs/bscinterp.h>
#include <gbs/curveonsurface.h>
#include <gbs/nsided.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <functional>
#include <initializer_list>
#include <map>
#include <memory>
#include <optional>
#include <ranges>
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
            Corner, ///< fillets of all the edges at a vertex: the sphere of the corner
            Lune,   ///< a kink of the radius law between two pieces of a fillet: the lune of the ball
            Mitre,  ///< blends of edges meeting at a vertex not tangentially: cut along the curves where they meet
            VertexBlend ///< fillets of all the edges at a vertex, no exact sphere: the face of the vertex blend
        };

        /**
         * A piece of an end made of several curves (`End::pieces`): its edge, its pcurve on the
         * blend, its sense (from the side of the contact 0 of the blend toward its contact 1), its
         * ends in that order, and the face it cuts (none: a mitre curve between two blends) with
         * its pcurve there.
         */
        template <std::floating_point T>
        struct EndPiece
        {
            EdgeId edge;
            std::shared_ptr<Curve<T, 2>> on_blend;
            bool forward{true};
            std::array<VertexId, 2> v;
            FaceId face;
            std::shared_ptr<Curve<T, 2>> on_face;
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
            /// an end made of several curves, from the contact 0 to the contact 1 (a mitre, a free
            /// end at a vertex of more than three edges); `edge` its first
            std::vector<EndPiece<T>> pieces;
        };

        template <std::floating_point T>
        struct Blend
        {
            EdgeBlend<T> g;
            std::array<End<T>, 2> ends; ///< at a (lower parameter) and at b
            std::array<EdgeId, 2> contacts;
            /// a piece of a fillet in pieces (`EdgeBlend::pieces`): `g` the piece alone; its ends
            /// at a vertex of its edge, or at a kink of the law (`kink`)
            bool part{};
            std::array<bool, 2> kink{};
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
            v.creases.clear();
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

        /// A curve across a blend where it meets another surface (`cross_curve`).
        template <std::floating_point T>
        struct CrossCurve
        {
            std::shared_ptr<BSCurve<T, 3>> c3;
            std::shared_ptr<BSCurve<T, 2>> on_blend, on_other; ///< on the blend, on the other surface
            T dev{};                                           ///< largest gap found between the points
            T t0{}, t1{};                                      ///< the parameters on the blend at its ends
        };

        /**
         * The curve across the blend g from s0 to s1 (`across`) where it meets the surface s3,
         * from (t0 on g, uv0 on s3) at s0: at each s the point of the blend on s3 (`meet_across`,
         * from the point before), n + 1 points interpolated in cubic over [0, 1], n doubled from 8
         * until the curve and its pcurves are within tol / 10 of the points solved between them
         * (512 at most); the best kept. Refused (`why`) where a point cannot be solved. `last`:
         * the point at s1 given (t on g, uv on s3), not solved (a double point of the intersection,
         * where g and s3 touch).
         */
        template <std::floating_point T>
        auto cross_curve(const EdgeBlend<T> &g, const Surface<T, 3> &s3, T s0, T s1, T t0, const point<T, 2> &uv0, T tol, const std::string &op,
                         const std::string &why, std::vector<ShapeId> involved, const std::optional<std::pair<T, point<T, 2>>> &last = std::nullopt)
            -> CrossCurve<T>
        {
            struct Pt
            {
                T t;
                point<T, 2> uv3;
                point<T, 3> p;
            };
            auto solve = [&](T tau, const Pt &from) -> Pt {
                const T s = s0 + tau * (s1 - s0);
                if (last && tau == T(1))
                {
                    const auto w = across(g, last->first, s);
                    return {last->first, last->second, g.surface->value(w[0], w[1])};
                }
                const auto x = meet_across(g, s3, s, from.t, from.uv3, tol);
                if (!x)
                    refuse<T>(op, why, involved);
                const auto w = across(g, x->first, s);
                return {x->first, x->second, g.surface->value(w[0], w[1])};
            };
            CrossCurve<T> r;
            r.dev = std::numeric_limits<T>::infinity();
            // n + 1 points across, n doubled from 8: the points at the middles, solved for the check,
            // are the new points
            std::size_t n = 8;
            std::vector<Pt> pts{solve(T(0), Pt{t0, uv0, {}})};
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
                    qb.push_back(across(g, pts[i].t, s0 + s[i] * (s1 - s0)));
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
                    d = std::max({d, norm(k3->value(sm) - x.p), norm(g.surface->value(wb[0], wb[1]) - x.p), norm(s3.value(wf[0], wf[1]) - x.p)});
                }
                if (d < r.dev)
                {
                    r.dev = d, r.c3 = k3, r.on_blend = kb, r.on_other = kf;
                    r.t0 = pts.front().t, r.t1 = pts.back().t;
                }
                if (r.dev <= tol / T(10) || 2 * n > 512)
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
            return r;
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
            const auto cc = cross_curve(g, *s3, T(0), T(1), end.tc[0], uv3[0], tol, op, "the blend does not meet the face at an end of its edge", {g.edge, f3});
            if (!(cc.dev <= tol))
                refuse<T>(op, "the end curve of a blend cannot be approximated within the tolerance", {g.edge, f3});
            end.edge = add_edge(m, std::shared_ptr<Curve<T, 3>>{cc.c3}, T(0), T(1), end.v[0], end.v[1], std::max(tol, cc.dev),
                                {{g.surface.get(), cc.on_blend}, {m.face(f3).surface.get(), cc.on_other}});
            end.on_blend = cc.on_blend;
            end.on_f3 = cc.on_other;
            end.forward = true;
        }

        /**
         * The edges and faces round the vertex v in turn, from the edge e0: e[i] between f[i - 1]
         * and f[i] (f[n - 1] before e[0]). None when an edge at v is not shared by two faces, a
         * face has not two edges at v, or they do not close round v.
         */
        template <std::floating_point T>
        auto ring(const Model<T> &m, const SolidView<T> &sv, VertexId v, EdgeId e0, FaceId f0) -> std::optional<std::pair<std::vector<EdgeId>, std::vector<FaceId>>>
        {
            const auto at = sv.edges_at(v);
            std::vector<EdgeId> es{e0};
            std::vector<FaceId> fs;
            FaceId f = f0;
            EdgeId e = e0;
            for (std::size_t i = 0; i < at.size(); ++i)
            {
                // the other edge of f at v
                std::vector<EdgeId> mine;
                for (auto w : m.face(f).wires)
                    for (const auto &ce : m.wire(w).coedges)
                        if (std::ranges::contains(at, ce.edge) && !std::ranges::contains(mine, ce.edge))
                            mine.push_back(ce.edge);
                if (mine.size() != 2 || !std::ranges::contains(mine, e))
                    return std::nullopt;
                fs.push_back(f);
                e = mine[0] == e ? mine[1] : mine[0];
                if (e == e0)
                    break;
                es.push_back(e);
                const auto us = sv.users(e);
                if (us.size() != 2 || us[0] == us[1])
                    return std::nullopt;
                f = us[0] == f ? us[1] : us[0];
            }
            if (e != e0 || es.size() != at.size() || fs.size() != at.size())
                return std::nullopt;
            return std::pair{es, fs};
        }

        /// Where the curves c1 and c2 cross (Gauss-Newton on their parameters from t1, t2).
        template <std::floating_point T>
        auto crossing(const Curve<T, 3> &c1, const Curve<T, 3> &c2, T t1, T t2, T tol) -> std::optional<std::array<T, 2>>
        {
            Vec<T> x(2);
            x << t1, t2;
            const T res = thick_detail::gauss_newton<T>(x, [&](const Vec<T> &y, Vec<T> &f, Mat<T> &J) {
                const auto d = c1.value(y[0]) - c2.value(y[1]);
                const auto d1 = c1.value(y[0], 1), d2 = c2.value(y[1], 1);
                f.resize(3);
                J.resize(3, 2);
                for (int i = 0; i < 3; ++i)
                {
                    f[i] = d[std::size_t(i)];
                    J.row(i) << d1[std::size_t(i)], -d2[std::size_t(i)];
                }
            });
            if (!(res <= tol / T(10)))
                return std::nullopt;
            return std::array<T, 2>{x[0], x[1]};
        }

        /// The point p on the blend g: its (t, s) (`across`), by Gauss-Newton from (t0, s0).
        template <std::floating_point T>
        auto locate(const EdgeBlend<T> &g, const point<T, 3> &p, T t0, T s0, T tol) -> std::optional<std::array<T, 2>>
        {
            Vec<T> x(2);
            x << t0, s0;
            const auto &b = *g.surface;
            const T res = thick_detail::gauss_newton<T>(x, [&](const Vec<T> &y, Vec<T> &f, Mat<T> &J) {
                const auto w = across(g, y[0], y[1]), dw = across_dt(g, y[0], y[1]);
                const auto ds = g.contacts[1].on_blend->value(y[0]) - g.contacts[0].on_blend->value(y[0]);
                const auto q = b.value(w[0], w[1]), bu = b.value(w[0], w[1], 1, 0), bv = b.value(w[0], w[1], 0, 1);
                const auto qt = dw[0] * bu + dw[1] * bv, qs = ds[0] * bu + ds[1] * bv;
                f.resize(3);
                J.resize(3, 2);
                for (int i = 0; i < 3; ++i)
                {
                    f[i] = q[std::size_t(i)] - p[std::size_t(i)];
                    J.row(i) << qt[std::size_t(i)], qs[std::size_t(i)];
                }
            });
            if (!(res <= tol / T(10)))
                return std::nullopt;
            return std::array<T, 2>{x[0], x[1]};
        }

        /// The surface of the face f continued past its bounds by about four widths of the blend g (as `third_face`).
        template <std::floating_point T>
        auto continued(const Model<T> &m, FaceId f, const EdgeBlend<T> &g, T tol) -> std::shared_ptr<Surface<T, 3>>
        {
            const T reach = T(4) * (norm(g.contacts[0].curve->value(g.a) - g.contacts[1].curve->value(g.a)) +
                                    norm(g.contacts[0].curve->value(g.b) - g.contacts[1].curve->value(g.b))) +
                            T(10) * tol;
            return thick_detail::extended_surface(m.face(f).surface, classify_detail::face_uv_box(m, f), reach, tol).surface;
        }

        /// The edge e split at the point p, the vertex vx put there; refused at a vertex of e or off it.
        template <std::floating_point T>
        void split_at(Model<T> &m, const SolidView<T> &sv, EdgeId e, VertexId vx, T tol, const std::string &op, const std::vector<ShapeId> &involved)
        {
            const auto &ed = m.edge(e);
            const auto [u, d] = fillet_detail::edge_point(ed, m.vertex(vx).pnt);
            const T eps = T(1e-9) * (ed.u2 - ed.u1);
            if (!(d <= T(100) * tol))
                refuse<T>(op, "blends at a vertex wider than a face next to it", involved);
            if (!(u > ed.u1 + eps && u < ed.u2 - eps))
                refuse<T>(op, "blends at a vertex meeting an edge next to it at its other vertex (not handled)", involved);
            cover(m, vx, ed.curve->value(u), std::max(tol, ed.tol));
            if (!split_detail::split_at_vertex(m, e, u, vx))
                refuse<T>(op, "an edge that cannot be split where blends meet it", involved);
            sv.invalidate();
        }

        /// The parameter t within the blend g (continued).
        template <std::floating_point T>
        bool within(const EdgeBlend<T> &g, T t)
        {
            return t >= g.a - T(1e-9) * (g.b - g.a) && t <= g.b + T(1e-9) * (g.b - g.a);
        }

        /// Where the contact k of the blend g meets the face f (continued) from its end at v, at t on
        /// g, on the edge e of f at v; nullopt where it does not within g.
        template <std::floating_point T>
        auto meets_face(const Model<T> &m, const SolidView<T> &sv, const EdgeBlend<T> &g, std::size_t k, FaceId f, EdgeId e, VertexId v, T t, T tol)
            -> std::optional<std::pair<T, point<T, 2>>>
        {
            const auto &ed = m.edge(e);
            const auto pc = sv.pcurve(f, e);
            if (!pc)
                return std::nullopt;
            const auto x = meet(*g.contacts[k].curve, *continued(m, f, g, tol), t, pc->value(ed.v1 == v ? ed.u1 : ed.u2), tol);
            if (!x || !within(g, x->first))
                return std::nullopt;
            return x;
        }

        /**
         * Whether the end at v (t on g) of the blend g reaches what cuts it there, before the
         * blends are cut (the blend continued further otherwise, as at a free end): in a mitre with
         * the blend g2 (a vertex of three edges, the edges not tangent), its contact not on their
         * common face meets the face across the sharp edge; at a free end at a vertex of more than
         * three edges, its contacts meet the faces next to theirs round v. True where the
         * configuration is not one of these (refused or handled later).
         */
        template <std::floating_point T>
        bool reaches(const Model<T> &m, const SolidView<T> &sv, const EdgeBlend<T> &g, const EdgeBlend<T> *g2, VertexId v, T t, T tol)
        {
            if (g2)
            {
                std::optional<std::size_t> ko; // the contact of g not on the common face
                for (std::size_t k1 = 0; k1 < 2; ++k1)
                    for (std::size_t k2 = 0; k2 < 2; ++k2)
                        if (g.contacts[k1].face == g2->contacts[k2].face)
                            ko = 1 - k1;
                std::optional<EdgeId> c;
                for (auto e : sv.edges_at(v))
                    if (e != g.edge && e != g2->edge)
                        c = e;
                if (!ko || !c)
                    return true;
                for (auto f : sv.users(*c))
                    if (f != g.contacts[*ko].face)
                        return meets_face(m, sv, g, *ko, f, *c, v, t, tol).has_value();
                return true;
            }
            const auto rg = ring(m, sv, v, g.edge, g.contacts[0].face);
            if (!rg)
                return true;
            const auto &[es, fs] = *rg;
            const std::size_t n = es.size();
            if (n < 4 || fs[n - 1] != g.contacts[1].face)
                return true;
            return meets_face(m, sv, g, 0, fs[1], es[1], v, t, tol) && meets_face(m, sv, g, 1, fs[n - 2], es[n - 1], v, t, tol);
        }

        /**
         * A mitre of two blends at a vertex v of three edges, the third edge c sharp (file header):
         * the blends meet along a curve from P to X, where their contacts cross on the face between
         * them (both touch that face there: a double point of the curve, given rather than solved). The narrower blend S is cut along that curve across its whole width: P is where its
         * other contact meets the wider blend B. B is cut along the same curve, then by the face
         * of that contact (beyond c for B, as at a free end) from P to where its own other contact
         * meets c. When both other contacts meet c at the same point (same widths), P is there and
         * the mitre curve crosses both. The end curve of B on that face goes with the contacts of
         * that face (`face_cuts`). Refused: a convex and a concave blend, no face between the
         * blends, contacts that do not cross, the other contact of the wider (of both, where they
         * meet c at one point) not meeting c within its blend, a mitre curve that cannot be solved.
         */
        template <std::floating_point T>
        void mitre(Model<T> &m, const SolidView<T> &sv, std::vector<Blend<T>> &bls, std::pair<std::size_t, std::size_t> pa, std::pair<std::size_t, std::size_t> pb,
                   T tol, T approx_tol, const std::string &op, std::vector<std::pair<FaceId, Cut<T>>> &face_cuts)
        {
            const VertexId v = bls[pa.first].ends[pa.second].vertex;
            const auto &ga = bls[pa.first].g, &gb = bls[pb.first].g;
            const std::vector<ShapeId> involved{v, ga.edge, gb.edge};
            const T limit = T(10) * tol + (ga.exact ? T(0) : approx_tol) + (gb.exact ? T(0) : approx_tol);
            if (ga.convex != gb.convex)
                refuse<T>(op, "a mitre of convex and concave blended edges (not handled)", involved);
            // the face between them
            std::optional<std::array<std::size_t, 2>> common;
            for (std::size_t k1 = 0; k1 < 2; ++k1)
                for (std::size_t k2 = 0; k2 < 2; ++k2)
                    if (ga.contacts[k1].face == gb.contacts[k2].face)
                        common = std::array{k1, k2};
            if (!common)
                refuse<T>(op, "blended edges meeting at a vertex with no face between them (not handled)", involved);
            // the sharp edge
            std::optional<EdgeId> c;
            for (auto e : sv.edges_at(v))
                if (e != ga.edge && e != gb.edge)
                    c = e;
            const auto ka = (*common)[0], kb = (*common)[1];
            if (!c || !std::ranges::contains(sv.users(*c), ga.contacts[1 - ka].face) || !std::ranges::contains(sv.users(*c), gb.contacts[1 - kb].face))
                refuse<T>(op, "blended edges at a vertex whose third edge does not lie between their other faces", involved);
            // X: the contacts crossing on the face between the blends
            const auto x = crossing(*ga.contacts[ka].curve, *gb.contacts[kb].curve, bls[pa.first].ends[pa.second].t, bls[pb.first].ends[pb.second].t, tol);
            if (!x || !within(ga, (*x)[0]) || !within(gb, (*x)[1]))
                refuse<T>(op, "the contacts of two blends at a vertex do not cross on the face between them within the blends continued", involved);
            // where each other contact meets the face of the other blend (on c)
            struct Side
            {
                std::size_t i, j, k;     // the blend, its end at v, the index of its contact on the common face
                FaceId face;             // the face of its other contact
                std::shared_ptr<Surface<T, 3>> surface; // that face continued
                T tx{};                  // its parameter at X
                std::optional<std::pair<T, point<T, 2>>> on_c; // its other contact meeting the face of the other one (none beyond its blend)
            };
            std::array<Side, 2> sd{Side{pa.first, pa.second, ka, ga.contacts[1 - ka].face, nullptr, (*x)[0], {}},
                                   Side{pb.first, pb.second, kb, gb.contacts[1 - kb].face, nullptr, (*x)[1], {}}};
            for (auto &q : sd)
                q.surface = continued(m, q.face, bls[q.i].g, tol);
            for (std::size_t q = 0; q < 2; ++q)
            {
                const auto &g = bls[sd[q].i].g;
                const auto &other = sd[1 - q];
                const auto pc = sv.pcurve(other.face, *c);
                const auto &ec = m.edge(*c);
                const auto uv = pc->value(ec.v1 == v ? ec.u1 : ec.u2);
                sd[q].on_c = meet(*g.contacts[1 - sd[q].k].curve, *other.surface, bls[sd[q].i].ends[sd[q].j].t, uv, tol);
                if (sd[q].on_c && !within(g, sd[q].on_c->first))
                    sd[q].on_c.reset();
            }
            const auto off_c = [&] { refuse<T>(op, "a contact of a blend at a vertex does not meet the face beyond the sharp edge within the blend continued", involved); };
            // the narrower blend: its other contact meets the wider blend before the sharp edge (where
            // it meets c, needed for the wider only, may lie past its blend: from its end there)
            auto cut_short = [&](std::size_t q) -> std::optional<std::pair<T, point<T, 2>>> {
                const auto &g = bls[sd[q].i].g, &o = bls[sd[1 - q].i].g;
                if (!sd[1 - q].on_c)
                    return std::nullopt;
                const T t0 = sd[q].on_c ? sd[q].on_c->first : bls[sd[q].i].ends[sd[q].j].t;
                const T tc = sd[q].on_c ? sd[q].on_c->first : sd[q].j == 1 ? g.b : g.a, to = sd[1 - q].on_c->first;
                const auto hit = meet(*g.contacts[1 - sd[q].k].curve, *o.surface, t0, o.contacts[1 - sd[1 - q].k].on_blend->value(to), tol);
                if (!hit || !within(g, hit->first))
                    return std::nullopt;
                // on the blend itself, not on the rest of its surface (the back of a cylinder)
                const auto on = locate(o, g.contacts[1 - sd[q].k].curve->value(hit->first), to, T(0.5), tol);
                if (!on || !within(o, (*on)[0]) || !((*on)[1] > -T(1e-9) && (*on)[1] < T(1) + T(1e-9)))
                    return std::nullopt;
                // past the sharp edge from v: farther along the contact than where it meets c
                const T away = sd[q].j == 1 ? T(-1) : T(1);
                const T gap = norm(g.contacts[1 - sd[q].k].curve->value(hit->first) - g.contacts[1 - sd[q].k].curve->value(tc));
                if (!((hit->first - tc) * away > T(0)) || !(gap > limit))
                    return std::nullopt;
                return hit;
            };
            std::size_t sq = 0;
            std::optional<std::pair<T, point<T, 2>>> hit;
            bool same{};
            if ((hit = cut_short(0)))
                sq = 0;
            else if ((hit = cut_short(1)))
                sq = 1;
            else
            {
                same = true;
                if (!sd[0].on_c || !sd[1].on_c)
                    off_c();
                const auto pa3 = ga.contacts[1 - ka].curve->value(sd[0].on_c->first), pb3 = gb.contacts[1 - kb].curve->value(sd[1].on_c->first);
                if (norm(pa3 - pb3) > limit)
                    refuse<T>(op, "two blends at a vertex that neither cut one another nor meet on the sharp edge", involved);
            }
            const auto &S = sd[sq], &B = sd[1 - sq];
            auto &bs = bls[S.i], &bb = bls[B.i];
            const auto &gs = bs.g, &gw = bb.g;
            // the mitre curve across S, from P (its other contact on the wider blend, or on c) to X:
            // both blends touch the common face at X, a double point of the curve, given
            const T tps = same ? S.on_c->first : hit->first;
            const auto uvp = same ? gw.contacts[1 - B.k].on_blend->value(B.on_c->first) : hit->second;
            const auto mc = cross_curve(gs, *gw.surface, T(1 - S.k), T(S.k), tps, uvp, tol, op, "two blends at a vertex do not meet along a curve across the narrower",
                                        involved, std::optional<std::pair<T, point<T, 2>>>{std::pair<T, point<T, 2>>{S.tx, gw.contacts[B.k].on_blend->value(B.tx)}});
            if (!(mc.dev <= tol))
                refuse<T>(op, "the mitre curve of two blends cannot be approximated within the tolerance", involved);
            const auto px = gs.contacts[S.k].curve->value(S.tx);
            const auto vx = m.add(Vertex<T>{px, tol});
            cover(m, vx, gw.contacts[B.k].curve->value(B.tx), tol);
            const auto pp = mc.c3->value(T(0));
            const auto vp = m.add(Vertex<T>{pp, tol});
            cover(m, vp, gs.contacts[1 - S.k].curve->value(mc.t0), tol);
            const auto ec = add_edge(m, std::shared_ptr<Curve<T, 3>>{mc.c3}, T(0), T(1), vp, vx, std::max(tol, mc.dev),
                                     {{gs.surface.get(), mc.on_blend}, {gw.surface.get(), mc.on_other}});
            // S: its end the mitre curve (running from its contact 1 - S.k to S.k)
            {
                auto &end = bs.ends[S.j];
                end.kind = EndKind::Mitre;
                end.tc[S.k] = S.tx;
                end.tc[1 - S.k] = mc.t0;
                end.v[S.k] = vx;
                end.v[1 - S.k] = vp;
                end.edge = ec;
                end.on_blend = mc.on_blend;
                end.forward = S.k == 1;
            }
            auto &end = bb.ends[B.j];
            end.kind = EndKind::Mitre;
            end.tc[B.k] = B.tx;
            end.v[B.k] = vx;
            if (same)
            {
                // both other contacts end on c at P
                cover(m, vp, gw.contacts[1 - B.k].curve->value(B.on_c->first), tol);
                split_at(m, sv, *c, vp, tol, op, involved);
                end.tc[1 - B.k] = B.on_c->first;
                end.v[1 - B.k] = vp;
                end.edge = ec;
                end.on_blend = mc.on_other;
                end.forward = B.k == 1;
                return;
            }
            // B past P: cut by the face of the other contact of S, from P to its own other contact (on c)
            const auto lp = locate(gw, pp, B.on_c->first, T(1 - B.k), tol);
            if (!lp || !within(gw, (*lp)[0]))
                refuse<T>(op, "the mitre curve of two blends does not end on the wider", involved);
            const auto uvf = gs.contacts[1 - S.k].pcurve->value(mc.t0);
            const auto fc = cross_curve(gw, *S.surface, (*lp)[1], T(1 - B.k), (*lp)[0], uvf, tol, op,
                                        "the wider of two blends at a vertex does not meet the face beyond the sharp edge", involved);
            if (!(fc.dev <= tol))
                refuse<T>(op, "the end curve of a blend cannot be approximated within the tolerance", involved);
            const auto pa_ = gw.contacts[1 - B.k].curve->value(fc.t1);
            const auto va = m.add(Vertex<T>{pa_, tol});
            cover(m, va, fc.c3->value(T(1)), tol);
            split_at(m, sv, *c, va, tol, op, involved);
            const auto ef = add_edge(m, std::shared_ptr<Curve<T, 3>>{fc.c3}, T(0), T(1), vp, va, std::max(tol, fc.dev),
                                     {{gw.surface.get(), fc.on_blend}, {m.face(S.face).surface.get(), fc.on_other}});
            end.tc[1 - B.k] = fc.t1;
            end.v[1 - B.k] = va;
            // its pieces from the contact 0: the mitre curve (it runs P -> X) and the end curve on the face (P -> A)
            if (B.k == 0)
                end.pieces = {EndPiece<T>{ec, mc.on_other, false, {vx, vp}, {}, nullptr}, EndPiece<T>{ef, fc.on_blend, true, {vp, va}, S.face, fc.on_other}};
            else
                end.pieces = {EndPiece<T>{ef, fc.on_blend, false, {va, vp}, S.face, fc.on_other}, EndPiece<T>{ec, mc.on_other, true, {vp, vx}, {}, nullptr}};
            end.edge = end.pieces.front().edge;
            end.on_blend = end.pieces.front().on_blend;
            end.forward = end.pieces.front().forward;
            face_cuts.emplace_back(S.face, Cut<T>{ef, fc.on_other, gs.edge});
        }

        /**
         * Three chamfers at a vertex of three edges: a mitre (file header). Each two meet along a
         * line from where their contacts cross on the face between them to the point Q common to
         * the three (Gauss-Newton on the (t, s) of the three), which must lie within each.
         */
        template <std::floating_point T>
        void mitre3(Model<T> &m, const SolidView<T> &sv, std::vector<Blend<T>> &bls, const std::array<std::pair<std::size_t, std::size_t>, 3> &at, T tol,
                    T approx_tol, const std::string &op)
        {
            const VertexId v = bls[at[0].first].ends[at[0].second].vertex;
            std::vector<ShapeId> involved{v};
            T limit = T(10) * tol;
            for (const auto &[i, j] : at)
            {
                involved.push_back(bls[i].g.edge);
                limit += bls[i].g.exact ? T(0) : approx_tol;
            }
            // the pairs and the faces between them: (q1, k1, q2, k2)
            struct Pair
            {
                std::size_t q1, k1, q2, k2;
                T t1{}, t2{};
                VertexId vx;
            };
            std::vector<Pair> pairs;
            for (std::size_t q1 = 0; q1 < 3; ++q1)
                for (std::size_t q2 = q1 + 1; q2 < 3; ++q2)
                    for (std::size_t k1 = 0; k1 < 2; ++k1)
                        for (std::size_t k2 = 0; k2 < 2; ++k2)
                            if (bls[at[q1].first].g.contacts[k1].face == bls[at[q2].first].g.contacts[k2].face)
                                pairs.push_back({q1, k1, q2, k2});
            if (pairs.size() != 3)
                refuse<T>(op, "blends at a corner whose faces are not shared two by two", involved);
            point<T, 3> mid{};
            for (auto &pr : pairs)
            {
                const auto &g1 = bls[at[pr.q1].first].g, &g2 = bls[at[pr.q2].first].g;
                const auto x = crossing(*g1.contacts[pr.k1].curve, *g2.contacts[pr.k2].curve, bls[at[pr.q1].first].ends[at[pr.q1].second].t,
                                        bls[at[pr.q2].first].ends[at[pr.q2].second].t, tol);
                if (!x || !within(g1, (*x)[0]) || !within(g2, (*x)[1]))
                    refuse<T>(op, "the contacts of two blends at a corner do not cross on the face between them within the blends continued", involved);
                pr.t1 = (*x)[0], pr.t2 = (*x)[1];
                const auto p = g1.contacts[pr.k1].curve->value(pr.t1);
                pr.vx = m.add(Vertex<T>{p, tol});
                cover(m, pr.vx, g2.contacts[pr.k2].curve->value(pr.t2), tol);
                mid = mid + p / T(3);
            }
            // Q: on the three blends
            Vec<T> y(6);
            for (std::size_t q = 0; q < 3; ++q)
            {
                const auto &bl = bls[at[q].first];
                const auto l = locate(bl.g, mid, bl.ends[at[q].second].t, T(0.5), std::numeric_limits<T>::infinity());
                y[Eigen::Index(2 * q)] = l ? (*l)[0] : bl.ends[at[q].second].t;
                y[Eigen::Index(2 * q + 1)] = l ? (*l)[1] : T(0.5);
            }
            auto point_on = [&](std::size_t q, T t, T s, point<T, 3> &dt, point<T, 3> &ds) {
                const auto &g = bls[at[q].first].g;
                const auto w = across(g, t, s), dw = across_dt(g, t, s);
                const auto dd = g.contacts[1].on_blend->value(t) - g.contacts[0].on_blend->value(t);
                const auto bu = g.surface->value(w[0], w[1], 1, 0), bv = g.surface->value(w[0], w[1], 0, 1);
                dt = dw[0] * bu + dw[1] * bv;
                ds = dd[0] * bu + dd[1] * bv;
                return g.surface->value(w[0], w[1]);
            };
            const T res = thick_detail::gauss_newton<T>(y, [&](const Vec<T> &z, Vec<T> &f, Mat<T> &J) {
                f.resize(6);
                J = Mat<T>::Zero(6, 6);
                std::array<point<T, 3>, 3> p, dt, ds;
                for (std::size_t q = 0; q < 3; ++q)
                    p[q] = point_on(q, z[Eigen::Index(2 * q)], z[Eigen::Index(2 * q + 1)], dt[q], ds[q]);
                for (std::size_t q = 0; q < 2; ++q)
                    for (int c = 0; c < 3; ++c)
                    {
                        const auto r = Eigen::Index(3 * q) + c;
                        f[r] = p[q][std::size_t(c)] - p[q + 1][std::size_t(c)];
                        J(r, Eigen::Index(2 * q)) = dt[q][std::size_t(c)];
                        J(r, Eigen::Index(2 * q + 1)) = ds[q][std::size_t(c)];
                        J(r, Eigen::Index(2 * q + 2)) = -dt[q + 1][std::size_t(c)];
                        J(r, Eigen::Index(2 * q + 3)) = -ds[q + 1][std::size_t(c)];
                    }
            });
            if (!(res <= tol / T(10)))
                refuse<T>(op, "three blends at a corner that do not meet at a point", involved);
            point<T, 3> pq{}, dt, ds;
            for (std::size_t q = 0; q < 3; ++q)
            {
                const T t = y[Eigen::Index(2 * q)], s = y[Eigen::Index(2 * q + 1)];
                if (!within(bls[at[q].first].g, t) || !(s > T(0) && s < T(1)))
                    refuse<T>(op, "three blends at a corner whose common point is not within each (not handled)", involved);
                pq = pq + point_on(q, t, s, dt, ds) / T(3);
            }
            const auto vq = m.add(Vertex<T>{pq, tol});
            // each mitre line, from X to Q, across the first blend of its pair
            std::vector<EdgeId> lines;
            std::vector<CrossCurve<T>> curves;
            for (const auto &pr : pairs)
            {
                const auto &g1 = bls[at[pr.q1].first].g, &g2 = bls[at[pr.q2].first].g;
                const auto cc = cross_curve(g1, *g2.surface, T(pr.k1), y[Eigen::Index(2 * pr.q1 + 1)], pr.t1, g2.contacts[pr.k2].on_blend->value(pr.t2), tol, op,
                                            "two blends at a corner do not meet along a curve to their common point", involved);
                if (!(cc.dev <= tol) || norm(cc.c3->value(T(1)) - pq) > limit)
                    refuse<T>(op, "the mitre curves of three blends at a corner do not meet at their common point", involved);
                // (the vertices as loose as the ends of the line accepted)
                cover(m, pr.vx, cc.c3->value(T(0)), tol);
                cover(m, vq, cc.c3->value(T(1)), tol);
                lines.push_back(add_edge(m, std::shared_ptr<Curve<T, 3>>{cc.c3}, T(0), T(1), pr.vx, vq, std::max(tol, cc.dev),
                                         {{g1.surface.get(), cc.on_blend}, {g2.surface.get(), cc.on_other}}));
                curves.push_back(cc);
            }
            // each blend: from its contact 0 (X of the pair on that face) to Q, then to its contact 1
            for (std::size_t q = 0; q < 3; ++q)
            {
                auto &bl = bls[at[q].first];
                auto &end = bl.ends[at[q].second];
                end.kind = EndKind::Mitre;
                std::array<EndPiece<T>, 2> ps;
                for (std::size_t pi = 0; pi < 3; ++pi)
                {
                    const auto &pr = pairs[pi];
                    if (pr.q1 != q && pr.q2 != q)
                        continue;
                    const std::size_t k = pr.q1 == q ? pr.k1 : pr.k2;
                    end.tc[k] = pr.q1 == q ? pr.t1 : pr.t2;
                    end.v[k] = pr.vx;
                    const auto on = pr.q1 == q ? curves[pi].on_blend : curves[pi].on_other;
                    // the line runs X -> Q: along from the contact 0, against toward the contact 1
                    ps[k] = k == 0 ? EndPiece<T>{lines[pi], on, true, {pr.vx, vq}, {}, nullptr} : EndPiece<T>{lines[pi], on, false, {vq, pr.vx}, {}, nullptr};
                }
                end.pieces = {ps[0], ps[1]};
                end.edge = ps[0].edge;
                end.on_blend = ps[0].on_blend;
                end.forward = true;
            }
        }

        /**
         * The free end of the blend at its end j, at a vertex v of n > 3 edges: cut by the n - 2
         * faces round v that do not bound its edge, in turn (file header). From the contact on one
         * face of the edge, where it meets the next face; across each face to the next edge round
         * v (where it meets the blend), split there; to the other contact, where it meets the face
         * before it. Refused: a point not found within the blend, or out of turn across it.
         */
        template <std::floating_point T>
        void free_end_ring(Model<T> &m, const SolidView<T> &sv, Blend<T> &bl, std::size_t j, T tol, const std::string &op)
        {
            auto &end = bl.ends[j];
            const auto &g = bl.g;
            const VertexId v = end.vertex;
            const std::vector<ShapeId> involved{g.edge, v};
            const auto rg = ring(m, sv, v, g.edge, g.contacts[0].face);
            if (!rg)
                refuse<T>(op, "a blend ending at a vertex whose edges and faces do not close round it", involved);
            const auto [es, fs] = *rg;
            const std::size_t n = es.size();
            if (fs[n - 1] != g.contacts[1].face)
                refuse<T>(op, "a blend ending at a vertex whose edges and faces do not close round it", involved);
            // the points in turn from the contact 0: on es[1], ..., es[n - 1]; (t, s) on the blend and (u, v) on the faces after and before
            struct Pt
            {
                T t, s;
                point<T, 3> p;
                VertexId vx;
            };
            std::vector<Pt> pts;
            std::vector<std::shared_ptr<Surface<T, 3>>> sf(n);
            for (std::size_t i = 1; i + 1 < n; ++i)
                sf[i] = continued(m, fs[i], g, tol);
            auto uv_at_v = [&](FaceId f, EdgeId e) {
                const auto &ed = m.edge(e);
                return sv.pcurve(f, e)->value(ed.v1 == v ? ed.u1 : ed.u2);
            };
            for (std::size_t k : {std::size_t{0}, std::size_t{1}})
            {
                // the contact k meets the face next to its own past the edge es[1] (k = 0) or es[n - 1]
                const FaceId f = k == 0 ? fs[1] : fs[n - 2];
                const EdgeId e = k == 0 ? es[1] : es[n - 1];
                const auto x = meet(*g.contacts[k].curve, *sf[k == 0 ? 1 : n - 2], end.t, uv_at_v(f, e), tol);
                if (!x || !within(g, x->first))
                    refuse<T>(op, "a contact does not meet the face at an end of its edge within the blend continued", involved);
                const auto p = g.contacts[k].curve->value(x->first);
                if (k == 0)
                    pts.push_back({x->first, T(0), p, {}});
                else
                {
                    // the edges in between first
                    for (std::size_t i = 2; i + 1 < n; ++i)
                    {
                        const auto &ed = m.edge(es[i]);
                        const auto w0 = across(g, end.t, T(0.5));
                        const auto hit = meet(*ed.curve, *g.surface, ed.v1 == v ? ed.u1 : ed.u2, w0, tol);
                        if (!hit)
                            refuse<T>(op, "a blend at a vertex of more than three edges does not meet an edge round it", involved);
                        const auto q = ed.curve->value(hit->first);
                        const auto l = locate(g, q, end.t, T(0.5), tol);
                        if (!l || !within(g, (*l)[0]) || !((*l)[1] > pts.back().s))
                            refuse<T>(op, "a blend at a vertex of more than three edges meets the edges round it out of turn (not handled)", involved);
                        pts.push_back({(*l)[0], (*l)[1], q, {}});
                    }
                    if (!(T(1) > pts.back().s))
                        refuse<T>(op, "a blend at a vertex of more than three edges meets the edges round it out of turn (not handled)", involved);
                    pts.push_back({x->first, T(1), p, {}});
                }
            }
            // (u, v) of each point on the face after it (es[i + 1] bounds fs[i + 1]), before the splits
            std::vector<point<T, 2>> start;
            for (std::size_t i = 0; i + 1 < pts.size(); ++i)
            {
                const auto &ed = m.edge(es[i + 1]);
                start.push_back(sv.pcurve(fs[i + 1], es[i + 1])->value(fillet_detail::edge_point(ed, pts[i].p)[0]));
            }
            for (std::size_t i = 0; i < pts.size(); ++i)
            {
                pts[i].vx = m.add(Vertex<T>{pts[i].p, tol});
                split_at(m, sv, es[i + 1], pts[i].vx, tol, op, involved);
            }
            end.f3 = fs[1];
            end.tc = {pts.front().t, pts.back().t};
            end.v = {pts.front().vx, pts.back().vx};
            end.pieces.clear();
            for (std::size_t i = 0; i + 1 < pts.size(); ++i)
            {
                const FaceId f = fs[i + 1];
                const auto cc = cross_curve(g, *sf[i + 1], pts[i].s, pts[i + 1].s, pts[i].t, start[i], tol, op, "the blend does not meet a face at an end of its edge", involved);
                if (!(cc.dev <= tol))
                    refuse<T>(op, "the end curve of a blend cannot be approximated within the tolerance", involved);
                const auto eid = add_edge(m, std::shared_ptr<Curve<T, 3>>{cc.c3}, T(0), T(1), pts[i].vx, pts[i + 1].vx, std::max(tol, cc.dev),
                                          {{g.surface.get(), cc.on_blend}, {m.face(f).surface.get(), cc.on_other}});
                end.pieces.push_back(EndPiece<T>{eid, cc.on_blend, true, {pts[i].vx, pts[i + 1].vx}, f, cc.on_other});
            }
            end.edge = end.pieces.front().edge;
            end.on_blend = end.pieces.front().on_blend;
            end.forward = true;
        }

        /// Two blends joined at the vertex v: the arc across both there, their contacts meeting at its ends.
        template <std::floating_point T>
        auto lune_face(Model<T> &m, const SolidView<T> &sv, Blend<T> &b1, std::size_t j1, Blend<T> &b2, std::size_t j2,
                       const std::shared_ptr<SphericalSurface<T>> &sphere, T tol, const std::string &op) -> LuneFace;

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
            std::shared_ptr<Curve<T, 2>> pc2 = fillet_detail::segment(across(g2, e2.t, T(k2[0])), across(g2, e2.t, T(k2[1])), T(0), T(1));
            T arc_gap{};
            for (int i = 1; i < 8; ++i)
            {
                const T s = T(i) / T(8);
                const auto w1 = pc1->value(s), w2 = pc2->value(s);
                arc_gap = std::max(arc_gap, norm(g1.surface->value(w1[0], w1[1]) - g2.surface->value(w2[0], w2[1])));
            }
            // the same arc parametrized otherwise across each blend (its rational segments split
            // differently: a bridge through a pocket, `brep_swallowed_faces.md`), each point of
            // the arc of b1 on that of b2 (not two arcs of a ball in two planes: the lune below)
            auto same_arc = [&] {
                const CurveOnSurface<T, 3> a2(pc2, g2.surface);
                for (int i = 1; i < 16; ++i)
                {
                    const T s = T(i) / T(16);
                    const auto w1 = pc1->value(s);
                    if (extrema_curve_point(a2, g1.surface->value(w1[0], w1[1]), s, T(1e-12))[1] > limit)
                        return false;
                }
                return true;
            };
            if (arc_gap > limit && k2[0] != k2[1] && gap[0] <= limit && gap[1] <= limit && same_arc())
            {
                // the pcurve on b2 projected from the arc of b1, as the edge between them is
                Edge<T> arc;
                arc.curve = std::make_shared<CurveOnSurface<T, 3>>(pc1, g1.surface);
                arc.u1 = T(0), arc.u2 = T(1);
                MakeFaceOptions<T> mo;
                mo.pcurve_tol = tol / T(10);
                if (const auto fit = project_pcurve(arc, *g2.surface, surface_closure(*g2.surface, tol), std::optional<point<T, 2>>{}, mo); fit && fit->deviation <= limit)
                    pc2 = fit->pcurve, arc_gap = fit->deviation;
            }
            if (k2[0] == k2[1] || gap[0] > limit || gap[1] > limit)
                refuse<T>(op, "blended edges meeting at a vertex but not tangentially, or blends that do not join there (not handled)", {g1.edge, g2.edge, e1.vertex});
            std::optional<BlendCrease<T>> crease;
            if (arc_gap > limit)
            {
                // the arcs of the blends differ at the vertex: only fillets whose radius follows a
                // law (a kink of the law there, or a radius varying where the curvature of the edges
                // jumps: the centre line has a corner); the balls next to the vertex cover the
                // middle of the lune of its ball between the arcs (the envelopes on either side meet
                // along a crease, both blends continued past the vertex by 1.5 radius for it), or
                // not (the lune is part of the fillet)
                if (!g1.variable && !g2.variable)
                    refuse<T>(op, "blended edges meeting at a vertex but not tangentially, or blends that do not join there (not handled)", {g1.edge, g2.edge, e1.vertex});
                // (the lune of the ball meets arcs of circle only: as within an edge, `fillet_law`)
                if (!g1.circular || !g2.circular)
                    refuse<T>(op, "conic fillets of a chain whose radius law kinks at their vertex (the lune of the ball meets arcs of circle only: not handled)",
                              {g1.edge, g2.edge, e1.vertex});
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
                    constexpr T pi = std::numbers::pi_v<T>;
                    const auto x = fillet_detail::unit((a1 + a2) / T(2) - c), z = fillet_detail::unit(cross(x, p2 - p1));
                    const auto sphere = std::make_shared<SphericalSurface<T>>(r, ax2<T, 3>{c, z, x}, std::array<T, 4>{-pi, pi, -pi / T(2), pi / T(2)});
                    for (std::size_t k = 0; k < 2; ++k)
                    {
                        e1.v[k] = m.add(Vertex<T>{g1.contacts[k].curve->value(e1.t), std::max(tol, gap[k])});
                        cover(m, e1.v[k], g2.contacts[k2[k]].curve->value(e2.t), tol);
                        e2.v[k2[k]] = e1.v[k];
                    }
                    lunes.push_back(lune_face(m, sv, b1, j1, b2, j2, sphere, tol, op));
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
            // (its tolerance: the gap of the two arcs, projected or not, or of the crease)
            const auto q1 = crease ? crease->on_pieces[0] : pc1, q2 = crease ? crease->on_pieces[1] : pc2;
            const auto eid = add_edge(m, std::shared_ptr<Curve<T, 3>>{std::make_shared<CurveOnSurface<T, 3>>(q1, g1.surface)}, T(0), T(1), e1.v[0], e1.v[1],
                                      std::max(tol, (crease ? crease->deviation : arc_gap) * (T(1) + T(1e-6))), {{g1.surface.get(), q1}, {g2.surface.get(), q2}});
            e1.kind = e2.kind = EndKind::Joint;
            e1.edge = e2.edge = eid;
            e1.on_blend = q1;
            e2.on_blend = q2;
            e1.forward = true;
            e2.forward = k2[0] == 0;
        }

        /**
         * The frame at o of the sphere of a corner, its poles and its seam clear of the patch of
         * the corner: the spherical polygon of the feet f (directions from o, in turn round it; its
         * sides the arcs, great circles), `dirs` the directions of points along the arcs. The pole z: of 2048
         * directions spread over the sphere (Fibonacci), the farthest from the arcs with neither
         * it nor its opposite in the patch; the seam (u = +-pi) in the middle of the widest gap
         * the arcs leave round z, x opposite. (x toward the sum of the feet, z across, held for a
         * box but nothing kept the pole out of the patch nor the seam off the arcs.) At a sharp
         * corner the patch nearly fills a hemisphere and every pole is a few degrees from the
         * arcs (4 at a spike of faces 160 degrees apart): their pcurves are interpolated finer
         * there (`corner`). None when no pole is more than 3 degrees from the arcs.
         */
        template <std::floating_point T>
        auto corner_frame(const point<T, 3> &o, const std::vector<point<T, 3>> &f, const std::vector<point<T, 3>> &dirs) -> std::optional<ax2<T, 3>>
        {
            constexpr T pi = std::numbers::pi_v<T>;
            // within the spherical polygon of the feet (in turn round it): on the side of its
            // middle of each side
            point<T, 3> mid{};
            for (const auto &p : f)
                mid = mid + p;
            auto inside = [&](const point<T, 3> &q) {
                for (std::size_t i = 0; i < f.size(); ++i)
                {
                    const auto n = cross(f[i], f[(i + 1) % f.size()]);
                    if (!((q * n) * (mid * n) > T(0)))
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
         * The n fillets at a vertex of n edges (the corner of a box, the top of a pyramid), at the
         * ends j of the blends i of `at`: their centre lines meet at the centre c of a sphere of
         * their radius, at r from the n faces (Gauss-Newton on the parameters of the n lines). Each
         * fillet ends on its arc through c (a great circle of the sphere); the contacts of two
         * fillets on their common face meet at the foot of c on it. The face of the sphere,
         * bounded by the n arcs (`corner_frame`: its poles and its seam clear of them), the
         * solid inside it at a convex corner and outside at a concave one, is returned with the
         * shell it goes in. Refused: centre lines that do not meet within the blends continued
         * (the centre past a vertex by more than `FilletOptions::extension`), contacts that do not
         * meet on a face, a face not shared by two of the fillets, edges convex and concave at the
         * corner.
         */
        template <std::floating_point T>
        auto corner(Model<T> &m, const SolidView<T> &sv, std::vector<Blend<T>> &bls, const std::vector<std::pair<std::size_t, std::size_t>> &at, T tol,
                    T approx_tol, const std::string &op) -> std::pair<ShellId, FaceUse>
        {
            const std::size_t n = at.size();
            const auto nx = Eigen::Index(n);
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
            Vec<T> x(nx);
            for (std::size_t q = 0; q < n; ++q)
                x[Eigen::Index(q)] = bls[at[q].first].ends[at[q].second].t;
            auto spine = [&](std::size_t q, T t) {
                const auto &g = bls[at[q].first].g;
                return g.spine->value(std::clamp(t, g.a, g.b));
            };
            const T res = thick_detail::gauss_newton<T>(x, [&](const Vec<T> &y, Vec<T> &f, Mat<T> &J) {
                f.resize(3 * (nx - 1));
                J = Mat<T>::Zero(3 * (nx - 1), nx);
                for (std::size_t q = 0; q + 1 < n; ++q)
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
            for (std::size_t q = 0; q < n; ++q)
                c = c + spine(q, x[Eigen::Index(q)]) / T(n);
            for (std::size_t q = 0; q < n; ++q)
            {
                const auto &g = bls[at[q].first].g;
                const T t = x[Eigen::Index(q)];
                if (!(res <= limit) || !(t >= g.a && t <= g.b) || norm(spine(q, t) - c) > limit)
                    refuse<T>(op, "fillets at a corner whose centre lines do not meet within the blends continued", involved);
            }
            // the feet of the centre on the faces: where the contacts of two fillets meet
            std::map<FaceId, std::vector<point<T, 3>>> feet;
            T r{};
            for (std::size_t q = 0; q < n; ++q)
            {
                const auto &g = bls[at[q].first].g;
                for (std::size_t k = 0; k < 2; ++k)
                {
                    const auto p = g.contacts[k].curve->value(x[Eigen::Index(q)]);
                    feet[g.contacts[k].face].push_back(p);
                    r += norm(p - c) / T(2 * n);
                }
            }
            if (feet.size() != n || std::ranges::any_of(feet, [](const auto &f) { return f.second.size() != 2; }))
                refuse<T>(op, "fillets at a corner whose faces are not shared two by two", involved);
            // the faces in turn round the vertex (the polygon of the feet on the sphere)
            const auto rg = ring(m, sv, vx, bls[at[0].first].g.edge, bls[at[0].first].g.contacts[0].face);
            if (!rg)
                refuse<T>(op, "fillets at a corner whose edges and faces do not close round it", involved);
            std::map<FaceId, VertexId> foot;
            std::vector<point<T, 3>> fdir; // the directions of the feet from c
            for (auto f : rg->second)
            {
                const auto it = feet.find(f);
                if (it == feet.end())
                    refuse<T>(op, "fillets at a corner whose faces are not shared two by two", involved);
                const auto &ps = it->second;
                const T gap = norm(ps[0] - ps[1]);
                if (gap > limit || std::abs(norm(ps[0] - c) - r) > limit || std::abs(norm(ps[1] - c) - r) > limit)
                    refuse<T>(op, "fillets at a corner whose contacts do not meet on a face", involved);
                foot[f] = m.add(Vertex<T>{(ps[0] + ps[1]) / T(2), std::max(tol, gap)});
                const auto e = (ps[0] + ps[1]) / T(2) - c;
                fdir.push_back(e / norm(e));
            }
            // the arcs: the ends of the fillets on their surfaces, and the directions of their points
            std::vector<std::shared_ptr<CurveOnSurface<T, 3>>> curves(n);
            std::vector<point<T, 3>> dirs;
            for (std::size_t q = 0; q < n; ++q)
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
            for (std::size_t q = 0; q < n; ++q)
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
            // the face of the sphere: the n arcs chained
            std::vector<CoEdge<T>> cs{arcs[0]};
            std::vector<bool> used(n, false);
            used[0] = true;
            for (std::size_t a = 0; a + 1 < n; ++a)
            {
                const auto last = ends_of(m, cs.back())[1];
                bool found = false;
                for (std::size_t q = 0; q < n && !found; ++q)
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

        /// A curve of space given by a function on [0, 1], its derivatives by differences: the cross directions of the ribbons of a vertex blend.
        template <std::floating_point T>
        class FieldCurve : public Curve<T, 3>
        {
        public:
            explicit FieldCurve(std::function<point<T, 3>(T)> f) : f_{std::move(f)} {}
            auto value(T u, size_t d = 0) const -> point<T, 3> override
            {
                if (d == 0)
                    return f_(u);
                const T e = T(1e-6);
                if (d == 1)
                    return (f_(u + e) - f_(u - e)) / (T(2) * e);
                return (f_(u + e) - T(2) * f_(u) + f_(u - e)) / (e * e);
            }
            auto bounds() const -> std::array<T, 2> override { return {T(0), T(1)}; }

        private:
            std::function<point<T, 3>(T)> f_;
        };

        /// The station of the blend g at `distance` along its edge from the vertex v (arc length; the parameters of the blend at the vertices those of the edge, `EdgeBlend::t1`, `t2`).
        template <std::floating_point T>
        auto setback_station(const Model<T> &m, const EdgeBlend<T> &g, VertexId v, T distance) -> std::optional<T>
        {
            const auto &ed = m.edge(g.edge);
            const bool at1 = ed.v1 == v;
            const T u0 = at1 ? ed.u1 : ed.u2, u1 = at1 ? ed.u2 : ed.u1;
            if (!(length(*ed.curve, std::min(u0, u1), std::max(u0, u1)) > distance))
                return std::nullopt;
            T lo = u0, hi = u1;
            for (int i = 0; i < 60; ++i)
            {
                const T mid = (lo + hi) / T(2);
                (length(*ed.curve, std::min(u0, mid), std::max(u0, mid)) < distance ? lo : hi) = mid;
            }
            const T u = (lo + hi) / T(2);
            const T t0 = at1 ? g.t1 : g.t2, tt = at1 ? g.t2 : g.t1;
            return t0 + (u - u0) * (tt - t0) / (u1 - u0);
        }

        /**
         * The sphere the face of a vertex blend leans toward (stage 6, § 4.5; note
         * `brep6_pr03_vertex_blends.md`): the sphere nearest the sides of the hole (least squares on
         * |p|^2 = 2 c . p + k at their points), weighted as far as the sides lie on it: by
         * (1 - g)^2, g the larger of its gap to the points over 5 % of its radius and of the sine of
         * the angle of the ribbons' normals to its own over 0.2 (none beyond either). Fillets of
         * radii nearly equal at a corner make a hole bounded by arcs near its sphere: the face near
         * that sphere, and the same as it in the limit.
         */
        template <std::floating_point T>
        auto guide_sphere(const std::vector<Ribbon<T>> &ribbons) -> std::optional<GuideSphere<T>>
        {
            constexpr T weight = T(1e4), gap_limit = T(0.05), angle_limit = T(0.2);
            std::vector<point<T, 3>> pts, nrm;
            for (const auto &rb : ribbons)
                for (int i = 0; i < 16; ++i)
                {
                    const T u = rb.t0 + (rb.t1 - rb.t0) * T(i) / T(16);
                    pts.push_back(rb.curve->value(u));
                    nrm.push_back(fillet_detail::unit(cross(rb.curve->value(u, 1), rb.cross->value(u))));
                }
            Mat<T> A(Eigen::Index(pts.size()), 4);
            Vec<T> rhs(Eigen::Index(pts.size()));
            for (std::size_t i = 0; i < pts.size(); ++i)
            {
                const auto &p = pts[i];
                A.row(Eigen::Index(i)) << T(2) * p[0], T(2) * p[1], T(2) * p[2], T(1);
                rhs[Eigen::Index(i)] = p * p;
            }
            const Vec<T> x = A.colPivHouseholderQr().solve(rhs);
            const point<T, 3> c{x[0], x[1], x[2]};
            const T r2 = x[3] + c * c;
            if (!(r2 > T(0)) || !std::isfinite(r2))
                return std::nullopt;
            const T r = std::sqrt(r2);
            T g{};
            for (std::size_t i = 0; i < pts.size(); ++i)
            {
                const auto d = pts[i] - c;
                g = std::max({g, std::abs(norm(d) - r) / (gap_limit * r), norm(cross(nrm[i], d / norm(d))) / angle_limit});
            }
            if (!(g < T(1)))
                return std::nullopt;
            return GuideSphere<T>{{c[0], c[1], c[2]}, r, weight * (T(1) - g) * (T(1) - g)};
        }

        /**
         * The vertex blend of the n fillets at a vertex of n edges whose sphere is not exact (stage
         * 6, § 4.2 to § 4.5; note `brep6_pr03_vertex_blends.md`), at the ends j of the blends i of
         * `at`. Round the vertex, the fillets and the faces between them in turn (`ring`); on each
         * face, the contacts of the two fillets cross at X (Gauss-Newton).
         * - **Setbacks.** Each fillet ends on its section at the farther from the vertex of the two
         *   stations where one of its contacts passes a crossing, or farther as asked (`setbacks`,
         *   the station at that distance along the edge).
         * - **Close ends.** On a face where the contacts of two fillets then end within 7e-3 of the
         *   size of the hole of one another: where the contacts cross obliquely (cosine above 0.05)
         *   between fillets whose radii differ by more than that, both fillets are set back by a
         *   quarter of their mean radius more (the sections, curves of curvature of the fillets,
         *   would ask for curvatures that do not agree at their common point); otherwise (radii
         *   nearly equal) both contacts end at their crossing, the ends there slightly oblique: no
         *   side shorter than that, and the hole of arcs on a sphere where the radii agree. Not
         *   next to a setback asked, which is kept: refused where it leaves the contacts apart by
         *   less than that, or by less than 5e-2 of the hole on both faces of its fillet.
         * - **Short joining curves.** The end of a fillet between two joining curves shorter than
         *   a tenth of the hole (one alone is filled; two round an end, the fill folds there or
         *   misses the tangency): the fillet set back further, until both are that long.
         * - **Joining curves.** On each face where the contacts of the two fillets do not end at one
         *   point, a curve from one to the other, tangent to both: a cubic of Hermite in the
         *   parameters of the face, its tangents of the length of its chord.
         * - **Face.** The hole bounded by the ends and the joining curves filled (`fill_n_sided`, up
         *   to 128 spans): ribbons along the ends, the derivative of the fillet along its edge toward
         *   the vertex; along the joining curves, in the tangent plane of the face, normal to the
         *   curve, into the hole; leaning toward the sphere of the sides as far as they lie on one
         *   (`guide_sphere`). The B-spline trimmed by the polygon of its domain; the ends and the
         *   curves keep their exact curves, the gap to it in their tolerance (`settle`).
         *
         * The face goes into `lunes`, to be oriented against the face of the first fillet along its
         * end; the joining curves into `face_cuts`, with the contacts of their faces. Refused: the
         * fillets and faces do not close round the vertex, contacts that do not cross within the
         * blends continued, a setback shorter than the automatic one, past the other end of its
         * edge or so close to the automatic one that its joining curves would be too short, a hole
         * that cannot be filled within the tolerances.
         */
        template <std::floating_point T>
        void vertex_blend(Model<T> &m, const SolidView<T> &sv, std::vector<Blend<T>> &bls, const std::vector<std::pair<std::size_t, std::size_t>> &at,
                          const std::map<std::pair<VertexId, EdgeId>, T> &setbacks, T tol, const FilletOptions<T> &o, const std::string &op,
                          std::vector<std::pair<FaceId, Cut<T>>> &face_cuts, std::vector<LuneFace> &lunes)
        {
            const std::size_t n = at.size();
            const VertexId vx = bls[at[0].first].ends[at[0].second].vertex;
            std::vector<ShapeId> involved{vx};
            T limit = T(10) * tol;
            for (const auto &[i, j] : at)
            {
                involved.push_back(bls[i].g.edge);
                limit += bls[i].g.exact ? T(0) : o.approx_tol;
            }
            // the fillets in turn round the vertex: b[q] on the edge es[q], between the faces fs[q - 1] and fs[q]
            const auto rg = ring(m, sv, vx, bls[at[0].first].g.edge, bls[at[0].first].g.contacts[0].face);
            if (!rg)
                refuse<T>(op, "fillets at a vertex whose edges and faces do not close round it", involved);
            const auto &[es, fs] = *rg;
            std::vector<std::pair<std::size_t, std::size_t>> b;
            for (auto e : es)
            {
                const auto it = std::ranges::find_if(at, [&](const auto &p) { return bls[p.first].g.edge == e; });
                if (it == at.end())
                    refuse<T>(op, "fillets at a vertex whose edges and faces do not close round it", involved);
                b.push_back(*it);
            }
            std::vector<std::size_t> kin(n), kout(n); // the contacts of b[q] on fs[q - 1] and on fs[q]
            std::vector<T> toward(n);                 // the sense of the parameter of b[q] toward the vertex
            for (std::size_t q = 0; q < n; ++q)
            {
                const auto &g = bls[b[q].first].g;
                kout[q] = g.contacts[0].face == fs[q] ? 0 : 1;
                kin[q] = 1 - kout[q];
                if (g.contacts[kout[q]].face != fs[q] || g.contacts[kin[q]].face != fs[(q + n - 1) % n])
                    refuse<T>(op, "fillets at a vertex whose faces are not shared two by two", involved);
                toward[q] = b[q].second == 1 ? T(1) : T(-1);
            }
            auto end_of = [&](std::size_t q) -> End<T> & { return bls[b[q].first].ends[b[q].second]; };
            auto blend_of = [&](std::size_t q) -> const EdgeBlend<T> & { return bls[b[q].first].g; };
            // the crossings: on fs[q], the parameters of b[q] (its contact kout) and of b[q + 1] (its contact kin)
            std::vector<std::array<T, 2>> xs(n);
            for (std::size_t q = 0; q < n; ++q)
            {
                const std::size_t q1 = (q + 1) % n;
                const auto &g1 = blend_of(q), &g2 = blend_of(q1);
                const auto x = crossing(*g1.contacts[kout[q]].curve, *g2.contacts[kin[q1]].curve, end_of(q).t, end_of(q1).t, tol);
                if (!x || !within(g1, (*x)[0]) || !within(g2, (*x)[1]))
                    refuse<T>(op, "the contacts of two fillets at a vertex do not cross on the face between them within the blends continued", involved);
                xs[q] = *x;
            }
            // the setbacks: each fillet ends on its section at the farther crossing from the vertex
            // (or farther as asked); on a face where the contacts of two fillets would then end
            // close to one another (a joining curve shorter than `snap` of the hole: radii nearly
            // equal), both end at their crossing instead, the ends there slightly oblique
            std::vector<T> sec(n);
            for (std::size_t q = 0; q < n; ++q)
            {
                const auto &g = blend_of(q);
                const T ta = xs[(q + n - 1) % n][1], tb = xs[q][0];
                sec[q] = toward[q] > T(0) ? std::min(ta, tb) : std::max(ta, tb);
                if (const auto it = setbacks.find({vx, g.edge}); it != setbacks.end())
                {
                    const auto at_d = setback_station(m, g, vx, it->second);
                    if (!at_d)
                        refuse<T>(op, "a setback longer than its edge", {vx, g.edge});
                    if ((*at_d - sec[q]) * toward[q] > T(1e-9) * (g.b - g.a))
                        refuse<T>(op, "a setback shorter than where the contacts of the fillets cross at the vertex", {vx, g.edge});
                    sec[q] = *at_d;
                }
            }
            // the radius of b[q] at its section (from its centre line to its contacts)
            auto radius_at = [&](std::size_t q, T t) {
                const auto &g = blend_of(q);
                return (norm(g.contacts[0].curve->value(t) - g.spine->value(t)) + norm(g.contacts[1].curve->value(t) - g.spine->value(t))) / T(2);
            };
            // the gap on fs[q] between the ends of the contacts of b[q] and b[q + 1] at their sections
            auto gap_on = [&](std::size_t q) {
                const std::size_t q1 = (q + 1) % n;
                return norm(blend_of(q).contacts[kout[q]].curve->value(sec[q]) - blend_of(q1).contacts[kin[q1]].curve->value(sec[q1]));
            };
            T size{};
            for (std::size_t q = 0; q < n; ++q)
                for (std::size_t k = 0; k < 2; ++k)
                    for (std::size_t q2 = 0; q2 < n; ++q2)
                        size = std::max(size, norm(blend_of(q).contacts[k].curve->value(sec[q]) - blend_of(q2).contacts[0].curve->value(sec[q2])));
            constexpr T snap = T(7e-3), short_side = T(0.1);
            // the setbacks asked are kept (no end at a crossing next to them): refused where the
            // contacts of their fillet and of a neighbour would then end apart but closer than
            // `snap` of the hole (a joining curve that short is not filled), or closer than half
            // `short_side` on both faces of their fillet (two short joining curves round its end)
            std::vector<bool> asked(n);
            for (std::size_t q = 0; q < n; ++q)
                asked[q] = setbacks.contains({vx, blend_of(q).edge});
            auto apart = [&](std::size_t q, T fraction) {
                const T gap = gap_on(q);
                return gap > limit && gap < fraction * size;
            };
            for (std::size_t q = 0; q < n; ++q)
            {
                const std::size_t before = (q + n - 1) % n;
                const bool tiny = apart(before, snap) || apart(q, snap), two = apart(before, short_side / T(2)) && apart(q, short_side / T(2));
                if (asked[q] && (tiny || two))
                    refuse<T>(op,
                              std::format("a setback asked so close to where the contacts of the fillets cross at the vertex that they would be joined "
                                          "by curves too short (under {:.2g} of the hole): lengthen it, or leave it out",
                                          tiny ? snap : short_side / T(2)),
                              {vx, blend_of(q).edge});
            }
            // contacts ending close to one another on a face where the sections of the two fillets
            // would ask for curvatures that do not agree at their common point (contacts crossing
            // obliquely, radii that differ): both fillets set back by a quarter of their radii
            // more, a joining curve between them
            std::vector<T> push(n, T(0));
            for (std::size_t q = 0; q < n; ++q)
            {
                const std::size_t q1 = (q + 1) % n;
                if (!(gap_on(q) < snap * size))
                    continue;
                const auto d1 = fillet_detail::unit(blend_of(q).contacts[kout[q]].curve->value(xs[q][0], 1));
                const auto d2 = fillet_detail::unit(blend_of(q1).contacts[kin[q1]].curve->value(xs[q][1], 1));
                const T r1 = radius_at(q, sec[q]), r2 = radius_at(q1, sec[q1]);
                if (std::abs(d1 * d2) < T(0.05) || std::abs(r1 - r2) < snap * size)
                    continue;
                const T d = (r1 + r2) / T(4);
                push[q] = std::max(push[q], d), push[q1] = std::max(push[q1], d);
            }
            for (std::size_t q = 0; q < n; ++q)
                if (push[q] > T(0))
                {
                    const auto &g = blend_of(q);
                    sec[q] -= toward[q] * push[q] / norm(g.spine->value(sec[q], 1));
                }
            // on fs[q], the contacts end at their crossing: close (radii nearly equal), no setback
            // asked there
            auto at_crossing = [&](std::size_t q) { return gap_on(q) < snap * size && !asked[q] && !asked[(q + 1) % n]; };
            // the end of a fillet between two joining curves shorter than `short_side` of the hole
            // (one alone is filled; two round an end: the fill folds there, or misses the
            // tangency, seen up to 3 % of the hole): the fillet set back further, until both are
            // that long (none set back as asked here: refused above)
            auto short_on = [&](std::size_t q) { return !at_crossing(q) && apart(q, short_side); };
            for (std::size_t p = 0; p < n; ++p)
            {
                const std::size_t before = (p + n - 1) % n;
                if (asked[p] || !short_on(before) || !short_on(p))
                    continue;
                const T s0 = sec[p], speed = norm(blend_of(p).spine->value(s0, 1));
                auto long_enough = [&](T d) {
                    sec[p] = s0 - toward[p] * d / speed;
                    return gap_on(before) >= short_side * size && gap_on(p) >= short_side * size;
                };
                T lo{}, hi = short_side * size;
                for (int i = 0; i < 8 && !long_enough(hi); ++i)
                    lo = hi, hi *= T(2);
                if (!long_enough(hi))
                    refuse<T>(op, "fillets at a vertex whose joining curves stay short however far they are set back", involved);
                for (int i = 0; i < 40; ++i)
                {
                    const T mid = (lo + hi) / T(2);
                    (long_enough(mid) ? hi : lo) = mid;
                }
                long_enough(hi);
            }
            std::vector<std::array<T, 2>> st(n);
            for (std::size_t q = 0; q < n; ++q)
                st[q] = {sec[q], sec[q]};
            // contacts ending close to one another elsewhere (radii nearly equal, no setback
            // asked): both end at their crossing, the ends there slightly oblique
            for (std::size_t q = 0; q < n; ++q)
                if (at_crossing(q))
                    st[q][1] = xs[q][0], st[(q + 1) % n][0] = xs[q][1];
            for (std::size_t q = 0; q < n; ++q)
            {
                // before the other end of the edge
                const auto &g = blend_of(q);
                const T other = m.edge(g.edge).v1 == vx ? g.t2 : g.t1;
                for (const T t : st[q])
                    if (!((t - other) * toward[q] > T(0)))
                        refuse<T>(op, "fillets at a vertex whose setback passes the other end of an edge", {vx, g.edge});
            }
            // the ends: their contacts at one point on a face where they end at a crossing
            for (std::size_t q = 0; q < n; ++q)
            {
                auto &end = end_of(q);
                const auto &g = blend_of(q);
                end.kind = EndKind::VertexBlend;
                end.tc[kin[q]] = st[q][0];
                end.tc[kout[q]] = st[q][1];
                end.t = (st[q][0] + st[q][1]) / T(2);
                for (std::size_t k = 0; k < 2; ++k)
                    end.v[k] = m.add(Vertex<T>{g.contacts[k].curve->value(end.tc[k]), tol});
            }
            std::vector<bool> joined(n); // on fs[q], the contacts end at one point
            for (std::size_t q = 0; q < n; ++q)
            {
                const std::size_t q1 = (q + 1) % n;
                auto &e1 = end_of(q), &e2 = end_of(q1);
                const auto pa = blend_of(q).contacts[kout[q]].curve->value(e1.tc[kout[q]]), pc = blend_of(q1).contacts[kin[q1]].curve->value(e2.tc[kin[q1]]);
                if (norm(pa - pc) <= limit)
                {
                    joined[q] = true;
                    const auto old = e2.v[kin[q1]];
                    e2.v[kin[q1]] = e1.v[kout[q]];
                    cover(m, e1.v[kout[q]], pc, tol);
                    m.erase(old);
                }
            }
            // the end of each fillet: from its contact 0 to its contact 1 in its (u, v) (a section,
            // or oblique from one crossing to the other)
            for (std::size_t q = 0; q < n; ++q)
            {
                auto &end = end_of(q);
                const auto &g = blend_of(q);
                end.on_blend = fillet_detail::segment(g.contacts[0].on_blend->value(end.tc[0]), g.contacts[1].on_blend->value(end.tc[1]), T(0), T(1));
                end.edge = add_edge(m, std::shared_ptr<Curve<T, 3>>{std::make_shared<CurveOnSurface<T, 3>>(end.on_blend, g.surface)}, T(0), T(1), end.v[0],
                                    end.v[1], tol, {{g.surface.get(), end.on_blend}});
                end.forward = true;
            }
            // the sides of the hole in turn: the end of b[q] from fs[q - 1] to fs[q], then the joining curve on fs[q]
            struct SideOf
            {
                EdgeId edge;
                bool along; // the side runs along its edge
            };
            std::vector<Ribbon<T>> ribbons;
            std::vector<SideOf> sides;
            for (std::size_t q = 0; q < n; ++q)
            {
                const std::size_t q1 = (q + 1) % n;
                const auto &g = blend_of(q);
                const auto &end = end_of(q);
                // along the end from its contact kin (t0, s0) to its contact kout (t1, s1); the cross
                // direction the derivative of the fillet along its edge, toward the vertex
                const T t0 = end.tc[kin[q]], t1 = end.tc[kout[q]], s0 = kin[q] == 0 ? T(0) : T(1), s1 = T(1) - s0, sign = toward[q];
                Ribbon<T> arc;
                arc.curve = std::make_shared<CurveOnSurface<T, 3>>(
                    fillet_detail::segment(g.contacts[kin[q]].on_blend->value(t0), g.contacts[kout[q]].on_blend->value(t1), T(0), T(1)), g.surface);
                arc.cross = std::make_shared<FieldCurve<T>>([g, t0, t1, s0, s1, sign](T u) {
                    const T t = t0 + u * (t1 - t0), s = s0 + u * (s1 - s0);
                    const auto w = across(g, t, s), dw = across_dt(g, t, s);
                    return sign * (dw[0] * g.surface->value(w[0], w[1], 1, 0) + dw[1] * g.surface->value(w[0], w[1], 0, 1));
                });
                ribbons.push_back(arc);
                sides.push_back({end.edge, kin[q] == 0});
                if (joined[q])
                    continue;
                // the joining curve on fs[q], from the contact of b[q] to that of b[q + 1]
                const auto &g2 = blend_of(q1);
                const auto &ca = g.contacts[kout[q]], &cc = g2.contacts[kin[q1]];
                const T ta = end.tc[kout[q]], tc = end_of(q1).tc[kin[q1]];
                const auto uva = ca.pcurve->value(ta), uvc = cc.pcurve->value(tc);
                auto da = toward[q] * ca.pcurve->value(ta, 1), dc = -toward[q1] * cc.pcurve->value(tc, 1);
                const T chord = norm(uvc - uva);
                if (!(norm(da) > T(0)) || !(norm(dc) > T(0)) || !(chord > T(0)))
                    refuse<T>(op, "a joining curve of a vertex blend reduced to a point", involved);
                da = da * (chord / norm(da)), dc = dc * (chord / norm(dc));
                const auto hermite = std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{uva, uva + da / T(3), uvc - dc / T(3), uvc},
                                                                    std::vector<T>{T(0), T(0), T(0), T(0), T(1), T(1), T(1), T(1)}, 3);
                const auto support = ca.support;
                const auto c3 = std::make_shared<CurveOnSurface<T, 3>>(hermite, support);
                const auto ej = add_edge(m, std::shared_ptr<Curve<T, 3>>{c3}, T(0), T(1), end.v[kout[q]], end_of(q1).v[kin[q1]], tol,
                                         {{m.face(fs[q]).surface.get(), hermite}});
                face_cuts.emplace_back(fs[q], Cut<T>{ej, hermite, g.edge});
                // its cross direction: in the tangent plane of the face, normal to the curve, into
                // the hole (the side of the normal of the end before it, P' x X, at their corner)
                const auto pe = arc.curve->value(T(1), 1), xe = arc.cross->value(T(1));
                const auto nh = cross(pe, xe);
                const auto n0 = normal_direction(*support, uva[0], uva[1]);
                const T side = n0 * nh > T(0) ? T(1) : T(-1), len = norm(c3->value(T(1)) - c3->value(T(0)));
                Ribbon<T> jr;
                jr.curve = c3;
                jr.cross = std::make_shared<FieldCurve<T>>([c3, hermite, support, side, len](T u) {
                    const auto uv = hermite->value(u);
                    const auto x = cross(normal_direction(*support, uv[0], uv[1]), c3->value(u, 1));
                    return (side * len / norm(x)) * x;
                });
                ribbons.push_back(jr);
                sides.push_back({ej, true});
            }
            // the face: the hole filled, trimmed by the polygon of its domain, leaning toward the
            // sphere of its sides as far as they lie on one
            NSidedOptions<T> no;
            no.approx_tol = o.approx_tol;
            no.angular_tol = o.angular_tol;
            no.sphere = guide_sphere(ribbons);
            // (a hole of six sides or more round a corner set back may need 128 spans to settle)
            no.max_spans = 128;
            NSidedSurface<T> fill;
            try
            {
                fill = fill_n_sided(ribbons, no);
            }
            catch (const std::exception &e)
            {
                refuse<T>(op, std::string("the vertex blend of fillets cannot be filled: ") + e.what(), involved);
            }
            const std::size_t ms = sides.size();
            std::vector<CoEdge<T>> cs;
            for (std::size_t i = 0; i < ms; ++i)
            {
                const point<T, 2> v0{fill.polygon[i][0], fill.polygon[i][1]}, v1{fill.polygon[(i + 1) % ms][0], fill.polygon[(i + 1) % ms][1]};
                const auto pc = sides[i].along ? fillet_detail::segment(v0, v1, T(0), T(1)) : fillet_detail::segment(v1, v0, T(0), T(1));
                cs.push_back(CoEdge<T>{sides[i].edge, sides[i].along ? Orientation::Forward : Orientation::Reversed, pc});
                const std::vector<OnSurface<T>> on{{fill.surface.get(), pc}};
                settle(m, sides[i].edge, std::span<const OnSurface<T>>{on});
            }
            if (uv_signed_area(m, std::span<const CoEdge<T>>{cs}) < T(0))
                cs = reversed(cs);
            const auto arc0 = end_of(0).edge;
            const auto in_wire = std::ranges::find(cs, arc0, &CoEdge<T>::edge)->orient;
            Wire<T> wire;
            wire.coedges = std::move(cs);
            wire.closed = true;
            const auto wid = m.add(std::move(wire));
            const auto fid = m.add(Face<T>{fill.surface, {wid}, tol, false});
            lunes.push_back(LuneFace{sv.shell_of(blend_of(0).contacts[0].face), fid, arc0, in_wire, b[0].first});
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
        auto lune_face(Model<T> &m, const SolidView<T> &sv, Blend<T> &b1, std::size_t j1, Blend<T> &b2, std::size_t j2,
                       const std::shared_ptr<SphericalSurface<T>> &sphere, T tol, const std::string &op) -> LuneFace
        {
            std::array<CoEdge<T>, 2> arcs;
            for (std::size_t q = 0; q < 2; ++q)
            {
                auto &bl = q == 0 ? b1 : b2;
                auto &end = bl.ends[q == 0 ? j1 : j2];
                const auto &g = bl.g;
                end.kind = EndKind::Lune;
                end.tc = {end.t, end.t};
                end.on_blend = fillet_detail::segment(across(g, end.t, T(0)), across(g, end.t, T(1)), T(0), T(1));
                const auto curve = std::make_shared<CurveOnSurface<T, 3>>(end.on_blend, g.surface);
                const auto [on, gap] = fillet_detail::sphere_pcurve(*sphere, *curve, tol);
                end.on_f3 = on;
                end.edge = add_edge(m, std::shared_ptr<Curve<T, 3>>{curve}, T(0), T(1), end.v[0], end.v[1], tol, {{g.surface.get(), end.on_blend}, {sphere.get(), on}});
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
            const auto fid = m.add(Face<T>{sphere, {wid}, tol, false});
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
            // an end of several pieces: from the contact 0 to the contact 1 at the high end, back at the low one
            if (!hi.pieces.empty())
                for (const auto &p : hi.pieces)
                    cs.push_back(CoEdge<T>{p.edge, p.forward ? Orientation::Forward : Orientation::Reversed, p.on_blend});
            else if (hi.edge.valid())
                cs.push_back(CoEdge<T>{hi.edge, hi.forward ? Orientation::Forward : Orientation::Reversed, hi.on_blend});
            cs.push_back(CoEdge<T>{bl.contacts[1], Orientation::Reversed, g.contacts[1].on_blend});
            if (!lo.pieces.empty())
                for (const auto &p : std::views::reverse(lo.pieces))
                    cs.push_back(CoEdge<T>{p.edge, p.forward ? Orientation::Reversed : Orientation::Forward, p.on_blend});
            else if (lo.edge.valid())
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

        /// Whether the edges e1 and e2 at v go on one into the other there (a chain): their
        /// directions leaving v opposite, within 1e-6 of their cosine.
        template <std::floating_point T>
        bool tangent_at(const Model<T> &m, EdgeId e1, EdgeId e2, VertexId v)
        {
            if (e1 == e2)
                return true;
            auto out = [&](EdgeId e) {
                const auto &ed = m.edge(e);
                const auto d = ed.v1 == v ? ed.curve->value(ed.u1, 1) : T(-1) * ed.curve->value(ed.u2, 1);
                return fillet_detail::unit<T>(d);
            };
            return out(e1) * out(e2) < -(T(1) - T(1e-6));
        }

        /// Whether two of the edges at v go on one into the other there (`tangent_at`): a chain
        /// through v among them, whatever their order.
        template <std::floating_point T>
        bool chain_through(const Model<T> &m, std::span<const EdgeId> edges, VertexId v)
        {
            for (std::size_t i = 0; i < edges.size(); ++i)
                for (std::size_t j = i + 1; j < edges.size(); ++j)
                    if (tangent_at(m, edges[i], edges[j], v))
                        return true;
            return false;
        }

        /**
         * The blends of the edges of the solid s (their geometry from `geometry`, on a copy):
         * the new solid (file header). `fillets`: the blends are fillets, the n of them at a
         * vertex of n edges meeting on a sphere (`corner`); chamfers otherwise, three at a vertex
         * of three edges meeting at a point (`mitre3`). Two at a vertex of three edges, not
         * tangent, meet in a mitre either way (`mitre`).
         */
        template <std::floating_point T, typename Geometry>
        auto blend_on_copy(Model<T> &m, detail::ShapeCopier<T> &copier, SolidId s, std::span<const EdgeId> edges, Geometry &&geometry,
                           const FilletOptions<T> &o, const std::string &op, bool fillets) -> SolidId
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
            // the setbacks asked, on the copy: at a vertex of an edge blended
            std::map<std::pair<VertexId, EdgeId>, T> setbacks;
            for (const auto &sb : o.setbacks)
            {
                if (!m.alive(sb.vertex) || !m.alive(sb.edge) || !std::ranges::contains(edges, sb.edge) ||
                    !(m.edge(sb.edge).v1 == sb.vertex || m.edge(sb.edge).v2 == sb.vertex))
                    fillet_detail::refuse<T>(BuildErrc::InvalidId, op, "a setback at a vertex of an edge not blended", {sb.vertex, sb.edge});
                if (!fillets)
                    fillet_detail::refuse<T>(BuildErrc::InvalidBounds, op, "a setback of a chamfer (fillets only)", {sb.vertex, sb.edge});
                setbacks[{copier.copy(sb.vertex), copier.copy(sb.edge)}] = sb.distance;
            }
            std::set<std::pair<VertexId, EdgeId>> used_setbacks;
            // the edges of two blends going on one into the other at v (a chain)
            auto tangent = [&](const Blend<T> &b1, const Blend<T> &b2, VertexId v) { return tangent_at(m, b1.g.edge, b2.g.edge, v); };
            // an end met past the blend continued (an edge oblique to the face at a free end; at a
            // mitre, the face across the sharp edge; at a free end at a vertex of more than three
            // edges, the faces round it): the blend continued further, by 2, 4 and 8 times as much,
            // while its geometry can be
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
                        std::vector<const Blend<T> *> others;
                        for (const auto &b2 : bls)
                            if (const auto &e2 = m.edge(b2.g.edge); b2.g.edge != g.edge && (e2.v1 == v || e2.v2 == v))
                                others.push_back(&b2);
                        const bool vanishing = (t == g.a || t == g.b) && norm(g.contacts[0].curve->value(t) - g.contacts[1].curve->value(t)) <= T(10) * tol;
                        const std::size_t nv = sv.edges_at(v).size();
                        if (vanishing)
                            continue;
                        const auto gv = view_at(g, t);
                        if (others.size() == 1 && nv == 3 && !tangent(bl, *others.front(), v))
                            short_end = short_end || !reaches(m, sv, gv, &others.front()->g, v, t, tol);
                        else if (others.empty() && nv > 3)
                            short_end = short_end || !reaches(m, sv, gv, static_cast<const EdgeBlend<T> *>(nullptr), v, t, tol);
                        else if (others.empty() && nv == 3)
                        {
                            const auto tf = third_face(m, sv, gv, v, tol, op);
                            for (std::size_t k = 0; k < 2; ++k)
                                short_end = short_end || !meet_contact(gv, k, tf, t, tol);
                        }
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
                        // (the cuts lie within the edge: the first piece and the last run on past its vertices)
                        const auto &pc = g.pieces[i];
                        Blend<T> b;
                        b.g = piece_view(g, i);
                        b.part = true;
                        b.kink = g.closed ? std::array<bool, 2>{true, true} : std::array<bool, 2>{pc.a > lo + eps, pc.b < hi - eps};
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
                    lunes.push_back(lune_face(m, sv, b1, 1, b2, 0, bw.lune->sphere, tol, op));
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
            // joints, mitres and corners (all the blends at once), then the vanishing and free ends
            std::vector<std::pair<FaceId, Cut<T>>> face_cuts;
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
                    if (others.empty())
                        continue;
                    const std::size_t nv = sv.edges_at(end.vertex).size();
                    if (others.size() == 1)
                    {
                        const auto [i2, j2] = others.front();
                        if (tangent(bls[i], bls[i2], end.vertex))
                            joint(m, sv, bls[i], j, bls[i2], j2, tol, o.approx_tol, op, lunes, i);
                        else if (nv == 3)
                            mitre(m, sv, bls, {i, j}, others.front(), tol, o.approx_tol, op, face_cuts);
                        else
                            refuse<T>(op, "two blended edges meeting at a vertex of more than three edges, not tangentially (not handled)", {end.vertex});
                        continue;
                    }
                    std::vector<std::pair<std::size_t, std::size_t>> all{{i, j}};
                    all.insert(all.end(), others.begin(), others.end());
                    // (any two of them, whichever blend comes first)
                    std::vector<EdgeId> at_v;
                    for (const auto &[i2, j2] : all)
                        at_v.push_back(bls[i2].g.edge);
                    if (chain_through(m, std::span<const EdgeId>{at_v}, end.vertex))
                        refuse<T>(op, "a chain of tangent edges through a vertex with other blended edges (not handled)", {end.vertex});
                    if (fillets && all.size() == nv)
                    {
                        // the sphere where it is exact and no setback is asked; the vertex blend otherwise
                        const bool asked = std::ranges::any_of(all, [&](const auto &p) { return setbacks.contains({end.vertex, bls[p.first].g.edge}); });
                        bool sphere = false;
                        // (a sphere only where every section is an arc of circle)
                        if (!asked && std::ranges::all_of(all, [&](const auto &p) { return bls[p.first].g.circular; }))
                        {
                            const auto size = detail::arena_size(m);
                            std::vector<End<T>> kept;
                            for (const auto &[i2, j2] : all)
                                kept.push_back(bls[i2].ends[j2]);
                            try
                            {
                                corners.push_back(corner(m, sv, bls, all, tol, o.approx_tol, op));
                                sphere = true;
                            }
                            // (any failure: a refusal, or a curve evaluated beyond its bounds)
                            catch (const std::exception &)
                            {
                                detail::rollback(m, size);
                                for (std::size_t q = 0; q < all.size(); ++q)
                                    bls[all[q].first].ends[all[q].second] = kept[q];
                                sv.invalidate();
                            }
                        }
                        if (!sphere)
                        {
                            vertex_blend(m, sv, bls, all, setbacks, tol, o, op, face_cuts, lunes);
                            for (const auto &p : all)
                                used_setbacks.insert({end.vertex, bls[p.first].g.edge});
                        }
                    }
                    else if (!fillets && all.size() == 3 && nv == 3)
                        mitre3(m, sv, bls, {all[0], all[1], all[2]}, tol, o.approx_tol, op);
                    else
                        refuse<T>(op, "blended edges at a vertex where some edges are not blended (not handled)", {end.vertex});
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
                    if (end.kind == EndKind::Joint || end.kind == EndKind::Corner || end.kind == EndKind::Lune || end.kind == EndKind::Loop ||
                        end.kind == EndKind::Mitre || end.kind == EndKind::VertexBlend)
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
                    if (sv.edges_at(end.vertex).size() > 3)
                        free_end_ring(m, sv, bl, j, tol, op);
                    else
                        free_end(m, sv, bl, j, tol, op);
                }
            }

            for (const auto &[key, d] : setbacks)
                if (!used_setbacks.contains(key))
                    refuse<T>(op, "a setback at a vertex without vertex blend (fillets of all its edges, no exact sphere)", {key.first, key.second});

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
            // the end curves of mitres on faces go on from contacts there
            for (const auto &[f, c] : face_cuts)
                cuts[f].push_back(c);
            for (const auto &[f, cs] : cuts)
                cut_face(m, f, cs, tol, removed, op);
            // the third faces at the free ends: the run through the vertex replaced by the end curve
            for (const auto &bl : bls)
                for (const auto &end : bl.ends)
                    if (end.kind == EndKind::Free && !end.pieces.empty())
                    {
                        // a free end across several faces: each the run through the vertex replaced by its piece
                        for (const auto &p : end.pieces)
                        {
                            std::set<EdgeId> doomed;
                            for (auto w : m.face(p.face).wires)
                                for (const auto &ce : m.wire(w).coedges)
                                    if (const auto &e = m.edge(ce.edge); !e.degenerate && (e.v1 == end.vertex || e.v2 == end.vertex))
                                        doomed.insert(ce.edge);
                            const std::vector<CoEdge<T>> chain{CoEdge<T>{p.edge, p.forward ? Orientation::Forward : Orientation::Reversed, p.on_face}};
                            if (!replace_run(m, p.face, chain, p.v[0], p.v[1], doomed, removed))
                                refuse<T>(op, "the end curve of a blend does not cut the corner of a face at the end", {bl.g.edge, p.face});
                        }
                    }
                    else if (end.kind == EndKind::Free)
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

        /// `blend_on_copy` on a copy of s: a failure names the entities given (the copies of
        /// them are erased with the attempt).
        template <std::floating_point T, typename Geometry>
        auto blend(Model<T> &m, SolidId s, std::span<const EdgeId> edges, Geometry &&geometry, const FilletOptions<T> &o, const std::string &op,
                   bool fillets) -> SolidId
        {
            detail::ShapeCopier<T> copier{m, Placement<T>{}};
            try
            {
                return blend_on_copy(m, copier, s, edges, std::forward<Geometry>(geometry), o, op, fillets);
            }
            catch (const BuildFailure &e)
            {
                throw BuildFailure(detail::to_originals(e.error(), {&copier}));
            }
        }

        // ---- pockets swallowed by a blend (`brep_swallowed_faces.md`)

        /**
         * A pocket (§ 3.1): between the ends v and w of two blended edges e and e2 bounded by the
         * same faces f[0] and f[1], the wires of both faces run from v to w through other faces
         * only, which they enclose. The edges of the runs, the faces enclosed and their edges and
         * vertices are the pocket. (The two ends of one edge are not looked for: e2 != e.)
         */
        template <std::floating_point T>
        struct Pocket
        {
            std::array<FaceId, 2> f;                  ///< the faces of the blended edges
            EdgeId e, e2;                             ///< the blended edges at v and at w
            VertexId v, w;                            ///< where the runs start and end
            std::array<std::set<EdgeId>, 2> runs;     ///< the edges of f[k] from v to w through the pocket
            std::set<FaceId> faces;                   ///< the faces enclosed
            std::set<EdgeId> edges;                   ///< the edges of those faces (the runs included)
            std::vector<point<T, 3>> samples;         ///< points of the pocket, all to be swallowed
        };

        /**
         * The pockets between the edges `edges` of the solid (§ 3.1): at each end v of an open edge
         * e of two faces where no edge given goes on tangentially, the runs of both faces from v,
         * away from e, reach the same end w of an edge given of the same two faces; no edge given
         * at v or w but e, e2 and edges of the pocket; the faces across the runs, and those met
         * from them, bounded by the runs where they touch f[0] or f[1]. Candidates: the same
         * pocket may come from both of its ends, and the faces round the other way (the rest of
         * the solid) form one too (`bridge` tells them apart).
         */
        template <std::floating_point T>
        auto pockets(const Model<T> &m, const SolidView<T> &sv, std::span<const EdgeId> edges) -> std::vector<Pocket<T>>
        {
            auto faces_of = [&](EdgeId x) -> std::optional<std::array<FaceId, 2>> {
                const auto us = sv.users(x);
                if (us.size() != 2 || us[0] == us[1])
                    return std::nullopt;
                return std::array{us[0], us[1]};
            };
            auto same_faces = [&](EdgeId x, const std::array<FaceId, 2> &f) {
                const auto g = faces_of(x);
                return g && ((*g)[0] == f[0] && (*g)[1] == f[1] || (*g)[0] == f[1] && (*g)[1] == f[0]);
            };
            // the run of the face f from v away from e: its edges, and the edge given of the faces
            // f2 met at its end, with that end
            auto run = [&](FaceId f, EdgeId e, VertexId v, const std::array<FaceId, 2> &f2)
                -> std::optional<std::tuple<std::set<EdgeId>, EdgeId, VertexId>> {
                for (auto wid : m.face(f).wires)
                {
                    const auto &cs = m.wire(wid).coedges;
                    const std::size_t n = cs.size();
                    for (std::size_t i = 0; i < n; ++i)
                    {
                        if (cs[i].edge != e)
                            continue;
                        const auto ends = ends_of(m, cs[i]);
                        if (ends[0] == ends[1])
                            return std::nullopt;
                        const bool fwd = ends[1] == v;
                        std::set<EdgeId> r;
                        for (std::size_t k = 1; k < n; ++k)
                        {
                            const auto &ce = cs[fwd ? (i + k) % n : (i + n - k) % n];
                            if (std::ranges::contains(edges, ce.edge) && same_faces(ce.edge, f2))
                            {
                                if (r.empty())
                                    return std::nullopt;
                                return std::tuple{r, ce.edge, ends_of(m, ce)[fwd ? 0 : 1]};
                            }
                            if (!m.edge(ce.edge).degenerate)
                                r.insert(ce.edge);
                        }
                        return std::nullopt;
                    }
                }
                return std::nullopt;
            };
            std::vector<Pocket<T>> r;
            for (auto e : edges)
            {
                const auto f = faces_of(e);
                const auto &ed = m.edge(e);
                if (!f || ed.degenerate || ed.v1 == ed.v2)
                    continue;
                for (auto v : {ed.v1, ed.v2})
                {
                    // (the edges given at v only: `tangent_at` of an edge away from v compares
                    // its end, and two opposite edges anywhere would make a chain through v)
                    std::vector<EdgeId> at_v;
                    for (auto x : edges)
                        if (m.edge(x).v1 == v || m.edge(x).v2 == v)
                            at_v.push_back(x);
                    if (chain_through(m, std::span<const EdgeId>{at_v}, v))
                        continue;
                    const auto r0 = run((*f)[0], e, v, *f), r1 = run((*f)[1], e, v, *f);
                    if (!r0 || !r1 || std::get<1>(*r0) != std::get<1>(*r1) || std::get<2>(*r0) != std::get<2>(*r1))
                        continue;
                    Pocket<T> p{*f, e, std::get<1>(*r0), v, std::get<2>(*r0), {std::get<0>(*r0), std::get<0>(*r1)}, {}, {}, {}};
                    // the faces across the runs, and those met from them
                    bool ok = true;
                    std::vector<FaceId> todo;
                    for (std::size_t k = 0; k < 2 && ok; ++k)
                        for (auto x : p.runs[k])
                        {
                            const auto us = sv.users(x);
                            if (us.size() != 2 || std::ranges::count(us, p.f[k]) != 1 || std::ranges::contains(us, p.f[1 - k]))
                                ok = false;
                            else
                                todo.push_back(us[0] == p.f[k] ? us[1] : us[0]);
                        }
                    while (ok && !todo.empty())
                    {
                        const auto g = todo.back();
                        todo.pop_back();
                        if (!p.faces.insert(g).second)
                            continue;
                        for (auto wid : m.face(g).wires)
                            for (const auto &ce : m.wire(wid).coedges)
                            {
                                p.edges.insert(ce.edge);
                                for (auto u : sv.users(ce.edge))
                                    if (u == p.f[0] || u == p.f[1])
                                        ok = ok && (p.runs[0].contains(ce.edge) || p.runs[1].contains(ce.edge));
                                    else if (u != g && !p.faces.contains(u))
                                        todo.push_back(u);
                            }
                    }
                    if (!ok || p.faces.empty() || std::ranges::any_of(p.faces, [&](FaceId g) { return sv.shell_of(g) != sv.shell_of(p.f[0]); }))
                        continue;
                    // no other edge given at v or w
                    for (auto x : {v, p.w})
                        for (auto y : sv.edges_at(x))
                            if (y != p.e && y != p.e2 && !p.edges.contains(y) && std::ranges::contains(edges, y))
                                ok = false;
                    if (!ok)
                        continue;
                    r.push_back(std::move(p));
                }
            }
            return r;
        }

        /**
         * The points of the pocket p that must all be swallowed (§ 3.3): inside each face, on each
         * edge (at 8 parameters), at each vertex but v and w (those of the blended edges).
         */
        template <std::floating_point T>
        void sample(const Model<T> &m, Pocket<T> &p, T tol)
        {
            for (auto g : p.faces)
                for (const auto &uv : interior_points(m, g, tol))
                    p.samples.push_back(m.face(g).surface->value(uv[0], uv[1]));
            for (auto x : p.edges)
            {
                const auto &ed = m.edge(x);
                if (ed.degenerate || !ed.curve)
                    continue;
                for (int i = 1; i < 8; ++i)
                    p.samples.push_back(ed.curve->value(ed.u1 + (ed.u2 - ed.u1) * T(i) / T(8)));
                for (auto y : {ed.v1, ed.v2})
                    if (y != p.v && y != p.w)
                        p.samples.push_back(m.vertex(y).pnt);
            }
        }

        /// The edge to come in place of the runs of a pocket: on the intersection of the surfaces of its faces continued.
        template <std::floating_point T>
        struct Bridge
        {
            std::shared_ptr<Curve<T, 3>> curve;
            T u1{}, u2{};
            VertexId a, b;                                    ///< at u1 and u2 (v and w in some order)
            std::array<std::shared_ptr<Curve<T, 2>>, 2> pc;   ///< on the surfaces continued, same parameter
            std::array<std::shared_ptr<Surface<T, 3>>, 2> surface; ///< those of f[0] and f[1] continued (same parameters)
            T deviation{};
        };

        /**
         * The bridge of the pocket p (§ 3.2): the surfaces of f[0] and f[1] continued over it (by
         * twice its reach from v and w), the branch of their intersection through v and w, from
         * one to the other on the side of the pocket (of the two arcs of a closed curve, the one
         * whose middle is nearer the middle of the pocket, across the origin of its parameter or
         * not, then the curve and its pcurves on that arc alone; a closed B-spline only within
         * its bounds). None where there is no such branch, or where it runs along an edge of both
         * faces (the faces round the other way, the rest of the solid, are no pocket).
         */
        template <std::floating_point T>
        auto bridge(const Model<T> &m, const SolidView<T> &sv, const Pocket<T> &p, const FilletOptions<T> &o) -> std::optional<Bridge<T>>
        {
            const T tol = o.tol;
            if (p.samples.empty())
                return std::nullopt;
            const auto pv = m.vertex(p.v).pnt, pw = m.vertex(p.w).pnt;
            T reach = norm(pw - pv);
            point<T, 3> mid{};
            for (const auto &q : p.samples)
            {
                reach = std::max({reach, norm(q - pv), norm(q - pw)});
                mid = mid + q;
            }
            mid = (T(1) / T(p.samples.size())) * mid;
            std::array<thick_detail::Extended<T>, 2> ext;
            for (std::size_t k = 0; k < 2; ++k)
                ext[k] = thick_detail::extended_surface(m.face(p.f[k]).surface, classify_detail::face_uv_box(m, p.f[k]), T(2) * reach, tol);
            IntersectOptions<T> io = o.intersect;
            io.tol = std::max(io.tol, tol / T(10));
            const auto res = intersect(*ext[0].surface, ext[0].box, *ext[1].surface, ext[1].box, io);
            const T near = T(10) * std::max({tol, m.vertex(p.v).tol, m.vertex(p.w).tol});
            // the edges of both faces
            std::vector<EdgeId> common;
            for (auto wid : m.face(p.f[0]).wires)
                for (const auto &ce : m.wire(wid).coedges)
                    if (std::ranges::contains(sv.users(ce.edge), p.f[1]) && !m.edge(ce.edge).degenerate)
                        common.push_back(ce.edge);
            for (const auto &c : res.curves)
            {
                const T eps = T(1e-12) * (T(1) + std::abs(c.b - c.a));
                auto proj = [&](const point<T, 3> &q) {
                    T best = c.a, gap = std::numeric_limits<T>::infinity();
                    for (int i = 0; i <= 16; ++i)
                    {
                        const auto x = extrema_curve_point(*c.curve, q, c.a + (c.b - c.a) * T(i) / T(16), eps);
                        if (x[1] < gap)
                            gap = x[1], best = x[0];
                    }
                    return std::pair{best, gap};
                };
                const auto [tv, dv] = proj(pv);
                const auto [tw, dw] = proj(pw);
                if (!(dv <= near && dw <= near))
                    continue;
                const bool elementary = dynamic_cast<const ElementaryCurve<T> *>(c.curve.get()) != nullptr;
                std::vector<std::array<T, 2>> arcs{{std::min(tv, tw), std::max(tv, tw)}};
                if (c.closed && elementary)
                    arcs.push_back({std::max(tv, tw), std::min(tv, tw) + (c.b - c.a)});
                std::ranges::sort(arcs, {}, [&](const std::array<T, 2> &x) { return norm(c.curve->value((x[0] + x[1]) / T(2)) - mid); });
                const auto [u1, u2] = arcs.front();
                for (int i = 1; i < 16; ++i)
                {
                    const auto q = c.curve->value(u1 + (u2 - u1) * T(i) / T(16));
                    for (auto x : common)
                        if (fillet_detail::edge_point(m.edge(x), q)[1] <= near)
                            return std::nullopt;
                }
                Bridge<T> b{c.curve, u1, u2, u1 == tv ? p.v : p.w, u1 == tv ? p.w : p.v, {c.pc1, c.pc2}, {ext[0].surface, ext[1].surface}, c.deviation};
                if (c.closed && elementary)
                {
                    // the arc alone, its pcurves from its start (the arc across the origin of the
                    // curve runs past the bounds of the pcurves of the whole turn, which may also
                    // have crossed the seam of a surface before reaching it)
                    b.curve = intersect_detail::restricted(dynamic_cast<const ElementaryCurve<T> &>(*c.curve), u1, u2);
                    b.deviation = T(0);
                    for (std::size_t k = 0; k < 2; ++k)
                    {
                        const auto *es = dynamic_cast<const ElementarySurface<T> *>(ext[k].surface.get());
                        if (!es)
                            return std::nullopt;
                        const auto fit = intersect_detail::primitive_pcurve(*b.curve, u1, u2, *es, ext[k].box, io.tol);
                        b.pc[k] = fit.pcurve;
                        b.deviation = std::max(b.deviation, fit.deviation);
                    }
                }
                return b;
            }
            return std::nullopt;
        }

        /**
         * The pocket p filled (§ 4) along its bridge b: the surfaces of f[0] and f[1] continued,
         * the edge of b in place of the run of each, the faces of the pocket out of the shell and
         * erased with its edges and vertices. The edge of the bridge, or none where a run cannot
         * be replaced (the model then partly edited: rolled back by the caller).
         */
        template <std::floating_point T>
        auto fill(Model<T> &m, const SolidView<T> &sv, const Pocket<T> &p, const Bridge<T> &b, T tol) -> std::optional<EdgeId>
        {
            for (std::size_t k = 0; k < 2; ++k)
                m.face(p.f[k]).surface = b.surface[k];
            const auto x = add_edge(m, b.curve, b.u1, b.u2, b.a, b.b, std::max(tol, b.deviation),
                                    {OnSurface<T>{b.surface[0].get(), b.pc[0]}, OnSurface<T>{b.surface[1].get(), b.pc[1]}});
            std::set<EdgeId> removed;
            for (std::size_t k = 0; k < 2; ++k)
                if (!replace_run(m, p.f[k], std::vector<CoEdge<T>>{CoEdge<T>{x, Orientation::Forward, b.pc[k]}}, b.a, b.b, p.runs[k], removed))
                    return std::nullopt;
            const auto sh = sv.shell_of(p.f[0]);
            std::erase_if(m.shell(sh).faces, [&](const FaceUse &fu) { return p.faces.contains(fu.face); });
            for (auto g : p.faces)
            {
                for (auto wid : m.face(g).wires)
                    m.erase(wid);
                m.erase(g);
            }
            for (auto y : p.edges)
            {
                for (auto z : {m.edge(y).v1, m.edge(y).v2})
                    if (z != p.v && z != p.w && m.alive(z))
                        m.erase(z);
                m.erase(y);
            }
            return x;
        }

        /**
         * The copy of a solid whose pockets are filled (`through_pockets`), for the blend of its
         * edges: the copies of the edges of the solid, the edges given that a pocket swallows
         * (their copies erased with it), the options on the copy.
         */
        template <std::floating_point T>
        struct Filled
        {
            std::map<EdgeId, EdgeId> copies; ///< of the edges of the solid
            std::set<EdgeId> swallowed;      ///< the edges given in a pocket
            FilletOptions<T> options;        ///< those given, the setbacks on the copies (none on an edge swallowed: without effect, as the edge)
        };

        /**
         * The blend `run` of the edges given, through the pockets they enclose (§ 3, § 4): the
         * pockets with a bridge, found on the solid (a candidate that throws is no pocket); a copy
         * of the solid where each is filled, blended on the edges given but those of a pocket and
         * on the bridges; kept if every point of each pocket is on the side of the result where
         * the bridge is (out of the solid, in for a concave blend), not on it, the filled copy
         * then erased. Otherwise, or without pocket (then without copy), the blend of the solid
         * as it is; a blend found to roll over a pocket (overflow, § 5) refused as such when that
         * one fails too. `run(s, edges, filled)`: the solid, its edges, the filled copy (none for
         * the solid as it is).
         */
        template <std::floating_point T, typename Run>
        auto through_pockets(Model<T> &m, SolidId s, std::span<const EdgeId> edges, const FilletOptions<T> &o, const std::string &op, Run &&run)
            -> SolidId
        {
            if (!m.alive(s) || std::ranges::any_of(edges, [&](EdgeId e) { return !m.alive(e); }))
                return run(s, edges, nullptr);
            const T tol = o.tol;
            // the pockets with a bridge, each once, none sharing a face with another
            std::vector<Pocket<T>> ps;
            std::vector<Bridge<T>> bs;
            try
            {
                const SolidView<T> sv{m, s};
                std::set<FaceId> taken;
                for (auto &p : pockets(m, sv, edges))
                {
                    if (std::ranges::any_of(p.faces, [&](FaceId g) { return taken.contains(g); }))
                        continue;
                    try
                    {
                        sample(m, p, tol);
                        if (auto b = bridge(m, sv, p, o))
                        {
                            taken.insert(p.faces.begin(), p.faces.end());
                            ps.push_back(std::move(p));
                            bs.push_back(std::move(*b));
                        }
                    }
                    catch (const std::exception &)
                    {
                    }
                }
            }
            catch (const std::exception &)
            {
                ps.clear(), bs.clear();
            }
            std::vector<ShapeId> overflow;
            if (!ps.empty())
            {
                const auto size = detail::arena_size(m);
                try
                {
                    detail::ShapeCopier<T> copier{m, Placement<T>{}};
                    const auto ns = std::get<SolidId>(copier.copy(ShapeId{s}));
                    Filled<T> filled;
                    for (auto e : explore<EdgeId>(m, ShapeId{s}))
                        filled.copies[e] = copier.copy(e);
                    const SolidView<T> sv{m, ns};
                    std::set<EdgeId> swallowed; // (on the copy)
                    std::vector<EdgeId> bridges;
                    for (std::size_t i = 0; i < ps.size(); ++i)
                    {
                        // the pocket and its bridge on the copy
                        const auto &p = ps[i];
                        Pocket<T> q{{copier.copy(p.f[0]), copier.copy(p.f[1])}, copier.copy(p.e), copier.copy(p.e2), copier.copy(p.v), copier.copy(p.w), {}, {}, {}, p.samples};
                        for (std::size_t k = 0; k < 2; ++k)
                            for (auto x : p.runs[k])
                                q.runs[k].insert(copier.copy(x));
                        for (auto g : p.faces)
                            q.faces.insert(copier.copy(g));
                        for (auto x : p.edges)
                            q.edges.insert(copier.copy(x));
                        auto b = bs[i];
                        b.a = copier.copy(b.a), b.b = copier.copy(b.b);
                        const auto x = fill(m, sv, q, b, tol);
                        if (!x)
                            throw BuildFailure({BuildErrc::NonGeneric, op + ": a pocket that cannot be filled", {}});
                        sv.invalidate();
                        swallowed.insert(q.edges.begin(), q.edges.end());
                        bridges.push_back(*x);
                    }
                    if (!check(m, ShapeId{ns}).ok())
                        throw BuildFailure({BuildErrc::NonGeneric, op + ": a pocket filled into an invalid solid", {}});
                    std::vector<EdgeId> through;
                    for (auto e : edges)
                    {
                        const auto it = filled.copies.find(e);
                        const auto c = it != filled.copies.end() ? it->second : e;
                        if (swallowed.contains(c))
                            filled.swallowed.insert(e);
                        else
                            through.push_back(c);
                    }
                    through.insert(through.end(), bridges.begin(), bridges.end());
                    // the setbacks at the copies of their vertices and edges (one at a vertex not
                    // of its edge, or on an edge of another solid, left to `blend` to refuse)
                    filled.options = o;
                    filled.options.setbacks.clear();
                    for (auto sb : o.setbacks)
                    {
                        if (filled.swallowed.contains(sb.edge))
                            continue;
                        if (filled.copies.contains(sb.edge) && m.alive(sb.vertex) && (m.edge(sb.edge).v1 == sb.vertex || m.edge(sb.edge).v2 == sb.vertex))
                            sb.vertex = copier.copy(sb.vertex), sb.edge = filled.copies.at(sb.edge);
                        filled.options.setbacks.push_back(sb);
                    }
                    const auto r = run(ns, std::span<const EdgeId>{through}, &filled);
                    // every point of each pocket where the bridge is, no face of the pocket
                    // (of the solid given) meeting the result
                    bool swallows = true;
                    for (std::size_t i = 0; i < ps.size() && swallows; ++i)
                    {
                        const auto &b = m.edge(bridges[i]);
                        const auto side = classify(m, r, b.curve->value((b.u1 + b.u2) / T(2)), tol);
                        swallows = side != Location::On && std::ranges::all_of(ps[i].samples, [&](const point<T, 3> &q) { return classify(m, r, q, tol) == side; });
                        for (auto g : ps[i].faces)
                        {
                            if (!swallows)
                                break;
                            const auto before = detail::arena_size(m);
                            SectionOptions<T> so;
                            so.tol = tol;
                            const auto x = section(m, ShapeId{r}, ShapeId{g}, so);
                            swallows = x && x->complete && x->edges.empty() && x->points.empty() && x->boundary_edges.empty();
                            detail::rollback(m, before);
                        }
                        if (!swallows)
                            for (auto e : edges)
                                if (e == ps[i].e || e == ps[i].e2)
                                    overflow.push_back(e);
                    }
                    if (swallows)
                    {
                        // the filled copy was a step: erased
                        boolean_detail::erase_shape(m, ShapeId{ns});
                        boolean_detail::erase_unused(m, size[1], size[0]);
                        return r;
                    }
                }
                catch (const std::exception &)
                {
                }
                detail::rollback(m, size);
            }
            try
            {
                return run(s, edges, nullptr);
            }
            catch (const BuildFailure &e)
            {
                if (overflow.empty() || e.error().code != BuildErrc::NonGeneric)
                    throw;
                refuse<T>(op, "a blend that rolls over a sharp edge next to it (overflow, not handled)", overflow);
            }
        }

        /**
         * The edges given grouped into chains along which a radius law runs: open edges meeting two
         * by two at a vertex where no other edge given meets them, one going on there the way the
         * other comes in (tangent, not turning back along it); each chain in order
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
            // the edge going on from e at v: the only other one given there, open, leaving v the
            // way e comes in (not turning back along it, `tangent_at`)
            auto next = [&](EdgeId e, VertexId v) -> std::optional<EdgeId> {
                const auto &l = at[v];
                if (l.size() != 2 || l[0] == l[1])
                    return std::nullopt;
                const EdgeId o = l[0] == e ? l[1] : l[0];
                if (m.edge(o).v1 == m.edge(o).v2 || !tangent_at(m, e, o, v))
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

        /**
         * The law on each edge given, from the chains of the edges (`law_chains`, `chain_radius`):
         * its value where it is the same along the whole chain, its function on the edge
         * otherwise. Refused (`InvalidBounds`): a point of the law on an edge not given, the
         * refusals of `chain_radius` (`what`: the name of the value in the messages).
         */
        template <std::floating_point T>
        auto law_on_edges(const Model<T> &m, std::span<const EdgeId> edges, const BlendLaw<T> &law, T tol, const std::string &op, const std::string &what)
            -> std::map<EdgeId, std::variant<T, fillet_detail::EdgeRadius<T>>>
        {
            for (auto e : edges)
                if (m.alive(e) && (m.edge(e).degenerate || !m.edge(e).curve))
                    fillet_detail::refuse<T>(BuildErrc::InvalidId, op, "a degenerate edge (no curve)", {e});
            for (const auto &p : law.points)
                if (!std::ranges::contains(edges, p.edge))
                    fillet_detail::refuse<T>(BuildErrc::InvalidBounds, op, "a point of the law on an edge not " + std::string(op.starts_with("chamfer") ? "chamfered" : "filleted"),
                                             {p.edge});
            std::map<EdgeId, std::variant<T, fillet_detail::EdgeRadius<T>>> r;
            for (const auto &[chain, closed] : law_chains(m, edges))
            {
                const auto R = fillet_detail::chain_radius(m, chain, closed, law, tol, op, what);
                for (std::size_t i = 0; i < chain.size(); ++i)
                {
                    if (R)
                        r.emplace(chain[i].edge, (*R)[i]);
                    else
                        r.emplace(chain[i].edge, std::ranges::find_if(law.points, [&](const BlendPoint<T> &p) {
                                                     return std::ranges::any_of(chain, [&](const auto &ce) { return ce.edge == p.edge; });
                                                 })->value);
                }
            }
            return r;
        }

        /// The distance on an edge from `law_on_edges`.
        template <std::floating_point T>
        auto distance_of(const std::variant<T, fillet_detail::EdgeRadius<T>> &v) -> fillet_detail::EdgeDistance<T>
        {
            if (const auto *d = std::get_if<T>(&v))
                return fillet_detail::edge_distance<T>(*d);
            return fillet_detail::edge_distance(std::get<fillet_detail::EdgeRadius<T>>(v));
        }

        /**
         * The section of the fillet of each edge given (`FilletOptions::section`, stage 7, § 4):
         * the arc of circle, a conic of rho constant or along the chains of tangent edges
         * (`law_on_edges`), or the G2 quartic. Refused (`InvalidBounds`): rho out of ]0, 1[.
         */
        template <std::floating_point T>
        auto sections_on_edges(const Model<T> &m, std::span<const EdgeId> edges, const BlendSection<T> &sec, T tol, const std::string &op)
            -> std::map<EdgeId, fillet_detail::SectionOf<T>>
        {
            std::map<EdgeId, fillet_detail::SectionOf<T>> r;
            auto inside = [](T x) { return x > T(0) && x < T(1); };
            for (auto e : edges)
            {
                fillet_detail::SectionOf<T> s;
                s.g2 = sec.kind == SectionKind::G2;
                if (sec.kind == SectionKind::Conic && !sec.rho_law)
                {
                    if (!inside(sec.rho))
                        fillet_detail::refuse<T>(BuildErrc::InvalidBounds, op, "the rho of a conic section must be within ]0, 1[", {e});
                    s.conic = true;
                    s.rho = [rho = sec.rho](T) { return rho; };
                }
                r.emplace(e, std::move(s));
            }
            if (sec.kind == SectionKind::Conic && sec.rho_law)
            {
                for (const auto &p : sec.rho_law->points)
                    if (!inside(p.value))
                        fillet_detail::refuse<T>(BuildErrc::InvalidBounds, op, "the rho of a conic section must be within ]0, 1[", {p.edge});
                for (const auto &[e, v] : law_on_edges(m, edges, *sec.rho_law, std::min(tol, T(1e-9)), op, "rho"))
                {
                    auto &s = r.at(e);
                    s.conic = true;
                    if (const auto *c = std::get_if<T>(&v))
                        s.rho = [rho = *c](T) { return rho; };
                    else
                    {
                        const auto &R = std::get<fillet_detail::EdgeRadius<T>>(v);
                        s.constant = false;
                        s.kinks = R.kinks;
                        // (continued past the vertices with its slope: kept within ]0, 1[)
                        s.rho = [f = R.r](T t) { return std::clamp(f(t), T(1e-3), T(1) - T(1e-3)); };
                    }
                }
            }
            return r;
        }

        /**
         * The law `law` on the filled copy (`through_pockets`; none: as it is): its points on the
         * copies of their edges, those on an edge swallowed left out (without effect, as the edge).
         */
        template <std::floating_point T>
        auto law_on_copies(BlendLaw<T> law, const Filled<T> *filled) -> BlendLaw<T>
        {
            if (filled)
            {
                std::erase_if(law.points, [&](const BlendPoint<T> &p) { return filled->swallowed.contains(p.edge); });
                for (auto &p : law.points)
                    if (const auto it = filled->copies.find(p.edge); it != filled->copies.end())
                        p.edge = it->second;
            }
            return law;
        }

        /// The section `sec` with its law of rho on the filled copy (`law_on_copies`).
        template <std::floating_point T>
        auto section_on_copies(BlendSection<T> sec, const Filled<T> *filled) -> BlendSection<T>
        {
            if (sec.rho_law)
                sec.rho_law = law_on_copies(*sec.rho_law, filled);
            return sec;
        }

        /// The face of an edge of the solid ns in its shells (the first use), for a chamfer on any face.
        template <std::floating_point T>
        auto any_face(const Model<T> &m, SolidId ns, EdgeId e) -> FaceId
        {
            const auto uses = shell_edge_uses(m, m.solid(ns).outer);
            if (const auto it = uses.find(e); it != uses.end() && !it->second.empty())
                return it->second.front().face;
            for (auto v : m.solid(ns).voids)
                if (const auto u = shell_edge_uses(m, v); u.contains(e) && !u.at(e).empty())
                    return u.at(e).front().face;
            return FaceId{};
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
                return std::unexpected(detail::rolled_back(size, e.error()));
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
     * common vertices), all those at a vertex (a corner: a sphere, or a vertex blend with the
     * setbacks of `o`), or two at a vertex of three edges (a mitre); on failure the model is left
     * as it was.
     */
    template <std::floating_point T>
    [[nodiscard]] auto fillet(Model<T> &m, SolidId s, std::span<const EdgeId> edges, std::type_identity_t<T> r, const FilletOptions<T> &o = {})
        -> BuildResult<SolidId>
    {
        const std::string op = "fillet";
        return blend_detail::guarded(m, op, [&] {
            return blend_detail::through_pockets(m, s, edges, o, op, [&](SolidId s0, std::span<const EdgeId> e0, const blend_detail::Filled<T> *filled) {
                return blend_detail::blend(
                    m, s0, e0,
                    [&, sections = blend_detail::sections_on_edges(m, e0, blend_detail::section_on_copies(o.section, filled), o.tol, op),
                     original = std::map<EdgeId, EdgeId>{}](SolidId ns, EdgeId e, detail::ShapeCopier<T> &c, const FilletOptions<T> &oe) mutable {
                        if (original.empty())
                            for (auto x : e0)
                                original[c.copy(x)] = x;
                        return fillet_detail::fillet(m, ns, e, r, oe, sections.at(original.at(e)));
                    },
                    filled ? filled->options : o, op, true);
            });
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
     * ends, chains, loops and corners of `fillet` (a sphere where the radii and their arcs agree
     * on it, a vertex blend otherwise). On failure the model is left as it was.
     */
    template <std::floating_point T>
    [[nodiscard]] auto fillet(Model<T> &m, SolidId s, std::span<const EdgeId> edges, const RadiusLaw<T> &law, const FilletOptions<T> &o = {})
        -> BuildResult<SolidId>
    {
        const std::string op = "fillet";
        return blend_detail::guarded(m, op, [&] {
            return blend_detail::through_pockets(m, s, edges, o, op, [&](SolidId s0, std::span<const EdgeId> e0, const blend_detail::Filled<T> *filled) {
                // the laws on the copies of their edges, the bridges in the chains
                const auto radius = blend_detail::law_on_edges(m, e0, blend_detail::law_on_copies(law, filled), o.tol, op, "radius");
                const auto sections = blend_detail::sections_on_edges(m, e0, blend_detail::section_on_copies(o.section, filled), o.tol, op);
                std::map<EdgeId, EdgeId> original; // the copies of the edges
                return blend_detail::blend(
                    m, s0, e0,
                    [&](SolidId ns, EdgeId e, detail::ShapeCopier<T> &c, const FilletOptions<T> &oe) {
                        if (original.empty())
                            for (auto x : e0)
                                original[c.copy(x)] = x;
                        const auto &r = radius.at(original.at(e));
                        if (const auto *v = std::get_if<T>(&r))
                            return fillet_detail::fillet(m, ns, e, *v, oe, sections.at(original.at(e)));
                        return fillet_detail::fillet_law(m, ns, e, std::get<fillet_detail::EdgeRadius<T>>(r), oe, T(0), sections.at(original.at(e)));
                    },
                    filled ? filled->options : o, op, true);
            });
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

    /**
     * @brief The solid `s` with its edges `edges` chamfered at distances following laws (stage 7,
     * `docs/sources/design/brep_palier7_blend_sections.md`, § 3): `d1` on the face `f` (every edge
     * a boundary of `f`), `d2` on the other face of each edge, each along the chains of tangent
     * edges given (`BlendLaw`, as the radius of `fillet`). Where the slope of a law jumps, the
     * chamfer comes in ruled pieces with a crease between them. The ends, chains, mitres of
     * `chamfer`; on failure the model is left as it was.
     */
    template <std::floating_point T>
    [[nodiscard]] auto chamfer(Model<T> &m, SolidId s, std::span<const EdgeId> edges, FaceId f, const BlendLaw<T> &d1, const BlendLaw<T> &d2,
                               const FilletOptions<T> &o = {}) -> BuildResult<SolidId>
    {
        const std::string op = "chamfer";
        return blend_detail::guarded(m, op, [&] {
            if (!m.alive(f))
                fillet_detail::refuse<T>(BuildErrc::InvalidId, op, "dead or invalid face", {f});
            const auto l1 = blend_detail::law_on_edges(m, edges, d1, o.tol, op, "distance"), l2 = blend_detail::law_on_edges(m, edges, d2, o.tol, op, "distance");
            std::map<EdgeId, EdgeId> original;
            return blend_detail::blend(
                m, s, edges,
                [&](SolidId ns, EdgeId e, detail::ShapeCopier<T> &c, const FilletOptions<T> &oe) {
                    if (original.empty())
                        for (auto x : edges)
                            original[c.copy(x)] = x;
                    const auto e0 = original.at(e);
                    return fillet_detail::chamfer_law(m, ns, e, c.copy(f), blend_detail::distance_of(l1.at(e0)),
                                                      std::optional{blend_detail::distance_of(l2.at(e0))}, T(0), oe);
                },
                o, op, false);
        });
    }

    template <std::floating_point T>
    [[nodiscard]] auto chamfer(Model<T> &m, SolidId s, const std::vector<EdgeId> &edges, FaceId f, const BlendLaw<T> &d1, const BlendLaw<T> &d2,
                               const FilletOptions<T> &o = {}) -> BuildResult<SolidId>
    {
        return chamfer(m, s, std::span<const EdgeId>{edges}, f, d1, d2, o);
    }

    /**
     * @brief The solid `s` with its edges `edges` chamfered at the same distance on both faces,
     * following the law `d` (stage 7, § 3), as `chamfer` with two laws.
     */
    template <std::floating_point T>
    [[nodiscard]] auto chamfer(Model<T> &m, SolidId s, std::span<const EdgeId> edges, const BlendLaw<T> &d, const FilletOptions<T> &o = {})
        -> BuildResult<SolidId>
    {
        const std::string op = "chamfer";
        return blend_detail::guarded(m, op, [&] {
            const auto l = blend_detail::law_on_edges(m, edges, d, o.tol, op, "distance");
            std::map<EdgeId, EdgeId> original;
            return blend_detail::blend(
                m, s, edges,
                [&](SolidId ns, EdgeId e, detail::ShapeCopier<T> &c, const FilletOptions<T> &oe) {
                    if (original.empty())
                        for (auto x : edges)
                            original[c.copy(x)] = x;
                    const auto D = blend_detail::distance_of(l.at(original.at(e)));
                    return fillet_detail::chamfer_law(m, ns, e, blend_detail::any_face(m, ns, e), D, std::optional{D}, T(0), oe);
                },
                o, op, false);
        });
    }

    template <std::floating_point T>
    [[nodiscard]] auto chamfer(Model<T> &m, SolidId s, const std::vector<EdgeId> &edges, const BlendLaw<T> &d, const FilletOptions<T> &o = {})
        -> BuildResult<SolidId>
    {
        return chamfer(m, s, std::span<const EdgeId>{edges}, d, o);
    }

    /**
     * @brief The solid `s` with its edges `edges` chamfered at the distance `d1` (a law) on the
     * face `f` and at the angle `angle` (radians, within ]0, pi/2[) from it (stage 7, § 3.2): the
     * other contact where the line from the first at that angle, in the plane normal to the
     * edge, meets the other face (between two planes at right angles, d2 = d1 tan(angle)). As
     * `chamfer` otherwise.
     */
    template <std::floating_point T>
    [[nodiscard]] auto chamfer_angle(Model<T> &m, SolidId s, std::span<const EdgeId> edges, FaceId f, const BlendLaw<T> &d1, std::type_identity_t<T> angle,
                                     const FilletOptions<T> &o = {}) -> BuildResult<SolidId>
    {
        const std::string op = "chamfer_angle";
        return blend_detail::guarded(m, op, [&] {
            if (!m.alive(f))
                fillet_detail::refuse<T>(BuildErrc::InvalidId, op, "dead or invalid face", {f});
            const auto l = blend_detail::law_on_edges(m, edges, d1, o.tol, op, "distance");
            std::map<EdgeId, EdgeId> original;
            return blend_detail::blend(
                m, s, edges,
                [&](SolidId ns, EdgeId e, detail::ShapeCopier<T> &c, const FilletOptions<T> &oe) {
                    if (original.empty())
                        for (auto x : edges)
                            original[c.copy(x)] = x;
                    return fillet_detail::chamfer_law(m, ns, e, c.copy(f), blend_detail::distance_of(l.at(original.at(e))), std::nullopt, angle, oe);
                },
                o, op, false);
        });
    }

    /// `chamfer_angle` at the same distance d1 everywhere.
    template <std::floating_point T>
    [[nodiscard]] auto chamfer_angle(Model<T> &m, SolidId s, std::span<const EdgeId> edges, FaceId f, std::type_identity_t<T> d1, std::type_identity_t<T> angle,
                                     const FilletOptions<T> &o = {}) -> BuildResult<SolidId>
    {
        BlendLaw<T> law;
        for (auto e : edges)
            law.points.push_back({e, T(0), d1});
        return chamfer_angle(m, s, edges, f, law, angle, o);
    }

    template <std::floating_point T>
    [[nodiscard]] auto chamfer_angle(Model<T> &m, SolidId s, const std::vector<EdgeId> &edges, FaceId f, std::type_identity_t<T> d1, std::type_identity_t<T> angle,
                                     const FilletOptions<T> &o = {}) -> BuildResult<SolidId>
    {
        return chamfer_angle(m, s, std::span<const EdgeId>{edges}, f, d1, angle, o);
    }

    template <std::floating_point T>
    [[nodiscard]] auto chamfer_angle(Model<T> &m, SolidId s, const std::vector<EdgeId> &edges, FaceId f, const BlendLaw<T> &d1, std::type_identity_t<T> angle,
                                     const FilletOptions<T> &o = {}) -> BuildResult<SolidId>
    {
        return chamfer_angle(m, s, std::span<const EdgeId>{edges}, f, d1, angle, o);
    }

    /**
     * @brief The solid `s` with its edges `edges` filleted at the chord c (stage 7, § 4.3): at
     * every station the contacts of the fillet c apart, its radius following from the angle of
     * the faces (between two planes at a constant angle: a cylinder), the section of `o`. The
     * chains, ends, mitres and vertex blends of `fillet` (a corner of radii that differ: a vertex
     * blend); on failure the model is left as it was.
     */
    template <std::floating_point T>
    [[nodiscard]] auto fillet_chordal(Model<T> &m, SolidId s, std::span<const EdgeId> edges, std::type_identity_t<T> c, const FilletOptions<T> &o = {})
        -> BuildResult<SolidId>
    {
        const std::string op = "fillet_chordal";
        return blend_detail::guarded(m, op, [&] {
            return blend_detail::through_pockets(m, s, edges, o, op, [&](SolidId s0, std::span<const EdgeId> e0, const blend_detail::Filled<T> *filled) {
                const auto sections = blend_detail::sections_on_edges(m, e0, blend_detail::section_on_copies(o.section, filled), o.tol, op);
                std::map<EdgeId, EdgeId> original;
                std::vector<EdgeId> copies;
                return blend_detail::blend(
                    m, s0, e0,
                    [&](SolidId ns, EdgeId e, detail::ShapeCopier<T> &cp, const FilletOptions<T> &oe) {
                        if (original.empty())
                            for (auto x : e0)
                                original[cp.copy(x)] = x, copies.push_back(cp.copy(x));
                        auto radius = fillet_detail::chord_radius(m, ns, e, c, oe);
                        // between two planes: the exact fillet of constant radius
                        if (const auto *r = std::get_if<T>(&radius))
                            return fillet_detail::fillet(m, ns, e, *r, oe, sections.at(original.at(e)));
                        auto &R = std::get<fillet_detail::EdgeRadius<T>>(radius);
                        // a chain going on at a vertex: the blend continued there for the joint
                        const auto &ed = m.edge(e);
                        for (std::size_t k = 0; k < 2; ++k)
                        {
                            const auto v = k == 0 ? ed.v1 : ed.v2;
                            R.joined[k] = std::ranges::any_of(copies, [&](EdgeId x) {
                                return x != e && (m.edge(x).v1 == v || m.edge(x).v2 == v) && blend_detail::tangent_at(m, e, x, v);
                            });
                        }
                        return fillet_detail::fillet_law(m, ns, e, R, oe, T(0), sections.at(original.at(e)));
                    },
                    filled ? filled->options : o, op, true);
            });
        });
    }

    template <std::floating_point T>
    [[nodiscard]] auto fillet_chordal(Model<T> &m, SolidId s, const std::vector<EdgeId> &edges, std::type_identity_t<T> c, const FilletOptions<T> &o = {})
        -> BuildResult<SolidId>
    {
        return fillet_chordal(m, s, std::span<const EdgeId>{edges}, c, o);
    }
}
