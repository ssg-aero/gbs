#pragma once
/**
 * @file boolean.h
 * @brief Boolean operations on solids: fuse, cut, common.
 *
 * A simplified general fuse, on copies of the operands (which are not modified):
 * 1. Section of the copies (`section.h`): the section edges are shared by the faces of both.
 * 2. Split of the faces of both copies along the section (`split.h`).
 * 3. Classification of every face against the other operand, by an interior point
 *    (`classify.h`): inside or outside.
 * 4. Selection: fuse keeps A outside B and B outside A; common keeps A inside B and B inside A;
 *    cut keeps A outside B and B inside A, reversed.
 * 5. Assembly: the kept faces already share their edges; shells are their connected components,
 *    solids their closed shells, cavities the shells of negative volume, given to the smallest
 *    solid around them.
 *
 * Coincident faces (on the same surface within the confusion radius, `brep_confusion_factor` tol,
 * as their vertices; primitive or not) are imprinted on each other first (`coincident.h`); their
 * coincident pieces are kept once or dropped according to the sides of the solids.
 *
 * Contacts: a point where the operands touch splits nothing; a line where faces are tangent is a
 * section edge like the others. Where more than two kept faces meet at an edge (solids touching
 * along it), the faces are paired around the edge, each with its neighbour across its material,
 * and every pair but the first gets its own copy of the edge: the shells stay manifold.
 *
 * An edge of one operand lying on a face of the other is imprinted on it (`imprint_edge`): split
 * where it crosses the boundary of the face, merged with the boundary edges it runs along, a cut
 * of the face inside it. A singular intersection (a tangency along a curve found by marching),
 * or an edge that cannot be given a pcurve on the face it lies on, give `BuildErrc::NonGeneric`.
 * The edges and vertices created by the operation and not used by the result are erased.
 */
#include <gbs-brep/builders.h>
#include <gbs-brep/classify.h>
#include <gbs-brep/coincident.h>
#include <gbs-brep/edit.h>
#include <gbs-brep/explore.h>
#include <gbs-brep/partial.h>
#include <gbs-brep/section.h>
#include <gbs-brep/split.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <optional>
#include <ranges>
#include <set>
#include <tuple>
#include <vector>

namespace gbs::brep
{
    enum class BooleanOp
    {
        Fuse,
        Cut,
        Common
    };

    template <std::floating_point T>
    struct BooleanOptions
    {
        T tol{brep_default_tolerance<T>};
        IntersectOptions<T> intersect{};
    };

    /// Operand a face of the result comes from.
    enum class FaceOrigin
    {
        A,
        B
    };

    template <std::floating_point T>
    struct BooleanResult
    {
        ShapeId shape;                                   ///< a solid, or a compound of solids (empty compound for an empty result)
        std::vector<SolidId> solids;
        std::vector<std::pair<FaceId, FaceOrigin>> faces; ///< faces of the result and their operand
    };

    namespace boolean_detail
    {
        template <std::floating_point T>
        auto fail_non_generic(const std::string &why) -> std::unexpected<BuildError>
        {
            return detail::fail(BuildErrc::NonGeneric, "boolean: " + why);
        }

        /// Erases the faces, wires, shells, solids and compounds reachable from a shape.
        template <std::floating_point T>
        void erase_shape(Model<T> &m, const ShapeId &s, const std::set<std::uint32_t> &keep_faces = {})
        {
            // (listed before anything is erased: the exploration goes through the containers)
            const auto compounds = explore<CompoundId>(m, s);
            const auto solids = explore<SolidId>(m, s);
            const auto shells = explore<ShellId>(m, s);
            const auto faces = explore<FaceId>(m, s);
            for (auto c : compounds)
                m.erase(c);
            for (auto so : solids)
                m.erase(so);
            for (auto sh : shells)
                m.erase(sh);
            for (auto f : faces)
                if (!keep_faces.contains(f.index) && m.alive(f))
                {
                    for (auto w : m.face(f).wires)
                        m.erase(w);
                    m.erase(f);
                }
        }

        /**
         * The piece `fu` lies on a face of `other` (one of `only`, when given): are the outward
         * normals of both on the same side? The face is the nearest one carrying the interior point
         * of the piece, within twice the confusion radius (coincident faces are decided within it,
         * their axes or normals and their gap each within it); nullopt when there is none.
         */
        template <std::floating_point T>
        auto same_side(const Model<T> &m, const FaceUse &fu, const ShapeId &other, T tol, const std::set<std::uint32_t> *only = nullptr)
            -> std::optional<bool>
        {
            const auto uv = interior_point(m, fu.face, tol);
            if (!uv)
                return std::nullopt;
            const auto &s = *m.face(fu.face).surface;
            const auto p = s.value((*uv)[0], (*uv)[1]);
            auto n = cross(s.value((*uv)[0], (*uv)[1], 1, 0), s.value((*uv)[0], (*uv)[1], 0, 1));
            if (fu.orient == Orientation::Reversed)
                n = T(-1) * n;
            std::optional<bool> same;
            T nearest = T(2) * brep_confusion_factor<T> * tol;
            for (auto sh : explore<ShellId>(m, other))
                for (const auto &gu : m.shell(sh).faces)
                {
                    if (only && !only->contains(gu.face.index))
                        continue;
                    const auto &gs = *m.face(gu.face).surface;
                    const auto box = classify_detail::face_uv_box(m, gu.face);
                    // parameters and distance of p on the surface
                    std::array<T, 3> q;
                    if (const auto *es = dynamic_cast<const ElementarySurface<T> *>(&gs))
                    {
                        const auto uv = intersect_detail::primitive_uv(*es, box, p);
                        q = {uv[0], uv[1], norm(es->value(uv[0], uv[1]) - p)};
                    }
                    else
                        q = intersect_detail::project_on_surface(gs, box, p);
                    if (q[2] > nearest || classify(m, gu.face, point<T, 2>{q[0], q[1]}, tol) == Location::Out)
                        continue;
                    auto g = cross(gs.value(q[0], q[1], 1, 0), gs.value(q[0], q[1], 0, 1));
                    if (gu.orient == Orientation::Reversed)
                        g = T(-1) * g;
                    nearest = q[2];
                    same = n * g > T(0);
                }
            return same;
        }

        /// A use of an edge by a kept face: the co-edge `index` of the wire `wire` of `face`.
        struct RadialUse
        {
            FaceId face;
            Orientation face_orient;
            WireId wire;
            std::uint32_t index;
        };

        /**
         * Pairs the uses of the edge `e` around it, each face with its neighbour across its
         * material. The faces are ordered by the angle, around the edge, of a point of each face
         * at a small distance from the middle of the edge: faces tangent along the edge are
         * ordered by their curvature. nullopt when the sides of the materials do not match.
         */
        template <std::floating_point T>
        auto pair_radially(const Model<T> &m, EdgeId e, const std::vector<RadialUse> &uses)
            -> std::optional<std::vector<std::pair<std::size_t, std::size_t>>>
        {
            const auto &ed = m.edge(e);
            const T tm = (ed.u1 + ed.u2) / T(2);
            const auto p = ed.curve->value(tm);
            auto a = ed.curve->value(tm, 1);
            a = a / norm(a);
            point<T, 3> e1 = std::abs(a[0]) < T(0.9) ? point<T, 3>{1, 0, 0} : point<T, 3>{0, 1, 0};
            e1 = e1 - (e1 * a) * a;
            e1 = e1 / norm(e1);
            const auto e2 = cross(a, e1);
            const T h = T(1e-2) * norm(ed.curve->value(ed.u2) - ed.curve->value(ed.u1)) + T(1e-6);
            struct Side
            {
                T angle;
                bool ccw; ///< the material lies counter-clockwise from the face
            };
            std::vector<Side> sides;
            for (const auto &u : uses)
            {
                const auto &ce = m.wire(u.wire).coedges[u.index];
                const auto &s = *m.face(u.face).surface;
                const auto uv = ce.pcurve->value(tm);
                const auto su = s.value(uv[0], uv[1], 1, 0), sv = s.value(uv[0], uv[1], 0, 1);
                auto n = cross(su, sv);
                if (u.face_orient == Orientation::Reversed)
                    n = T(-1) * n;
                // the face lies on the left of its boundary, seen from the outside
                auto d = ed.curve->value(tm, 1);
                if (compose(ce.orient, u.face_orient) == Orientation::Reversed)
                    d = T(-1) * d;
                auto w = cross(n, d);
                w = (h / norm(w)) * w;
                // the step in (u, v) closest to w
                const T a11 = su * su, a12 = su * sv, a22 = sv * sv, b1 = su * w, b2 = sv * w;
                const T det = a11 * a22 - a12 * a12;
                if (!(std::abs(det) > T(0)))
                    return std::nullopt;
                const auto q = s.value(uv[0] + (b1 * a22 - b2 * a12) / det, uv[1] + (a11 * b2 - a12 * b1) / det);
                auto r = q - p;
                r = r - (r * a) * a;
                const T angle = std::atan2(r * e2, r * e1);
                // the material is behind the outward normal
                sides.push_back({angle, T(-1) * (n * cross(a, r)) > T(0)});
            }
            auto order = std::views::iota(std::size_t{0}, uses.size()) | std::ranges::to<std::vector>();
            std::ranges::sort(order, {}, [&](std::size_t i) { return sides[i].angle; });
            const auto k = order.size();
            std::vector<std::size_t> partner(k);
            for (std::size_t i = 0; i < k; ++i)
                partner[order[i]] = sides[order[i]].ccw ? order[(i + 1) % k] : order[(i + k - 1) % k];
            std::vector<std::pair<std::size_t, std::size_t>> pairs;
            for (std::size_t i = 0; i < k; ++i)
            {
                if (partner[partner[i]] != i || partner[i] == i)
                    return std::nullopt;
                if (i < partner[i])
                    pairs.emplace_back(i, partner[i]);
            }
            return pairs;
        }

        /**
         * Gives every pair of uses of an edge used by more than two kept faces, but the first,
         * its own copy of the edge. False when the faces cannot be paired.
         */
        template <std::floating_point T>
        bool separate_non_manifold_edges(Model<T> &m, const std::vector<std::pair<FaceId, Orientation>> &faces)
        {
            std::map<std::uint32_t, std::vector<RadialUse>> by_edge;
            for (const auto &[f, o] : faces)
                for (auto w : m.face(f).wires)
                    for (std::uint32_t i = 0; i < m.wire(w).coedges.size(); ++i)
                    {
                        const auto e = m.wire(w).coedges[i].edge;
                        if (!m.edge(e).degenerate)
                            by_edge[e.index].push_back({f, o, w, i});
                    }
            for (const auto &[index, uses] : by_edge)
            {
                if (uses.size() <= 2)
                    continue;
                const EdgeId e{index};
                const auto pairs = pair_radially(m, e, uses);
                if (!pairs)
                    return false;
                for (std::size_t k = 1; k < pairs->size(); ++k)
                {
                    const auto copy = m.add(Edge<T>{m.edge(e)});
                    for (auto i : {(*pairs)[k].first, (*pairs)[k].second})
                        m.wire(uses[i].wire).coedges[uses[i].index].edge = copy;
                }
            }
            return true;
        }

        /// Erases the edges and vertices from the indices `first_edge` and `first_vertex` that no wire or edge uses.
        template <std::floating_point T>
        void erase_unused(Model<T> &m, std::size_t first_edge, std::size_t first_vertex)
        {
            std::set<std::uint32_t> edges, vertices;
            for (auto w : m.template ids<WireId>())
                for (const auto &ce : m.wire(w).coedges)
                    edges.insert(ce.edge.index);
            for (auto c : m.template ids<CompoundId>())
                for (const auto &s : explore<EdgeId>(m, ShapeId{c}))
                    edges.insert(s.index);
            for (auto e : m.template ids<EdgeId>())
                if (e.index >= first_edge && !edges.contains(e.index))
                    m.erase(e);
            for (auto e : m.template ids<EdgeId>())
                vertices.insert({m.edge(e).v1.index, m.edge(e).v2.index});
            for (auto c : m.template ids<CompoundId>())
                for (const auto &v : explore<VertexId>(m, ShapeId{c}))
                    vertices.insert(v.index);
            for (auto v : m.template ids<VertexId>())
                if (v.index >= first_vertex && !vertices.contains(v.index))
                    m.erase(v);
        }

        /// Position of a face against the solids of an operand: inside if inside one of them.
        template <std::floating_point T>
        auto locate(const Model<T> &m, FaceId f, const std::vector<SolidId> &solids, T tol) -> std::optional<Location>
        {
            Location r = Location::Out;
            for (auto s : solids)
            {
                const auto l = classify(m, f, s, tol);
                if (!l)
                    return std::nullopt;
                if (*l == Location::On)
                    return Location::On;
                if (*l == Location::In)
                    r = Location::In;
            }
            return r;
        }
    }

    namespace boolean_detail
    {
        /// The operation on copies of the operands; on failure, `boolean` rolls the model back.
        template <std::floating_point T>
        auto run(Model<T> &m, BooleanOp op, const ShapeId &a, const ShapeId &b, const BooleanOptions<T> &opts, std::size_t first_edge,
                 std::size_t first_vertex) -> BuildResult<BooleanResult<T>>
        {
            const T tol = opts.tol;
            // ---- copies, so that the operands are not modified
            const auto ca = copy_shape(m, a, Placement<T>{}), cb = copy_shape(m, b, Placement<T>{});
            const auto solids_a = explore<SolidId>(m, ca), solids_b = explore<SolidId>(m, cb);

            // ---- 0. faces coincident on a part of their overlap: the B-spline face split along the
            //         border of the common region, its pieces off the other surface intersected with
            //         it only on the patches where they cross it
            SectionOptions<T> so{tol, opts.intersect};
            {
                auto r = split_partial_coincidences(m, ca, cb, tol);
                if (!r)
                    return std::unexpected(r.error());
                so.restricted = std::move(*r);
            }

            // ---- 1. section
            auto sec = section(m, ca, cb, so);
            if (!sec)
                return std::unexpected(sec.error());
            // (contact points split nothing; edges of tangency are section edges like the others)
            if (sec->singular || !sec->complete)
                return fail_non_generic<T>("singular or incomplete intersection");

            // ---- 2. coincident faces imprinted on each other, then split of the faces of both
            //         copies by the same edges
            std::map<std::uint32_t, std::vector<FaceCut<T>>> cuts;
            for (const auto &se : sec->edges)
            {
                cuts[se.face_a.index].push_back({se.edge, se.pcurve_a});
                cuts[se.face_b.index].push_back({se.edge, se.pcurve_b});
            }
            // the coincident faces and those each one coincides with
            std::map<std::uint32_t, std::vector<std::uint32_t>> partners;
            section_detail::VertexPool<T> pool{m, {}, tol};
            for (const auto &s : {ca, cb})
                pool.ids.append_range(explore<VertexId>(m, s));
            for (const auto &se : sec->edges)
                pool.ids.append_range(std::array{m.edge(se.edge).v1, m.edge(se.edge).v2});
            auto vertex_at = [&](const point<T, 3> &p) { return pool.at(p).first; };
            for (const auto &[fa, fb] : sec->coincident)
            {
                auto c = imprint_coincident(m, fa, fb, tol, vertex_at);
                if (!c)
                    return fail_non_generic<T>("coincident faces: an edge of one is not on the other");
                cuts[fa.index].append_range(c->first);
                cuts[fb.index].append_range(c->second);
                partners[fa.index].push_back(fb.index);
                partners[fb.index].push_back(fa.index);
            }
            // the edges lying on a face of the other operand (away from coincident faces, which the
            // imprint above has handled) are imprinted on it
            for (const auto &[e, f] : sec->edges_on_faces)
            {
                if (!m.alive(e) || !m.alive(f) || partners.contains(f.index))
                    continue;
                auto c = imprint_edge(m, e, f, tol, vertex_at);
                if (!c)
                    return fail_non_generic<T>("an edge lying on a face: no pcurve on it");
                cuts[f.index].append_range(*c);
            }
            std::map<std::uint32_t, std::vector<std::uint32_t>> pieces; // of the coincident faces
            for (auto &[index, c] : cuts)
            {
                // (an edge shared by two coincident faces of the other operand is imprinted by both)
                std::vector<FaceCut<T>> unique;
                for (auto &fc : c)
                    if (!std::ranges::contains(unique, fc.edge, &FaceCut<T>::edge))
                        unique.push_back(std::move(fc));
                c = std::move(unique);
                auto s = split_face(m, FaceId{index}, c, tol);
                if (!s)
                    return std::unexpected(s.error());
                if (partners.contains(index))
                    pieces[index] = s->faces | std::views::transform(&FaceId::index) | std::ranges::to<std::vector>();
                // a face touched along a line inside it is not split there
                for (auto e : s->dangling)
                    if (std::ranges::none_of(sec->edges, [&](const SectionEdge<T> &se) { return se.edge == e && se.tangent; }))
                        return fail_non_generic<T>("a section edge does not split a face");
            }
            // the faces a piece of a coincident face may lie on: the pieces of the faces it coincides
            // with (a face without cuts is its own piece)
            auto pieces_of = [&](std::uint32_t i) {
                const auto it = pieces.find(i);
                return it != pieces.end() ? it->second : std::vector{i};
            };
            std::map<std::uint32_t, std::set<std::uint32_t>> carriers;
            for (const auto &[i, js] : partners)
                for (auto p : pieces_of(i))
                    for (auto j : js)
                        carriers[p].insert_range(pieces_of(j));
            // the section wires and compound are not part of the result (their edges are)
            for (auto w : sec->wires)
                m.erase(w);
            m.erase(sec->compound);

            // ---- 3. classification and 4. selection
            struct Kept
            {
                FaceId face;
                Orientation orient;
                FaceOrigin origin;
            };
            std::vector<Kept> kept;
            for (auto [shape, other, other_shape, origin] :
                 {std::tuple{ca, &solids_b, cb, FaceOrigin::A}, std::tuple{cb, &solids_a, ca, FaceOrigin::B}})
                for (auto sh : explore<ShellId>(m, shape))
                    for (const auto &fu : m.shell(sh).faces)
                    {
                        // a piece of a coincident face lies on a piece of a face it coincides with, within
                        // the confusion radius, beyond the tolerance of the classification
                        const auto c = carriers.find(fu.face.index);
                        const auto on_other = c != carriers.end() ? same_side(m, fu, other_shape, tol, &c->second) : std::nullopt;
                        const auto where = on_other ? std::optional{Location::On} : locate(m, fu.face, *other, tol);
                        if (!where)
                            return fail_non_generic<T>("a face without inside");
                        bool keep{}, flip{};
                        if (*where == Location::On)
                        {
                            // a piece lying on a face of the other operand: one copy of the coincident
                            // pieces is kept (A's) when the solids are on the same side, none otherwise;
                            // for a cut, A's piece stays where B touches it from outside
                            const auto same = on_other ? on_other : same_side(m, fu, other_shape, tol);
                            if (!same)
                                return fail_non_generic<T>("a face lies on the boundary of the other operand");
                            if (origin == FaceOrigin::A)
                                keep = op == BooleanOp::Cut ? !*same : *same;
                        }
                        else
                        {
                            const bool in = *where == Location::In;
                            switch (op)
                            {
                            case BooleanOp::Fuse: keep = !in; break;
                            case BooleanOp::Common: keep = in; break;
                            case BooleanOp::Cut:
                                keep = origin == FaceOrigin::A ? !in : in;
                                flip = origin == FaceOrigin::B;
                                break;
                            }
                        }
                        if (keep)
                            kept.push_back({fu.face, flip ? reverse(fu.orient) : fu.orient, origin});
                    }

            // ---- 5. assembly: faces paired around the edges where solids touch, shells by shared
            //         edges, solids by closed shells
            {
                std::vector<std::pair<FaceId, Orientation>> faces;
                for (const auto &k : kept)
                    faces.emplace_back(k.face, k.orient);
                if (!separate_non_manifold_edges(m, faces))
                    return fail_non_generic<T>("faces meeting at an edge cannot be paired");
            }
            std::map<std::uint32_t, std::vector<std::size_t>> by_edge;
            for (std::size_t k = 0; k < kept.size(); ++k)
                for (auto e : explore<EdgeId>(m, ShapeId{kept[k].face}))
                    by_edge[e.index].push_back(k);
            auto parent = std::views::iota(std::size_t{0}, kept.size()) | std::ranges::to<std::vector>();
            auto root = [&](this auto &self, std::size_t x) -> std::size_t { return parent[x] == x ? x : parent[x] = self(parent[x]); };
            for (const auto &[e, ks] : by_edge)
                for (std::size_t i = 1; i < ks.size(); ++i)
                    parent[root(ks[i])] = root(ks[0]);
            std::map<std::size_t, Shell> shells;
            for (std::size_t k = 0; k < kept.size(); ++k)
                shells[root(k)].faces.push_back(FaceUse{kept[k].face, kept[k].orient});
            // the shells keep the sense of their faces: outwards for a solid (positive volume), into
            // a cavity (negative); computed once, the solids are made from them directly
            std::vector<std::pair<ShellId, T>> outers, voids;
            for (auto &[r, sh] : shells)
            {
                const auto id = m.add(sh);
                if (!is_closed(m, id))
                    return detail::fail(BuildErrc::ShellNotClosed, "boolean: a shell of the result is open", {id});
                if (!is_orientable(m, id))
                    return detail::fail(BuildErrc::ShellNotOrientable, "boolean: a shell of the result is not orientable", {id});
                m.shell(id).closed = true;
                const T v = signed_volume(m, id, 64);
                if (!(std::abs(v) > T(0)) || !std::isfinite(v))
                    return detail::fail(BuildErrc::ZeroVolume, "boolean: a shell of the result encloses no volume", {id});
                (v > T(0) ? outers : voids).emplace_back(id, v);
            }
            BooleanResult<T> res;
            for (const auto &[oid, ov] : outers)
                res.solids.push_back(m.add(Solid{oid, {}}));
            // each cavity in the smallest solid holding a point of it (the solids without their
            // cavities while they are located)
            std::vector<std::pair<std::size_t, ShellId>> cavities;
            for (const auto &[vid, vv] : voids)
            {
                const auto f = m.shell(vid).faces.front().face;
                const auto uv = interior_point(m, f, tol);
                if (!uv)
                    return fail_non_generic<T>("a face of a cavity without inside");
                const auto p = m.face(f).surface->value((*uv)[0], (*uv)[1]);
                std::optional<std::size_t> owner;
                for (std::size_t k = 0; k < outers.size(); ++k)
                {
                    const auto l = classify_detail::in_solid(m, res.solids[k], p, tol);
                    if (!l)
                        return fail_non_generic<T>("a cavity cannot be located");
                    if (*l == Location::In && (!owner || outers[k].second < outers[*owner].second))
                        owner = k;
                }
                if (!owner)
                    return detail::fail(BuildErrc::HoleOutsideOuter, "boolean: a cavity of the result is outside every solid", {vid});
                cavities.emplace_back(*owner, vid);
            }
            for (const auto &[k, vid] : cavities)
                m.solid(res.solids[k]).voids.push_back(vid);
            // the copies give their faces to the result
            std::set<std::uint32_t> kept_faces;
            for (const auto &k : kept)
                kept_faces.insert(k.face.index);
            erase_shape(m, ca, kept_faces);
            erase_shape(m, cb, kept_faces);
            erase_unused(m, first_edge, first_vertex);

            for (const auto &k : kept)
                res.faces.emplace_back(k.face, k.origin);
            if (res.solids.size() == 1)
                res.shape = res.solids.front();
            else
                res.shape = m.add(Compound{std::vector<ShapeId>(res.solids.begin(), res.solids.end())});
            return res;
        }
    }

    /**
     * @brief Boolean operation `op` between the solids of `a` and those of `b` (a solid or a
     * compound of solids). The operands are not modified; the result is built from copies of
     * their faces. On failure the model is left as it was: everything the operation added is
     * erased, also when a step throws `BRepError`: a `BuildFailure` keeps its reason and its
     * entities (a dead id is `InvalidId`, a classification whose rays are all capped `NonGeneric`),
     * any other gives `NonGeneric`.
     */
    template <std::floating_point T>
    [[nodiscard]] auto boolean(Model<T> &m, BooleanOp op, const ShapeId &a, const ShapeId &b, std::type_identity_t<BooleanOptions<T>> opts = {})
        -> BuildResult<BooleanResult<T>>
    {
        using namespace boolean_detail;
        if (!m.alive(a) || !m.alive(b))
            return detail::fail(BuildErrc::InvalidId, "boolean: dead or invalid shape");
        if (explore<SolidId>(m, a).empty() || explore<SolidId>(m, b).empty())
            return detail::fail(BuildErrc::InvalidId, "boolean: both operands must contain solids");
        const auto size = detail::arena_size(m);
        try
        {
            auto r = run(m, op, a, b, opts, size[1], size[0]);
            if (!r)
                detail::rollback(m, size);
            return r;
        }
        catch (const BuildFailure &e)
        {
            detail::rollback(m, size);
            auto err = e.error();
            err.message = "boolean: " + err.message;
            return std::unexpected(std::move(err));
        }
        catch (const BRepError &e)
        {
            detail::rollback(m, size);
            return fail_non_generic<T>(e.what());
        }
    }

    template <std::floating_point T>
    [[nodiscard]] auto fuse(Model<T> &m, const ShapeId &a, const ShapeId &b, std::type_identity_t<BooleanOptions<T>> opts = {})
    {
        return boolean(m, BooleanOp::Fuse, a, b, opts);
    }

    template <std::floating_point T>
    [[nodiscard]] auto cut(Model<T> &m, const ShapeId &a, const ShapeId &b, std::type_identity_t<BooleanOptions<T>> opts = {})
    {
        return boolean(m, BooleanOp::Cut, a, b, opts);
    }

    template <std::floating_point T>
    [[nodiscard]] auto common(Model<T> &m, const ShapeId &a, const ShapeId &b, std::type_identity_t<BooleanOptions<T>> opts = {})
    {
        return boolean(m, BooleanOp::Common, a, b, opts);
    }
}
