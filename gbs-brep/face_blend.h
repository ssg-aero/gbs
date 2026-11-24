#pragma once
/**
 * @file face_blend.h
 * @brief Fillets between two faces, or two sets of faces joined tangentially, that have no edge
 * in common (stage 8, `docs/sources/design/brep_palier8_face_blends.md`, PR 1:
 * `docs/sources/design/brep8_pr01_face_blends.md`).
 *
 * The faces between them (the strip) are taken out of a copy of the solid, and the two sides
 * continued until they meet: a virtual sharp edge, which the fillets of `blend.h` round off. The
 * result is kept only if the strip lies wholly on the side the fillet takes off (or adds), its
 * faces out of reach of the fillet: the pockets of `brep_swallowed_faces.md`, whose strip ends at
 * end faces here rather than between two edges given.
 *
 * Building the copy (§ 2.2 of the design note):
 * - the **spokes**: the edges that end on the strip, between two faces kept (two faces of the
 *   sides, a face of a side and an end face, a seam), continued past the strip until they meet the
 *   face across (of the other side, or of the side the end face does not touch): their ends there
 *   are the vertices of the virtual edge;
 * - the **pieces** of the virtual edge: the intersection of a face of each side, continued,
 *   between two of those vertices, along the strip;
 * - in each face kept, the run of its wire along the strip (the old spokes, the edges with the
 *   strip) replaced by the spokes continued and the pieces; the strip erased.
 */
#include <gbs-brep/blend.h>
#include <gbs/curveonsurface.h>

#include <algorithm>
#include <array>
#include <map>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <ranges>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace gbs::brep
{
    namespace face_blend_detail
    {
        /**
         * The faces between the two sides f[0] and f[1] (§ 2.1): a connected set of faces of
         * neither side, every face of both sides touching it, each face next to it of a side,
         * of it, or an end face (one that touches both sides).
         */
        template <std::floating_point T>
        struct Strip
        {
            std::array<std::set<FaceId>, 2> f;     ///< the two sides
            std::set<FaceId> faces;                ///< the faces swallowed
            std::set<FaceId> ends;                 ///< the end faces
            std::set<EdgeId> edges;                ///< the edges of the faces swallowed
            std::set<VertexId> vertices;           ///< and their vertices
            std::vector<point<T, 3>> samples;      ///< points of the strip, all to be swallowed
        };

        /// The faces next to the face g across its edges (g itself across a seam).
        template <std::floating_point T>
        auto neighbours(const Model<T> &m, const blend_detail::SolidView<T> &sv, FaceId g) -> std::set<FaceId>
        {
            std::set<FaceId> r;
            for (auto wid : m.face(g).wires)
                for (const auto &ce : m.wire(wid).coedges)
                    if (!m.edge(ce.edge).degenerate)
                        for (auto u : sv.users(ce.edge))
                            if (u != g)
                                r.insert(u);
            return r;
        }

        /**
         * The candidate strips between the sides f, smallest first. First pass: a face next to
         * the strip that touches both sides is an end face. Second pass (a strip of several
         * faces round a closed side, a chamfer round the foot of a boss in a stadium): such a
         * face is part of the strip, which then has no end face.
         */
        template <std::floating_point T>
        auto strips(const Model<T> &m, const blend_detail::SolidView<T> &sv, const std::array<std::set<FaceId>, 2> &f) -> std::vector<Strip<T>>
        {
            auto in_sides = [&](FaceId g) { return f[0].contains(g) || f[1].contains(g); };
            std::map<FaceId, std::set<FaceId>> nb;
            auto next = [&](FaceId g) -> const std::set<FaceId> & {
                auto it = nb.find(g);
                if (it == nb.end())
                    it = nb.emplace(g, neighbours(m, sv, g)).first;
                return it->second;
            };
            auto touches = [&](FaceId g, std::size_t k) { return std::ranges::any_of(next(g), [&](FaceId h) { return f[k].contains(h); }); };
            std::vector<Strip<T>> r;
            std::set<std::set<FaceId>> seen;
            for (const bool second : {false, true})
                for (auto g0 : f[0])
                    for (auto seed : next(g0))
                    {
                        if (in_sides(seed) || sv.shell_of(seed) != sv.shell_of(g0))
                            continue;
                        Strip<T> s;
                        s.f = f;
                        s.faces.insert(seed);
                        std::vector<FaceId> todo{seed};
                        while (!todo.empty())
                        {
                            const auto g = todo.back();
                            todo.pop_back();
                            for (auto h : next(g))
                            {
                                if (in_sides(h) || s.faces.contains(h) || s.ends.contains(h))
                                    continue;
                                if (!second && touches(h, 0) && touches(h, 1))
                                    s.ends.insert(h);
                                else
                                    s.faces.insert(h), todo.push_back(h);
                            }
                        }
                        // every face of both sides along the strip
                        bool ok = true;
                        for (std::size_t k = 0; k < 2 && ok; ++k)
                            for (auto g : f[k])
                                ok = ok && std::ranges::any_of(next(g), [&](FaceId h) { return s.faces.contains(h); });
                        if (!ok || !seen.insert(s.faces).second)
                            continue;
                        for (auto g : s.faces)
                            for (auto wid : m.face(g).wires)
                                for (const auto &ce : m.wire(wid).coedges)
                                {
                                    s.edges.insert(ce.edge);
                                    s.vertices.insert(m.edge(ce.edge).v1);
                                    s.vertices.insert(m.edge(ce.edge).v2);
                                }
                        r.push_back(std::move(s));
                    }
            std::ranges::stable_sort(r, {}, [](const Strip<T> &s) { return s.faces.size(); });
            return r;
        }

        /**
         * Whether the strip st is the rest of the solid: its faces, its end faces and the sides
         * all those of the shell (no end face touching both sides bounds it, as at a chamfer that
         * ends on another chamfer).
         */
        template <std::floating_point T>
        bool rest_of_solid(const Model<T> &m, const blend_detail::SolidView<T> &sv, const Strip<T> &st)
        {
            return std::ranges::all_of(m.shell(sv.shell_of(*st.f[0].begin())).faces, [&](const FaceUse &fu) {
                return st.faces.contains(fu.face) || st.ends.contains(fu.face) || st.f[0].contains(fu.face) || st.f[1].contains(fu.face);
            });
        }

        /// The points of the strip that must all be swallowed: inside each face, on each edge (at 8 parameters), at each vertex.
        template <std::floating_point T>
        void sample(const Model<T> &m, Strip<T> &s, T tol)
        {
            for (auto g : s.faces)
                for (const auto &uv : interior_points(m, g, tol))
                    s.samples.push_back(m.face(g).surface->value(uv[0], uv[1]));
            for (auto x : s.edges)
            {
                const auto &ed = m.edge(x);
                if (ed.degenerate || !ed.curve)
                    continue;
                for (int i = 1; i < 8; ++i)
                    s.samples.push_back(ed.curve->value(ed.u1 + (ed.u2 - ed.u1) * T(i) / T(8)));
            }
            for (auto v : s.vertices)
                s.samples.push_back(m.vertex(v).pnt);
        }

        /// The distance from p to the nearest sample of the strip.
        template <std::floating_point T>
        auto distance_to(const Strip<T> &s, const point<T, 3> &p) -> T
        {
            T d = std::numeric_limits<T>::infinity();
            for (const auto &q : s.samples)
                d = std::min(d, norm(q - p));
            return d;
        }

        /**
         * The curve of the edge continued by about d past each end, with its parameter range:
         * a line or a circle on a wider range (an arc short of a turn), a B-spline by the
         * polynomials of its end spans; none for any other curve.
         */
        template <std::floating_point T>
        auto continued(const Edge<T> &ed, T d, const std::vector<std::pair<std::shared_ptr<Surface<T, 3>>, std::shared_ptr<Surface<T, 3>>>> &supports)
            -> std::optional<std::tuple<std::shared_ptr<Curve<T, 3>>, T, T>>
        {
            T d1 = d / std::max(norm(ed.curve->value(ed.u1, 1)), T(1e-300)), d2 = d / std::max(norm(ed.curve->value(ed.u2, 1)), T(1e-300));
            if (const auto *ec = dynamic_cast<const ElementaryCurve<T> *>(ed.curve.get()))
            {
                if (dynamic_cast<const EllipseCurve<T> *>(ec))
                {
                    const T room = T(2) * std::numbers::pi_v<T> - (ed.u2 - ed.u1);
                    if (!(room > T(0)))
                        return std::nullopt;
                    const T k = std::min(T(1), room / (T(2) * (d1 + d2)));
                    d1 *= k, d2 *= k;
                }
                return std::tuple{intersect_detail::restricted(*ec, ed.u1 - d1, ed.u2 + d2), ed.u1 - d1, ed.u2 + d2};
            }
            auto go_on = [&]<bool rational>(const BSCurveGeneral<T, 3, rational> &b) -> std::optional<std::tuple<std::shared_ptr<Curve<T, 3>>, T, T>> {
                const auto bb = b.bounds();
                const std::array<T, 2> dt{std::max(T(0), bb[0] - (ed.u1 - d1)), std::max(T(0), (ed.u2 + d2) - bb[1])};
                try
                {
                    auto x = std::make_shared<std::conditional_t<rational, BSCurveRational<T, 3>, BSCurve<T, 3>>>(extrapolated(b, dt));
                    const auto xb = x->bounds();
                    return std::tuple{std::shared_ptr<Curve<T, 3>>{x}, xb[0], xb[1]};
                }
                catch (const std::invalid_argument &)
                {
                    return std::nullopt;
                }
            };
            if (const auto *b = dynamic_cast<const BSCurve<T, 3> *>(ed.curve.get()))
                return go_on(*b);
            if (const auto *b = dynamic_cast<const BSCurveRational<T, 3> *>(ed.curve.get()))
                return go_on(*b);
            // a curve on the surface of one of its faces: its pcurve continued, on that surface continued
            if (const auto *cs = dynamic_cast<const CurveOnSurface<T, 3> *>(ed.curve.get()))
                for (const auto &[old, wider] : supports)
                {
                    if (old.get() != cs->p_basisSurface().get())
                        continue;
                    const auto base = cs->p_basisCurve();
                    auto go_on2 = [&]<bool rational>(const BSCurveGeneral<T, 2, rational> &b) -> std::shared_ptr<Curve<T, 2>> {
                        const auto bb = b.bounds();
                        const std::array<T, 2> dt{std::max(T(0), bb[0] - (ed.u1 - d1)), std::max(T(0), (ed.u2 + d2) - bb[1])};
                        try
                        {
                            return std::make_shared<std::conditional_t<rational, BSCurveRational<T, 2>, BSCurve<T, 2>>>(extrapolated(b, dt));
                        }
                        catch (const std::invalid_argument &)
                        {
                            return nullptr;
                        }
                    };
                    std::shared_ptr<Curve<T, 2>> pc;
                    if (const auto *b = dynamic_cast<const BSCurve<T, 2> *>(base.get()))
                        pc = go_on2(*b);
                    else if (const auto *b = dynamic_cast<const BSCurveRational<T, 2> *>(base.get()))
                        pc = go_on2(*b);
                    if (!pc)
                        return std::nullopt;
                    const auto pb = pc->bounds();
                    return std::tuple{std::shared_ptr<Curve<T, 3>>{std::make_shared<CurveOnSurface<T, 3>>(pc, wider)}, pb[0], pb[1]};
                }
            return std::nullopt;
        }

        /// An edge ending on the strip between two faces kept (§ 2.2), continued to the face across.
        template <std::floating_point T>
        struct Spoke
        {
            EdgeId edge;                         ///< the edge (on the copy)
            std::array<FaceId, 2> faces;         ///< its faces (a seam: the same twice)
            VertexId far, near;                  ///< its end away from the strip, on it
            FaceId across;                       ///< the face it is continued to
            std::shared_ptr<Curve<T, 3>> curve;  ///< continued
            T t_far{}, t_x{};                    ///< the parameters of its ends, continued
            point<T, 3> x{};                     ///< where it meets the face across
            std::size_t vertex{};                ///< the index of that point among the vertices of the virtual edge
            EdgeId continued;                    ///< the new edge, from far to the virtual edge
        };

        /// A vertex of the virtual edge: where spokes meet the faces across, the faces it lies on.
        template <std::floating_point T>
        struct Node
        {
            point<T, 3> p{};
            std::set<FaceId> on;
            VertexId v;
        };

        /// A piece of the virtual edge: on a face of each side, between two nodes, its pcurves on both.
        template <std::floating_point T>
        struct Piece
        {
            EdgeId edge;
            std::array<FaceId, 2> faces;
            std::array<std::size_t, 2> nodes; ///< at its start and end
            std::map<FaceId, std::shared_ptr<Curve<T, 2>>> pc;
        };

        /// The pcurve of the edge e on the surface s by projection (`hint`: where it starts when it lies on a seam).
        template <std::floating_point T>
        auto projected(const Model<T> &m, EdgeId e, const Surface<T, 3> &s, const std::optional<point<T, 2>> &hint, T tol) -> std::shared_ptr<Curve<T, 2>>
        {
            MakeFaceOptions<T> mo;
            mo.tol = tol;
            mo.pcurve_tol = std::max(tol / T(10), brep_pcurve_approx_tol<T>);
            const auto fit = project_pcurve(m.edge(e), s, surface_closure(s, tol), hint, mo);
            if (!fit)
                throw BuildFailure({BuildErrc::NonGeneric, "a curve of the virtual edge off its face", {e}});
            return fit->pcurve;
        }

        /**
         * The copy of the solid with the strip taken out and its sides continued until they meet
         * (§ 2.2), built on the copy ns of the solid (the strip given on the copy): the pieces of
         * the virtual edge, in no particular order. Throws BuildFailure where it cannot be built.
         */
        template <std::floating_point T>
        auto close_over(Model<T> &m, SolidId ns, const Strip<T> &s, const FilletOptions<T> &o, const std::string &op) -> std::vector<EdgeId>
        {
            const T tol = o.tol;
            const blend_detail::SolidView<T> sv{m, ns};
            auto side_of = [&](FaceId g) -> int { return s.f[0].contains(g) ? 0 : s.f[1].contains(g) ? 1 : s.ends.contains(g) ? -1 : -2; };
            // the faces kept along the strip, continued past it by twice its size
            T reach = T(0);
            for (const auto &p : s.samples)
                for (const auto &q : {s.samples.front(), s.samples.back()})
                    reach = std::max(reach, norm(p - q));
            std::map<FaceId, thick_detail::Extended<T>> ext;
            for (const auto &set : {s.f[0], s.f[1], s.ends})
                for (auto g : set)
                    ext[g] = thick_detail::extended_surface(m.face(g).surface, classify_detail::face_uv_box(m, g), T(2) * reach, tol);
            // the spokes
            std::vector<Spoke<T>> spokes;
            std::set<EdgeId> seen;
            for (auto g : std::views::keys(ext))
                for (auto wid : m.face(g).wires)
                    for (const auto &ce : m.wire(wid).coedges)
                    {
                        const auto &ed = m.edge(ce.edge);
                        if (ed.degenerate || s.edges.contains(ce.edge) || !seen.insert(ce.edge).second)
                            continue;
                        const bool on1 = s.vertices.contains(ed.v1), on2 = s.vertices.contains(ed.v2);
                        if (!on1 && !on2)
                            continue;
                        if (on1 && on2)
                            blend_detail::refuse<T>(op, "an edge between the faces kept that runs from the strip to the strip (not handled)", {ce.edge});
                        const auto us = sv.users(ce.edge);
                        if (us.size() != 2 || side_of(us[0]) < -1 || side_of(us[1]) < -1)
                            blend_detail::refuse<T>(op, "an edge ending on the faces between that is not between the faces given or their end faces", {ce.edge});
                        Spoke<T> sp{ce.edge, {us[0], us[1]}, on1 ? ed.v2 : ed.v1, on1 ? ed.v1 : ed.v2, {}, {}, {}, {}, {}, {}, {}};
                        const int a = side_of(us[0]), b = side_of(us[1]);
                        if (a == -1 && b == -1)
                            blend_detail::refuse<T>(op, "an edge between two end faces ending on the strip (not handled)", {ce.edge});
                        if (a + b == 1 && a >= 0 && b >= 0)
                            blend_detail::refuse<T>(op, "faces given that already meet next to the faces between them (not handled)", {ce.edge});
                        // the side across: the other one than that of its faces
                        const std::size_t k = (a == 0 || b == 0) ? 1 : 0;
                        std::vector<std::pair<std::shared_ptr<Surface<T, 3>>, std::shared_ptr<Surface<T, 3>>>> supports;
                        for (auto u : us)
                            supports.push_back({m.face(u).surface, ext.at(u).surface});
                        const auto cont = continued(ed, T(2) * reach, supports);
                        if (!cont)
                            blend_detail::refuse<T>(op, "an edge ending on the strip that cannot be continued", {ce.edge});
                        const auto &[c, lo, hi] = *cont;
                        sp.curve = c;
                        sp.t_far = sp.far == ed.v1 ? ed.u1 : ed.u2;
                        const T t_near = sp.near == ed.v1 ? ed.u1 : ed.u2;
                        const T out = t_near > sp.t_far ? T(1) : T(-1); // toward the strip and past it
                        const auto pn = m.vertex(sp.near).pnt;
                        T best = std::numeric_limits<T>::infinity();
                        for (auto r : s.f[k])
                        {
                            const auto &x = ext.at(r);
                            const auto uv0 = extrema_surf_pnt(*x.surface, pn, tol / T(100), std::array<std::size_t, 2>{20, 20});
                            const auto hit = blend_detail::meet(*c, *x.surface, t_near, point<T, 2>{uv0[0], uv0[1]}, tol);
                            if (!hit || (hit->first - t_near) * out < -tol || hit->first < lo - tol || hit->first > hi + tol)
                                continue;
                            const T gap = std::abs(hit->first - t_near);
                            if (gap < best)
                                best = gap, sp.across = r, sp.t_x = hit->first, sp.x = c->value(hit->first);
                        }
                        if (!(best < std::numeric_limits<T>::infinity()))
                            blend_detail::refuse<T>(op, "faces given that do not meet once continued over the faces between them (parallel faces: a full round fillet), or that no strip of faces joins (not handled)", {ce.edge});
                        spokes.push_back(std::move(sp));
                    }

            // the nodes: the points where spokes meet the faces across, merged
            const T near = T(10) * tol;
            std::vector<Node<T>> nodes;
            for (auto &sp : spokes)
            {
                std::size_t i = 0;
                while (i < nodes.size() && norm(nodes[i].p - sp.x) > near)
                    ++i;
                if (i == nodes.size())
                    nodes.push_back({sp.x, {}, {}});
                nodes[i].on.insert(sp.faces.begin(), sp.faces.end());
                nodes[i].on.insert(sp.across);
                sp.vertex = i;
            }
            IntersectOptions<T> io = o.intersect;
            io.tol = std::max(io.tol, tol / T(10));
            std::map<std::pair<FaceId, FaceId>, SurfaceSurfaceResult<T>> meets; // (of the sides continued)
            auto meet_of = [&](FaceId p0, FaceId p1) -> const SurfaceSurfaceResult<T> & {
                auto it = meets.find({p0, p1});
                if (it == meets.end())
                {
                    const auto &x0 = ext.at(p0), &x1 = ext.at(p1);
                    it = meets.emplace(std::pair{p0, p1}, intersect(*x0.surface, x0.box, *x1.surface, x1.box, io)).first;
                }
                return it->second;
            };
            // (a face round a period has a seam, which ends on the strip: `make_face` refuses one without)
            if (spokes.empty())
                blend_detail::refuse<T>(op, "a strip without an edge ending on it between the faces kept (not handled)");
            for (auto &nd : nodes)
                nd.v = m.add(Vertex<T>{nd.p, tol});
            // the spokes continued: new edges from their far ends to the nodes, a pcurve on each face (each use of a seam)
            std::map<std::pair<EdgeId, FaceId>, std::vector<std::pair<point<T, 2>, std::shared_ptr<Curve<T, 2>>>>> spoke_pc; // (old edge, face) -> (hint, pcurve)
            for (auto &sp : spokes)
            {
                const bool fwd = sp.t_x > sp.t_far;
                Edge<T> e;
                e.curve = sp.curve;
                e.u1 = fwd ? sp.t_far : sp.t_x, e.u2 = fwd ? sp.t_x : sp.t_far;
                e.v1 = fwd ? sp.far : nodes[sp.vertex].v, e.v2 = fwd ? nodes[sp.vertex].v : sp.far;
                e.tol = std::max(tol, m.edge(sp.edge).tol);
                sp.continued = m.add(std::move(e));
                std::vector<blend_detail::OnSurface<T>> pcs;
                for (auto g : std::set<FaceId>{sp.faces.begin(), sp.faces.end()})
                    for (auto wid : m.face(g).wires)
                        for (const auto &ce : m.wire(wid).coedges)
                            if (ce.edge == sp.edge)
                            {
                                const auto &old = m.edge(sp.edge);
                                const auto hint = coedge_uv(m, CoEdge<T>{sp.edge, Orientation::Forward, ce.pcurve}, sp.far == old.v1 ? old.u1 : old.u2);
                                const auto pc = projected(m, sp.continued, *ext.at(g).surface, std::optional<point<T, 2>>{hint}, tol);
                                spoke_pc[{sp.edge, g}].push_back({hint, pc});
                                pcs.push_back({ext.at(g).surface.get(), pc});
                            }
                blend_detail::settle(m, sp.continued, std::span<const blend_detail::OnSurface<T>>{pcs});
            }
            // the pieces: on a face of each side, between two nodes on both, along the strip
            T dmax = T(0);
            for (const auto &nd : nodes)
                dmax = std::max(dmax, distance_to(s, nd.p));
            const T along = T(4) * dmax + near;
            std::vector<Piece<T>> pieces;
            for (auto p0 : s.f[0])
                for (auto p1 : s.f[1])
                {
                    std::vector<std::size_t> at;
                    for (std::size_t i = 0; i < nodes.size(); ++i)
                        if (nodes[i].on.contains(p0) && nodes[i].on.contains(p1))
                            at.push_back(i);
                    if (at.empty())
                        continue;
                    const auto &x0 = ext.at(p0), &x1 = ext.at(p1);
                    for (const auto &c : meet_of(p0, p1).curves)
                    {
                        const T eps = T(1e-12) * (T(1) + std::abs(c.b - c.a));
                        std::vector<std::pair<T, std::size_t>> ts;
                        for (auto i : at)
                        {
                            T bt = c.a, gap = std::numeric_limits<T>::infinity();
                            for (int j = 0; j <= 16; ++j)
                            {
                                const auto x = extrema_curve_point(*c.curve, nodes[i].p, c.a + (c.b - c.a) * T(j) / T(16), eps);
                                if (x[1] < gap)
                                    gap = x[1], bt = x[0];
                            }
                            if (gap <= near)
                                ts.push_back({bt, i});
                        }
                        if (ts.empty())
                            continue;
                        std::ranges::sort(ts);
                        const T period = c.b - c.a;
                        std::vector<std::array<std::size_t, 2>> arcs; // indices in ts
                        for (std::size_t j = 0; j + 1 < ts.size(); ++j)
                            arcs.push_back({j, j + 1});
                        if (c.closed)
                            arcs.push_back({ts.size() - 1, 0});
                        for (const auto &[ia, ib] : arcs)
                        {
                            const T u1 = ts[ia].first, u2 = ib > ia ? ts[ib].first : ts[ib].first + period;
                            if (!(u2 - u1 > eps))
                                continue;
                            const bool across = u2 > c.b + eps;
                            const auto *ec = dynamic_cast<const ElementaryCurve<T> *>(c.curve.get());
                            if (across && !ec)
                                continue; // (a closed B-spline across its origin: not handled)
                            const auto curve = across ? intersect_detail::restricted(*ec, u1, u2) : c.curve;
                            if (distance_to(s, curve->value((u1 + u2) / T(2))) > along)
                                continue;
                            const auto va = nodes[ts[ia].second].v, vb = nodes[ts[ib].second].v;
                            Edge<T> e;
                            e.curve = curve;
                            e.u1 = u1, e.u2 = u2, e.v1 = va, e.v2 = vb, e.tol = std::max(tol, c.deviation);
                            const auto id = m.add(std::move(e));
                            Piece<T> pc{id, {p0, p1}, {ts[ia].second, ts[ib].second}, {}};
                            pc.pc[p0] = across ? projected(m, id, *x0.surface, std::optional<point<T, 2>>{}, tol) : c.pc1;
                            pc.pc[p1] = across ? projected(m, id, *x1.surface, std::optional<point<T, 2>>{}, tol) : c.pc2;
                            const std::vector<blend_detail::OnSurface<T>> pcs{{x0.surface.get(), pc.pc[p0]}, {x1.surface.get(), pc.pc[p1]}};
                            blend_detail::settle(m, id, std::span<const blend_detail::OnSurface<T>>{pcs});
                            pieces.push_back(std::move(pc));
                        }
                    }
                }
            if (pieces.empty())
                blend_detail::refuse<T>(op, "faces given that do not meet under the faces between them once continued (parallel faces: a full round fillet?)");
            // the faces kept continued
            for (const auto &[g, x] : ext)
                m.face(g).surface = x.surface;
            // in each face kept, the run of its wire along the strip replaced
            std::map<EdgeId, const Spoke<T> *> spoke_of;
            for (const auto &sp : spokes)
                spoke_of[sp.edge] = &sp;
            std::set<std::pair<std::size_t, FaceId>> used; // (a piece, a face it has been put in)
            std::set<EdgeId> removed;
            for (const auto &[g, x] : ext)
            {
                auto doomed = [&](EdgeId e) { return s.edges.contains(e) || (spoke_of.contains(e) && std::ranges::contains(spoke_of.at(e)->faces, g)); };
                const auto wires = m.face(g).wires;
                for (auto wid : wires)
                {
                    const auto cs = m.wire(wid).coedges;
                    const std::size_t n = cs.size();
                    std::size_t nd = 0;
                    for (const auto &ce : cs)
                        nd += doomed(ce.edge) ? 1 : 0;
                    if (nd == 0)
                        continue;
                    // the direction of the strip along the face, from its first edge with it
                    // (the edge of the run with the strip nearest the point p)
                    auto strip_dir = [&](const std::vector<CoEdge<T>> &run, const point<T, 3> &p) -> std::pair<point<T, 3>, point<T, 3>> {
                        std::pair<point<T, 3>, point<T, 3>> best{};
                        T gap = std::numeric_limits<T>::infinity();
                        for (const auto &ce : run)
                            if (s.edges.contains(ce.edge))
                            {
                                const auto &ed = m.edge(ce.edge);
                                for (int i = 1; i < 8; ++i)
                                {
                                    const T u = ed.u1 + (ed.u2 - ed.u1) * T(i) / T(8);
                                    const auto q = ed.curve->value(u);
                                    if (norm(q - p) < gap)
                                    {
                                        const auto d = ed.curve->value(u, 1);
                                        gap = norm(q - p), best = {q, ce.orient == Orientation::Forward ? d : T(-1) * d};
                                    }
                                }
                            }
                        return best;
                    };
                    // the pieces on g from node a to node b (a == b: round a loop, or nothing), as co-edges
                    auto path = [&](std::size_t a, std::size_t b, const std::vector<CoEdge<T>> &along_run) {
                        std::vector<CoEdge<T>> r;
                        std::size_t cur = a;
                        bool moved = false;
                        while (!(cur == b && moved))
                        {
                            std::optional<std::size_t> pick;
                            for (std::size_t i = 0; i < pieces.size() && !pick; ++i)
                                if (!used.contains({i, g}) && std::ranges::contains(pieces[i].faces, g) && (pieces[i].nodes[0] == cur || pieces[i].nodes[1] == cur))
                                    pick = i;
                            if (!pick)
                                break;
                            const auto &pc = pieces[*pick];
                            bool fwd = pc.nodes[0] == cur;
                            if (pc.nodes[0] == pc.nodes[1])
                            {
                                // a loop: the way the strip runs along the face
                                const auto &ed = m.edge(pc.edge);
                                const auto dir = strip_dir(along_run, ed.curve->value((ed.u1 + ed.u2) / T(2)));
                                const auto x = extrema_curve_point(*ed.curve, dir.first, (ed.u1 + ed.u2) / T(2), T(1e-12));
                                fwd = ed.curve->value(x[0], 1) * dir.second >= T(0);
                            }
                            used.insert({*pick, g});
                            r.push_back(CoEdge<T>{pc.edge, fwd ? Orientation::Forward : Orientation::Reversed, pc.pc.at(g)});
                            cur = pc.nodes[fwd ? 1 : 0];
                            moved = true;
                        }
                        if (cur != b)
                            blend_detail::refuse<T>(op, "the virtual edge does not run along a face given from end to end", {g});
                        return r;
                    };
                    if (nd == n && std::ranges::none_of(cs, [&](const CoEdge<T> &ce) { return spoke_of.contains(ce.edge); }))
                    {
                        // a whole wire along the strip (the foot of a boss on the plate): the loop of the pieces on g
                        std::optional<std::size_t> start;
                        for (std::size_t i = 0; i < pieces.size() && !start; ++i)
                            if (!used.contains({i, g}) && std::ranges::contains(pieces[i].faces, g))
                                start = pieces[i].nodes[0];
                        if (!start)
                            blend_detail::refuse<T>(op, "a face given along the strip without the virtual edge", {g});
                        auto loop = path(*start, *start, cs);
                        // the loop the way the strip ran: its first piece compared to the strip
                        {
                            const auto &ed = m.edge(loop.front().edge);
                            const auto dir = strip_dir(cs, ed.curve->value((ed.u1 + ed.u2) / T(2)));
                            const auto xx = extrema_curve_point(*ed.curve, dir.first, (ed.u1 + ed.u2) / T(2), T(1e-12));
                            const auto d = ed.curve->value(xx[0], 1);
                            if ((loop.front().orient == Orientation::Forward ? d : T(-1) * d) * dir.second < T(0))
                                loop = blend_detail::reversed(loop);
                        }
                        for (const auto &ce : cs)
                            removed.insert(ce.edge);
                        m.wire(wid).coedges = loop;
                        continue;
                    }
                    // the run of doomed co-edges (one per wire), from a co-edge not doomed
                    std::size_t i0 = 0;
                    while (i0 < n && doomed(cs[i0].edge))
                        ++i0;
                    if (i0 == n)
                        blend_detail::refuse<T>(op, "a face kept whose wire runs only along the strip", {g});
                    std::vector<CoEdge<T>> run;
                    for (std::size_t k = 1; k <= n; ++k)
                    {
                        const auto &ce = cs[(i0 + k) % n];
                        if (doomed(ce.edge))
                        {
                            if (!run.empty() && !doomed(cs[(i0 + k - 1) % n].edge))
                                blend_detail::refuse<T>(op, "a face kept along the strip in two places (not handled)", {g});
                            run.push_back(ce);
                        }
                    }
                    std::vector<CoEdge<T>> chain;
                    std::optional<std::size_t> at; // the node reached
                    std::set<EdgeId> dead;
                    std::map<EdgeId, std::size_t> use; // (the uses of a seam, in order)
                    for (const auto &ce : run)
                    {
                        dead.insert(ce.edge);
                        if (!spoke_of.contains(ce.edge))
                            continue;
                        const auto &sp = *spoke_of.at(ce.edge);
                        const auto ends = blend_detail::ends_of(m, ce);
                        const auto &pcs = spoke_pc.at({sp.edge, g});
                        // the pcurve of this use: the one started where this co-edge has its far end
                        const auto &old = m.edge(sp.edge);
                        const auto hint = coedge_uv(m, CoEdge<T>{sp.edge, Orientation::Forward, ce.pcurve}, sp.far == old.v1 ? old.u1 : old.u2);
                        std::size_t best = 0;
                        for (std::size_t j = 1; j < pcs.size(); ++j)
                            if (norm(pcs[j].first - hint) < norm(pcs[best].first - hint))
                                best = j;
                        const auto &ne = m.edge(sp.continued);
                        if (ends[0] == sp.far)
                        {
                            // toward the strip: far -> node
                            chain.push_back(CoEdge<T>{sp.continued, ne.v1 == sp.far ? Orientation::Forward : Orientation::Reversed, pcs[best].second});
                            at = sp.vertex;
                        }
                        else
                        {
                            // away from it: the pieces from the node reached, then node -> far
                            if (at)
                            {
                                auto p = path(*at, sp.vertex, run);
                                chain.insert(chain.end(), p.begin(), p.end());
                            }
                            chain.push_back(CoEdge<T>{sp.continued, ne.v2 == sp.far ? Orientation::Forward : Orientation::Reversed, pcs[best].second});
                            at.reset();
                        }
                    }
                    if (chain.empty())
                        blend_detail::refuse<T>(op, "a face kept along the strip without an edge ending on it", {g});
                    const auto a = blend_detail::ends_of(m, run.front())[0], b = blend_detail::ends_of(m, run.back())[1];
                    if (!blend_detail::replace_run(m, g, chain, a, b, dead, removed))
                        blend_detail::refuse<T>(op, "the virtual edge cannot take the place of the strip in a face kept", {g});
                }
            }
            if (used.size() != 2 * pieces.size())
                blend_detail::refuse<T>(op, "pieces of the virtual edge left over (the faces given meet elsewhere under the strip)");
            // the strip out of the shell, erased with its edges, the old spokes and their vertices on it
            const auto sh = sv.shell_of(*s.f[0].begin());
            std::erase_if(m.shell(sh).faces, [&](const FaceUse &fu) { return s.faces.contains(fu.face); });
            for (auto g : s.faces)
            {
                for (auto wid : m.face(g).wires)
                    m.erase(wid);
                m.erase(g);
            }
            std::set<EdgeId> gone{s.edges.begin(), s.edges.end()};
            for (const auto &sp : spokes)
                gone.insert(sp.edge);
            for (auto e : gone)
                if (m.alive(e))
                    m.erase(e);
            for (auto v : s.vertices)
                if (m.alive(v))
                    m.erase(v);
            std::vector<EdgeId> r;
            for (const auto &pc : pieces)
                r.push_back(pc.edge);
            return r;
        }

        /**
         * The fillet `run` of the virtual edge of the faces f[0] and f[1] of the solid s (§ 2):
         * the strips between them, smallest first; for each, a copy of the solid closed over it
         * (`close_over`), its virtual edge blended (`run(copy, edges)`); kept if every point of
         * the strip is where the virtual edge is (out of the result, in for a concave fillet)
         * and no face of the strip meets the result, the copy then erased (each candidate under
         * try: one that throws is dropped). Otherwise the next strip; none: refused, as spilling
         * over the strip if a blend was built, as no strip if every candidate is the rest of the
         * solid (no end face touching both sides, as at a chamfer ending on another chamfer).
         */
        template <std::floating_point T, typename Run>
        auto through_strip(Model<T> &m, SolidId s, const std::array<std::set<FaceId>, 2> &f, const FilletOptions<T> &o, const std::string &op, Run &&run)
            -> SolidId
        {
            const T tol = o.tol;
            std::vector<Strip<T>> candidates;
            // whether every candidate is the rest of the solid: no strip bounded by end faces that touch both sides
            bool round = true;
            {
                const blend_detail::SolidView<T> sv{m, s};
                candidates = strips(m, sv, f);
                round = std::ranges::all_of(candidates, [&](const Strip<T> &st) { return rest_of_solid(m, sv, st); });
            }
            if (candidates.empty())
                blend_detail::refuse<T>(op, "faces that no strip of faces joins (each face given next to the faces between them, those bounded by faces that touch both sides)");
            std::optional<BuildError> last;
            bool spilled = false;
            for (auto &st : candidates)
            {
                const auto size = detail::arena_size(m);
                try
                {
                    sample(m, st, tol);
                    detail::ShapeCopier<T> copier{m, Placement<T>{}};
                    const auto ns = std::get<SolidId>(copier.copy(ShapeId{s}));
                    Strip<T> q;
                    for (std::size_t k = 0; k < 2; ++k)
                        for (auto g : st.f[k])
                            q.f[k].insert(copier.copy(g));
                    for (auto g : st.faces)
                        q.faces.insert(copier.copy(g));
                    for (auto g : st.ends)
                        q.ends.insert(copier.copy(g));
                    for (auto e : st.edges)
                        q.edges.insert(copier.copy(e));
                    for (auto v : st.vertices)
                        q.vertices.insert(copier.copy(v));
                    q.samples = st.samples;
                    const auto virt = close_over(m, ns, q, o, op);
                    if (const auto report = check(m, ShapeId{ns}); !report.ok())
                        blend_detail::refuse<T>(op, std::string("the faces given continued over the strip into an invalid solid: ") + to_string(report.entries.front().issue) + " " +
                                                         report.entries.front().detail);
                    const auto r = run(ns, std::span<const EdgeId>{virt});
                    // every point of the strip where the virtual edge is, no face of the strip
                    // meeting the result but on an end face (a strip filled in: its edges with the
                    // end faces then on the result)
                    auto on_end = [&](const point<T, 3> &p) {
                        return std::ranges::any_of(st.ends, [&](FaceId g) {
                            const auto &srf = *m.face(g).surface;
                            return extrema_surf_pnt(srf, p, tol / T(100), std::array<std::size_t, 2>{10, 10})[2] <= T(10) * tol;
                        });
                    };
                    auto on_ends = [&](EdgeId e) {
                        const auto &ed = m.edge(e);
                        return ed.degenerate || std::ranges::all_of(std::array{T(0), T(0.25), T(0.5), T(0.75), T(1)},
                                                                    [&](T a) { return on_end(ed.curve->value(ed.u1 + (ed.u2 - ed.u1) * a)); });
                    };
                    const auto &b = m.edge(virt.front());
                    const auto side = classify(m, r, b.curve->value((b.u1 + b.u2) / T(2)), tol);
                    bool swallows = side != Location::On && std::ranges::all_of(st.samples, [&](const point<T, 3> &p) {
                        const auto at = classify(m, r, p, tol);
                        return at == side || (at == Location::On && on_end(p));
                    });
                    for (auto g : st.faces)
                    {
                        if (!swallows)
                            break;
                        const auto before = detail::arena_size(m);
                        SectionOptions<T> so;
                        so.tol = tol;
                        const auto x = section(m, ShapeId{r}, ShapeId{g}, so);
                        swallows = x && x->complete && std::ranges::all_of(x->edges, [&](const SectionEdge<T> &se) { return on_ends(se.edge); }) &&
                                   std::ranges::all_of(x->points, [&](VertexId v) { return on_end(m.vertex(v).pnt); }) &&
                                   std::ranges::all_of(x->boundary_edges, [&](EdgeId e) { return on_ends(e); });
                        detail::rollback(m, before);
                    }
                    if (swallows)
                    {
                        boolean_detail::erase_shape(m, ShapeId{ns});
                        boolean_detail::erase_unused(m, size[1], size[0]);
                        return r;
                    }
                    spilled = true;
                }
                catch (const BuildFailure &e)
                {
                    if (!last)
                        last = e.error();
                }
                catch (const std::exception &e)
                {
                    if (!last)
                        last = BuildError{BuildErrc::NonGeneric, op + ": " + e.what(), {}};
                }
                detail::rollback(m, size);
            }
            if (spilled)
                blend_detail::refuse<T>(op, "the fillet spills over the faces between the faces given (overflow, not handled): a larger radius, or fillet their edges");
            if (round)
                blend_detail::refuse<T>(op, "faces given with no strip between them bounded by end faces that touch both: the faces next to them run round the rest of the "
                                            "solid (at an end, two faces that each touch one side, as a chamfer ending on another chamfer: not handled)");
            throw BuildFailure(*last);
        }

        /// Whether the faces g and h are tangent along their edge e: their outward normals within `angular_tol` of each other at 9 points.
        template <std::floating_point T>
        bool smooth(const Model<T> &m, const blend_detail::SolidView<T> &sv, EdgeId e, FaceId g, FaceId h, T angular_tol)
        {
            const auto &ed = m.edge(e);
            const auto pg = sv.pcurve(g, e), ph = sv.pcurve(h, e);
            if (!pg || !ph || !ed.curve)
                return false;
            auto normal = [&](FaceId f, const Curve<T, 2> &pc, T t) {
                const auto uv = pc.value(t);
                const auto n = normal_direction(*m.face(f).surface, uv[0], uv[1]);
                return sv.orient(f) == Orientation::Forward ? n : T(-1) * n;
            };
            try
            {
                for (int i = 0; i < 9; ++i)
                {
                    const T t = ed.u1 + (ed.u2 - ed.u1) * (T(i) + T(0.5)) / T(9);
                    const auto a = normal(g, *pg, t), b = normal(h, *ph, t);
                    if (!(a * b > T(0) && norm(cross(a, b)) <= angular_tol))
                        return false;
                }
            }
            catch (const std::domain_error &)
            {
                return false; // (no normal there: not shown tangent)
            }
            return true;
        }

        /**
         * The faces given as two sides, each once (§ 2.1): refused (`InvalidId`) where a face is
         * not one of the solid s, where a face is on both sides, and where the faces of a side are
         * not joined along edges where they are tangent (within `angular_tol`).
         */
        template <std::floating_point T>
        auto sides(const Model<T> &m, SolidId s, std::span<const FaceId> f1, std::span<const FaceId> f2, const FilletOptions<T> &o, const std::string &op)
            -> std::array<std::set<FaceId>, 2>
        {
            const std::array<std::set<FaceId>, 2> f{f1 | std::ranges::to<std::set>(), f2 | std::ranges::to<std::set>()};
            const blend_detail::SolidView<T> sv{m, s};
            for (const auto &side : f)
                for (auto g : side)
                    if (!sv.uses.contains(g))
                        fillet_detail::refuse<T>(BuildErrc::InvalidId, op, "a face given that is not a face of the solid", {g});
            for (auto g : f[0])
                if (f[1].contains(g))
                    fillet_detail::refuse<T>(BuildErrc::InvalidId, op, "a face given on both sides", {g});
            for (const auto &side : f)
            {
                // the faces of the side reached from its first one across tangent edges
                std::set<FaceId> reached{*side.begin()};
                std::vector<FaceId> todo{*side.begin()};
                while (!todo.empty())
                {
                    const auto g = todo.back();
                    todo.pop_back();
                    for (auto wid : m.face(g).wires)
                        for (const auto &ce : m.wire(wid).coedges)
                            if (!m.edge(ce.edge).degenerate)
                                for (auto h : sv.users(ce.edge))
                                    if (h != g && side.contains(h) && !reached.contains(h) && smooth(m, sv, ce.edge, g, h, o.angular_tol))
                                        reached.insert(h), todo.push_back(h);
                }
                if (reached.size() != side.size())
                {
                    std::vector<ShapeId> apart;
                    for (auto g : side)
                        if (!reached.contains(g))
                            apart.push_back(g);
                    fillet_detail::refuse<T>(BuildErrc::InvalidId, op,
                                             "faces of a side that are not joined tangentially (a side: a face, or faces joined along edges where they are tangent)",
                                             std::move(apart));
                }
            }
            return f;
        }
    }

    /**
     * @brief The solid `s` with a fillet of radius r between the faces `faces1` and `faces2`,
     * which have no edge in common (stage 8, `docs/sources/design/brep_palier8_face_blends.md`):
     * a new solid, `s` unchanged. Each side is a face, or faces joined tangentially along
     * their common edges. The faces between them (the strip, bounded by end faces that touch
     * both sides) are swallowed: the fillet is that of the sharp edge where both sides meet once
     * continued over the strip (`fillet`, with its ends, chains and options), kept only if the
     * strip lies wholly on the side it takes off, or adds (a concave fillet). Faces of both
     * sides that all have edges in common: the fillet of those edges. Refused (`InvalidId`):
     * a side empty, a face dead, not of the solid or on both sides, the faces of a side not
     * joined tangentially (within `angular_tol`). Refused (`NonGeneric`): faces that meet along
     * edges in places only, no strip joining the sides (or only the rest of the solid, as at a
     * chamfer ending on another chamfer), sides that do not meet under it (parallel faces), a
     * fillet that would spill over the strip. On failure the model is left as it was.
     */
    template <std::floating_point T>
    [[nodiscard]] auto fillet_faces(Model<T> &m, SolidId s, std::span<const FaceId> faces1, std::span<const FaceId> faces2, std::type_identity_t<T> r,
                                    const FilletOptions<T> &o = {}) -> BuildResult<SolidId>
    {
        const std::string op = "fillet_faces";
        return blend_detail::guarded(m, op, [&] {
            if (!m.alive(s))
                fillet_detail::refuse<T>(BuildErrc::InvalidId, op, "dead or invalid solid", {s});
            for (const auto &fs : {faces1, faces2})
            {
                if (fs.empty())
                    fillet_detail::refuse<T>(BuildErrc::InvalidId, op, "no face on a side");
                for (auto g : fs)
                    if (!m.alive(g))
                        fillet_detail::refuse<T>(BuildErrc::InvalidId, op, "dead or invalid face", {g});
            }
            const auto f = face_blend_detail::sides(m, s, faces1, faces2, o, op);
            // faces with edges in common: the fillet of those edges, when every face of both sides
            // has one and no strip joins them elsewhere (not edge to edge in places, across a strip
            // in others)
            {
                const blend_detail::SolidView<T> sv{m, s};
                std::set<EdgeId> common;
                std::set<FaceId> meeting;
                for (auto g : f[0])
                    for (auto wid : m.face(g).wires)
                        for (const auto &ce : m.wire(wid).coedges)
                            if (!m.edge(ce.edge).degenerate)
                                for (auto u : sv.users(ce.edge))
                                    if (f[1].contains(u))
                                        common.insert(ce.edge), meeting.insert(g), meeting.insert(u);
                if (!common.empty())
                {
                    const bool strip = std::ranges::any_of(face_blend_detail::strips(m, sv, f),
                                                           [&](const face_blend_detail::Strip<T> &st) { return !face_blend_detail::rest_of_solid(m, sv, st); });
                    if (strip || meeting.size() != f[0].size() + f[1].size())
                    {
                        std::vector<ShapeId> apart;
                        for (const auto &side : f)
                            for (auto g : side)
                                if (!meeting.contains(g))
                                    apart.push_back(g);
                        blend_detail::refuse<T>(op, "faces given that meet along edges in places only, with faces between them elsewhere (not handled)", std::move(apart));
                    }
                    const std::vector<EdgeId> edges{common.begin(), common.end()};
                    auto x = fillet(m, s, std::span<const EdgeId>{edges}, r, o);
                    if (!x)
                        throw BuildFailure(x.error());
                    return *x;
                }
            }
            return face_blend_detail::through_strip(m, s, f, o, op, [&](SolidId s0, std::span<const EdgeId> e0) {
                return blend_detail::blend(
                    m, s0, e0,
                    [&, sections = blend_detail::sections_on_edges(m, e0, o.section, o.tol, op), original = std::map<EdgeId, EdgeId>{}](
                        SolidId ns, EdgeId e, detail::ShapeCopier<T> &c, const FilletOptions<T> &oe) mutable {
                        if (original.empty())
                            for (auto x : e0)
                                original[c.copy(x)] = x;
                        return fillet_detail::fillet(m, ns, e, r, oe, sections.at(original.at(e)));
                    },
                    o, op, true);
            });
        });
    }

    template <std::floating_point T>
    [[nodiscard]] auto fillet_faces(Model<T> &m, SolidId s, const std::vector<FaceId> &faces1, const std::vector<FaceId> &faces2, std::type_identity_t<T> r,
                                    const FilletOptions<T> &o = {}) -> BuildResult<SolidId>
    {
        return fillet_faces(m, s, std::span<const FaceId>{faces1}, std::span<const FaceId>{faces2}, r, o);
    }
}
