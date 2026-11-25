#pragma once
/**
 * @file full_round.h
 * @brief Full round fillets (stage 8, `docs/sources/design/brep_palier8_face_blends.md`, § 3;
 * PR 3: `docs/sources/design/brep8_pr03_full_round.md`, PR 4: `brep8_pr04_full_round_loops.md`):
 * the middle face M between two faces F1
 * and F2 replaced by the fillet tangent to the three, whose radius follows the geometry (the
 * "Tritangent Fillet" of CATIA, the "Full round fillet" of SolidWorks): the top of a rib rounded
 * into a half cylinder, the blunt trailing edge of a blade rounded, the rim of a disc rounded
 * round a loop.
 *
 * **Geometry.** A ball touching the three faces rolls along M:
 * 1. its **stations** are on the **guide**, the edge between M and F1: in the plane normal to the
 *    guide at its parameter t, the centre c, the radius r and the three contacts are the solution
 *    of c = S_k(u_k, v_k) -+ r N_k(u_k, v_k) for the three faces and (c - g(t)) . g'(t) = 0
 *    (Gauss-Newton on 10 unknowns, from the station next to it; the first from the circle
 *    inscribed in the tangents of the faces in that plane). There is no virtual edge: between
 *    parallel flanks F1 and F2 never meet. The march (PR 4) is that of `fillet_law`: 512 steps
 *    over the guide from its middle, each halved where Gauss-Newton fails or the ball jumps to
 *    another root (its contacts moved further than 4 times what the guide carries, measured on
 *    the step halved: past a thin nose the plane of the station meets M again); faces closing
 *    across a seam without being periodic seen past it (`SeamlessSurface`);
 * 2. the **section** is the arc of the characteristic circle of the ball, the circle through its
 *    three contacts (the ball of a varying radius touches its envelope along that circle, not
 *    along a great circle), from the contact on F1 to the contact on F2 by the contact on M: a
 *    half turn between parallel flanks, more where they come nearer away from M (`arc_through`);
 * 3. the **surface** is read on the stations (PR 4, § 3.7 of the stage), whatever the faces:
 *    the centres on a line and the radius constant, a cylinder; on a line and the radius linear
 *    along it, a cone (its apex where the radius vanishes); on a circle and the radius constant,
 *    a torus. Kept only if, at the stations and between them, the exact arcs lie on it within
 *    `approx_tol` / 2 and its normal within `angular_tol` / 2, the contacts on it within tol;
 *    the pcurves of the contacts interpolated through the stations (in pieces between the knots
 *    of the guide below C3, or where a guide on a surface crosses its knots), the centre line
 *    too. Otherwise the rational B-spline skinned through the sections (`skinned`), within
 *    `approx_tol` of them and its normal within `angular_tol`. The fillet is continued past the
 *    ends of the guide by `extension` times the radius, and by as far as the other edge of M
 *    runs past them (an oblique end).
 * The ball is on the side of the material where M is convex along both its edges with the sides
 * (the fillet takes material off), outside where it is concave along both (the bottom of a slot:
 * the fillet adds material).
 *
 * **Folding** (PR 4, question 18), which `check` does not see, refused at the contacts of 65
 * stations: a face more curved toward the ball than the ball, (a) across, along the arc of the
 * section (the arc would leave the face), (b) along the contact (the centre line would fold
 * back); (b) along the edges too where the march stops (no ball where a face folds).
 *
 * **Topology.** On a copy of the solid, M is taken out and collapsed onto the guide: the edge
 * between M and F2 replaced by the guide in F2 (its pcurve there the old one of that edge,
 * reparametrized on the guide: `guide_on_other`), the edges of M with its end faces reduced to
 * their vertex. M round a loop (PR 4, § 3.5): its two edges with the sides closed, M a ring of
 * two wires or closed by seams; no end face, the blend of a closed edge. The copy has the topology of the solid with a sharp edge in place of M (not its
 * geometry: F2 is a width of M away from that edge), which is what the blends of `blend.h` need:
 * F1 and F2 cut along the contacts, each end face cutting the fillet where it meets it (a free
 * end, `free_end`), the result checked (`check`). The collapsed copy is erased.
 *
 * **One face on both sides** (F1 = F2: a blade of one patch lofted from open sections, a U swept
 * from one curve): the face meets M along two edges, the guide the first. On the copy it is split
 * along the isoparametric line halfway between them (`split_between`: both edges isoparametric
 * in one direction), from where that line crosses the contour of the face to where it crosses it
 * again (an end cut obliquely, the root of a blade fused to its hub), a piece on each side; the
 * pieces are merged again on the result and the edges the cut split joined back at its ends
 * (`merge_pieces`): the result has the faces and the edges of the solid given. Refused: a face
 * closing on itself across a seam, a line crossing the contour other than twice.
 *
 * Refused (`InvalidId`): a face dead or not of the solid, a face both a side and the middle, a
 * middle face that has no edge with a side (or not two with a face on both sides). Refused
 * (`InvalidBounds`): a section other than the arc of circle, setbacks, tolerances not positive.
 * Refused (`NonGeneric`): sides or a middle of several faces, a middle face with holes, meeting a
 * side along several edges, narrowing to a point, with several faces at an end, or ending on a
 * face tangent to a side (a rounded end: a quarter of a sphere, PR 5), all not handled; a loop
 * that does not close, an edge of M closed and the other not; folding (a) or (b);
 * one face on both sides whose edges with M are not isoparametric in one direction (split it);
 * a middle face convex along a side and concave along the other; no ball touching the three
 * faces; a contact out of its face; balls that nest (the middle face narrows faster than the ball
 * moves along it); an end face that does not cut the fillet, or a middle face too narrow at its
 * end for the fillet to meet it. On failure the model is left as it was.
 *
 * Sources: Vida, Martin and Várady, "A survey of blending methods that use parametric surfaces"
 * (Computer-Aided Design 26(5), 1994), § 5 (blends of three surfaces); the envelope of a family
 * of spheres (canal surfaces) for the characteristic circle; the behaviour of CATIA and of
 * SolidWorks.
 */
#include <gbs-brep/face_blend.h>
#include <gbs-brep/split.h>
#include <gbs-brep/unify.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <numbers>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <vector>

namespace gbs::brep
{
    namespace full_round_detail
    {
        /// The three faces of a full round and the edges of the middle one.
        template <std::floating_point T>
        struct Setup
        {
            std::array<FaceId, 3> f;                    ///< F1, M, F2
            std::array<EdgeId, 2> guide;                ///< M with F1 (the guide), M with F2
            std::array<fillet_detail::Side<T>, 3> side; ///< F1 and M along the guide, F2 along the other edge
            bool convex{};                              ///< M convex along both edges: the ball in the material
            std::array<std::vector<EdgeId>, 2> ends;    ///< the edges of M at the end v1 of the guide, at its end v2
            bool same{};                                ///< the other edge runs as the guide (its v1 at the end v1 of the guide)
            bool one{};                                 ///< one face on both sides (F1 = F2: a blade of one patch)
            bool closed{};                              ///< M round a loop (the rim of a disc): its edges with the sides closed, no end face
        };

        /// The ball of a station: its centre, its radius, its contacts on F1, M and F2 and their parameters.
        template <std::floating_point T>
        struct Ball
        {
            point<T, 3> c;
            T r{};
            std::array<point<T, 3>, 3> p;
            std::array<point<T, 2>, 3> uv;
        };

        /**
         * Whether the faces g and h, which meet at p, are tangent there: their outward normals at
         * the feet of p (within 10 tol of it on both) within 1e-6 of each other, the same way.
         */
        template <std::floating_point T>
        bool tangent_at(const Model<T> &m, const blend_detail::SolidView<T> &sv, FaceId g, FaceId h, const point<T, 3> &p, T tol)
        {
            std::array<point<T, 3>, 2> n;
            for (std::size_t k = 0; k < 2; ++k)
            {
                const FaceId f = k == 0 ? g : h;
                const auto &srf = *m.face(f).surface;
                const auto y = extrema_surf_pnt(srf, p, tol / T(100), std::array<std::size_t, 2>{10, 10});
                if (!(y[2] <= T(10) * tol))
                    return false;
                try
                {
                    n[k] = sv.outward(f, srf, point<T, 2>{y[0], y[1]});
                }
                catch (const std::exception &)
                {
                    return false;
                }
            }
            return n[0] * n[1] > T(0) && norm(cross(n[0], n[1])) <= T(1e-6);
        }

        /**
         * The faces given checked, the guide and the ends of the middle face (file header for
         * the refusals).
         */
        template <std::floating_point T>
        auto setup(const Model<T> &m, SolidId s, std::span<const FaceId> faces1, std::span<const FaceId> middle, std::span<const FaceId> faces2,
                   const FilletOptions<T> &o, const std::string &op) -> Setup<T>
        {
            using fillet_detail::refuse;
            if (!m.alive(s))
                refuse<T>(BuildErrc::InvalidId, op, "dead or invalid solid", {s});
            for (const auto &fs : {faces1, middle, faces2})
            {
                if (fs.empty())
                    refuse<T>(BuildErrc::InvalidId, op, "no face on a side or in the middle");
                for (auto g : fs)
                    if (!m.alive(g))
                        refuse<T>(BuildErrc::InvalidId, op, "dead or invalid face", {g});
            }
            fillet_detail::check_options(o, op, EdgeId{});
            if (o.section.kind != SectionKind::Circular)
                refuse<T>(BuildErrc::InvalidBounds, op, "a conic or G2 section (a full round has the arc of circle through its three contacts)");
            if (!o.setbacks.empty())
                refuse<T>(BuildErrc::InvalidBounds, op, "setbacks (a full round ends on the end faces of its middle face)");
            const blend_detail::SolidView<T> sv{m, s};
            for (const auto &fs : {faces1, middle, faces2})
                for (auto g : fs)
                    if (!sv.uses.contains(g))
                        refuse<T>(BuildErrc::InvalidId, op, "a face given that is not a face of the solid", {g});
            if (faces1.size() != 1 || middle.size() != 1 || faces2.size() != 1)
                blend_detail::refuse<T>(op, "sides or a middle of several faces (not handled yet): one face on each side and one in the middle");
            Setup<T> q;
            q.f = {faces1.front(), middle.front(), faces2.front()};
            if (q.f[0] == q.f[1] || q.f[1] == q.f[2])
                refuse<T>(BuildErrc::InvalidId, op, "a face given both as a side and as the middle face", {q.f[1]});
            // one face on both sides (a blade of one patch, a U swept from one curve): its two
            // edges with the middle face, the guide the first of them
            q.one = q.f[0] == q.f[2];
            if (q.one)
            {
                const auto common = blend_detail::common_edges(m, sv, std::array<std::set<FaceId>, 2>{std::set<FaceId>{q.f[1]}, std::set<FaceId>{q.f[0]}});
                if (common.size() != 2)
                    refuse<T>(BuildErrc::InvalidId, op,
                              "the same face on both sides, which does not meet the middle face along two edges (one on each side of it)", {q.f[0], q.f[1]});
                q.guide = {*common.begin(), *std::next(common.begin())};
            }
            else
                for (std::size_t k = 0; k < 2; ++k)
                {
                    const FaceId side = q.f[2 * k];
                    const auto common = blend_detail::common_edges(m, sv, std::array<std::set<FaceId>, 2>{std::set<FaceId>{q.f[1]}, std::set<FaceId>{side}});
                    if (common.empty())
                        refuse<T>(BuildErrc::InvalidId, op, "the middle face has no edge with a side (the middle: the face between the two sides, next to both)",
                                  {q.f[1], side});
                    if (common.size() != 1)
                        blend_detail::refuse<T>(op, "the middle face meets a side along several edges (not handled)", {q.f[1], side});
                    q.guide[k] = *common.begin();
                }
            // M round a loop (§ 3.5 of the stage, PR 4): both its edges with the sides closed; M a
            // ring of two wires (the top of a tube), or one wire closed by seams (the rim of a disc)
            const std::array<bool, 2> loop{m.edge(q.guide[0]).v1 == m.edge(q.guide[0]).v2, m.edge(q.guide[1]).v1 == m.edge(q.guide[1]).v2};
            if (loop[0] != loop[1])
                blend_detail::refuse<T>(op, "a middle face whose edge with one side is closed, with the other not (not handled)", {q.f[1]});
            q.closed = loop[0];
            const auto &mw = m.face(q.f[1]).wires;
            const bool ring = q.closed && mw.size() == 2 && std::ranges::all_of(mw, [&](WireId w) {
                const auto &cs = m.wire(w).coedges;
                return cs.size() == 1 && (cs.front().edge == q.guide[0] || cs.front().edge == q.guide[1]);
            });
            if (mw.size() != 1 && !ring)
                blend_detail::refuse<T>(op, "a middle face with holes (not handled: split it)", {q.f[1]});
            // the sides along their edges with M, and whether M is convex there
            std::array<bool, 2> convex{};
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto st = fillet_detail::setup(m, s, q.guide[k], op);
                convex[k] = st.convex;
                for (const auto &sd : st.sides)
                    if (sd.face != q.f[1])
                        q.side[2 * k] = sd;
                    else if (k == 0)
                        q.side[1] = sd;
            }
            if (convex[0] != convex[1])
                blend_detail::refuse<T>(op, "a middle face convex along a side and concave along the other (no ball touches the three faces on one side of them)",
                                        {q.f[1]});
            q.convex = convex[0];
            // the uses of the guide and of the other edge by M
            std::array<Orientation, 2> use{};
            for (auto w : mw)
                for (const auto &ce : m.wire(w).coedges)
                    for (std::size_t k = 0; k < 2; ++k)
                        if (ce.edge == q.guide[k])
                            use[k] = ce.orient;
            // M runs along the other edge against the guide: as edges, the same way when their uses by M differ
            q.same = use[0] != use[1];
            if (ring)
                return q;
            // the wire of M: the guide, the edges at its end, the other edge, the edges at its start
            const auto &cs = m.wire(mw.front()).coedges;
            const std::size_t n = cs.size();
            std::size_t i0 = n, i1 = n;
            for (std::size_t i = 0; i < n; ++i)
            {
                if (cs[i].edge == q.guide[0])
                    i0 = i;
                if (cs[i].edge == q.guide[1])
                    i1 = i;
            }
            const bool fwd = cs[i0].orient == Orientation::Forward;
            for (std::size_t i = (i0 + 1) % n; i != i1; i = (i + 1) % n)
                q.ends[fwd ? 1 : 0].push_back(cs[i].edge);
            for (std::size_t i = (i1 + 1) % n; i != i0; i = (i + 1) % n)
                q.ends[fwd ? 0 : 1].push_back(cs[i].edge);
            if (q.closed)
            {
                // (between the two closed edges, the seams of M only: no end face)
                for (const auto &run : q.ends)
                    for (auto e : run)
                        if (std::ranges::count_if(cs, [&](const CoEdge<T> &x) { return x.edge == e; }) != 2)
                            blend_detail::refuse<T>(op, "a middle face round a loop with an edge between its two closed edges that is not a seam (not handled)",
                                                    {q.f[1], e});
                return q;
            }
            for (const auto &run : q.ends)
            {
                if (std::ranges::all_of(run, [&](EdgeId e) { return m.edge(e).degenerate; }))
                    blend_detail::refuse<T>(op, "a middle face narrowing to a point (its sides meet at an end of it: not handled; fillet_faces takes such a face)",
                                            {q.f[1]});
                if (run.size() != 1)
                    blend_detail::refuse<T>(op, "several faces at an end of the middle face (a rib ending in a prow: not handled)", {q.f[1]});
            }
            // an end face tangent to a side where it meets it at the end of M (the rounded end of a
            // rib, a half cylinder tangent to its flanks): the full round would end in a quarter of
            // a sphere, which is not built
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto e = q.ends[k].front();
                const auto us = sv.users(e);
                const FaceId end = us.size() == 2 ? (us[0] == q.f[1] ? us[1] : us[0]) : FaceId{};
                if (!end.valid())
                    continue;
                for (const VertexId v : {m.edge(e).v1, m.edge(e).v2})
                    for (const FaceId side : {q.f[0], q.f[2]})
                        if (tangent_at(m, sv, end, side, m.vertex(v).pnt, o.tol))
                            blend_detail::refuse<T>(op,
                                                    "an end face tangent to a side (the rounded end of a rib): not handled, the full round would end "
                                                    "there in a quarter of a sphere; end the middle face on a face across the sides",
                                                    {q.f[1], end, side});
            }
            return q;
        }

        /// The direction of the side sd away from its edge at t, in the plane across the unit tangent te of the edge there.
        template <std::floating_point T>
        auto away(const fillet_detail::Side<T> &sd, T t, const point<T, 3> &te) -> point<T, 3>
        {
            const auto uv = sd.pcurve->value(t);
            const auto su = sd.surface->value(uv[0], uv[1], 1, 0), sv = sd.surface->value(uv[0], uv[1], 0, 1);
            // the face is on the left of its co-edge in (u, v)
            auto d = sd.pcurve->value(t, 1);
            if (!sd.forward)
                d = T(-1) * d;
            auto b = (-d[1]) * su + d[0] * sv;
            b = b - (b * te) * te;
            return fillet_detail::unit(b);
        }

        /// The step of the parameters of the surface s at uv that moves its point by d (least squares in its tangent plane).
        template <std::floating_point T>
        auto moved(const Surface<T, 3> &s, const point<T, 2> &uv, const point<T, 3> &d) -> point<T, 2>
        {
            const auto su = s.value(uv[0], uv[1], 1, 0), sv = s.value(uv[0], uv[1], 0, 1);
            const T a = su * su, b = su * sv, c = sv * sv, det = a * c - b * b;
            if (!(det > T(1e-24) * a * c))
                return uv;
            const T x = su * d, y = sv * d;
            return {uv[0] + (c * x - b * y) / det, uv[1] + (a * y - b * x) / det};
        }

        /**
         * The normal curvature of the surface s at uv in the direction d (in space, projected on
         * its tangent plane), positive where it bends toward `toward` (a unit vector off the
         * surface: toward the ball). Zero where the direction or the frame degenerates.
         */
        template <std::floating_point T>
        auto bending(const Surface<T, 3> &s, const point<T, 2> &uv, const point<T, 3> &d, const point<T, 3> &toward) -> T
        {
            const auto su = s.value(uv[0], uv[1], 1, 0), sv = s.value(uv[0], uv[1], 0, 1);
            const T A = su * su, B = su * sv, C = sv * sv, det = A * C - B * B;
            if (!(det > T(1e-24) * A * C))
                return T(0);
            const T x = su * d, y = sv * d;
            const T w0 = (C * x - B * y) / det, w1 = (A * y - B * x) / det;
            const auto dd = w0 * su + w1 * sv;
            if (!(dd * dd > T(0)))
                return T(0);
            const auto suu = s.value(uv[0], uv[1], 2, 0), suv = s.value(uv[0], uv[1], 1, 1), svv = s.value(uv[0], uv[1], 0, 2);
            return ((w0 * w0) * suu + (T(2) * w0 * w1) * suv + (w1 * w1) * svv) * toward / (dd * dd);
        }

        /**
         * The geometry of the full round of the setup q (file header): an `EdgeBlend` on the
         * parameter of the guide, its contacts on F1 and F2, its spine the centres of the balls.
         */
        template <std::floating_point T>
        auto geometry(const Model<T> &m, const Setup<T> &q, const FilletOptions<T> &o, const std::string &op) -> EdgeBlend<T>
        {
            using fillet_detail::unit;
            const EdgeId eid = q.guide[0];
            const auto &ed = m.edge(eid), &ed2 = m.edge(q.guide[1]);
            const T tol = o.tol;
            const T sgn = q.convex ? T(-1) : T(1);
            std::array<T, 3> sign{};
            for (std::size_t k = 0; k < 3; ++k)
                sign[k] = sgn * (q.side[k].orient == Orientation::Forward ? T(1) : T(-1));
            auto no_ball = [&] [[noreturn]] () {
                blend_detail::refuse<T>(op,
                                        "no ball touches the three faces at a station of the middle face (sides that part from it, or a middle face too wide "
                                        "for them): fillet its two edges with the sides instead",
                                        {q.f[1]});
            };

            // ---- the circle inscribed in the tangents of the three faces in the plane across
            // the guide at t (within the edge): the start of the first station
            auto seed = [&](T t, const std::array<const Surface<T, 3> *, 3> &srf) -> std::array<T, 10> {
                t = std::clamp(t, ed.u1, ed.u2);
                const auto p = ed.curve->value(t);
                const auto te = unit(ed.curve->value(t, 1));
                // the point of the other edge in that plane: from its nearest of 33 points
                T t2 = ed2.u1, gap = std::numeric_limits<T>::infinity();
                for (int i = 0; i <= 32; ++i)
                {
                    const T u = ed2.u1 + (ed2.u2 - ed2.u1) * T(i) / T(32);
                    if (const T d = norm(ed2.curve->value(u) - p); d < gap)
                        gap = d, t2 = u;
                }
                for (int i = 0; i < 30; ++i)
                {
                    const T h = (ed2.curve->value(t2) - p) * te, dh = ed2.curve->value(t2, 1) * te;
                    if (!(std::abs(dh) > T(1e-12) * norm(ed2.curve->value(t2, 1))))
                        break;
                    const T next = std::clamp(t2 - h / dh, ed2.u1, ed2.u2);
                    if (std::abs(next - t2) <= T(1e-14) * (T(1) + std::abs(t2)))
                        break;
                    t2 = next;
                }
                const auto pq = ed2.curve->value(t2) - p;
                auto d = pq - (pq * te) * te;
                const T w = norm(d);
                if (!(w > tol))
                    no_ball();
                d = d / w;
                const auto b1 = away(q.side[0], t, te), b2 = away(q.side[2], t2, te);
                const T t1h = std::tan(std::acos(std::clamp(d * b1, T(-1), T(1))) / T(2)), t2h = std::tan(std::acos(std::clamp(T(-1) * (d * b2), T(-1), T(1))) / T(2));
                const T x = w * t2h / (t1h + t2h), r = w * t1h * t2h / (t1h + t2h);
                // (the centre on the side of M where the sides go)
                auto nin = b1 - (b1 * d) * d;
                if (!(norm(nin) > T(1e-9)) || !(r > tol) || !std::isfinite(r))
                    no_ball();
                nin = unit(nin);
                const auto c = p + x * d + r * nin;
                const auto uv1 = moved(*srf[0], q.side[0].pcurve->value(t), x * b1);
                const auto uvm = moved(*srf[1], q.side[1].pcurve->value(t), x * d);
                const auto uv2 = moved(*srf[2], q.side[2].pcurve->value(t2), (w - x) * b2);
                return {uv1[0], uv1[1], uvm[0], uvm[1], uv2[0], uv2[1], c[0], c[1], c[2], r};
            };
            const T tm = (ed.u1 + ed.u2) / T(2);
            const auto x0 = seed(tm, {q.side[0].surface.get(), q.side[1].surface.get(), q.side[2].surface.get()});
            const T r0 = x0[9];
            const T width = norm(q.side[1].surface->value(x0[2], x0[3]) - ed.curve->value(tm)) * T(2);

            // ---- the supports and the guide, continued past the ends of the middle face: by the
            // extension, and as far as the other edge of M runs past the ends of the guide (an
            // oblique end, the guide on its short side)
            T overhang{};
            for (const bool start : {true, false})
            {
                const T u = start ? ed.u1 : ed.u2;
                const auto out = (start ? T(-1) : T(1)) * unit(ed.curve->value(u, 1));
                const VertexId w = (start == q.same) ? ed2.v1 : ed2.v2;
                overhang = std::max(overhang, (m.vertex(w).pnt - ed.curve->value(u)) * out);
            }
            // (a loop: none, no end face to cut)
            const T ext = q.closed ? T(0) : o.extension * r0 + overhang;
            const T reach = T(2) * (r0 + width) + ext;
            std::array<fillet_detail::Support<T>, 3> sup;
            for (std::size_t k = 0; k < 3; ++k)
            {
                sup[k] = fillet_detail::support(m, q.side[k], k == 2 ? ed2 : ed, reach, tol);
                // a face closing across a seam without being periodic (a closed B-spline: the
                // profile of a fin) seen past its seam, as the supports of the fillets of edges
                if (!dynamic_cast<const ElementarySurface<T> *>(sup[k].surface.get()))
                    if (const auto per = fillet_detail::unrepeated(*sup[k].surface, fillet_detail::periods_of(*sup[k].surface, tol), tol); per[0] > T(0) || per[1] > T(0))
                        sup[k].surface = std::make_shared<SeamlessSurface<T>>(sup[k].surface, per);
            }
            T lo{}, hi{};
            const auto guide = fillet_detail::continued_edge(ed, ext, std::array<fillet_detail::Side<T>, 2>{q.side[0], q.side[1]}, {sup[0].surface, sup[1].surface}, tol, lo, hi);
            const T v1 = norm(ed.curve->value(ed.u1, 1)), v2 = norm(ed.curve->value(ed.u2, 1));
            T a = std::max(lo, ed.u1 - ext / v1), b = std::min(hi, ed.u2 + ext / v2);

            // ---- the stations
            const T scale = T(1) + norm(ed.curve->value(tm)) + r0;
            auto solve = [&](T t, const std::array<T, 10> &start) -> std::optional<std::array<T, 10>> {
                const auto e = guide->value(t);
                const auto te = unit(guide->value(std::clamp(t, lo + T(1e-9) * (hi - lo), hi - T(1e-9) * (hi - lo)), 1));
                thick_detail::Vec<T> x(10);
                for (int i = 0; i < 10; ++i)
                    x[i] = start[std::size_t(i)];
                const T res = thick_detail::gauss_newton<T>(x, [&](const thick_detail::Vec<T> &y, thick_detail::Vec<T> &f, thick_detail::Mat<T> &J) {
                    f.resize(10);
                    J = thick_detail::Mat<T>::Zero(10, 10);
                    const point<T, 3> c{y[6], y[7], y[8]};
                    const T r = y[9];
                    for (int k = 0; k < 3; ++k)
                    {
                        const auto &srf = sup[std::size_t(k)].surface;
                        const T sg = sign[std::size_t(k)];
                        const OffsetSurface<T> off{srf, sg * r};
                        std::array<point<T, 3>, 4> D;
                        // (derivatives: [0] value, [1] d/dv, [2] d/du, [3] d2/dudv)
                        off.derivatives(y[2 * k], y[2 * k + 1], 1, 1, D.data());
                        const auto N = normal_direction(*srf, y[2 * k], y[2 * k + 1]);
                        for (int i = 0; i < 3; ++i)
                        {
                            const auto row = 3 * k + i;
                            f[row] = D[0][std::size_t(i)] - c[std::size_t(i)];
                            J(row, 2 * k) = D[2][std::size_t(i)];
                            J(row, 2 * k + 1) = D[1][std::size_t(i)];
                            J(row, 6 + i) = T(-1);
                            J(row, 9) = sg * N[std::size_t(i)];
                        }
                    }
                    f[9] = (c - e) * te;
                    for (int i = 0; i < 3; ++i)
                        J(9, 6 + i) = te[std::size_t(i)];
                });
                if (!(res <= T(1e-12) * scale) || !(x[9] > tol))
                    return std::nullopt;
                std::array<T, 10> r;
                for (int i = 0; i < 10; ++i)
                    r[std::size_t(i)] = x[i];
                return r;
            };
            auto ball_of = [&](const std::array<T, 10> &x) {
                Ball<T> bl;
                bl.c = {x[6], x[7], x[8]};
                bl.r = x[9];
                for (std::size_t k = 0; k < 3; ++k)
                {
                    bl.uv[k] = {x[2 * k], x[2 * k + 1]};
                    bl.p[k] = sup[k].surface->value(x[2 * k], x[2 * k + 1]);
                }
                return bl;
            };
            // a station from one solved, the step halved where Gauss-Newton fails or the ball jumps
            // to another root (its contacts moved further than the guide carries them: 4 times its
            // speed, on the step divided, not growing with the step; past a thin nose the plane of
            // the station meets M again, and a ball tangent to the three faces there is a solution
            // too, PR 4, § 3.5)
            std::map<T, std::array<T, 10>> solved;
            auto jumped = [&](const std::array<T, 10> &x0, T t0, const std::array<T, 10> &x1, T t1) {
                const auto b0 = ball_of(x0), b1 = ball_of(x1);
                const T speed = std::max(norm(guide->value(t0, 1)), norm(guide->value(t1, 1)));
                T far{};
                for (std::size_t k = 0; k < 3; ++k)
                    far = std::max(far, norm(b1.p[k] - b0.p[k]));
                return far > T(4) * speed * std::abs(t1 - t0) + T(10) * tol;
            };
            auto advance = [&](auto &&self, T t0, T t1, int depth) -> bool {
                const auto &x0 = solved.at(t0);
                const auto x = solve(t1, x0);
                if (x && !jumped(x0, t0, *x, t1))
                {
                    solved.insert_or_assign(t1, *x);
                    return true;
                }
                if (depth == 0)
                    return false;
                const T mid = (t0 + t1) / T(2);
                return self(self, t0, mid, depth - 1) && self(self, mid, t1, depth - 1);
            };
            // (a station asked by the skin or a check: from the nearest station of the march)
            auto find = [&](T t) -> const std::array<T, 10> * {
                if (const auto it = solved.find(t); it != solved.end())
                    return &it->second;
                if (solved.empty())
                    return nullptr;
                auto it = solved.lower_bound(t);
                if (it == solved.end() || (it != solved.begin() && t - std::prev(it)->first < it->first - t))
                    --it;
                if (!advance(advance, it->first, t, 8))
                    return nullptr;
                return &solved.at(t);
            };
            auto raw = [&](T t) -> const std::array<T, 10> & {
                const auto *x = find(t);
                if (!x)
                    no_ball();
                return *x;
            };
            // (the march stopped, no ball: first the folding of a face at the radius of the first
            // ball, which leaves no ball where it folds (the fin thicker than the nose of its edge
            // is round, question 18 (b)); otherwise no ball)
            const char *face_names[3] = {"the first side", "the middle face", "the second side"};
            auto stopped = [&] [[noreturn]] () {
                for (std::size_t k = 0; k < 3; ++k)
                {
                    const auto &sd = q.side[k];
                    const auto &e = m.edge(k == 2 ? q.guide[1] : q.guide[0]);
                    for (int i = 0; i <= 2048; ++i)
                    {
                        const T t = e.u1 + (e.u2 - e.u1) * T(i) / T(2048);
                        const auto uv = sd.pcurve->value(t);
                        const auto toward = sign[k] * normal_direction(*sd.surface, uv[0], uv[1]);
                        if (bending(*sd.surface, uv, unit(e.curve->value(t, 1)), toward) > (T(1) + T(1e-6)) / r0)
                            blend_detail::refuse<T>(op,
                                                    std::string("the radius of the full round (") + std::to_string(r0) + ") exceeds a radius of curvature of " +
                                                        face_names[k] + " along the guide: its centre line would fold back",
                                                    {q.f[k]});
                    }
                }
                no_ball();
            };
            // the march of `fillet_law` (stage 5): 512 steps over the edge from its middle, each
            // divided down to 1/1024 where it must; then past the vertices of the guide as far as
            // there is a ball (a middle face narrowing toward its end: none past the point where
            // its sides meet): the blend ends at the last station solved. A loop: round it, its
            // stations at both ends the same
            {
                const auto x = solve(tm, seed(tm, {sup[0].surface.get(), sup[1].surface.get(), sup[2].surface.get()}));
                if (!x)
                    no_ball();
                solved.emplace(tm, *x);
            }
            constexpr int nm = 256;
            for (const bool up : {true, false})
            {
                const T end = up ? ed.u2 : ed.u1, limit = up ? b : a;
                T prev = tm;
                for (int i = 1; i <= nm; ++i)
                {
                    const T t = i == nm ? end : tm + (end - tm) * T(i) / T(nm);
                    if (!advance(advance, prev, t, 10))
                        stopped();
                    prev = t;
                }
                T last = end;
                for (int i = 1; i <= 32 && limit != end; ++i)
                {
                    const T t = i == 32 ? limit : end + (limit - end) * T(i) / T(32);
                    if (!advance(advance, last, t, 6))
                        break;
                    last = t;
                }
                (up ? b : a) = last;
            }
            if (q.closed)
            {
                const auto b0 = ball_of(solved.at(a)), b1 = ball_of(solved.at(b));
                T gap = norm(b1.c - b0.c);
                for (std::size_t k = 0; k < 3; ++k)
                    gap = std::max(gap, norm(b1.p[k] - b0.p[k]));
                if (!(gap <= T(10) * tol))
                    blend_detail::refuse<T>(op, "the full round round a loop does not close (its balls at the start and at the end of the loop differ)", {q.f[1]});
            }
            constexpr int ns = 64;
            auto at = [&](int i) { return i == 0 ? a : i == ns ? b : a + (b - a) * T(i) / T(ns); };
            // the balls nest where the radius varies faster than the centre moves: no envelope there
            for (int i = 0; i < ns; ++i)
            {
                const auto b0 = ball_of(raw(at(i))), b1 = ball_of(raw(at(i + 1)));
                if (!(std::abs(b1.r - b0.r) < norm(b1.c - b0.c) * (T(1) - T(1e-6))))
                    blend_detail::refuse<T>(op,
                                            "the balls of the full round nest (the middle face narrows or widens faster than the ball moves along it): "
                                            "fillet its two edges with the sides instead",
                                            {q.f[1]});
            }
            auto station = [&](T t) {
                const auto bl = ball_of(raw(t));
                // the characteristic circle: through the three contacts
                const auto circle = fillet_detail::circle3(bl.p[0], bl.p[1], bl.p[2]);
                if (!circle)
                    blend_detail::refuse<T>(op, "the contacts of a ball of the full round are aligned", {q.f[1]});
                fillet_detail::Station<T> st;
                st.c = bl.c;
                st.p1 = bl.p[0], st.p2 = bl.p[2];
                st.uv1 = bl.uv[0], st.uv2 = bl.uv[2];
                st.m = std::get<0>(*circle);
                st.via = bl.p[1];
                return st;
            };

            EdgeBlend<T> rb;
            rb.edge = eid;
            rb.convex = q.convex;
            rb.a = a, rb.b = b;
            rb.t1 = ed.u1, rb.t2 = ed.u2;
            rb.closed = q.closed;
            rb.circular = true;
            rb.contacts[0].face = q.f[0], rb.contacts[1].face = q.f[2];
            rb.contacts[0].support = sup[0].surface, rb.contacts[1].support = sup[2].surface;

            // ---- folding, which `check` does not see (question 18), at the contacts of 65
            // stations: (a) across, a face more curved toward the ball than the ball along the arc
            // of the section (the arc leaves the face); (b) along, a face more curved toward the
            // ball than the ball along the contact (the centre line folds back, the surface cuts
            // itself; the offset of `fillet_law` folds there). The normal curvatures at the
            // contacts, not the offset of the whole face: a U whose walls are the sides folds at its
            // bottom, far from the fillet
            {
                const auto &names = face_names;
                for (int i = 0; i <= ns; ++i)
                {
                    const auto bl = ball_of(raw(at(i)));
                    const auto prev = ball_of(raw(at(std::max(0, i - 1)))), next = ball_of(raw(at(std::min(ns, i + 1))));
                    const auto circle = fillet_detail::circle3(bl.p[0], bl.p[1], bl.p[2]);
                    for (std::size_t k = 0; k < 3; ++k)
                    {
                        const auto &srf = *sup[k].surface;
                        const auto toward = unit(bl.c - bl.p[k]);
                        // (a) across: along the arc of the section at the contact
                        if (circle)
                            if (const auto d = cross(std::get<1>(*circle), bl.p[k] - std::get<0>(*circle)); norm(d) > T(0))
                                if (const T kappa = bending(srf, bl.uv[k], unit(d), toward); kappa > (T(1) + T(1e-6)) / bl.r)
                                    blend_detail::refuse<T>(op,
                                                            std::string("the ball of the full round is less curved than ") + names[k] + " across, at its contact (radius of "
                                                                "curvature " + std::to_string(T(1) / kappa) + " along the arc of the section, the ball " + std::to_string(bl.r) +
                                                                "): the arc would leave the face",
                                                            {q.f[k]});
                        // (b) along: along the contact, from station to station
                        if (const auto d = next.p[k] - prev.p[k]; norm(d) > T(0))
                            if (const T kappa = bending(srf, bl.uv[k], unit(d), toward); kappa > (T(1) + T(1e-6)) / bl.r)
                                blend_detail::refuse<T>(op,
                                                        std::string("the radius of the full round (") + std::to_string(bl.r) + ") exceeds a radius of curvature of " +
                                                            names[k] + " along the guide (" + std::to_string(T(1) / kappa) + "): its centre line would fold back",
                                                        {q.f[k]});
                    }
                }
            }

            // ---- an exact surface from the stations (§ 3.7 of the stage): the centres on a line or
            // on a circle, the radius constant or linear along the line: a cylinder, a cone, a
            // torus, whatever the faces; kept if, at the stations and between them, the exact arcs
            // lie on it and its normal is that of the envelope, and the contacts drawn on it lie on
            // their faces. The pcurves of the contacts on their faces by interpolation of the
            // stations (projected: the parameters of the contacts solved), the centre line too
            constexpr int nc = 16;
            std::vector<Ball<T>> balls;
            for (int i = 0; i <= nc; ++i)
                balls.push_back(ball_of(raw(i == nc ? b : a + (b - a) * T(i) / T(nc))));
            const auto &ba = balls.front(), &bb = balls.back();
            const T r_lo = std::ranges::min(balls, {}, &Ball<T>::r).r, r_hi = std::ranges::max(balls, {}, &Ball<T>::r).r;
            const bool constant = r_hi - r_lo <= tol;
            // (the point of a curve f through the stations, its parameters continuous across the periods per)
            // (in pieces between the knots of the guide below C3, joined: the speed of a guide from
            // a revolution, a rational circle, kinks there, the stations no smoother)
            std::vector<T> cuts{a};
            {
                std::vector<T> ks = fillet_detail::edge_weak_knots(*guide, 3);
                // (a curve on a surface, the edge of a revolution: where its pcurve crosses a knot
                // of the surface below C3, the arcs of a rational circle)
                if (const auto *cs = dynamic_cast<const CurveOnSurface<T, 3> *>(guide.get()))
                    for (std::size_t d = 0; d < 2; ++d)
                        for (const T x : fillet_detail::weak_knots(*cs->p_basisSurface(), d))
                            for (const T t : fillet_detail::crossings(*cs->p_basisCurve(), a, b, d, x, fillet_detail::view_periods(*cs->p_basisSurface())[d]))
                                ks.push_back(t);
                for (const T k : ks)
                    if (k > a + T(1e-9) * (b - a) && k < b - T(1e-9) * (b - a) && std::ranges::none_of(cuts, [&](T c) { return std::abs(c - k) <= T(1e-9) * (b - a); }))
                        cuts.push_back(k);
            }
            std::ranges::sort(cuts);
            cuts.push_back(b);
            auto fit = [&]<std::size_t dim>(auto &&f, const std::array<T, 2> &per, T within) -> std::shared_ptr<BSCurve<T, dim>> {
                std::vector<BSCurve<T, dim>> parts;
                std::optional<point<T, dim>> last;
                for (std::size_t piece = 0; piece + 1 < cuts.size(); ++piece)
                {
                    const T lo = cuts[piece], hi = cuts[piece + 1];
                    std::optional<BSCurve<T, dim>> got;
                    for (std::size_t n = 8; n <= 1024 && !got; n *= 2)
                    {
                        points_vector<T, dim> q;
                        std::vector<T> ts;
                        for (std::size_t i = 0; i <= n; ++i)
                        {
                            const T t = i == n ? hi : lo + (hi - lo) * T(i) / T(n);
                            auto p = f(t);
                            if constexpr (dim == 2)
                            {
                                const auto ref = q.empty() ? last : std::optional<point<T, dim>>{q.back()};
                                if (ref)
                                    for (std::size_t d = 0; d < 2; ++d)
                                        if (per[d] > T(0))
                                            p[d] += per[d] * std::round(((*ref)[d] - p[d]) / per[d]);
                            }
                            q.push_back(p);
                            ts.push_back(t);
                        }
                        BSCurve<T, dim> c = interpolate(q, ts, 3);
                        T dev{};
                        for (std::size_t i = 0; i < n; ++i)
                        {
                            const T t = (ts[i] + ts[i + 1]) / T(2);
                            auto p = f(t);
                            const auto cv = c.value(t);
                            if constexpr (dim == 2)
                                for (std::size_t d = 0; d < 2; ++d)
                                    if (per[d] > T(0))
                                        p[d] += per[d] * std::round((cv[d] - p[d]) / per[d]);
                            dev = std::max(dev, norm(cv - p));
                        }
                        if (dev <= within)
                            got = c;
                    }
                    if (!got)
                        return nullptr;
                    last = got->value(hi);
                    parts.push_back(*got);
                }
                return std::make_shared<BSCurve<T, dim>>(parts.size() == 1 ? parts.front() : fillet_detail::joined(parts));
            };
            // the distance from p to the primitive and its normal there (unit), for the checks
            struct Primitive
            {
                int kind{}; // 1 cylinder, 2 cone, 3 torus
                point<T, 3> o, z;
                T R{}, r{}, alpha{};
            };
            auto measure = [&](const Primitive &pr, const point<T, 3> &p) -> std::pair<T, point<T, 3>> {
                const auto q = p - pr.o;
                const auto along = q * pr.z;
                auto rad = q - along * pr.z;
                const T rr = norm(rad);
                const auto e = rr > T(0) ? rad / rr : point<T, 3>{};
                if (pr.kind == 1)
                    return {std::abs(rr - pr.r), e};
                if (pr.kind == 2)
                {
                    const auto n = std::cos(pr.alpha) * e - std::sin(pr.alpha) * pr.z;
                    return {std::abs(q * n), n};
                }
                const auto ring = pr.o + pr.R * e;
                const auto w = p - ring;
                return {std::abs(norm(w) - pr.r), unit(w)};
            };
            // the line or the circle of the centres
            std::optional<Primitive> prim;
            {
                const auto ca = ba.c, cb = bb.c;
                const auto axis = cb - ca;
                bool on_line = norm(axis) > T(1e-9) * scale;
                for (const auto &x : balls)
                    on_line = on_line && norm(cross(x.c - ca, unit(axis))) <= tol;
                if (on_line && constant)
                    prim = Primitive{1, ca, unit(axis), T(0), (r_lo + r_hi) / T(2), T(0)};
                else if (on_line)
                {
                    // the radius linear along the line: a cone of apex where it vanishes
                    const T ds = norm(axis), dr = bb.r - ba.r, sin_a = std::abs(dr) / ds;
                    bool lin = sin_a < T(1) - T(1e-9);
                    for (const auto &x : balls)
                        lin = lin && std::abs(x.r - (ba.r + dr * ((x.c - ca) * unit(axis)) / ds)) <= tol;
                    if (lin)
                    {
                        const auto z = dr > T(0) ? unit(axis) : T(-1) * unit(axis);
                        const auto &near = dr > T(0) ? ba : bb;
                        prim = Primitive{2, near.c - (near.r / sin_a) * z, z, T(0), T(0), std::asin(sin_a)};
                    }
                }
                else if (constant)
                {
                    const auto &cm = balls[nc / 2];
                    if (const auto cc = fillet_detail::circle3(ca, cm.c, q.closed ? balls[nc / 4].c : cb); cc)
                    {
                        auto [o3, z3, R3] = *cc;
                        // (the axis such that the centres turn positively round it along the guide)
                        if (cross(ca - o3, balls[1].c - o3) * z3 < T(0))
                            z3 = T(-1) * z3;
                        bool on_circle = true;
                        for (const auto &x : balls)
                            on_circle = on_circle && std::abs((x.c - o3) * z3) <= tol && std::abs(norm(x.c - o3) - R3) <= tol;
                        if (on_circle)
                            prim = Primitive{3, o3, z3, R3, (r_lo + r_hi) / T(2), T(0)};
                    }
                }
            }
            // the checks between the stations (and at them): the arcs exact on the primitive
            bool exact = prim.has_value();
            for (int i = 0; exact && i <= 2 * ns; ++i)
            {
                const auto st = station(a + (b - a) * T(i) / T(2 * ns));
                const auto arc = fillet_detail::arc_through(*st.m, st.p1, st.p2, *st.via);
                for (int j = 0; j <= 8 && exact; ++j)
                {
                    const auto p = arc.value(T(j) / T(8));
                    const auto [d, n] = measure(*prim, p);
                    const auto ne = unit(p - st.c);
                    exact = d <= o.approx_tol / T(2) && norm(cross(n, ne)) <= o.angular_tol / T(2);
                }
            }
            if (exact)
            {
                // the contacts: their pcurves on their faces through the stations, on the primitive within tol
                std::array<std::shared_ptr<BSCurve<T, 2>>, 2> pc;
                for (std::size_t k = 0; k < 2 && exact; ++k)
                {
                    // (the parameters within 1e-9 of those of the stations: the contacts on their faces exactly)
                    pc[k] = fit.template operator()<2>([&](T t) { return ball_of(raw(t)).uv[2 * k]; }, fillet_detail::view_periods(*sup[2 * k].surface),
                                                       T(1e-9) * (T(1) + norm(ba.uv[2 * k])));
                    exact = pc[k] != nullptr;
                    for (int i = 0; exact && i <= 4 * ns; ++i)
                    {
                        const T t = a + (b - a) * T(i) / T(4 * ns);
                        const auto uv = pc[k]->value(t);
                        exact = measure(*prim, sup[2 * k].surface->value(uv[0], uv[1])).first <= tol;
                    }
                }
                // (the centre line within tol: it locates the balls, the surface is the primitive; the
                // speed of a guide from a revolution may kink, its centres then no smoother)
                const auto spine = exact ? fit.template operator()<3>([&](T t) { return ball_of(raw(t)).c; }, std::array<T, 2>{}, tol) : nullptr;
                exact = exact && spine;
                if (exact)
                {
                    const auto &pr = *prim;
                    // the frame: x toward the contact with M at the first station (the seam opposite it)
                    auto dx = pr.kind == 3 ? ba.c - pr.o : ba.p[1] - (pr.kind == 2 ? pr.o : ba.c);
                    dx = unit(dx - (dx * pr.z) * pr.z);
                    const Frame3<T> f{ax2<T, 3>{pr.o, pr.z, dx}};
                    constexpr T pi = std::numbers::pi_v<T>;
                    T u1 = std::numeric_limits<T>::infinity(), u2 = -u1, v1 = u1, v2 = -u1;
                    if (pr.kind != 3)
                    {
                        // the angle round the axis of the contacts (both sides of the contact with M), v along the axis
                        for (const auto &x : balls)
                            for (std::size_t k : {std::size_t{0}, std::size_t{2}})
                            {
                                const auto w = x.p[k] - pr.o;
                                const T u = std::atan2(w * f.y, w * f.x);
                                u1 = std::min(u1, u), u2 = std::max(u2, u);
                                const T v = w * pr.z;
                                v1 = std::min(v1, v), v2 = std::max(v2, v);
                            }
                        if (!(u1 > -pi + T(1e-6) && u2 < pi - T(1e-6)))
                            exact = false;
                        if (exact && pr.kind == 1)
                            rb.surface = std::make_shared<CylindricalSurface<T>>(pr.r, Frame3<T>{ax2<T, 3>{ba.c, pr.z, dx}},
                                                                                 std::array<T, 4>{u1, u2, v1 - (ba.c - pr.o) * pr.z, v2 - (ba.c - pr.o) * pr.z});
                        else if (exact)
                            rb.surface = std::make_shared<ConicalSurface<T>>(T(0), pr.alpha, f, std::array<T, 4>{u1, u2, v1, v2});
                    }
                    else
                    {
                        // u the angle of the centres round the axis (continuous along the guide: the
                        // axis turned to make it grow), v the angle of the contacts round the circle of
                        // the centres, on the side of the contact with M
                        auto angle = [&](const point<T, 3> &c) { const auto w = c - pr.o; return std::atan2(w * f.y, w * f.x); };
                        std::vector<T> us;
                        for (const auto &x : balls)
                        {
                            T u = angle(x.c);
                            if (!us.empty())
                                u += T(2) * pi * std::round((us.back() - u) / (T(2) * pi));
                            us.push_back(u);
                        }
                        u1 = us.front(), u2 = q.closed ? us.front() + T(2) * pi : us.back();
                        T vm_first{};
                        for (std::size_t i = 0; i < balls.size() && exact; ++i)
                        {
                            const auto &x = balls[i];
                            const auto e = f.radial(us[i]);
                            auto minor = [&](const point<T, 3> &p) { const auto w = p - (pr.o + pr.R * e); return std::atan2(w * pr.z, w * e); };
                            // (the contact with M continuous from station to station: atan2 jumps at pi, the bottom of a groove)
                            T vm = minor(x.p[1]);
                            if (i > 0)
                                vm += T(2) * pi * std::round((vm_first - vm) / (T(2) * pi));
                            else
                                vm_first = vm;
                            for (std::size_t k : {std::size_t{0}, std::size_t{2}})
                            {
                                T v = minor(x.p[k]);
                                v += T(2) * pi * std::round((vm - v) / (T(2) * pi));
                                v1 = std::min({v1, v, vm}), v2 = std::max({v2, v, vm});
                            }
                        }
                        if (exact && !(v2 - v1 < T(2) * pi - T(1e-6)))
                            exact = false;
                        if (exact)
                            rb.surface = std::make_shared<ToroidalSurface<T>>(pr.R, pr.r, f, std::array<T, 4>{u1, u2, v1, v2});
                    }
                    if (exact)
                    {
                        rb.exact = true;
                        rb.variable = pr.kind == 2;
                        rb.spine = spine;
                        for (std::size_t k = 0; k < 2; ++k)
                        {
                            rb.contacts[k].pcurve = pc[k];
                            rb.contacts[k].curve = std::make_shared<CurveOnSurface<T, 3>>(pc[k], rb.contacts[k].support);
                        }
                    }
                    else
                        rb.surface = nullptr;
                }
            }
            if (!rb.surface)
            {
                rb.exact = false;
                // a station where the guide is below C3, and where a contact crosses a knot of its face below C3
                std::vector<T> breaks;
                for (const T k : fillet_detail::edge_weak_knots(*guide, 3))
                    if (k > a && k < b)
                        breaks.push_back(k);
                for (std::size_t k = 0; k < 3; ++k)
                    for (std::size_t d = 0; d < 2; ++d)
                        for (const T x : fillet_detail::weak_knots(*sup[k].surface, d))
                            for (int i = 0; i < ns; ++i)
                            {
                                auto g = [&](T t) { return ball_of(raw(t)).uv[k][d] - x; };
                                T tl = at(i), th = at(i + 1), gl = g(tl);
                                if (!(gl * g(th) < T(0)))
                                    continue;
                                for (int j = 0; j < 60 && th - tl > T(1e-13) * (T(1) + std::abs(a) + std::abs(b)); ++j)
                                {
                                    const T mdl = (tl + th) / T(2), gm = g(mdl);
                                    if (gm * gl > T(0))
                                        tl = mdl, gl = gm;
                                    else
                                        th = mdl;
                                }
                                breaks.push_back((tl + th) / T(2));
                            }
                rb.variable = true;
                auto sk = fillet_detail::with_edge(eid, [&] {
                    return fillet_detail::skinned(station, std::array<const Surface<T, 3> *, 2>{sup[0].surface.get(), sup[2].surface.get()}, a, b, o.approx_tol,
                                                  o.angular_tol, breaks);
                });
                rb.surface = sk.surface;
                rb.deviation = sk.deviation;
                rb.spine = sk.spine;
                for (std::size_t k = 0; k < 2; ++k)
                {
                    rb.contacts[k].pcurve = sk.pcurves[k];
                    rb.contacts[k].curve = std::make_shared<CurveOnSurface<T, 3>>(sk.pcurves[k], rb.contacts[k].support);
                }
            }

            // ---- the contacts within their faces: on the sides as those of a fillet of an edge;
            // on the middle face, in at the middle of the guide and out only across an end of it
            fillet_detail::check_on_face(m, q.side[0], *rb.contacts[0].pcurve, ed, rb.t1, rb.t2, tol, op, eid);
            // (on F2 between the stations at the ends of its edge with M, not those of the guide: an
            // oblique end puts them apart, the middle of the guide off F2)
            {
                auto foot = [&](VertexId v) {
                    const T t0 = v == (q.same ? ed2.v1 : ed2.v2) ? ed.u1 : ed.u2;
                    return std::clamp(extrema_curve_point(*guide, m.vertex(v).pnt, t0, T(1e-14) * (T(1) + std::abs(b - a)))[0], a, b);
                };
                if (q.closed)
                    fillet_detail::check_on_face(m, q.side[2], *rb.contacts[1].pcurve, ed2, rb.t1, rb.t2, tol, op, q.guide[1]);
                else
                    fillet_detail::check_on_face(m, q.side[2], *rb.contacts[1].pcurve, ed2, foot(ed2.v1), foot(ed2.v2), tol, op, q.guide[1]);
            }
            {
                const auto fbox = classify_detail::face_uv_box(m, q.f[1]);
                const auto per = fillet_detail::periods_of(*q.side[1].surface, tol);
                for (int i = 0; i <= 16; ++i)
                {
                    const T t = rb.t1 + (rb.t2 - rb.t1) * T(i) / T(16);
                    const auto bl = ball_of(raw(t));
                    if (classify(m, q.f[1], fillet_detail::into_face(per, fbox, bl.uv[1]), tol) != Location::Out)
                        continue;
                    const auto across = fillet_detail::nearest_edge(m, q.f[1], EdgeId{}, bl.p[1]);
                    const bool at_end = across && (std::ranges::contains(q.ends[0], *across) || std::ranges::contains(q.ends[1], *across));
                    if (i == 8 || !at_end)
                        blend_detail::refuse<T>(op, "the ball of the full round leaves the middle face (it touches it beyond one of its sides)", {q.f[1]});
                }
            }
            fillet_detail::on_blend(rb, std::max(o.intersect.tol, tol / T(10)), true);
            return rb;
        }

        /**
         * The pcurve on F2 of the guide put in the place of the other edge of M (its pcurve pc
         * there, on the parameter of that edge): at the parameter t of the guide, the point of
         * the old pcurve where the plane normal to the guide at t meets the other edge (the
         * correspondence of the stations), interpolated (cubic) through 2^k + 1 points, doubled
         * from 16 until within 1e-9 of the old pcurve between them. Not the segment between the
         * ends of the old pcurve (PR 3): for a closed edge they are the same point, and a curved
         * pcurve (an arc of circle on a plane) is not its chord.
         */
        template <std::floating_point T>
        auto guide_on_other(const Model<T> &m, const Setup<T> &q, const std::shared_ptr<Curve<T, 2>> &pc) -> std::shared_ptr<Curve<T, 2>>
        {
            const auto &g = m.edge(q.guide[0]), &g2 = m.edge(q.guide[1]);
            const T period = g2.u2 - g2.u1;
            // the parameter on the other edge of the point in the plane normal to the guide at t, from s0
            auto across = [&](T t, T s0) {
                const auto p = g.curve->value(t), te = fillet_detail::unit(g.curve->value(t, 1));
                T s = s0;
                for (int i = 0; i < 40; ++i)
                {
                    auto sv = s;
                    if (q.closed)
                        sv = g2.u1 + std::fmod(std::fmod(s - g2.u1, period) + period, period);
                    const T h = (g2.curve->value(sv) - p) * te, dh = g2.curve->value(sv, 1) * te;
                    if (!(std::abs(dh) > T(1e-14)))
                        break;
                    const T step = h / dh;
                    s -= step;
                    if (!q.closed)
                        s = std::clamp(s, g2.u1, g2.u2);
                    if (std::abs(step) <= T(1e-15) * (T(1) + std::abs(s)))
                        break;
                }
                return s;
            };
            auto on_old = [&](T s) {
                if (q.closed)
                    s = g2.u1 + std::fmod(std::fmod(s - g2.u1, period) + period, period);
                return pc->value(s);
            };
            // the start: the end of the other edge at the start of the guide
            const T s_start = q.closed ? across(g.u1, (q.same ? g2.u1 : g2.u2)) : (q.same ? g2.u1 : g2.u2);
            std::shared_ptr<Curve<T, 2>> best;
            for (std::size_t n = 16; n <= 1024; n *= 2)
            {
                points_vector<T, 2> pts;
                std::vector<T> ts, ss;
                T s = s_start;
                for (std::size_t i = 0; i <= n; ++i)
                {
                    const T t = i == n ? g.u2 : g.u1 + (g.u2 - g.u1) * T(i) / T(n);
                    s = i == 0 ? s_start : across(t, s);
                    ts.push_back(t), ss.push_back(s), pts.push_back(on_old(s));
                }
                if (!q.closed)
                    pts.front() = on_old(q.same ? g2.u1 : g2.u2), pts.back() = on_old(q.same ? g2.u2 : g2.u1);
                auto c = std::make_shared<BSCurve<T, 2>>(interpolate(pts, ts, 3));
                best = c;
                T dev{};
                for (std::size_t i = 0; i < n; ++i)
                {
                    const T t = (ts[i] + ts[i + 1]) / T(2);
                    dev = std::max(dev, norm(c->value(t) - on_old(across(t, (ss[i] + ss[i + 1]) / T(2)))));
                }
                if (dev <= T(1e-9))
                    break;
            }
            return best;
        }

        /**
         * The middle face of the setup q (on a copy ns of the solid) collapsed onto the guide
         * (file header): M out of its shell and erased, its edge with F2 replaced by the guide in
         * F2, its edges at each end reduced to the vertex of the guide there (out of the wires of
         * the end faces), the vertices of the other edge replaced by those of the guide.
         */
        template <std::floating_point T>
        void collapse(Model<T> &m, SolidId ns, const Setup<T> &q)
        {
            const blend_detail::SolidView<T> sv{m, ns};
            const auto &g = m.edge(q.guide[0]);
            const auto g2 = m.edge(q.guide[1]);
            // the use of the guide by M: that of F2 once M is gone
            Orientation om{};
            for (auto wid : m.face(q.f[1]).wires)
                for (const auto &ce : m.wire(wid).coedges)
                    if (ce.edge == q.guide[0])
                        om = ce.orient;
            for (auto wid : m.face(q.f[2]).wires)
                for (auto &ce : m.wire(wid).coedges)
                    if (ce.edge == q.guide[1])
                        ce = CoEdge<T>{q.guide[0], om, guide_on_other(m, q, ce.pcurve)};
            // the vertices of the ends of M: those of the guide
            std::map<VertexId, VertexId> to;
            for (std::size_t k = 0; k < 2; ++k)
            {
                const VertexId v = k == 0 ? g.v1 : g.v2;
                for (auto e : q.ends[k])
                    for (auto x : {m.edge(e).v1, m.edge(e).v2})
                        if (x != v)
                            to[x] = v;
            }
            std::set<EdgeId> gone{q.guide[1]};
            for (const auto &run : q.ends)
                gone.insert(run.begin(), run.end());
            const auto sh = sv.shell_of(q.f[1]);
            std::erase_if(m.shell(sh).faces, [&](const FaceUse &fu) { return fu.face == q.f[1]; });
            for (auto wid : m.face(q.f[1]).wires)
                m.erase(wid);
            m.erase(q.f[1]);
            for (const auto &fu : sv.faces())
            {
                if (fu.face == q.f[1])
                    continue;
                for (auto wid : m.face(fu.face).wires)
                {
                    std::erase_if(m.wire(wid).coedges, [&](const CoEdge<T> &ce) { return gone.contains(ce.edge); });
                    for (const auto &ce : m.wire(wid).coedges)
                    {
                        auto &e = m.edge(ce.edge);
                        if (const auto it = to.find(e.v1); it != to.end())
                            e.v1 = it->second;
                        if (const auto it = to.find(e.v2); it != to.end())
                            e.v2 = it->second;
                    }
                }
            }
            for (auto e : gone)
                m.erase(e);
            for (const auto &[x, v] : to)
                m.erase(x);
        }
        /**
         * One face on both sides (the setup q on the copy ns): split along the isoparametric
         * line halfway between its two edges with M, which must be isoparametric in the same
         * direction (a blade of one patch lofted from open sections, a U swept from one curve).
         * The pieces along the guide and along the other edge. Refused otherwise, saying to split
         * the face.
         */
        /// The pieces of the face on both sides (along the guide, along the other edge) and the ends of the cut between them.
        template <std::floating_point T>
        struct Split
        {
            std::array<FaceId, 2> pieces;
            std::array<point<T, 3>, 2> ends;
        };

        template <std::floating_point T>
        auto split_between(Model<T> &m, SolidId ns, const Setup<T> &q, T tol, const std::string &op) -> Split<T>
        {
            const blend_detail::SolidView<T> sv{m, ns};
            const FaceId f = q.f[0];
            auto refuse_split = [&] [[noreturn]] (const std::string &why) {
                blend_detail::refuse<T>(op, "the same face on both sides " + why + ": split it between its edges with the middle face (split_face) and give each piece as a side",
                                        {f});
            };
            // (a face closing across a seam: its pieces would have to be cut twice)
            for (auto w : m.face(f).wires)
                for (const auto &ce : m.wire(w).coedges)
                    if (!m.edge(ce.edge).degenerate && std::ranges::count_if(m.wire(w).coedges, [&](const CoEdge<T> &x) { return x.edge == ce.edge; }) > 1)
                        refuse_split("closing on itself across a seam (not handled)");
            const auto box = classify_detail::face_uv_box(m, f);
            std::optional<std::size_t> dir;
            std::array<T, 2> at{};
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto pc = sv.pcurve(f, q.guide[k]);
                const auto &ed = m.edge(q.guide[k]);
                if (!pc)
                    refuse_split("without a pcurve along the middle face");
                std::optional<std::size_t> d;
                for (std::size_t c = 0; c < 2 && !d; ++c)
                {
                    // (within tol in space: the step in the parameter c times the speed of the
                    // surface along it, at 17 points of the edge)
                    const T x0 = pc->value(ed.u1)[c];
                    const auto &srf = *m.face(f).surface;
                    bool iso = true;
                    for (int i = 1; i <= 16 && iso; ++i)
                    {
                        const auto uv = pc->value(ed.u1 + (ed.u2 - ed.u1) * T(i) / T(16));
                        const T speed = norm(srf.value(uv[0], uv[1], c == 0 ? 1 : 0, c == 1 ? 1 : 0));
                        iso = std::abs(uv[c] - x0) * speed <= tol;
                    }
                    if (iso)
                        d = c, at[k] = x0;
                }
                if (!d || (dir && *dir != *d))
                    refuse_split("whose edges with the middle face are not isoparametric in one direction");
                dir = d;
            }
            const std::size_t c = *dir, o = 1 - c;
            const T x = (at[0] + at[1]) / T(2);
            // the isoparametric line x across the face: from where it crosses the contour of the
            // face to where it crosses it again (an end of the face cut obliquely, the root of a
            // blade fused to its hub: not the box of the face)
            std::vector<T> cross;
            const T span = box[2 * o + 1] - box[2 * o];
            for (auto w : m.face(f).wires)
                for (const auto &ce : m.wire(w).coedges)
                {
                    const auto &ed = m.edge(ce.edge);
                    if (ed.degenerate || !ce.pcurve)
                        continue;
                    auto g = [&](T t) { return ce.pcurve->value(t)[c] - x; };
                    constexpr int n = 64;
                    for (int i = 0; i < n; ++i)
                    {
                        T lo = ed.u1 + (ed.u2 - ed.u1) * T(i) / T(n), hi = ed.u1 + (ed.u2 - ed.u1) * T(i + 1) / T(n);
                        T glo = g(lo);
                        const T ghi = g(hi);
                        T t{};
                        if (glo == T(0))
                            t = lo;
                        else if (ghi == T(0))
                            t = hi;
                        else if (glo * ghi < T(0))
                        {
                            for (int j = 0; j < 80 && hi - lo > T(1e-15) * (T(1) + std::abs(ed.u2 - ed.u1)); ++j)
                            {
                                const T mid = (lo + hi) / T(2), gm = g(mid);
                                if (gm * glo > T(0))
                                    lo = mid, glo = gm;
                                else
                                    hi = mid;
                            }
                            t = (lo + hi) / T(2);
                        }
                        else
                            continue;
                        const T y = ce.pcurve->value(t)[o];
                        if (std::ranges::none_of(cross, [&](T z) { return std::abs(z - y) <= T(1e-9) * (T(1) + span); }))
                            cross.push_back(y);
                    }
                }
            std::ranges::sort(cross);
            if (cross.size() != 2)
                refuse_split("whose isoparametric line halfway between its edges with the middle face crosses its contour " + std::to_string(cross.size()) +
                             " times, not twice (holes, or a contour that turns back)");
            point<T, 2> p0{}, p1{};
            p0[c] = x, p1[c] = x;
            p0[o] = cross[0], p1[o] = cross[1];
            const auto pc = fillet_detail::segment<T>(p0, p1, T(0), T(1));
            const auto curve = std::make_shared<CurveOnSurface<T, 3>>(pc, m.face(f).surface);
            const auto e = make_edge(m, std::shared_ptr<Curve<T, 3>>{curve}, T(0), T(1), tol);
            if (!e)
                refuse_split("that cannot be cut halfway between its edges with the middle face");
            const auto ends = std::array{curve->value(T(0)), curve->value(T(1))};
            const std::array<FaceCut<T>, 1> cuts{FaceCut<T>{*e, pc}};
            const auto sp = split_face(m, f, std::span<const FaceCut<T>>{cuts}, tol);
            if (!sp || sp->faces.size() != 2)
                refuse_split("that cannot be cut halfway between its edges with the middle face");
            // the piece along the guide first
            auto holds = [&](FaceId g, EdgeId y) {
                for (auto w : m.face(g).wires)
                    for (const auto &ce : m.wire(w).coedges)
                        if (ce.edge == y)
                            return true;
                return false;
            };
            if (holds(sp->faces[0], q.guide[0]) && holds(sp->faces[1], q.guide[1]))
                return {{sp->faces[0], sp->faces[1]}, ends};
            if (holds(sp->faces[1], q.guide[0]) && holds(sp->faces[0], q.guide[1]))
                return {{sp->faces[1], sp->faces[0]}, ends};
            refuse_split("whose pieces, cut halfway, do not each hold an edge with the middle face");
        }

        /**
         * The two pieces of one face split by `split_between` merged again on the result r, the
         * edges split by the cut joined back. The pieces: the two faces of r on one surface object
         * (a copy shares the copy of a surface between its faces) that is the surface s of the
         * face given (the same point at the middle of the box b of that face).
         */
        template <std::floating_point T>
        void merge_pieces(Model<T> &m, SolidId r, const Surface<T, 3> &s, const std::array<T, 4> &b, const std::array<point<T, 3>, 2> &ends, T tol,
                          const std::string &op)
        {
            const auto sh = m.solid(r).outer;
            const T u = (b[0] + b[1]) / T(2), v = (b[2] + b[3]) / T(2);
            const auto p = s.value(u, v);
            std::map<const Surface<T, 3> *, std::vector<FaceUse>> by_surface;
            for (const auto &fu : m.shell(sh).faces)
                by_surface[m.face(fu.face).surface.get()].push_back(fu);
            std::vector<FaceUse> on;
            for (const auto &[x, fus] : by_surface)
                if (fus.size() == 2 && norm(x->value(u, v) - p) <= tol)
                    on = fus;
            if (on.size() != 2 || !merge_faces(m, on[0], on[1], tol))
                blend_detail::refuse<T>(op, "the pieces of the face on both sides cannot be merged again after the full round");
            // (the edges split by the cut only: at its ends, not elsewhere on the solid)
            std::vector<VertexId> at;
            for (auto v : explore<VertexId>(m, ShapeId{r}))
                if (std::ranges::any_of(ends, [&](const point<T, 3> &p) { return norm(m.vertex(v).pnt - p) <= T(10) * tol; }))
                    at.push_back(v);
            const auto faces = m.shell(sh).faces;
            join_edges(m, std::span<const FaceUse>{faces}, std::span<const VertexId>{at});
        }
    }

    /**
     * @brief The geometry of the full round fillet of the face `middle` of the solid `s` between
     * its neighbours `face1` and `face2` (file header): the surface, the centres of the balls
     * (`spine`; the radius at t is the distance from it to a contact) and the contacts on `face1`
     * and `face2`, every curve on the parameter of the edge between `middle` and `face1`
     * (`edge`). The model is not modified.
     */
    template <std::floating_point T>
    [[nodiscard]] auto full_round_geometry(const Model<T> &m, SolidId s, FaceId face1, FaceId middle, FaceId face2, const FilletOptions<T> &o = {})
        -> BuildResult<EdgeBlend<T>>
    {
        const std::string op = "full_round_geometry";
        return fillet_detail::guarded(op, [&] {
            const std::array<FaceId, 1> f1{face1}, fm{middle}, f2{face2};
            const auto q = full_round_detail::setup(m, s, std::span<const FaceId>{f1}, std::span<const FaceId>{fm}, std::span<const FaceId>{f2}, o, op);
            return full_round_detail::geometry(m, q, o, op);
        });
    }

    /**
     * @brief The solid `s` with its face `middle` replaced by the full round fillet tangent to it
     * and to its neighbours `faces1` and `faces2` (file header; stage 8, § 3 of
     * `docs/sources/design/brep_palier8_face_blends.md`, PR 3: `brep8_pr03_full_round.md`, PR 4:
     * `brep8_pr04_full_round_loops.md`): a new
     * solid, `s` unchanged. The radius is not given: it is that of the ball touching the three
     * faces, at each station along the middle face (half the width of a rib with parallel
     * flanks).
     *
     * The middle face is taken away; `faces1` and `faces2` are cut along the contacts of the
     * ball, the faces at the ends of the middle face (each next to it and to both sides) along
     * the curve where they meet the fillet. One face on each side and one in the middle (the
     * sets of faces of `fillet_faces` are not handled yet); the same face on both sides (a blade
     * of one patch) is split between its two edges with the middle face and merged again. A middle
     * face round a loop (the rim of a disc, the bottom of a groove round a shaft) has no end face.
     *
     * Refused (`InvalidId`): a side or the middle empty, a face dead or not of the solid, a face
     * both a side and the middle, a middle face with no edge with a side. `InvalidBounds`: a
     * conic or G2 section, setbacks, tolerances not positive. `NonGeneric`: no ball touching the
     * three faces (fillet the two edges of the middle face instead), a contact out of its face,
     * balls that nest, an end face that does not cut the fillet, a face more curved toward the
     * ball than the ball at its contact (across or along: folding), an end face tangent to a side
     * (a rounded end), a middle face with holes or narrowing to a point (not handled). On failure
     * the model is left as it was.
     */
    template <std::floating_point T>
    [[nodiscard]] auto fillet_full_round(Model<T> &m, SolidId s, std::span<const FaceId> faces1, std::span<const FaceId> middle, std::span<const FaceId> faces2,
                                         const FilletOptions<T> &o = {}) -> BuildResult<SolidId>
    {
        const std::string op = "fillet_full_round";
        return blend_detail::guarded(m, op, [&] {
            const auto q = full_round_detail::setup(m, s, faces1, middle, faces2, o, op);
            const auto size = detail::arena_size(m);
            detail::ShapeCopier<T> copier{m, Placement<T>{}};
            try
            {
                // the copy with the middle face collapsed onto the guide, then the blend of the guide
                const auto ns = std::get<SolidId>(copier.copy(ShapeId{s}));
                auto qc = q;
                for (auto &g : qc.f)
                    g = copier.copy(g);
                for (auto &e : qc.guide)
                    e = copier.copy(e);
                for (auto &run : qc.ends)
                    for (auto &e : run)
                        e = copier.copy(e);
                // one face on both sides: split on the copy between its two edges with M, a piece
                // on each side, merged again on the result
                std::array<point<T, 3>, 2> cut_ends{};
                if (q.one)
                {
                    const auto sp = full_round_detail::split_between(m, ns, qc, o.tol, op);
                    qc.f[0] = sp.pieces[0], qc.f[2] = sp.pieces[1];
                    cut_ends = sp.ends;
                }
                full_round_detail::collapse(m, ns, qc);
                const std::array<EdgeId, 1> edges{qc.guide[0]};
                auto r = blend_detail::blend(
                    m, ns, std::span<const EdgeId>{edges},
                    [&](SolidId, EdgeId e, detail::ShapeCopier<T> &c, const FilletOptions<T> &oe) {
                        // (the geometry on the solid given; its faces and its edge those of the copy blended)
                        auto g = full_round_detail::geometry(m, q, oe, op);
                        g.edge = e;
                        g.contacts[0].face = c.copy(qc.f[0]);
                        g.contacts[1].face = c.copy(qc.f[2]);
                        return g;
                    },
                    o, op, true);
                if (q.one)
                    full_round_detail::merge_pieces(m, r, *m.face(q.f[0]).surface, classify_detail::face_uv_box(m, q.f[0]), cut_ends, o.tol, op);
                boolean_detail::erase_shape(m, ShapeId{ns});
                boolean_detail::erase_unused(m, size[1], size[0]);
                return r;
            }
            catch (const BuildFailure &e)
            {
                throw BuildFailure(detail::to_originals(e.error(), {&copier}));
            }
        });
    }

    template <std::floating_point T>
    [[nodiscard]] auto fillet_full_round(Model<T> &m, SolidId s, const std::vector<FaceId> &faces1, const std::vector<FaceId> &middle,
                                         const std::vector<FaceId> &faces2, const FilletOptions<T> &o = {}) -> BuildResult<SolidId>
    {
        return fillet_full_round(m, s, std::span<const FaceId>{faces1}, std::span<const FaceId>{middle}, std::span<const FaceId>{faces2}, o);
    }
}
