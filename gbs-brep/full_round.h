#pragma once
/**
 * @file full_round.h
 * @brief Full round fillets (stage 8, `docs/sources/design/brep_palier8_face_blends.md`, § 3;
 * PR 3: `docs/sources/design/brep8_pr03_full_round.md`): the middle face M between two faces F1
 * and F2 replaced by the fillet tangent to the three, whose radius follows the geometry (the
 * "Tritangent Fillet" of CATIA, the "Full round fillet" of SolidWorks): the top of a rib rounded
 * into a half cylinder, the blunt trailing edge of a blade rounded.
 *
 * **Geometry.** A ball touching the three faces rolls along M:
 * 1. its **stations** are on the **guide**, the edge between M and F1: in the plane normal to the
 *    guide at its parameter t, the centre c, the radius r and the three contacts are the solution
 *    of c = S_k(u_k, v_k) -+ r N_k(u_k, v_k) for the three faces and (c - g(t)) . g'(t) = 0
 *    (Gauss-Newton on 10 unknowns, from the station next to it; the first from the circle
 *    inscribed in the tangents of the faces in that plane). There is no virtual edge: between
 *    parallel flanks F1 and F2 never meet;
 * 2. the **section** is the arc of the characteristic circle of the ball, the circle through its
 *    three contacts (the ball of a varying radius touches its envelope along that circle, not
 *    along a great circle), from the contact on F1 to the contact on F2 by the contact on M: a
 *    half turn between parallel flanks, more where they come nearer away from M (`arc_through`);
 * 3. the **surface** is a cylinder where the guide is a line, the three faces planes and the
 *    radius constant (a straight rib, with or without draft); otherwise the rational B-spline
 *    skinned through the sections (`skinned`), within `approx_tol` of them and its normal within
 *    `angular_tol`.
 * The ball is on the side of the material where M is convex along both its edges with the sides
 * (the fillet takes material off), outside where it is concave along both (the bottom of a slot:
 * the fillet adds material).
 *
 * **Topology.** On a copy of the solid, M is taken out and collapsed onto the guide: the edge
 * between M and F2 replaced by the guide in F2, the edges of M with its end faces reduced to
 * their vertex. The copy has the topology of the solid with a sharp edge in place of M (not its
 * geometry: F2 is a width of M away from that edge), which is what the blends of `blend.h` need:
 * F1 and F2 cut along the contacts, each end face cutting the fillet where it meets it (a free
 * end, `free_end`), the result checked (`check`). The collapsed copy is erased.
 *
 * Refused (`InvalidId`): a face dead, not of the solid or given twice, a middle face that has no
 * edge with a side. Refused (`InvalidBounds`): a section other than the arc of circle, tolerances
 * not positive. Refused (`NonGeneric`): sides or a middle of several faces, a middle face with
 * holes, meeting a side along several edges, closing on itself or narrowing to a point (not
 * handled), with several faces at an end; a middle face convex along a side and concave along the other; no ball touching the
 * three faces; a contact out of its face; balls that nest (the middle face narrows faster than
 * the ball moves along it); an end face that does not cut the fillet. On failure the model is
 * left as it was.
 *
 * Sources: Vida, Martin and Várady, "A survey of blending methods that use parametric surfaces"
 * (Computer-Aided Design 26(5), 1994), § 5 (blends of three surfaces); the envelope of a family
 * of spheres (canal surfaces) for the characteristic circle; the behaviour of CATIA and of
 * SolidWorks.
 */
#include <gbs-brep/face_blend.h>

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
            if (q.f[0] == q.f[1] || q.f[1] == q.f[2] || q.f[0] == q.f[2])
                refuse<T>(BuildErrc::InvalidId, op, "a face given twice (two sides and the face between them)", {q.f[0], q.f[1], q.f[2]});
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
            if (m.face(q.f[1]).wires.size() != 1)
                blend_detail::refuse<T>(op, "a middle face with holes (not handled)", {q.f[1]});
            for (auto e : q.guide)
                if (m.edge(e).v1 == m.edge(e).v2)
                    blend_detail::refuse<T>(op, "a middle face that closes on itself (the rim of a disc: not handled yet)", {q.f[1], e});
            // the sides along their edges with M, and whether M is convex there
            std::array<bool, 2> convex{};
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto st = fillet_detail::setup(m, s, q.guide[k], op);
                convex[k] = st.convex;
                for (const auto &sd : st.sides)
                    if (sd.face == q.f[2 * k])
                        q.side[2 * k] = sd;
                    else if (k == 0)
                        q.side[1] = sd;
            }
            if (convex[0] != convex[1])
                blend_detail::refuse<T>(op, "a middle face convex along a side and concave along the other (no ball touches the three faces on one side of them)",
                                        {q.f[1]});
            q.convex = convex[0];
            // the wire of M: the guide, the edges at its end, the other edge, the edges at its start
            const auto &cs = m.wire(m.face(q.f[1]).wires.front()).coedges;
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
            for (const auto &run : q.ends)
            {
                if (run.empty())
                    blend_detail::refuse<T>(op, "a middle face narrowing to a point (its sides meet at an end of it: not handled; fillet_faces takes such a face)",
                                            {q.f[1]});
                if (run.size() != 1)
                    blend_detail::refuse<T>(op, "several faces at an end of the middle face (a rib ending in a prow: not handled)", {q.f[1]});
                for (auto e : run)
                    if (m.edge(e).degenerate)
                        blend_detail::refuse<T>(op, "a middle face with a degenerate edge at an end (not handled)", {q.f[1], e});
            }
            // M runs along the other edge against the guide: as edges, the same way when their uses by M differ
            q.same = cs[i0].orient != cs[i1].orient;
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
            auto no_ball = [&]() {
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

            // ---- the supports and the guide, continued past the ends of the middle face
            const T ext = o.extension * r0;
            const T reach = T(2) * (r0 + width) + ext;
            std::array<fillet_detail::Support<T>, 3> sup;
            for (std::size_t k = 0; k < 3; ++k)
                sup[k] = fillet_detail::support(m, q.side[k], k == 2 ? ed2 : ed, reach, tol);
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
            // a station from the nearest one solved; within the edge, from the inscribed circle
            // where that fails or jumps to another root. None: no ball there
            std::map<T, std::array<T, 10>> solved;
            auto find = [&](T t) -> const std::array<T, 10> * {
                if (const auto it = solved.find(t); it != solved.end())
                    return &it->second;
                std::optional<std::array<T, 10>> x;
                if (!solved.empty())
                {
                    auto it = solved.lower_bound(t);
                    if (it == solved.end() || (it != solved.begin() && t - std::prev(it)->first < it->first - t))
                        --it;
                    x = solve(t, it->second);
                    // (a ball that jumped to another root: its contacts moved further than the guide could carry them)
                    if (x)
                    {
                        const auto b0 = ball_of(it->second), b1 = ball_of(*x);
                        const T dt = std::abs(t - it->first), speed = norm(guide->value(t, 1));
                        T far{};
                        for (std::size_t k = 0; k < 3; ++k)
                            far = std::max(far, norm(b1.p[k] - b0.p[k]));
                        if (far > T(20) * speed * dt + T(0.1) * b0.r)
                            x.reset();
                    }
                }
                if (!x && t >= ed.u1 && t <= ed.u2)
                    x = solve(t, seed(t, {sup[0].surface.get(), sup[1].surface.get(), sup[2].surface.get()}));
                return x ? &solved.emplace(t, *x).first->second : nullptr;
            };
            auto raw = [&](T t) -> const std::array<T, 10> & {
                const auto *x = find(t);
                if (!x)
                    no_ball();
                return *x;
            };
            // from the middle to the vertices of the guide, then past them as far as there is a
            // ball (a middle face narrowing toward its end: none past the point where its sides
            // meet): the blend ends at the last station solved
            (void)raw(tm);
            for (const bool up : {true, false})
            {
                const T end = up ? ed.u2 : ed.u1, limit = up ? b : a;
                for (int i = 1; i <= 32; ++i)
                    (void)raw(i == 32 ? end : tm + (end - tm) * T(i) / T(32));
                T last = end;
                for (int i = 1; i <= 8 && limit != end; ++i)
                {
                    const T t = i == 8 ? limit : end + (limit - end) * T(i) / T(8);
                    if (!find(t))
                        break;
                    last = t;
                }
                (up ? b : a) = last;
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
            rb.circular = true;
            rb.contacts[0].face = q.f[0], rb.contacts[1].face = q.f[2];
            rb.contacts[0].support = sup[0].surface, rb.contacts[1].support = sup[2].surface;

            // ---- a cylinder: three planes, a straight guide along them, one radius
            const auto ba = ball_of(raw(a)), bb = ball_of(raw(b));
            const bool planes = std::ranges::all_of(q.side, [](const auto &sd) { return dynamic_cast<const Plane<T> *>(sd.surface.get()) != nullptr; });
            const auto axis = bb.c - ba.c;
            // (a line run at a constant speed, whatever its curve: the contacts are then affine in its parameter)
            bool straight = true;
            for (int i = 1; i < 8 && straight; ++i)
            {
                const T w = T(i) / T(8);
                straight = norm(guide->value(a + (b - a) * w) - ((T(1) - w) * guide->value(a) + w * guide->value(b))) <= T(1e-12) * scale;
            }
            if (planes && straight && std::abs(bb.r - ba.r) <= T(1e-12) * scale && norm(cross(unit(axis), unit(guide->value(tm, 1)))) <= T(1e-12))
            {
                constexpr T pi = std::numbers::pi_v<T>;
                const auto z = unit(axis);
                // x toward the contact with the middle face: the seam of the cylinder opposite
                auto x = ba.p[1] - ba.c;
                x = unit(x - (x * z) * z);
                const Frame3<T> f{ax2<T, 3>{ba.c, z, x}};
                std::array<T, 2> u{};
                for (std::size_t k = 0; k < 2; ++k)
                {
                    const auto w = ba.p[2 * k] - ba.c;
                    u[k] = std::atan2(w * f.y, w * f.x);
                }
                if (!(std::abs(u[0]) < pi - T(1e-6) && std::abs(u[1]) < pi - T(1e-6) && u[0] * u[1] < T(0)))
                    blend_detail::refuse<T>(op, "the contacts of the full round with the sides are not on either side of its contact with the middle face", {q.f[1]});
                rb.surface = std::make_shared<CylindricalSurface<T>>(ba.r, f, std::array<T, 4>{std::min(u[0], u[1]), std::max(u[0], u[1]), T(0), norm(axis)});
                rb.exact = true;
                const auto dir = axis / (b - a);
                rb.spine = std::make_shared<LineCurve<T>>(ba.c - a * dir, dir, a, b);
                for (std::size_t k = 0; k < 2; ++k)
                {
                    auto &ct = rb.contacts[k];
                    ct.pcurve = fillet_detail::segment<T>(ba.uv[2 * k], bb.uv[2 * k], a, b);
                    ct.curve = std::make_shared<CurveOnSurface<T, 3>>(ct.pcurve, ct.support);
                }
            }
            else
            {
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
            fillet_detail::check_on_face(m, q.side[2], *rb.contacts[1].pcurve, ed2, q.same ? rb.t1 : rb.t2, q.same ? rb.t2 : rb.t1, tol, op, q.guide[1]);
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
            for (const auto &ce : m.wire(m.face(q.f[1]).wires.front()).coedges)
                if (ce.edge == q.guide[0])
                    om = ce.orient;
            for (auto wid : m.face(q.f[2]).wires)
                for (auto &ce : m.wire(wid).coedges)
                    if (ce.edge == q.guide[1])
                    {
                        // (a pcurve of the guide on F2 from end to end of its old edge: the place of the co-edge in (u, v))
                        const auto from = coedge_uv(m, ce, g2.u1), to = coedge_uv(m, ce, g2.u2);
                        const bool fwd = om == Orientation::Forward;
                        ce = CoEdge<T>{q.guide[0], om, fillet_detail::segment<T>(fwd ? from : to, fwd ? to : from, g.u1, g.u2)};
                    }
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
     * `docs/sources/design/brep_palier8_face_blends.md`, PR 3: `brep8_pr03_full_round.md`): a new
     * solid, `s` unchanged. The radius is not given: it is that of the ball touching the three
     * faces, at each station along the middle face (half the width of a rib with parallel
     * flanks).
     *
     * The middle face is taken away; `faces1` and `faces2` are cut along the contacts of the
     * ball, the faces at the ends of the middle face (each next to it and to both sides) along
     * the curve where they meet the fillet. One face on each side and one in the middle (the
     * sets of faces of `fillet_faces` are not handled yet).
     *
     * Refused (`InvalidId`): a side or the middle empty, a face dead, not of the solid or given
     * twice, a middle face with no edge with a side. `InvalidBounds`: a conic or G2 section,
     * setbacks, tolerances not positive. `NonGeneric`: no ball touching the three faces (fillet
     * the two edges of the middle face instead), a contact out of its face, balls that nest, an
     * end face that does not cut the fillet, a middle face closing on itself, with holes or
     * narrowing to a point (not handled). On failure the model is left as it was.
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
                full_round_detail::collapse(m, ns, qc);
                const std::array<EdgeId, 1> edges{qc.guide[0]};
                const auto r = blend_detail::blend(
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
