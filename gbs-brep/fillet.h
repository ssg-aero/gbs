#pragma once
/**
 * @file fillet.h
 * @brief Geometry of a constant-radius fillet (`fillet_geometry`) and of a chamfer
 * (`chamfer_geometry`) along an edge of a solid: the surface of the blend and its two contact
 * curves on the faces of the edge. The solid is not modified; the topology (the faces cut along
 * the contacts, the ends, the sewing) comes with `fillet` and `chamfer`.
 *
 * **Fillet.** A ball of radius r rolling between the two faces of the edge, inside the material
 * along a convex edge, outside along a concave one (the side read on the outward normals and on
 * the directions of the faces away from the edge):
 * 1. its centre runs along the **centre line**, the intersection of the surfaces of the faces
 *    offset by r toward the ball (`offset_surface`: exact for the primitives; `intersect`: closed
 *    forms between primitives, marching otherwise), the branch nearest the edge, between the
 *    points nearest its vertices and continued past them by `extension` r (none for a closed
 *    edge: its whole loop, from its vertex);
 * 2. it touches each face at the foot of the centre on the face, c - r n: the **contacts**, whose
 *    pcurves are those of the centre line on the offset surfaces (same parameters, but the shift
 *    of v of a cone);
 * 3. the **surface** is the envelope of the ball, between the contacts: around a straight centre
 *    line a cylinder of radius r, around a circular one a torus of minor radius r (whatever the
 *    faces); along any other, a rational B-spline: the arc of circle between the contacts, exact
 *    (quadratic, one or two segments), at stations of the centre line put back on the rolling
 *    ball by Gauss-Newton (the marched line only starts them: it loses its accuracy where the
 *    faces become nearly tangent), skinned by cubic Hermite, the stations doubled until it is
 *    within `approx_tol` of the exact arcs between them and its normal within `angular_tol`;
 *    the centre line and the contacts interpolated on the same stations. Where the faces
 *    become tangent at a vertex, the fillet vanishes there (its contacts meet) and ends there,
 *    not continued past it.
 *
 * **Chamfer.** The ruled surface between two contacts, at the distances d1 on the first face and
 * d2 on the second from the edge, measured in the plane across the edge: the point of the face
 * at the distance d from the point of the edge in its normal plane, on the side of the face
 * (Gauss-Newton from the stations next to it; a pcurve interpolated through the stations,
 * doubled until within `approx_tol`). Between two lines in a plane it is a plane; between two
 * coaxial circles a plane, a cylinder or a cone; a ruled B-spline otherwise.
 *
 * Every curve of the result (centre line or edge, contacts, their pcurves) has the same
 * parameter, on [a, b].
 *
 * Refused (`NonGeneric`): a radius beyond a radius of curvature of a face on the side of the ball
 * (its offset folds), a contact out of its face (the radius or the distance beyond the width of
 * the face), faces tangent along the edge, an edge convex in places and concave in others, no
 * centre line near the edge. `InvalidId`: an edge not of the solid, used by one face or by more
 * than two, or twice by one face (a seam). `InvalidBounds`: a radius or a distance not positive.
 *
 * Sources: Rossignac and Requicha, "Constant-radius blending in solid modelling" (Computers in
 * Mechanical Engineering, 1984); Choi and Ju, "Constant-radius blending in surface modelling"
 * (Computer-Aided Design 21(4), 1989); Vida, Martin and Várady, "A survey of blending methods
 * that use parametric surfaces" (Computer-Aided Design 26(5), 1994); Piegl and Tiller, *The NURBS
 * Book*, 2nd ed., § 7.5 (arcs of circle) and § 10.3 (skinning).
 */
#include <gbs-brep/classify.h>
#include <gbs-brep/thick.h>
#include <gbs/bscinterp.h>
#include <gbs/curveonsurface.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs::brep
{
    template <std::floating_point T>
    struct FilletOptions
    {
        T tol{brep_default_tolerance<T>}; ///< tolerance of the result
        T approx_tol{T(1e-5)};            ///< approximated surfaces and contacts: deviation from the exact blend
        T angular_tol{T(1e-5)};           ///< approximated fillets: sine of the angle of their normal to that of the exact blend
        IntersectOptions<T> intersect{};  ///< intersection of the offset surfaces (its tol is raised to `tol` / 10 at least)
        T extension{T(1)};                ///< length of the blend past each end of an open edge, in radius (fillet) or largest distance (chamfer)
    };

    /// A contact curve of a blend on a face of the edge.
    template <std::floating_point T>
    struct BlendContact
    {
        FaceId face;
        std::shared_ptr<Surface<T, 3>> support; ///< the surface of the face, continued past the ends of the edge (same parameters)
        std::shared_ptr<Curve<T, 2>> pcurve;    ///< on `support`
        std::shared_ptr<Curve<T, 3>> curve;     ///< in space: `pcurve` on `support`
        std::shared_ptr<Curve<T, 2>> on_blend;  ///< on the surface of the blend
        T deviation{};                          ///< largest gap between `curve` and the surface of the blend along `on_blend`
    };

    /// The geometry of a fillet or of a chamfer along an edge: every curve on the parameters [a, b].
    template <std::floating_point T>
    struct EdgeBlend
    {
        EdgeId edge;
        bool convex{};                          ///< the blend removes material (adds it along a concave edge)
        bool closed{};                          ///< a closed edge: the blend runs its whole loop
        std::shared_ptr<Surface<T, 3>> surface; ///< cylinder, torus, plane, cone, or B-spline (rational for a fillet)
        bool exact{};                           ///< a primitive; otherwise within `deviation` of the exact blend
        T deviation{};                          ///< largest distance of the approximated surface from the exact blend
        std::shared_ptr<Curve<T, 3>> spine;     ///< the centre line of a fillet, the edge (continued) of a chamfer
        T a{}, b{};                             ///< the blend, past the ends of the edge by the extension
        T t1{}, t2{};                           ///< the parameters at the vertices v1 and v2 of the edge (t1 > t2: the other way)
        std::array<BlendContact<T>, 2> contacts; ///< on the faces of the two uses of the edge
    };

    namespace fillet_detail
    {
        template <std::floating_point T>
        [[noreturn]] void refuse(BuildErrc c, const std::string &op, const std::string &why, std::vector<ShapeId> shapes = {})
        {
            throw BuildFailure(BuildError{c, op + ": " + why, std::move(shapes)});
        }

        template <std::floating_point T>
        point<T, 3> unit(const point<T, 3> &v)
        {
            return v / norm(v);
        }

        /// A face of the edge.
        template <std::floating_point T>
        struct Side
        {
            FaceId face;
            Orientation orient{};                   ///< of the face in its shell
            std::shared_ptr<Surface<T, 3>> surface; ///< of the face
            std::shared_ptr<Curve<T, 2>> pcurve;    ///< of the edge on the face
            bool forward{};                          ///< the co-edge runs along the edge
        };

        template <std::floating_point T>
        struct EdgeSetup
        {
            std::array<Side<T>, 2> sides;
            bool convex{};
            bool closed{};
            T min_angle{}; ///< smallest angle of the wedge between the faces where the blend goes
        };

        /**
         * The faces of the edge in the solid and the side of the blend: convex when the direction
         * of the first face away from the edge, b1, goes below the second face (b1 . n2 < 0, n2
         * its outward normal), at 9 points; the angle of the wedge between b1 and b2.
         */
        template <std::floating_point T>
        auto setup(const Model<T> &m, SolidId sid, EdgeId eid, const std::string &op) -> EdgeSetup<T>
        {
            if (!m.alive(sid) || !m.alive(eid))
                refuse<T>(BuildErrc::InvalidId, op, "dead or invalid solid or edge", {sid, eid});
            const auto &so = m.solid(sid);
            std::vector<ShellId> shells{so.outer};
            shells.insert(shells.end(), so.voids.begin(), so.voids.end());
            std::optional<std::vector<EdgeUse>> uses;
            std::map<std::uint32_t, Orientation> face_orient;
            for (auto sh : shells)
            {
                for (const auto &fu : m.shell(sh).faces)
                    face_orient[fu.face.index] = fu.orient;
                auto all = shell_edge_uses(m, sh);
                if (const auto it = all.find(eid); it != all.end())
                    uses = it->second;
            }
            const auto &ed = m.edge(eid);
            if (!uses || ed.degenerate)
                refuse<T>(BuildErrc::InvalidId, op, "not an edge of the solid", {eid});
            if (uses->size() != 2)
                refuse<T>(BuildErrc::InvalidId, op, "an edge not shared by two faces", {eid});
            if ((*uses)[0].face == (*uses)[1].face)
                refuse<T>(BuildErrc::InvalidId, op, "a seam edge (used twice by one face)", {eid});
            EdgeSetup<T> st;
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto &u = (*uses)[k];
                const auto &ce = m.wire(u.coedge.wire).coedges[u.coedge.index];
                if (!ce.pcurve)
                    refuse<T>(BuildErrc::MissingPCurve, op, "an edge without pcurve", {eid});
                st.sides[k] = Side<T>{u.face, face_orient.at(u.face.index), m.face(u.face).surface, ce.pcurve, ce.orient == Orientation::Forward};
            }
            st.closed = ed.v1 == ed.v2;
            // at 9 points: the outward normals and the directions of the faces away from the edge
            int convex = 0, concave = 0;
            st.min_angle = std::numbers::pi_v<T>;
            for (int i = 0; i < 9; ++i)
            {
                const T t = ed.u1 + (ed.u2 - ed.u1) * (T(i) + T(0.5)) / T(9);
                const auto te = unit(ed.curve->value(t, 1));
                std::array<point<T, 3>, 2> n, b;
                for (std::size_t k = 0; k < 2; ++k)
                {
                    const auto &sd = st.sides[k];
                    const auto uv = sd.pcurve->value(t);
                    const auto su = sd.surface->value(uv[0], uv[1], 1, 0), sv = sd.surface->value(uv[0], uv[1], 0, 1);
                    const auto nn = normal_direction(*sd.surface, uv[0], uv[1]);
                    n[k] = sd.orient == Orientation::Forward ? nn : T(-1) * nn;
                    // the face is on the left of its co-edge in (u, v) (outer wires counter-clockwise, holes clockwise)
                    auto d = sd.pcurve->value(t, 1);
                    if (!sd.forward)
                        d = T(-1) * d;
                    auto bk = (-d[1]) * su + d[0] * sv;
                    bk = bk - (bk * te) * te;
                    b[k] = unit(bk);
                }
                const T c = b[0] * n[1];
                if (!(std::abs(c) > T(1e-6)))
                    refuse<T>(BuildErrc::NonGeneric, op, "the faces of the edge are tangent along it", {eid});
                (c < T(0) ? convex : concave) += 1;
                st.min_angle = std::min(st.min_angle, std::acos(std::clamp(b[0] * b[1], T(-1), T(1))));
            }
            if (convex > 0 && concave > 0)
                refuse<T>(BuildErrc::NonGeneric, op, "an edge convex in places and concave in others", {eid});
            st.convex = convex > 0;
            return st;
        }

        /// The surface of a face continued around the edge, and the box where a contact may lie.
        template <std::floating_point T>
        struct Support
        {
            std::shared_ptr<Surface<T, 3>> surface;
            std::array<T, 4> box{};
        };

        /**
         * The surface of the side continued by `reach` past the face (`extended_surface`), on the
         * band of its parameters within `reach` of the edge.
         */
        template <std::floating_point T>
        auto support(const Model<T> &m, const Side<T> &sd, const Edge<T> &ed, T reach, T tol) -> Support<T>
        {
            const auto fbox = classify_detail::face_uv_box(m, sd.face);
            const auto ex = thick_detail::extended_surface(sd.surface, fbox, reach, tol);
            // the box of the pcurve, and the smallest speeds along it
            std::array<T, 4> band{std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity(),
                                  -std::numeric_limits<T>::infinity()};
            std::array<T, 2> speed{std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity()};
            T fastest{};
            std::vector<std::array<T, 2>> sp;
            for (int i = 0; i <= 32; ++i)
            {
                const auto uv = sd.pcurve->value(ed.u1 + (ed.u2 - ed.u1) * T(i) / T(32));
                band = {std::min(band[0], uv[0]), std::max(band[1], uv[0]), std::min(band[2], uv[1]), std::max(band[3], uv[1])};
                const T a = norm(sd.surface->value(uv[0], uv[1], 1, 0)), b = norm(sd.surface->value(uv[0], uv[1], 0, 1));
                sp.push_back({a, b});
                fastest = std::max({fastest, a, b});
            }
            for (const auto &[a, b] : sp)
            {
                if (a > T(1e-9) * fastest)
                    speed[0] = std::min(speed[0], a);
                if (b > T(1e-9) * fastest)
                    speed[1] = std::min(speed[1], b);
            }
            std::array<T, 4> box = ex.box;
            for (int k = 0; k < 2; ++k)
                if (std::isfinite(speed[std::size_t(k)]))
                {
                    const T w = reach / speed[std::size_t(k)];
                    box[2 * k] = std::max(box[2 * k], band[2 * k] - w);
                    box[2 * k + 1] = std::min(box[2 * k + 1], band[2 * k + 1] + w);
                }
            if (!(box[1] > box[0]) || !(box[3] > box[2]))
                box = ex.box;
            return {ex.surface, box};
        }

        /// u and v of a pcurve point brought into the box of the face by periods.
        template <std::floating_point T>
        auto into_face(const Surface<T, 3> &s, const std::array<T, 4> &fbox, point<T, 2> uv) -> point<T, 2>
        {
            if (const auto *e = dynamic_cast<const ElementarySurface<T> *>(&s))
            {
                const auto per = e->periods();
                for (std::size_t k = 0; k < 2; ++k)
                    if (per[k] > T(0))
                    {
                        const T lo = fbox[2 * k], eps = T(1e-9) * per[k];
                        while (uv[k] < lo - eps)
                            uv[k] += per[k];
                        while (uv[k] > fbox[2 * k + 1] + eps)
                            uv[k] -= per[k];
                    }
            }
            return uv;
        }

        /// The contact stays on its face between the vertices (15 points inside the edge).
        template <std::floating_point T>
        void check_on_face(const Model<T> &m, const Side<T> &sd, const Curve<T, 2> &pc, T t1, T t2, T tol, const std::string &op, EdgeId eid)
        {
            const auto fbox = classify_detail::face_uv_box(m, sd.face);
            for (int i = 1; i < 16; ++i)
            {
                const auto uv = into_face(*sd.surface, fbox, pc.value(t1 + (t2 - t1) * T(i) / T(16)));
                if (classify(m, sd.face, uv, tol) == Location::Out)
                    refuse<T>(BuildErrc::NonGeneric, op, "a contact runs off its face (beyond the width of the face)", {eid, sd.face});
            }
        }

        /// The pcurve on the support of a pcurve on its offset (parameters (u, v + dv)).
        template <std::floating_point T>
        auto unshifted(const std::shared_ptr<Curve<T, 2>> &pc, T dv) -> std::shared_ptr<Curve<T, 2>>
        {
            if (dv == T(0))
                return pc;
            if (auto s = shifted_pcurve(pc, point<T, 2>{T(0), -dv}))
                return s;
            throw BuildFailure(BuildError{BuildErrc::NonGeneric, "fillet: a pcurve on a cone that cannot be moved", {}});
        }

        /// The quadratic rational arc of circle from p1 to p2 around c (one or two segments, on [0, 1]).
        template <std::floating_point T>
        auto arc(const point<T, 3> &c, const point<T, 3> &p1, const point<T, 3> &p2, int segments) -> BSCurveRational<T, 3>
        {
            const auto e1 = p1 - c, e2 = p2 - c;
            const T rho = (norm(e1) + norm(e2)) / T(2);
            const auto u1 = e1 / norm(e1), u2 = e2 / norm(e2);
            const T phi = std::atan2(norm(cross(u1, u2)), u1 * u2);
            const auto s = u1 + u2;
            if (!(norm(s) > T(1e-6)))
                throw BuildFailure(BuildError{BuildErrc::NonGeneric, "fillet: contacts opposite on the ball", {}});
            const auto bis = s / norm(s);
            if (segments == 1)
            {
                const T w = std::cos(phi / T(2));
                return BSCurveRational<T, 3>(points_vector<T, 3>{p1, c + (rho / w) * bis, p2}, std::vector<T>{0, 0, 0, 1, 1, 1}, std::vector<T>{1, w, 1}, 2);
            }
            const T w = std::cos(phi / T(4));
            const auto q1 = u1 + bis, q2 = bis + u2;
            return BSCurveRational<T, 3>(points_vector<T, 3>{p1, c + (rho / w) * (q1 / norm(q1)), c + rho * bis, c + (rho / w) * (q2 / norm(q2)), p2},
                                         std::vector<T>{0, 0, 0, T(0.5), T(0.5), 1, 1, 1}, std::vector<T>{1, w, 1, w, 1}, 2);
        }

        /// A segment of (u, v) from (u1, v1) at t = a to (u2, v2) at t = b.
        template <std::floating_point T>
        auto segment(const point<T, 2> &p, const point<T, 2> &q, T a, T b) -> std::shared_ptr<Curve<T, 2>>
        {
            return std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{p, q}, std::vector<T>{a, a, b, b}, 1);
        }

        /// A station of the ball: its centre, its contacts and their parameters on the surfaces of the faces.
        template <std::floating_point T>
        struct Station
        {
            point<T, 3> c, p1, p2;
            point<T, 2> uv1, uv2;
        };

        /// The rational B-spline of a fillet and the curves fitted on its stations.
        template <std::floating_point T>
        struct Skinned
        {
            std::shared_ptr<Surface<T, 3>> surface;
            std::shared_ptr<Curve<T, 3>> spine;
            std::array<std::shared_ptr<Curve<T, 2>>, 2> pcurves; ///< of the contacts, on the surfaces of the faces
            T deviation{};
        };

        /**
         * The rational B-spline of the arcs of circle around the centres between the contacts of
         * the stations of `station` on [a, b] (u across, from p1 at 0 to p2 at 1; v along): the
         * homogeneous poles of the arcs at n + 1 stations and their derivatives along the centre
         * line (central differences), interpolated by cubic Hermite (C1, the tangent planes of the
         * arcs at the stations); the centres and the parameters of the contacts interpolated
         * (cubic) at the same stations. n is doubled from 8 until, at three points between each
         * pair of stations, the surface is within approx_tol of the exact arcs and its normal
         * within angular_tol of theirs (where the arc is not reduced to a point), and the centre
         * line and the contacts (their pcurves on `supports`) within approx_tol of the stations.
         */
        template <std::floating_point T, typename F>
        auto skinned(F &&station, const std::array<const Surface<T, 3> *, 2> &supports, T a, T b, T approx_tol, T angular_tol) -> Skinned<T>
        {
            T phi{};
            for (int i = 0; i <= 16; ++i)
            {
                const auto st = station(a + (b - a) * T(i) / T(16));
                const auto e1 = st.p1 - st.c, e2 = st.p2 - st.c;
                phi = std::max(phi, std::atan2(norm(cross(e1, e2)), e1 * e2));
            }
            const int segments = phi > T(2) * std::numbers::pi_v<T> / T(3) ? 2 : 1;
            auto section = [&](const Station<T> &st) { return arc(st.c, st.p1, st.p2, segments); };
            const auto ku = section(station(a)).knotsFlats();
            const std::size_t nu = ku.size() - 3;
            auto homogeneous = [&](T s) {
                const auto x = section(station(s));
                std::vector<std::array<T, 4>> h(nu);
                for (std::size_t j = 0; j < nu; ++j)
                    h[j] = x.poles()[j];
                return h;
            };
            const T hs = T(1e-6) * (b - a);
            for (std::size_t n = 8; n <= 2048; n *= 2)
            {
                std::vector<T> v(n + 1);
                std::vector<std::vector<constrType<T, 4, 2>>> rows(nu, std::vector<constrType<T, 4, 2>>(n + 1));
                points_vector<T, 3> centres(n + 1);
                std::array<std::vector<constrType<T, 2, 1>>, 2> uv{std::vector<constrType<T, 2, 1>>(n + 1), std::vector<constrType<T, 2, 1>>(n + 1)};
                for (std::size_t i = 0; i <= n; ++i)
                {
                    v[i] = i == n ? b : a + (b - a) * T(i) / T(n);
                    const auto st = station(v[i]);
                    centres[i] = st.c;
                    uv[0][i] = {st.uv1};
                    uv[1][i] = {st.uv2};
                    const T s0 = std::max(a, v[i] - hs), s1 = std::min(b, v[i] + hs);
                    const auto x = section(st);
                    const auto h0 = homogeneous(s0), h1 = homogeneous(s1);
                    for (std::size_t j = 0; j < nu; ++j)
                    {
                        rows[j][i][0] = x.poles()[j];
                        rows[j][i][1] = (h1[j] - h0[j]) / (s1 - s0);
                    }
                }
                // the parameters of the contacts continuous from station to station (a periodic
                // parameter brought back by its period across the seam of its surface)
                for (std::size_t k = 0; k < 2; ++k)
                    if (const auto *e = dynamic_cast<const ElementarySurface<T> *>(supports[k]))
                        for (std::size_t d = 0; d < 2; ++d)
                            if (const T per = e->periods()[d]; per > T(0))
                                for (std::size_t i = 1; i <= n; ++i)
                                    uv[k][i][0][d] += per * std::round((uv[k][i - 1][0][d] - uv[k][i][0][d]) / per);
                std::vector<BSCurve<T, 4>> along;
                for (const auto &r : rows)
                    along.push_back(interpolate<T, 4, 2>(r, v));
                const auto kv = along.front().knotsFlats();
                const std::size_t nv = along.front().poles().size();
                std::vector<std::array<T, 4>> poles(nu * nv);
                for (std::size_t i = 0; i < nv; ++i)
                    for (std::size_t j = 0; j < nu; ++j)
                        poles[j + i * nu] = along[j].poles()[i];
                Skinned<T> r;
                r.surface = std::make_shared<BSSurfaceRational<T, 3>>(poles, ku, kv, std::size_t{2}, std::size_t{3});
                r.spine = std::make_shared<BSCurve<T, 3>>(interpolate<T, 3>(centres, v, std::size_t{3}));
                for (std::size_t k = 0; k < 2; ++k)
                    r.pcurves[k] = std::make_shared<BSCurve<T, 2>>(interpolate<T, 2>(uv[k], v, std::size_t{3}));
                T dev{}, ang{};
                for (std::size_t i = 0; i < n && dev <= approx_tol && ang <= angular_tol; ++i)
                    for (int j = 1; j < 4; ++j)
                    {
                        const T s = v[i] + (v[i + 1] - v[i]) * T(j) / T(4);
                        const auto st = station(s);
                        const auto exact = section(st);
                        dev = std::max(dev, norm(r.spine->value(s) - st.c));
                        for (std::size_t k = 0; k < 2; ++k)
                        {
                            const auto q = r.pcurves[k]->value(s);
                            dev = std::max(dev, norm(supports[k]->value(q[0], q[1]) - (k == 0 ? st.p1 : st.p2)));
                        }
                        const T arc_length = norm(st.p1 - st.c) * std::atan2(norm(cross(st.p1 - st.c, st.p2 - st.c)), (st.p1 - st.c) * (st.p2 - st.c));
                        for (int k = 0; k <= 8; ++k)
                        {
                            const T u = T(k) / T(8);
                            const auto q = r.surface->value(u, s);
                            dev = std::max(dev, norm(q - exact.value(u)));
                            // the normal of the envelope, along the radius (undefined where the arc
                            // is reduced to a point: a fillet vanishing where its faces become tangent)
                            if (arc_length > T(1e3) * approx_tol)
                            {
                                const auto e = exact.value(u) - st.c;
                                ang = std::max(ang, norm(cross(normal_direction(*r.surface, u, s), e / norm(e))));
                            }
                        }
                    }
                if (dev <= approx_tol && ang <= angular_tol)
                {
                    r.deviation = dev;
                    return r;
                }
            }
            throw BuildFailure(BuildError{BuildErrc::NonGeneric, "fillet: the surface cannot be approximated within approx_tol", {}});
        }

        /**
         * The blend surface as a primitive around a straight (cylinder) or circular (torus) spine,
         * on the box of the contacts; nullptr for another spine.
         */
        template <std::floating_point T>
        auto canal(const Curve<T, 3> &spine, const Curve<T, 3> &p1, const Curve<T, 3> &p2, T a, T b, T r) -> std::shared_ptr<ElementarySurface<T>>
        {
            constexpr T pi = std::numbers::pi_v<T>;
            const T sm = (a + b) / T(2);
            if (const auto *l = dynamic_cast<const LineCurve<T> *>(&spine))
            {
                const auto o = spine.value(sm);
                auto z = unit(l->vector());
                const auto e1 = p1.value(sm) - o;
                const auto x = unit(e1 - (e1 * z) * z);
                if ((p2.value(sm) - o) * cross(z, x) < T(0))
                    z = T(-1) * z;
                const Frame3<T> f{ax2<T, 3>{o, z, x}};
                T u1 = std::numeric_limits<T>::infinity(), u2 = -u1;
                for (int i = 0; i <= 32; ++i)
                {
                    const T s = a + (b - a) * T(i) / T(32);
                    for (const auto *p : {&p1, &p2})
                    {
                        const auto q = p->value(s) - o;
                        const T u = std::atan2(q * f.y, q * f.x);
                        u1 = std::min(u1, u), u2 = std::max(u2, u);
                    }
                }
                const T v1 = (spine.value(a) - o) * z, v2 = (spine.value(b) - o) * z;
                return std::make_shared<CylindricalSurface<T>>(r, f, std::array<T, 4>{u1, u2, std::min(v1, v2), std::max(v1, v2)});
            }
            if (const auto *e = dynamic_cast<const EllipseCurve<T> *>(&spine))
            {
                if (!(std::abs(e->radius1() - e->radius2()) <= T(1e-12) * e->radius1()))
                    return nullptr;
                const auto &f = e->frame();
                auto minor = [&](const Curve<T, 3> &p, T s) {
                    const auto q = p.value(s) - spine.value(s);
                    return std::atan2(q * f.z, q * f.radial(s));
                };
                // the arc between the contacts shorter than a half turn, each contact continuous along
                const T w1 = minor(p1, sm);
                T w2 = minor(p2, sm);
                while (w2 > w1 + pi)
                    w2 -= T(2) * pi;
                while (w2 < w1 - pi)
                    w2 += T(2) * pi;
                T v1 = std::numeric_limits<T>::infinity(), v2 = -v1;
                for (int i = 0; i <= 32; ++i)
                {
                    const T s = a + (b - a) * T(i) / T(32);
                    for (const auto &[p, w] : {std::pair{&p1, w1}, std::pair{&p2, w2}})
                    {
                        T v = minor(*p, s);
                        v += T(2) * pi * std::round((w - v) / (T(2) * pi));
                        v1 = std::min(v1, v), v2 = std::max(v2, v);
                    }
                }
                return std::make_shared<ToroidalSurface<T>>(e->radius1(), r, f, std::array<T, 4>{a, b, v1, v2});
            }
            return nullptr;
        }

        /**
         * The contacts' pcurves on the blend: on a primitive fitted (`primitive_pcurve`), on a
         * skinned surface its sides u = 0 and u = 1; the gap of each contact to the surface.
         */
        template <std::floating_point T>
        void on_blend(EdgeBlend<T> &r, T tol)
        {
            if (const auto *e = dynamic_cast<const ElementarySurface<T> *>(r.surface.get()))
            {
                for (auto &c : r.contacts)
                {
                    const auto fit = intersect_detail::primitive_pcurve(*c.curve, r.a, r.b, *e, e->bounds(), tol);
                    c.on_blend = fit.pcurve;
                    c.deviation = fit.deviation;
                }
                return;
            }
            for (std::size_t k = 0; k < 2; ++k)
            {
                auto &c = r.contacts[k];
                const T u = T(k);
                c.on_blend = segment<T>({u, r.a}, {u, r.b}, r.a, r.b);
                T dev{};
                for (int i = 0; i <= 256; ++i)
                {
                    const T s = r.a + (r.b - r.a) * T(i) / T(256);
                    dev = std::max(dev, norm(r.surface->value(u, s) - c.curve->value(s)));
                }
                c.deviation = dev;
            }
        }

        template <std::floating_point T>
        auto fillet(const Model<T> &m, SolidId sid, EdgeId eid, T r, const FilletOptions<T> &o) -> EdgeBlend<T>
        {
            const std::string op = "fillet_geometry";
            if (!(r > o.tol) || !std::isfinite(r))
                refuse<T>(BuildErrc::InvalidBounds, op, "the radius must be larger than the tolerance", {eid});
            const auto st = setup(m, sid, eid, op);
            const auto &ed = m.edge(eid);
            const T tol = o.tol;
            const T ext = st.closed ? T(0) : o.extension * r;
            // the contacts within r / tan(a / 2) of the edge, a the angle of the wedge
            const T reach = r / std::tan(std::max(st.min_angle, T(1e-3)) / T(2)) + r + ext;
            const T side = st.convex ? T(-1) : T(1); // the ball along the outward normals
            std::array<Support<T>, 2> sup;
            std::array<ExactOffset<T>, 2> off;
            std::array<std::array<T, 4>, 2> obox;
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto &sd = st.sides[k];
                sup[k] = support(m, sd, ed, reach, tol);
                const T d = side * r * (sd.orient == Orientation::Forward ? T(1) : T(-1));
                const auto x = offset_surface<T>(sup[k].surface, d, tol, sup[k].box);
                if (!x)
                    refuse<T>(BuildErrc::NonGeneric, op, "the radius exceeds a radius of curvature of a face", {eid, sd.face});
                off[k] = *x;
                obox[k] = {sup[k].box[0], sup[k].box[1], sup[k].box[2] + x->dv, sup[k].box[3] + x->dv};
            }
            IntersectOptions<T> io = o.intersect;
            io.tol = std::max(io.tol, tol / T(10));
            const auto res = intersect(*off[0].surface, obox[0], *off[1].surface, obox[1], io);

            // ---- the branch of the centre line nearest the edge
            const auto pm = ed.curve->value((ed.u1 + ed.u2) / T(2));
            const auto q1 = m.vertex(ed.v1).pnt, q2 = m.vertex(ed.v2).pnt;
            const IntersectionCurve<T> *best = nullptr;
            T best_gap = std::numeric_limits<T>::infinity();
            for (const auto &c : res.curves)
            {
                const auto [s, d] = extrema_curve_point(*c.curve, pm, T(1e-12) * (T(1) + std::abs(c.b - c.a)));
                if (d < best_gap)
                    best_gap = d, best = &c;
            }
            if (!best)
                refuse<T>(BuildErrc::NonGeneric, op, "no centre line near the edge (the offsets of its faces do not meet)", {eid});
            const auto &c = *best;
            const T span = c.b - c.a, eps = T(1e-12) * (T(1) + std::abs(c.a) + std::abs(c.b));
            const bool loop = c.closed || norm(c.curve->value(c.a) - c.curve->value(c.b)) <= T(10) * std::max(tol, io.tol);
            auto proj = [&](const point<T, 3> &p) { return extrema_curve_point(*c.curve, p, T(1e-12) * (T(1) + std::abs(span)))[0]; };
            auto wrap = [&](T t) { return loop && t >= c.b - eps ? t - span : t; };
            T start{}, end{}, t1{}, t2{};
            if (st.closed)
            {
                if (!loop)
                    refuse<T>(BuildErrc::NonGeneric, op, "the centre line of a closed edge is not closed", {eid});
                start = t1 = wrap(proj(q1));
                end = t2 = start + span;
            }
            else
            {
                const T s1 = proj(q1), s2 = proj(q2);
                bool along = s1 <= s2; // the centre line runs from v1 to v2
                if (!loop)
                    start = std::min(s1, s2), end = std::max(s1, s2);
                else
                {
                    // the arc between the vertices whose middle is nearest that of the edge
                    const T w1 = wrap(s1), w2 = wrap(s2);
                    auto at = [&](T t) {
                        t = c.a + std::fmod(t - c.a, span);
                        return c.curve->value(t < c.a ? t + span : t);
                    };
                    const std::array<std::array<T, 2>, 2> arcs{{{w1, w2 > w1 ? w2 : w2 + span}, {w2, w1 > w2 ? w1 : w1 + span}}};
                    const std::size_t i = norm(at((arcs[0][0] + arcs[0][1]) / T(2)) - pm) <= norm(at((arcs[1][0] + arcs[1][1]) / T(2)) - pm) ? 0 : 1;
                    start = arcs[i][0], end = arcs[i][1];
                    along = i == 0;
                }
                t1 = along ? start : end;
                t2 = along ? end : start;
                // continued past the vertices
                const T room = loop ? std::max(T(0), span - (end - start)) / T(2) : std::numeric_limits<T>::infinity();
                const T da = std::min(room, ext / norm(c.curve->value(std::clamp(start, c.a, c.b), 1)));
                const T db = std::min(room, ext / norm(c.curve->value(std::clamp(end <= c.b ? end : end - span, c.a, c.b), 1)));
                start -= da, end += db;
                if (!loop && !dynamic_cast<const ElementaryCurve<T> *>(c.curve.get()))
                    start = std::max(start, c.a), end = std::min(end, c.b);
                if (loop && start < c.a - eps)
                    start += span, end += span, t1 += span, t2 += span;
            }
            if (!(end - start > T(0)))
                refuse<T>(BuildErrc::NonGeneric, op, "the centre line reduces to a point", {eid});

            // ---- the centre line on [start, end] and its pcurves on the offsets
            EdgeBlend<T> rb;
            rb.edge = eid;
            rb.convex = st.convex;
            rb.closed = st.closed;
            rb.a = start, rb.b = end, rb.t1 = t1, rb.t2 = t2;
            std::array<std::shared_ptr<Curve<T, 2>>, 2> pc;
            if (const auto *ec = dynamic_cast<const ElementaryCurve<T> *>(c.curve.get()))
            {
                rb.spine = intersect_detail::restricted(*ec, start, end);
                for (std::size_t k = 0; k < 2; ++k)
                {
                    const auto *e = dynamic_cast<const ElementarySurface<T> *>(off[k].surface.get());
                    if (!e)
                        refuse<T>(BuildErrc::NonGeneric, op, "a closed-form centre line on a surface that is not a primitive", {eid});
                    pc[k] = intersect_detail::primitive_pcurve(*rb.spine, start, end, *e, obox[k], io.tol).pcurve;
                }
            }
            else if (end > c.b + eps)
            {
                rb.spine = thick_detail::from_parameter(c.curve, c.a, c.b, start);
                pc[0] = thick_detail::from_parameter(c.pc1, c.a, c.b, start);
                pc[1] = thick_detail::from_parameter(c.pc2, c.a, c.b, start);
            }
            else
                rb.spine = c.curve, pc[0] = c.pc1, pc[1] = c.pc2;

            // ---- the contacts
            for (std::size_t k = 0; k < 2; ++k)
            {
                auto &ct = rb.contacts[k];
                ct.face = st.sides[k].face;
                ct.support = sup[k].surface;
                ct.pcurve = unshifted(pc[k], off[k].dv);
                ct.curve = std::make_shared<CurveOnSurface<T, 3>>(ct.pcurve, ct.support);
                check_on_face(m, st.sides[k], *ct.pcurve, std::min(t1, t2), std::max(t1, t2), tol, op, eid);
            }

            // a fillet vanishing at a vertex (its faces tangent there: its contacts meet) ends
            // there, not continued past it: no face is to cut it at that end
            if (!st.closed)
                for (const T t : {t1, t2})
                    if (norm(rb.contacts[0].curve->value(t) - rb.contacts[1].curve->value(t)) <= T(10) * tol)
                        (t == std::min(t1, t2) ? rb.a : rb.b) = t;

            // ---- the surface
            if (auto e = canal(*rb.spine, *rb.contacts[0].curve, *rb.contacts[1].curve, rb.a, rb.b, r))
            {
                rb.surface = e;
                rb.exact = true;
            }
            else
            {
                // the stations of the ball, exact: from the marched centre line and its pcurves on
                // the offsets, Gauss-Newton on O1(U1, V1) = O2(U2, V2) and (O1 - c) . c' = 0 (the
                // centre in the plane across the marched line at s; Choi and Ju). The marched
                // centre line is accurate to its tolerance divided by the sine of the angle of the
                // offsets, poor where the faces become tangent; its contacts would no longer be at
                // the same distance from it and the arcs no circles. The marched station is kept
                // where Gauss-Newton does not converge or goes further than r / 10 (another root)
                const T scale = T(1) + norm(rb.spine->value(rb.a)) + r;
                auto station = [&](T s) {
                    const auto c0 = rb.spine->value(s);
                    const auto t0 = unit(rb.spine->value(s, 1));
                    const auto w1 = pc[0]->value(s), w2 = pc[1]->value(s);
                    thick_detail::Vec<T> x(4);
                    x << w1[0], w1[1], w2[0], w2[1];
                    const T res = thick_detail::gauss_newton<T>(x, [&](const thick_detail::Vec<T> &y, thick_detail::Vec<T> &f, thick_detail::Mat<T> &J) {
                        std::array<point<T, 3>, 4> D1, D2;
                        off[0].surface->derivatives(y[0], y[1], 1, 1, D1.data());
                        off[1].surface->derivatives(y[2], y[3], 1, 1, D2.data());
                        // (derivatives: [0] value, [1] d/dv, [2] d/du, [3] d2/dudv)
                        const auto d = D1[0] - D2[0];
                        f.resize(4);
                        J.resize(4, 4);
                        f << d[0], d[1], d[2], (D1[0] - c0) * t0;
                        for (int i = 0; i < 3; ++i)
                            J.row(i) << D1[2][i], D1[1][i], -D2[2][i], -D2[1][i];
                        J.row(3) << D1[2] * t0, D1[1] * t0, T(0), T(0);
                    });
                    if (const T jump = norm(off[0].surface->value(x[0], x[1]) - c0); !(res <= T(1e-12) * scale) || jump > T(0.1) * r)
                        x << w1[0], w1[1], w2[0], w2[1];
                    Station<T> st;
                    st.c = off[0].surface->value(x[0], x[1]);
                    st.uv1 = {x[0], x[1] - off[0].dv};
                    st.uv2 = {x[2], x[3] - off[1].dv};
                    st.p1 = sup[0].surface->value(st.uv1[0], st.uv1[1]);
                    st.p2 = sup[1].surface->value(st.uv2[0], st.uv2[1]);
                    return st;
                };
                const auto sk = skinned(station, std::array<const Surface<T, 3> *, 2>{sup[0].surface.get(), sup[1].surface.get()}, rb.a, rb.b,
                                        o.approx_tol, o.angular_tol);
                rb.surface = sk.surface;
                rb.deviation = sk.deviation;
                rb.spine = sk.spine;
                for (std::size_t k = 0; k < 2; ++k)
                {
                    rb.contacts[k].pcurve = sk.pcurves[k];
                    rb.contacts[k].curve = std::make_shared<CurveOnSurface<T, 3>>(sk.pcurves[k], rb.contacts[k].support);
                }
            }
            on_blend(rb, io.tol);
            return rb;
        }

        // ---- chamfers

        /// The curve of the edge on [lo, hi] past its ends: a line or a B-spline continued; [u1, u2] otherwise.
        template <std::floating_point T>
        auto continued_edge(const Edge<T> &ed, T ext, T &lo, T &hi) -> std::shared_ptr<Curve<T, 3>>
        {
            lo = ed.u1, hi = ed.u2;
            if (!(ext > T(0)))
                return ed.curve;
            const T d1 = ext / norm(ed.curve->value(ed.u1, 1)), d2 = ext / norm(ed.curve->value(ed.u2, 1));
            if (const auto *l = dynamic_cast<const LineCurve<T> *>(ed.curve.get()))
            {
                lo -= d1, hi += d2;
                return std::make_shared<LineCurve<T>>(l->origin(), l->vector(), lo, hi);
            }
            auto spline = [&]<bool rational>(const BSCurveGeneral<T, 3, rational> &b) -> std::shared_ptr<Curve<T, 3>> {
                const auto bb = b.bounds();
                const std::array<T, 2> dt{std::max(T(0), bb[0] - (ed.u1 - d1)), std::max(T(0), (ed.u2 + d2) - bb[1])};
                try
                {
                    auto x = std::make_shared<std::conditional_t<rational, BSCurveRational<T, 3>, BSCurve<T, 3>>>(extrapolated(b, dt));
                    lo -= d1, hi += d2;
                    return x;
                }
                catch (const std::invalid_argument &)
                {
                    return ed.curve;
                }
            };
            if (const auto *b = dynamic_cast<const BSCurve<T, 3> *>(ed.curve.get()))
                return spline(*b);
            if (const auto *b = dynamic_cast<const BSCurveRational<T, 3> *>(ed.curve.get()))
                return spline(*b);
            return ed.curve;
        }

        /**
         * The point of `s` at the distance d from e(t) in the plane across the curve e at t,
         * from uv (Gauss-Newton on (|S - e|^2 - d^2) / 2d and (S - e) . e' / |e'|); nullopt
         * when it does not converge within tol / 100 or lands on the far side (b: the direction
         * of the face away from the edge).
         */
        template <std::floating_point T>
        auto foot(const Surface<T, 3> &s, const Curve<T, 3> &e, T t, T d, point<T, 2> uv, T tol) -> std::optional<point<T, 2>>
        {
            using thick_detail::Vec;
            using thick_detail::Mat;
            const auto p = e.value(t);
            const auto te = unit(e.value(t, 1));
            Vec<T> x(2);
            x << uv[0], uv[1];
            const T res = thick_detail::gauss_newton<T>(x, [&](const Vec<T> &y, Vec<T> &r, Mat<T> &J) {
                const auto q = s.value(y[0], y[1]) - p;
                const auto su = s.value(y[0], y[1], 1, 0), sv = s.value(y[0], y[1], 0, 1);
                r.resize(2);
                J.resize(2, 2);
                r << (q * q - d * d) / (T(2) * d), q * te;
                J << (q * su) / d, (q * sv) / d, su * te, sv * te;
            });
            if (!(res <= tol / T(100)))
                return std::nullopt;
            return point<T, 2>{x[0], x[1]};
        }

        /**
         * The contact at the distance d on a side, on [lo, hi]: the feet at n + 1 stations, each
         * from the one before (the first, at the middle of the edge, from its pcurve moved by d
         * toward the face), the move of the edge brought onto the tangent plane, the step halved
         * where it does not converge; interpolated (cubic); n doubled from 8 until the pcurve is
         * within approx_tol of the feet at the middles.
         */
        template <std::floating_point T>
        auto contact_pcurve(const Side<T> &sd, const Surface<T, 3> &s, const Curve<T, 3> &e, const Edge<T> &ed, T lo, T hi, T d, const FilletOptions<T> &o,
                            const std::string &op, EdgeId eid) -> std::shared_ptr<Curve<T, 2>>
        {
            // the seed: the pcurve of the edge at its middle moved by d along the face, away from it
            const T tm = (ed.u1 + ed.u2) / T(2);
            const auto uv0 = sd.pcurve->value(tm);
            const auto su = s.value(uv0[0], uv0[1], 1, 0), sv = s.value(uv0[0], uv0[1], 0, 1);
            auto dd = sd.pcurve->value(tm, 1);
            if (!sd.forward)
                dd = T(-1) * dd;
            const point<T, 2> left{-dd[1], dd[0]};
            auto b = left[0] * su + left[1] * sv;
            const auto te = unit(e.value(tm, 1));
            b = unit(b - (b * te) * te);
            // (u, v) of d b on the tangent plane, by the metric
            const T E = su * su, F = su * sv, G = sv * sv, det = E * G - F * F;
            const T bu = d * (G * (b * su) - F * (b * sv)) / det, bv = d * (E * (b * sv) - F * (b * su)) / det;
            const auto first = foot(s, e, tm, d, point<T, 2>{uv0[0] + bu, uv0[1] + bv}, o.tol);
            if (!first || (s.value((*first)[0], (*first)[1]) - e.value(tm)) * b < T(0.5) * d)
                refuse<T>(BuildErrc::NonGeneric, op, "no point of a face at the distance asked from the edge", {eid, sd.face});
            // the foot at t from the one at t0: the move of the edge brought onto the tangent plane at
            // the foot (first-order predictor), the step halved where Gauss-Newton does not converge
            auto predicted = [&](T t, T t0, const point<T, 2> &uv) {
                const auto a = s.value(uv[0], uv[1], 1, 0), c = s.value(uv[0], uv[1], 0, 1);
                const auto dp = e.value(t) - e.value(t0);
                const T E = a * a, F = a * c, G = c * c, det = E * G - F * F;
                if (!(std::abs(det) > T(0)))
                    return uv;
                return point<T, 2>{uv[0] + (G * (dp * a) - F * (dp * c)) / det, uv[1] + (E * (dp * c) - F * (dp * a)) / det};
            };
            auto march = [&](auto &&self, T t, T t0, const point<T, 2> &uv0, int depth) -> std::optional<point<T, 2>> {
                if (auto x = foot(s, e, t, d, predicted(t, t0, uv0), o.tol))
                    return x;
                if (depth == 0)
                    return std::nullopt;
                const T tm = (t + t0) / T(2);
                const auto m = self(self, tm, t0, uv0, depth - 1);
                return m ? self(self, t, tm, *m, depth - 1) : std::nullopt;
            };
            auto solve = [&](T t, T t0, const point<T, 2> &seed) {
                const auto x = march(march, t, t0, seed, 10);
                if (!x)
                    refuse<T>(BuildErrc::NonGeneric, op, "no point of a face at the distance asked from the edge", {eid, sd.face});
                return *x;
            };
            for (std::size_t n = 8; n <= 1024; n *= 2)
            {
                std::vector<T> t(n + 1);
                std::vector<point<T, 2>> uv(n + 1);
                for (std::size_t i = 0; i <= n; ++i)
                    t[i] = i == n ? hi : lo + (hi - lo) * T(i) / T(n);
                // from the station nearest the middle of the edge, both ways
                const std::size_t i0 = std::size_t(std::clamp<T>(std::round((tm - lo) / (hi - lo) * T(n)), T(0), T(n)));
                uv[i0] = solve(t[i0], tm, *first);
                for (std::size_t i = i0 + 1; i <= n; ++i)
                    uv[i] = solve(t[i], t[i - 1], uv[i - 1]);
                for (std::size_t i = i0; i-- > 0;)
                    uv[i] = solve(t[i], t[i + 1], uv[i + 1]);
                std::vector<constrType<T, 2, 1>> q(n + 1);
                for (std::size_t i = 0; i <= n; ++i)
                    q[i] = {uv[i]};
                auto pc = std::make_shared<BSCurve<T, 2>>(interpolate<T, 2>(q, t, std::size_t{3}));
                T dev{};
                for (std::size_t i = 0; i < n && dev <= o.approx_tol; ++i)
                {
                    const T tt = (t[i] + t[i + 1]) / T(2);
                    const auto y = pc->value(tt);
                    const auto x = solve(tt, t[i], uv[i]);
                    dev = std::max(dev, norm(s.value(y[0], y[1]) - s.value(x[0], x[1])));
                }
                if (dev <= o.approx_tol)
                    return pc;
            }
            refuse<T>(BuildErrc::NonGeneric, op, "a contact cannot be approximated within approx_tol", {eid, sd.face});
        }

        /// The circle through p, q, w: centre, unit normal, radius; nullopt for aligned points.
        template <std::floating_point T>
        auto circle3(const point<T, 3> &p, const point<T, 3> &q, const point<T, 3> &w) -> std::optional<std::tuple<point<T, 3>, point<T, 3>, T>>
        {
            const auto a = q - p, b = w - p;
            const auto n = cross(a, b);
            const T nn = n * n;
            if (!(nn > T(1e-24) * (a * a) * (b * b)))
                return std::nullopt;
            const auto c = p + ((b * b) * cross(n, a) + (a * a) * cross(b, n)) / (T(2) * nn);
            return std::tuple{c, n / std::sqrt(nn), norm(p - c)};
        }

        /**
         * The ruled surface between the contacts c1 and c2 on [a, b] as a primitive: a plane
         * between two lines in a plane, a plane, a cylinder or a cone between two coaxial circles;
         * nullptr otherwise.
         */
        template <std::floating_point T>
        auto ruled_primitive(const Curve<T, 3> &c1, const Curve<T, 3> &c2, T a, T b, T tol) -> std::shared_ptr<ElementarySurface<T>>
        {
            constexpr int n = 32;
            std::array<std::vector<point<T, 3>>, 2> P;
            for (int i = 0; i <= n; ++i)
            {
                const T s = a + (b - a) * T(i) / T(n);
                P[0].push_back(c1.value(s));
                P[1].push_back(c2.value(s));
            }
            auto plane_box = [&](const Frame3<T> &f) {
                std::array<T, 4> box{std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity(),
                                     -std::numeric_limits<T>::infinity()};
                for (const auto &Q : P)
                    for (const auto &p : Q)
                    {
                        const T u = (p - f.o) * f.x, v = (p - f.o) * f.y;
                        box = {std::min(box[0], u), std::max(box[1], u), std::min(box[2], v), std::max(box[3], v)};
                    }
                return box;
            };
            // two lines in a plane
            auto straight = [&](const std::vector<point<T, 3>> &Q) {
                const auto d = Q.back() - Q.front();
                if (!(norm(d) > tol))
                    return false;
                const auto u = unit(d);
                return std::ranges::all_of(Q, [&](const point<T, 3> &p) { return norm(cross(p - Q.front(), u)) <= tol; });
            };
            if (straight(P[0]) && straight(P[1]))
            {
                const auto x = unit(P[0].back() - P[0].front());
                const auto z = cross(x, P[1].front() - P[0].front()), z2 = cross(x, P[1].back() - P[0].front());
                const auto nz = norm(z) >= norm(z2) ? z : z2;
                if (!(norm(nz) > tol))
                    return nullptr;
                const Frame3<T> f{ax2<T, 3>{P[0].front(), unit(nz), x}};
                if (std::ranges::all_of(P[1], [&](const point<T, 3> &p) { return std::abs((p - f.o) * f.z) <= tol; }))
                    return std::make_shared<Plane<T>>(f, plane_box(f));
                return nullptr;
            }
            // two coaxial circles
            std::array<std::tuple<point<T, 3>, point<T, 3>, T>, 2> circ;
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto ck = circle3(P[k][0], P[k][n / 3], P[k][(2 * n) / 3]);
                if (!ck)
                    return nullptr;
                const auto &[cc, cn, cr] = *ck;
                for (const auto &p : P[k])
                    if (std::abs(norm(p - cc) - cr) > tol || std::abs((p - cc) * cn) > tol)
                        return nullptr;
                circ[k] = *ck;
            }
            const auto &[o1, n1, r1] = circ[0];
            const auto &[o2, n2, r2] = circ[1];
            if (norm(cross(n1, n2)) > T(1e-9) || norm(cross(o2 - o1, n1)) > tol)
                return nullptr;
            const T h = (o2 - o1) * n1;
            const auto z = h >= T(0) ? n1 : T(-1) * n1;
            const Frame3<T> f{ax2<T, 3>{o1, z, P[0].front() - o1}};
            if (std::abs(h) <= tol)
                return std::make_shared<Plane<T>>(f, plane_box(f));
            // the angles of the contacts, continuous along them
            T u1 = std::numeric_limits<T>::infinity(), u2 = -u1;
            for (const auto &Q : P)
            {
                T prev{};
                for (std::size_t i = 0; i < Q.size(); ++i)
                {
                    const auto q = Q[i] - f.o;
                    T u = std::atan2(q * f.y, q * f.x);
                    if (i > 0)
                        u += T(2) * std::numbers::pi_v<T> * std::round((prev - u) / (T(2) * std::numbers::pi_v<T>));
                    prev = u;
                    u1 = std::min(u1, u), u2 = std::max(u2, u);
                }
            }
            const std::array<T, 4> box{u1, u2, T(0), std::abs(h)};
            if (std::abs(r1 - r2) <= tol)
                return std::make_shared<CylindricalSurface<T>>(r1, f, box);
            return std::make_shared<ConicalSurface<T>>(r1, std::atan((r2 - r1) / std::abs(h)), f, box);
        }

        /// The ruled B-spline between c1 (v = 0) and c2 (v = 1), u along on [a, b]: both interpolated at the same stations.
        template <std::floating_point T>
        auto ruled_spline(const Curve<T, 3> &c1, const Curve<T, 3> &c2, T a, T b, T approx_tol, T &deviation) -> std::shared_ptr<Surface<T, 3>>
        {
            for (std::size_t n = 8; n <= 2048; n *= 2)
            {
                std::vector<T> t(n + 1);
                std::array<points_vector<T, 3>, 2> q;
                for (std::size_t i = 0; i <= n; ++i)
                {
                    t[i] = i == n ? b : a + (b - a) * T(i) / T(n);
                    q[0].push_back(c1.value(t[i]));
                    q[1].push_back(c2.value(t[i]));
                }
                const auto b1 = interpolate<T, 3>(q[0], t, std::size_t{3}), b2 = interpolate<T, 3>(q[1], t, std::size_t{3});
                T dev{};
                for (std::size_t i = 0; i < n && dev <= approx_tol; ++i)
                    for (int j = 1; j < 4; ++j)
                    {
                        const T s = t[i] + (t[i + 1] - t[i]) * T(j) / T(4);
                        dev = std::max({dev, norm(b1.value(s) - c1.value(s)), norm(b2.value(s) - c2.value(s))});
                    }
                if (dev <= approx_tol)
                {
                    deviation = dev;
                    auto poles = b1.poles();
                    const auto p2 = b2.poles();
                    poles.insert(poles.end(), p2.begin(), p2.end());
                    return std::make_shared<BSSurface<T, 3>>(poles, b1.knotsFlats(), std::vector<T>{0, 0, 1, 1}, b1.degree(), std::size_t{1});
                }
            }
            throw BuildFailure(BuildError{BuildErrc::NonGeneric, "chamfer: the surface cannot be approximated within approx_tol", {}});
        }

        template <std::floating_point T>
        auto chamfer(const Model<T> &m, SolidId sid, EdgeId eid, T d1, T d2, const FilletOptions<T> &o) -> EdgeBlend<T>
        {
            const std::string op = "chamfer_geometry";
            if (!(d1 > o.tol) || !(d2 > o.tol) || !std::isfinite(d1) || !std::isfinite(d2))
                refuse<T>(BuildErrc::InvalidBounds, op, "the distances must be larger than the tolerance", {eid});
            const auto st = setup(m, sid, eid, op);
            const auto &ed = m.edge(eid);
            const T tol = o.tol;
            const T ext = st.closed ? T(0) : o.extension * std::max(d1, d2);
            EdgeBlend<T> rb;
            rb.edge = eid;
            rb.convex = st.convex;
            rb.closed = st.closed;
            rb.spine = continued_edge(ed, ext, rb.a, rb.b);
            rb.t1 = ed.u1, rb.t2 = ed.u2;
            const std::array<T, 2> d{d1, d2};
            for (std::size_t k = 0; k < 2; ++k)
            {
                const auto &sd = st.sides[k];
                auto &ct = rb.contacts[k];
                ct.face = sd.face;
                ct.support = support(m, sd, ed, T(2) * d[k] + ext, tol).surface;
                ct.pcurve = contact_pcurve(sd, *ct.support, *rb.spine, ed, rb.a, rb.b, d[k], o, op, eid);
                ct.curve = std::make_shared<CurveOnSurface<T, 3>>(ct.pcurve, ct.support);
                check_on_face(m, sd, *ct.pcurve, ed.u1, ed.u2, tol, op, eid);
            }
            if (auto e = ruled_primitive(*rb.contacts[0].curve, *rb.contacts[1].curve, rb.a, rb.b, tol))
            {
                rb.surface = e;
                rb.exact = true;
            }
            else
                rb.surface = ruled_spline(*rb.contacts[0].curve, *rb.contacts[1].curve, rb.a, rb.b, o.approx_tol, rb.deviation);
            // the contacts on the ruled B-spline: its sides v = 0 and v = 1
            if (!rb.exact)
            {
                for (std::size_t k = 0; k < 2; ++k)
                {
                    auto &c = rb.contacts[k];
                    const T v = T(k);
                    c.on_blend = segment<T>({rb.a, v}, {rb.b, v}, rb.a, rb.b);
                    T dev{};
                    for (int i = 0; i <= 256; ++i)
                    {
                        const T s = rb.a + (rb.b - rb.a) * T(i) / T(256);
                        dev = std::max(dev, norm(rb.surface->value(s, v) - c.curve->value(s)));
                    }
                    c.deviation = dev;
                }
            }
            else
                on_blend(rb, std::max(o.intersect.tol, tol / T(10)));
            return rb;
        }

        template <typename F>
        auto guarded(F &&f) -> BuildResult<std::invoke_result_t<F>>
        {
            try
            {
                return f();
            }
            catch (const BuildFailure &e)
            {
                return std::unexpected(e.error());
            }
        }
    }

    /**
     * @brief The geometry of the fillet of radius r along the edge `e` of the solid `s` (file
     * header): its centre line, its contacts on the two faces of the edge and its surface,
     * continued by `extension` r past the vertices of an open edge. The model is not modified.
     */
    template <std::floating_point T>
    [[nodiscard]] auto fillet_geometry(const Model<T> &m, SolidId s, EdgeId e, std::type_identity_t<T> r, const FilletOptions<T> &o = {})
        -> BuildResult<EdgeBlend<T>>
    {
        return fillet_detail::guarded([&] { return fillet_detail::fillet(m, s, e, r, o); });
    }

    /**
     * @brief The geometry of the chamfer along the edge `e` of the solid `s`, at the distances d1
     * on the face of the first use of the edge and d2 on the other (file header): its contacts and
     * its ruled surface, continued by `extension` max(d1, d2) past the vertices of an open edge.
     * The model is not modified.
     */
    template <std::floating_point T>
    [[nodiscard]] auto chamfer_geometry(const Model<T> &m, SolidId s, EdgeId e, std::type_identity_t<T> d1, std::type_identity_t<T> d2,
                                        const FilletOptions<T> &o = {}) -> BuildResult<EdgeBlend<T>>
    {
        return fillet_detail::guarded([&] { return fillet_detail::chamfer(m, s, e, d1, d2, o); });
    }
}
