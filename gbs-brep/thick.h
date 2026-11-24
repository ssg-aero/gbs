#pragma once
/**
 * @file thick.h
 * @brief Offset of a solid with sharp joints (`offset_solid`), and thick solids, hollowed out to
 * a thickness and opened by removed faces (`make_thick`).
 *
 * Each face is carried by its offset surface (`offset_surface`: exact for the primitives, an
 * `OffsetSurface` otherwise) at its distance along the outward normal of the solid; the
 * topology is that of the solid:
 * 1. **Edges.** An edge where the offsets of its two faces meet along a common curve (the faces
 *    tangent along it at the same distance, a seam, two faces left in place) is that curve: its
 *    pcurve on the offset surface, continued past its ends when the edge gets longer. Any other
 *    edge is sharp: the intersection of the two offset surfaces (closed forms between primitives,
 *    marching otherwise), the branch through its two new vertices (on a loop, the arc nearest the
 *    old edge; a closed edge, the loop run from its vertex). The surfaces are first continued
 *    past the face by what the joint needs (|d2 - d1 cos a| / sin a on the first face, a the
 *    angle of the normals): a primitive on a wider box, a B-spline by the polynomials of its end
 *    spans, but across its seam or past a pole.
 * 2. **Vertices.** The point common to the offset surfaces around the vertex: on the offset of
 *    its tangent edge when it has one, met by the other faces; the point where every face meets
 *    when all its edges are sharp (Gauss-Newton, three faces or more agreeing; between two faces,
 *    a vertex between two arcs of a circle, the point of the intersection of their offsets in the
 *    plane across it). A pole of a face at a sharp vertex is left out of its wire.
 * 3. The faces, wires, shells and solid of the solid on the new entities, checked (`check`, the
 *    offset regular on each new face, and the sections of the pairs of faces sharing no vertex
 *    whose boxes overlap: none may meet).
 *
 * A thick solid is the solid minus its offset inward by the thickness, the removed faces offset
 * by 0: they lie on the faces of the solid, and the boolean of coincident faces opens it there.
 * A negative thickness grows outward.
 *
 * Refused (`NonGeneric`): an offset that folds (beyond a radius of curvature: the apex of a
 * cone), a face that vanishes (an edge turned back), an offset intersecting itself (two faces
 * sharing no vertex that meet), a vertex where the offsets do not meet at one point (more than
 * three sharp faces not agreeing, tangent edges not through the same point), a face whose offset
 * would cover the pole of its surface at a sharp vertex (a sphere cut through its poles, offset
 * outward), a joint too sharp to be continued, an edge partly tangent; the model is left as it
 * was.
 */
#include <gbs-brep/boolean.h>
#include <gbs-brep/check.h>
#include <gbs-brep/explore.h>
#include <gbs-brep/offset.h>
#include <gbs-brep/section.h>
#include <gbs/bsctools.h>
#include <gbs/bsstools.h>
#include <gbs/curvecomposite.h>
#include <gbs/curvetrimmed.h>
#include <gbs/extrema.h>
#include <gbs/intersect/surface_surface.h>

#include <Eigen/Dense>

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

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs::brep
{
    template <std::floating_point T>
    struct OffsetOptions
    {
        T tol{brep_default_tolerance<T>}; ///< tolerance of the result; edges whose offsets part by less are tangent
        IntersectOptions<T> intersect{};  ///< intersection of the offset surfaces (its tol is raised to `tol` / 10 at least)
        T max_extension{T(100)};          ///< largest continuation of a face past a sharp edge, in |d|
    };

    namespace thick_detail
    {
        template <std::floating_point T>
        [[noreturn]] void refuse(const std::string &why, std::vector<ShapeId> shapes = {})
        {
            throw BuildFailure(BuildError{BuildErrc::NonGeneric, "offset_solid: " + why, std::move(shapes)});
        }

        /// The clamped 2D B-spline (rational or not) continued by dt past both ends (`extrapolated`: the
        /// polynomials of its end spans); any other curve as it is.
        template <std::floating_point T>
        auto continued(const std::shared_ptr<Curve<T, 2>> &pc, T dt) -> std::shared_ptr<Curve<T, 2>>
        {
            if (!(dt > T(0)))
                return pc;
            auto go = [&]<bool rational>(const BSCurveGeneral<T, 2, rational> &b) -> std::shared_ptr<Curve<T, 2>> {
                try
                {
                    return std::make_shared<std::conditional_t<rational, BSCurveRational<T, 2>, BSCurve<T, 2>>>(extrapolated(b, std::array<T, 2>{dt, dt}));
                }
                catch (const std::invalid_argument &)
                {
                    return pc; // knots not clamped, an empty end span
                }
            };
            if (const auto *b = dynamic_cast<const BSCurve<T, 2> *>(pc.get()))
                return go(*b);
            if (const auto *b = dynamic_cast<const BSCurveRational<T, 2> *>(pc.get()))
                return go(*b);
            return pc;
        }

        /// The closed curve c of [a, b] run over a whole turn from t: on [t, t + b - a], its pieces
        /// [t, b] and [a, t] joined (as it is when t is at an end).
        template <std::floating_point T, std::size_t dim>
        auto from_parameter(const std::shared_ptr<Curve<T, dim>> &c, T a, T b, T t) -> std::shared_ptr<Curve<T, dim>>
        {
            const T eps = T(1e-12) * (T(1) + std::abs(a) + std::abs(b));
            if (!(t > a + eps) || !(t < b - eps))
                return c;
            return std::make_shared<CurveComposite<T, dim>>(
                std::vector<std::shared_ptr<Curve<T, dim>>>{std::make_shared<CurveTrimmed<T, dim>>(c, t, b), std::make_shared<CurveTrimmed<T, dim>>(c, a, t)});
        }

        /// The smallest speed |S_u| (k = 0) or |S_v| (k = 1) along the side `side` of the box, its
        /// points slower than 1e-9 of the fastest left out (a pole: its speed is the rounding of 0,
        /// 6e-17 r on a sphere); 0 when none is left.
        template <std::floating_point T>
        T side_speed(const Surface<T, 3> &s, const std::array<T, 4> &box, int k, T side)
        {
            std::array<T, 9> speed{};
            for (int i = 0; i <= 8; ++i)
            {
                const T w = (k == 0 ? box[2] + (box[3] - box[2]) * T(i) / T(8) : box[0] + (box[1] - box[0]) * T(i) / T(8));
                speed[std::size_t(i)] = norm(k == 0 ? s.value(side, w, 1, 0) : s.value(w, side, 0, 1));
            }
            const T largest = std::ranges::max(speed);
            T best = std::numeric_limits<T>::infinity();
            for (T v : speed)
                if (v > T(1e-9) * largest)
                    best = std::min(best, v);
            return std::isfinite(best) ? best : T(0);
        }

        template <std::floating_point T>
        struct Extended
        {
            std::shared_ptr<Surface<T, 3>> surface;
            std::array<T, 4> box{}; ///< the box of the face continued by the length asked, within the surface
            bool extended{};
        };

        /**
         * The surface of a face continued by the length L past each side of the box of the face:
         * a primitive on a wider box (a periodic parameter at most a period, a sphere up to its
         * poles, a cone short of its apex), a B-spline by the polynomials of its end spans, but
         * across a direction where it closes on itself (a NURBS cylinder: its seam) or past a
         * degenerate side (a pole; `surface_closure` within tol). Any other surface is left as it is.
         */
        template <std::floating_point T>
        auto extended_surface(const std::shared_ptr<Surface<T, 3>> &s, const std::array<T, 4> &fbox, T L, T tol) -> Extended<T>
        {
            if (!(L > T(0)))
                return {s, fbox, false};
            std::array<T, 4> want = fbox;
            for (int k = 0; k < 2; ++k)
                for (int e = 0; e < 2; ++e)
                {
                    const T speed = side_speed(*s, fbox, k, fbox[2 * k + e]);
                    if (speed > T(0))
                        want[2 * k + e] += (e == 0 ? -T(1) : T(1)) * L / speed;
                }
            const auto sb = s->bounds();
            if (auto e = std::dynamic_pointer_cast<ElementarySurface<T>>(s))
            {
                const auto periods = e->periods();
                std::array<T, 4> box{std::min(sb[0], want[0]), std::max(sb[1], want[1]), std::min(sb[2], want[2]), std::max(sb[3], want[3])};
                for (int k = 0; k < 2; ++k)
                    if (periods[k] > T(0) && box[2 * k + 1] - box[2 * k] > periods[k])
                    {
                        // at most a period, centred on the face
                        const T c = (fbox[2 * k] + fbox[2 * k + 1]) / T(2);
                        const T w = std::max(periods[k], fbox[2 * k + 1] - fbox[2 * k]);
                        box[2 * k] = std::min(sb[2 * k], c - w / T(2));
                        box[2 * k + 1] = std::max(sb[2 * k + 1], c + w / T(2));
                        if (box[2 * k + 1] - box[2 * k] > periods[k] + T(1e-12) * periods[k])
                            box[2 * k] = sb[2 * k], box[2 * k + 1] = sb[2 * k + 1];
                    }
                if (dynamic_cast<const SphericalSurface<T> *>(e.get()))
                {
                    box[2] = std::max(box[2], -std::numbers::pi_v<T> / T(2));
                    box[3] = std::min(box[3], std::numbers::pi_v<T> / T(2));
                }
                if (const auto *co = dynamic_cast<const ConicalSurface<T> *>(e.get()))
                {
                    // short of the apex v = -R / tan a
                    const T ta = std::tan(co->semi_angle());
                    if (std::abs(ta) > T(0))
                    {
                        const T apex = -co->radius() / ta, margin = T(1e-3) * (fbox[3] - fbox[2]);
                        if (apex < fbox[2])
                            box[2] = std::max(box[2], apex + margin);
                        else if (apex > fbox[3])
                            box[3] = std::min(box[3], apex - margin);
                    }
                }
                if (box == sb)
                    return {s, {std::max(want[0], sb[0]), std::min(want[1], sb[1]), std::max(want[2], sb[2]), std::min(want[3], sb[3])}, false};
                auto n = e->with_box(box);
                return {n, {std::max(want[0], box[0]), std::min(want[1], box[1]), std::max(want[2], box[2]), std::min(want[3], box[3])}, true};
            }
            auto within = [&](const std::array<T, 4> &b) {
                return std::array<T, 4>{std::max(want[0], b[0]), std::min(want[1], b[1]), std::max(want[2], b[2]), std::min(want[3], b[3])};
            };
            auto spline = [&]<bool rational>(const BSSurfaceGeneral<T, 3, rational> &b) -> std::optional<Extended<T>> {
                // past the bounds of the surface only, and neither across its seam nor past a pole
                std::array<T, 4> ext{std::max(T(0), sb[0] - want[0]), std::max(T(0), want[1] - sb[1]), std::max(T(0), sb[2] - want[2]),
                                     std::max(T(0), want[3] - sb[3])};
                const auto cl = surface_closure(*s, tol);
                if (cl.closed_u || cl.degenerate_u1)
                    ext[0] = T(0);
                if (cl.closed_u || cl.degenerate_u2)
                    ext[1] = T(0);
                if (cl.closed_v || cl.degenerate_v1)
                    ext[2] = T(0);
                if (cl.closed_v || cl.degenerate_v2)
                    ext[3] = T(0);
                if (ext == std::array<T, 4>{})
                    return Extended<T>{s, within(sb), false};
                try
                {
                    auto v = extrapolated_v(b, std::array<T, 2>{ext[2], ext[3]});
                    v.invertUV();
                    auto u = extrapolated_v(v, std::array<T, 2>{ext[0], ext[1]});
                    u.invertUV();
                    auto n = std::make_shared<std::remove_cvref_t<decltype(u)>>(std::move(u));
                    return Extended<T>{n, within(n->bounds()), true};
                }
                catch (const std::invalid_argument &)
                {
                    return std::nullopt;
                }
            };
            if (const auto *b = dynamic_cast<const BSSurface<T, 3> *>(s.get()))
                if (auto r = spline(*b))
                    return *r;
            if (const auto *b = dynamic_cast<const BSSurfaceRational<T, 3> *>(s.get()))
                if (auto r = spline(*b))
                    return *r;
            return {s, within(sb), false};
        }

        /// A face of the solid and its offset.
        template <std::floating_point T>
        struct OffsetFace
        {
            FaceId face;
            Orientation orient{};
            T d{};                                  ///< distance along the normal of the surface
            std::shared_ptr<Surface<T, 3>> surface; ///< the surface of the face, continued
            std::shared_ptr<Surface<T, 3>> offset;  ///< its offset, parameters (u, v + dv)
            T dv{};
            std::array<T, 4> box{};                 ///< on the parameters of the offset
            bool extended{};
            T length{};                              ///< continuation asked past its edges
        };

        enum class EdgeKind
        {
            Degenerate,
            Tangent, ///< the offsets of its faces meet along one curve
            Sharp
        };

        template <std::floating_point T>
        struct EdgeInfo
        {
            EdgeId edge;
            std::vector<EdgeUse> uses;
            EdgeKind kind{EdgeKind::Sharp};
        };

        /// A face around a vertex, and the parameters of the vertex on its offset.
        template <std::floating_point T>
        struct Corner
        {
            std::size_t face{}; ///< index in the offset faces
            std::array<T, 2> uv{};
        };

        template <std::floating_point T>
        struct NewVertex
        {
            point<T, 3> p{};
            std::map<std::pair<std::uint32_t, int>, T> t; ///< parameter on the continued tangent edges ending there (by edge index, end 0 or 1)
            bool sharp{};                                 ///< met by a sharp edge: found where the offsets meet
        };

        template <std::floating_point T>
        using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
        template <std::floating_point T>
        using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

        /**
         * Gauss-Newton on r(x) = 0, J its Jacobian (`f` fills both), the step solved by QR (least
         * squares when there are more equations), halved while it does not lower the residual; a
         * step where f throws (a parameter beyond the bounds of a curve) is halved as well. The
         * largest residual at the end (infinite if it cannot be evaluated at the start).
         */
        template <std::floating_point T, typename F>
        T gauss_newton(Vec<T> &x, F &&f, int iterations = 60)
        {
            auto eval = [&](const Vec<T> &y, Vec<T> &r, Mat<T> &J) {
                try
                {
                    f(y, r, J);
                    return r.allFinite() && J.allFinite();
                }
                catch (const std::exception &)
                {
                    return false;
                }
            };
            Vec<T> r;
            Mat<T> J;
            if (!eval(x, r, J))
                return std::numeric_limits<T>::infinity();
            for (int it = 0; it < iterations; ++it)
            {
                const Vec<T> dx = J.colPivHouseholderQr().solve(-r);
                if (!dx.allFinite())
                    break;
                T step = T(1);
                bool moved = false;
                Vec<T> ry;
                Mat<T> Jy;
                for (int h = 0; h < 12 && !moved; ++h, step /= T(2))
                    moved = eval(x + step * dx, ry, Jy) && ry.norm() <= r.norm() * (T(1) + T(1e-12));
                if (!moved)
                    break;
                step *= T(2); // the step taken
                x += step * dx;
                r = std::move(ry);
                J = std::move(Jy);
                if (dx.norm() * step <= T(1e-15) * (T(1) + x.norm()))
                    break;
            }
            return r.template lpNorm<Eigen::Infinity>();
        }

        /**
         * The parameters where Gauss-Newton starts on the offset s near the point x: uv as given
         * where the parametrization is regular; at a degenerate point (a pole: S_u x S_v vanishes,
         * the Jacobian is singular there), the parameters of x on s (closed form for a primitive,
         * the nearest point otherwise).
         */
        template <std::floating_point T>
        auto start_uv(const Surface<T, 3> &s, const std::array<T, 2> &uv, const point<T, 3> &x) -> std::array<T, 2>
        {
            const auto su = s.value(uv[0], uv[1], 1, 0), sv = s.value(uv[0], uv[1], 0, 1);
            const T l = std::max(norm(su), norm(sv));
            if (norm(cross(su, sv)) > T(1e-9) * l * l)
                return uv;
            if (const auto *e = dynamic_cast<const ElementarySurface<T> *>(&s))
                return e->parameters(x);
            const auto r = extrema_surf_pnt(s, x, T(1e-12));
            return {r[0], r[1]};
        }
    }

    namespace thick_detail
    {
        /**
         * The solid `s` with each face carried by its offset at distance(face) along the outward
         * normal (sharp joints). Throws `BuildFailure` on a refusal; the caller rolls back.
         */
        template <std::floating_point T, typename Distance>
        auto offset_faces(Model<T> &m, SolidId sid, Distance &&distance, const OffsetOptions<T> &o) -> SolidId
        {
            const T tol = o.tol;
            const auto &so = m.solid(sid);
            std::vector<ShellId> shells{so.outer};
            shells.insert(shells.end(), so.voids.begin(), so.voids.end());

            // ---- the faces, their distances along the normals of their surfaces
            std::vector<OffsetFace<T>> faces;
            std::map<std::uint32_t, std::size_t> face_index;
            T dmax{};
            for (auto sh : shells)
                for (const auto &fu : m.shell(sh).faces)
                {
                    if (face_index.contains(fu.face.index))
                        refuse<T>("a face used twice", {fu.face});
                    const T d = distance(fu.face);
                    OffsetFace<T> f;
                    f.face = fu.face;
                    f.orient = fu.orient;
                    f.d = fu.orient == Orientation::Forward ? d : -d;
                    f.surface = m.face(fu.face).surface;
                    dmax = std::max(dmax, std::abs(d));
                    face_index[fu.face.index] = faces.size();
                    faces.push_back(std::move(f));
                }
            auto outward = [&](const OffsetFace<T> &f, const std::array<T, 2> &uv) {
                const auto n = normal_direction(*f.surface, uv[0], uv[1]);
                return f.orient == Orientation::Forward ? n : point<T, 3>{-n[0], -n[1], -n[2]};
            };

            // ---- the edges: tangent or sharp, and how far their faces go on past them
            std::vector<EdgeInfo<T>> edges;
            std::map<std::uint32_t, std::size_t> edge_index;
            for (auto sh : shells)
                for (auto &[e, uses] : shell_edge_uses(m, sh))
                {
                    if (uses.size() != 2)
                        refuse<T>("an edge not shared by two faces", {e});
                    edge_index[e.index] = edges.size();
                    edges.push_back({e, uses, EdgeKind::Sharp});
                }
            auto pcurve_of = [&](const CoEdgeRef &r) { return m.wire(r.wire).coedges[r.index].pcurve; };
            for (auto &ei : edges)
            {
                const auto &e = m.edge(ei.edge);
                const auto &fa = faces[face_index.at(ei.uses[0].face.index)], &fb = faces[face_index.at(ei.uses[1].face.index)];
                const auto pa = pcurve_of(ei.uses[0].coedge), pb = pcurve_of(ei.uses[1].coedge);
                if (!pa || !pb)
                    refuse<T>("an edge without pcurves", {ei.edge});
                int tangent = 0;
                const int n = 9;
                T along_a{}, along_b{};
                for (int i = 0; i < n; ++i)
                {
                    const T t = e.u1 + (e.u2 - e.u1) * (T(i) + T(0.5)) / T(n);
                    const auto na = outward(fa, pa->value(t)), nb = outward(fb, pb->value(t));
                    const T da = fa.orient == Orientation::Forward ? fa.d : -fa.d, db = fb.orient == Orientation::Forward ? fb.d : -fb.d;
                    if (norm(da * na - db * nb) <= tol)
                    {
                        ++tangent;
                        continue;
                    }
                    const T c = std::clamp(na * nb, T(-1), T(1)), s = std::sqrt(std::max(T(0), T(1) - c * c));
                    if (!(s > T(1e-12)))
                        refuse<T>(c > T(0) ? "an edge between tangent faces at different distances" : "an edge whose faces meet back to back", {ei.edge});
                    along_a = std::max(along_a, std::abs(db - da * c) / s);
                    along_b = std::max(along_b, std::abs(da - db * c) / s);
                }
                if (tangent == n)
                    ei.kind = EdgeKind::Tangent;
                else if (tangent > 0)
                    refuse<T>("an edge partly tangent", {ei.edge});
                else
                {
                    if (std::max(along_a, along_b) > o.max_extension * std::max(dmax, tol))
                        refuse<T>("a joint too sharp to be continued", {ei.edge});
                    // a vertex of three faces goes further than its edges: twice
                    faces[face_index.at(ei.uses[0].face.index)].length = std::max(fa.length, T(2) * along_a + T(10) * tol);
                    faces[face_index.at(ei.uses[1].face.index)].length = std::max(fb.length, T(2) * along_b + T(10) * tol);
                }
            }
            // degenerate edges (poles), from the wires
            std::vector<std::pair<EdgeId, std::pair<std::size_t, CoEdgeRef>>> degenerate;
            for (std::size_t k = 0; k < faces.size(); ++k)
                for (auto w : m.face(faces[k].face).wires)
                    for (std::uint32_t i = 0; i < m.wire(w).coedges.size(); ++i)
                        if (m.edge(m.wire(w).coedges[i].edge).degenerate)
                            degenerate.push_back({m.wire(w).coedges[i].edge, {k, CoEdgeRef{w, i}}});

            // ---- the offset surfaces, continued
            for (auto &f : faces)
            {
                const auto box = classify_detail::face_uv_box(m, f.face);
                auto ext = extended_surface(f.surface, box, f.length, tol);
                f.extended = ext.extended;
                f.surface = ext.surface;
                if (f.d == T(0))
                {
                    f.offset = f.surface;
                    f.box = ext.box;
                    continue;
                }
                // the fold on the face (its continuation may fold where no new edge goes: the new
                // face is checked at the end)
                auto off = offset_surface(f.surface, f.d, tol / T(10), box);
                if (!off)
                    refuse<T>("the offset folds (beyond a radius of curvature)", {f.face});
                f.offset = off->surface;
                f.dv = off->dv;
                f.box = {ext.box[0], ext.box[1], ext.box[2] + f.dv, ext.box[3] + f.dv};
            }
            auto shifted = [&](const OffsetFace<T> &f, const std::shared_ptr<Curve<T, 2>> &pc) {
                auto s = offset_detail::shifted(pc, f.dv);
                if (!s)
                    refuse<T>("a pcurve that cannot be shifted", {f.face});
                return s;
            };
            // the pcurves of the tangent edges on the offsets, continued past their ends
            std::map<CoEdgeRef, std::shared_ptr<Curve<T, 2>>> tangent_pc;
            for (const auto &ei : edges)
                if (ei.kind == EdgeKind::Tangent)
                {
                    const auto &e = m.edge(ei.edge);
                    for (const auto &u : ei.uses)
                    {
                        const auto &f = faces[face_index.at(u.face.index)];
                        const auto pc = shifted(f, pcurve_of(u.coedge));
                        // f.length is a length: in the parameter of the edge, divided by its speed on
                        // the offset at its slower end
                        T speed = std::numeric_limits<T>::infinity();
                        for (T t : {e.u1, e.u2})
                        {
                            const auto q = pc->value(t), dq = pc->value(t, 1);
                            speed = std::min(speed, norm(f.offset->value(q[0], q[1], 1, 0) * dq[0] + f.offset->value(q[0], q[1], 0, 1) * dq[1]));
                        }
                        const T dt = (e.u2 - e.u1) + (speed > T(0) && std::isfinite(speed) ? f.length / speed : f.length);
                        tangent_pc[u.coedge] = continued(pc, dt);
                    }
                }

            // ---- the vertices
            std::map<std::uint32_t, std::vector<Corner<T>>> corners;           // by vertex: the faces around it
            std::map<std::uint32_t, std::vector<std::size_t>> vertex_edges;      // by vertex: its edges (index in edges)
            std::map<std::uint32_t, point<T, 3>> vertex_point;
            for (std::size_t k = 0; k < faces.size(); ++k)
                for (auto w : m.face(faces[k].face).wires)
                    for (const auto &ce : m.wire(w).coedges)
                    {
                        const auto &e = m.edge(ce.edge);
                        if (!ce.pcurve)
                            refuse<T>("a co-edge without pcurve", {ce.edge});
                        for (const auto &[v, t] : {std::pair{e.v1, e.u1}, std::pair{e.v2, e.u2}})
                        {
                            auto uv = ce.pcurve->value(t);
                            uv[1] += faces[k].dv;
                            auto &c = corners[v.index];
                            if (std::ranges::none_of(c, [&](const Corner<T> &x) { return x.face == k; }))
                                c.push_back({k, uv});
                            vertex_point[v.index] = m.vertex(v).pnt;
                            if (!e.degenerate)
                            {
                                auto &ve = vertex_edges[v.index];
                                const auto i = edge_index.at(ce.edge.index);
                                if (!std::ranges::contains(ve, i))
                                    ve.push_back(i);
                            }
                        }
                    }
            // by vertex: the faces with a pole there
            std::map<std::uint32_t, std::set<std::size_t>> pole_faces;
            for (const auto &[eid, use] : degenerate)
                pole_faces[m.edge(eid).v1.index].insert(use.first);
            std::map<std::uint32_t, NewVertex<T>> moved;
            for (const auto &[v, around] : corners)
            {
                NewVertex<T> nv;
                std::vector<std::size_t> tangents;
                bool sharp = false;
                for (auto i : vertex_edges[v])
                {
                    if (edges[i].kind == EdgeKind::Tangent)
                        tangents.push_back(i);
                    sharp = sharp || edges[i].kind == EdgeKind::Sharp;
                }
                const VertexId vid{v};
                if (!sharp)
                {
                    // the offsets of the faces meet there
                    const auto &c0 = around.front();
                    nv.p = faces[c0.face].offset->value(c0.uv[0], c0.uv[1]);
                    for (const auto &c : around)
                        if (norm(faces[c.face].offset->value(c.uv[0], c.uv[1]) - nv.p) > tol)
                            refuse<T>("the offsets of the faces around a vertex part", {vid});
                    for (auto i : tangents)
                    {
                        const auto &e = m.edge(edges[i].edge);
                        if (e.v1 == vid)
                            nv.t[{edges[i].edge.index, 0}] = e.u1;
                        if (e.v2 == vid)
                            nv.t[{edges[i].edge.index, 1}] = e.u2;
                    }
                }
                else if (!tangents.empty())
                {
                    // on the offset of a tangent edge, where the other faces meet it; then on the
                    // offsets of the other tangent edges
                    auto others_of = [&](const EdgeInfo<T> &ei) {
                        std::vector<std::size_t> r;
                        for (const auto &c : around)
                            if (c.face != face_index.at(ei.uses[0].face.index) && c.face != face_index.at(ei.uses[1].face.index))
                                r.push_back(c.face);
                        return r;
                    };
                    const auto lead = std::ranges::find_if(tangents, [&](std::size_t i) { return !others_of(edges[i]).empty(); });
                    if (lead == tangents.end())
                        refuse<T>("a vertex of a tangent edge without another face", {vid});
                    const auto &ei = edges[*lead];
                    const auto &e = m.edge(ei.edge);
                    if (e.v1 == e.v2)
                        refuse<T>("a closed tangent edge through a sharp vertex", {vid});
                    const auto &fa = faces[face_index.at(ei.uses[0].face.index)];
                    const auto pc = tangent_pc.at(ei.uses[0].coedge);
                    const auto others = others_of(ei);
                    const T t0 = e.v1 == vid ? e.u1 : e.u2;
                    const std::size_t n = 1 + 2 * others.size();
                    Vec<T> x(static_cast<Eigen::Index>(n));
                    x[0] = t0;
                    const auto q0 = pc->value(t0);
                    const auto p0 = fa.offset->value(q0[0], q0[1]);
                    for (std::size_t j = 0; j < others.size(); ++j)
                    {
                        const auto c = std::ranges::find(around, others[j], &Corner<T>::face);
                        const auto uv = start_uv(*faces[others[j]].offset, c->uv, p0);
                        x[Eigen::Index(1 + 2 * j)] = uv[0];
                        x[Eigen::Index(2 + 2 * j)] = uv[1];
                    }
                    auto f = [&](const Vec<T> &y, Vec<T> &r, Mat<T> &J) {
                        r.setZero(Eigen::Index(3 * others.size()));
                        J.setZero(Eigen::Index(3 * others.size()), Eigen::Index(n));
                        const auto q = pc->value(y[0]), dq = pc->value(y[0], 1);
                        const auto p = fa.offset->value(q[0], q[1]);
                        const auto dp = fa.offset->value(q[0], q[1], 1, 0) * dq[0] + fa.offset->value(q[0], q[1], 0, 1) * dq[1];
                        for (std::size_t j = 0; j < others.size(); ++j)
                        {
                            const auto &s = *faces[others[j]].offset;
                            const T u = y[Eigen::Index(1 + 2 * j)], w = y[Eigen::Index(2 + 2 * j)];
                            const auto g = s.value(u, w), gu = s.value(u, w, 1, 0), gv = s.value(u, w, 0, 1);
                            for (int a = 0; a < 3; ++a)
                            {
                                const auto row = Eigen::Index(3 * j + std::size_t(a));
                                r[row] = p[a] - g[a];
                                J(row, 0) = dp[a];
                                J(row, Eigen::Index(1 + 2 * j)) = -gu[a];
                                J(row, Eigen::Index(2 + 2 * j)) = -gv[a];
                            }
                        }
                    };
                    const T res = gauss_newton<T>(x, f);
                    if (!(res <= tol / T(10)))
                        refuse<T>("the offsets of the faces around a vertex do not meet", {vid});
                    const auto q = pc->value(x[0]);
                    nv.p = fa.offset->value(q[0], q[1]);
                    nv.sharp = true;
                    nv.t[{ei.edge.index, e.v1 == vid ? 0 : 1}] = x[0];
                    for (auto j : tangents)
                    {
                        if (j == *lead)
                            continue;
                        const auto &ej = edges[j];
                        const auto &ee = m.edge(ej.edge);
                        if (ee.v1 == ee.v2)
                            refuse<T>("a closed tangent edge through a sharp vertex", {vid});
                        const auto &fj = faces[face_index.at(ej.uses[0].face.index)];
                        const auto pj = tangent_pc.at(ej.uses[0].coedge);
                        Vec<T> y(1);
                        y[0] = ee.v1 == vid ? ee.u1 : ee.u2;
                        auto g = [&](const Vec<T> &z, Vec<T> &r, Mat<T> &J) {
                            const auto q = pj->value(z[0]), dq = pj->value(z[0], 1);
                            const auto p = fj.offset->value(q[0], q[1]);
                            const auto dp = fj.offset->value(q[0], q[1], 1, 0) * dq[0] + fj.offset->value(q[0], q[1], 0, 1) * dq[1];
                            r.resize(3);
                            J.resize(3, 1);
                            for (int a = 0; a < 3; ++a)
                                r[a] = p[a] - nv.p[a], J(a, 0) = dp[a];
                        };
                        if (!(gauss_newton<T>(y, g) <= tol / T(10)))
                            refuse<T>("the offsets of the tangent edges around a vertex part", {vid});
                        nv.t[{ej.edge.index, ee.v1 == vid ? 0 : 1}] = y[0];
                    }
                }
                else
                {
                    // where the offsets of every face around meet: x and (u, v) on each. Between two
                    // faces (a vertex between two arcs of a circle, at a pole on a cut), the point of
                    // the intersection of their offsets in the plane through the vertex across it:
                    // (x - p0) . (n1 x n2) = 0, the nearest one for a circle
                    const std::size_t k = around.size();
                    if (k < 2)
                        refuse<T>("a vertex of a single sharp face", {vid});
                    const auto p0 = vertex_point[v];
                    point<T, 3> across{};
                    if (k == 2)
                    {
                        const auto n1 = normal_direction(*faces[around[0].face].offset, around[0].uv[0], around[0].uv[1]);
                        const auto n2 = normal_direction(*faces[around[1].face].offset, around[1].uv[0], around[1].uv[1]);
                        across = cross(n1, n2);
                        if (!(norm(across) > T(1e-9)))
                            refuse<T>("a vertex of two sharp faces tangent there", {vid});
                        across = across / norm(across);
                    }
                    // first guess: the planes of the offsets, least squares, the smallest step (the
                    // point of the line of two planes nearest to the vertex)
                    Mat<T> A(Eigen::Index(k), 3);
                    Vec<T> b = Vec<T>::Zero(Eigen::Index(k));
                    for (std::size_t j = 0; j < k; ++j)
                    {
                        const auto &f = faces[around[j].face];
                        const auto nrm = normal_direction(*f.offset, around[j].uv[0], around[j].uv[1]);
                        const auto q = f.offset->value(around[j].uv[0], around[j].uv[1]);
                        for (int a = 0; a < 3; ++a)
                            A(Eigen::Index(j), a) = nrm[a];
                        b[Eigen::Index(j)] = nrm * (q - p0);
                    }
                    const Vec<T> dx = A.completeOrthogonalDecomposition().solve(b);
                    const point<T, 3> x0{p0[0] + dx[0], p0[1] + dx[1], p0[2] + dx[2]};
                    Vec<T> x(Eigen::Index(3 + 2 * k));
                    for (int a = 0; a < 3; ++a)
                        x[a] = x0[a];
                    for (std::size_t j = 0; j < k; ++j)
                    {
                        const auto uv = start_uv(*faces[around[j].face].offset, around[j].uv, x0);
                        x[Eigen::Index(3 + 2 * j)] = uv[0];
                        x[Eigen::Index(4 + 2 * j)] = uv[1];
                    }
                    const auto rows = Eigen::Index(3 * k + (k == 2 ? 1 : 0));
                    auto f = [&](const Vec<T> &y, Vec<T> &r, Mat<T> &J) {
                        r.setZero(rows);
                        J.setZero(rows, Eigen::Index(3 + 2 * k));
                        for (std::size_t j = 0; j < k; ++j)
                        {
                            const auto &s = *faces[around[j].face].offset;
                            const T u = y[Eigen::Index(3 + 2 * j)], w = y[Eigen::Index(4 + 2 * j)];
                            const auto g = s.value(u, w), gu = s.value(u, w, 1, 0), gv = s.value(u, w, 0, 1);
                            for (int a = 0; a < 3; ++a)
                            {
                                const auto row = Eigen::Index(3 * j + std::size_t(a));
                                r[row] = g[a] - y[a];
                                J(row, a) = T(-1);
                                J(row, Eigen::Index(3 + 2 * j)) = gu[a];
                                J(row, Eigen::Index(4 + 2 * j)) = gv[a];
                            }
                        }
                        if (k == 2)
                            for (int a = 0; a < 3; ++a)
                            {
                                r[rows - 1] += (y[a] - p0[a]) * across[a];
                                J(rows - 1, a) = across[a];
                            }
                    };
                    const T res = gauss_newton<T>(x, f);
                    if (!(res <= tol / T(10)))
                        refuse<T>("the offsets of the faces around a vertex do not meet at one point", {vid});
                    nv.p = {x[0], x[1], x[2]};
                    nv.sharp = true;
                    // a face with a pole there: the pole of its offset inside every other face (behind
                    // their offsets) would be inside the new face, which then needs another topology
                    // (a sphere cut through its poles, offset outward)
                    if (const auto pf = pole_faces.find(v); pf != pole_faces.end())
                        for (std::size_t j = 0; j < k; ++j)
                        {
                            if (!pf->second.contains(around[j].face))
                                continue;
                            const auto pole = faces[around[j].face].offset->value(around[j].uv[0], around[j].uv[1]);
                            bool inside = true;
                            for (std::size_t i = 0; i < k && inside; ++i)
                            {
                                if (i == j)
                                    continue;
                                const auto &fi = faces[around[i].face];
                                auto n = normal_direction(*fi.offset, x[Eigen::Index(3 + 2 * i)], x[Eigen::Index(4 + 2 * i)]);
                                if (fi.orient != Orientation::Forward)
                                    n = point<T, 3>{-n[0], -n[1], -n[2]};
                                inside = (pole - nv.p) * n < -tol;
                            }
                            if (inside)
                                refuse<T>("the offset of a face would cover its pole (a sphere cut through its poles, offset outward)", {vid});
                        }
                }
                moved[v] = std::move(nv);
            }

            // ---- the new entities
            std::map<std::uint32_t, VertexId> new_vertex;
            for (const auto &[v, nv] : moved)
                new_vertex[v] = m.add(Vertex<T>{nv.p, std::max(tol, m.vertex(VertexId{v}).tol)});
            std::map<CoEdgeRef, std::shared_ptr<Curve<T, 2>>> new_pc;
            std::map<std::uint32_t, std::pair<EdgeId, bool>> new_edge; // and reversed
            IntersectOptions<T> io = o.intersect;
            io.tol = std::max(io.tol, tol / T(10));
            for (const auto &ei : edges)
            {
                const auto &e = m.edge(ei.edge);
                const auto &fa = faces[face_index.at(ei.uses[0].face.index)], &fb = faces[face_index.at(ei.uses[1].face.index)];
                const auto &a = moved.at(e.v1.index), &b = moved.at(e.v2.index);
                if (ei.kind == EdgeKind::Tangent)
                {
                    const T t1 = a.t.contains({ei.edge.index, 0}) ? a.t.at({ei.edge.index, 0}) : e.u1;
                    const T t2 = b.t.contains({ei.edge.index, 1}) ? b.t.at({ei.edge.index, 1}) : e.u2;
                    if (!(t2 > t1))
                        refuse<T>("a face vanishes (an edge turned back)", {ei.edge});
                    const auto pa = tangent_pc.at(ei.uses[0].coedge), pb = tangent_pc.at(ei.uses[1].coedge);
                    for (const auto *pc : {pa.get(), pb.get()})
                    {
                        const auto pcb = pc->bounds();
                        if (t1 < pcb[0] - T(1e-12) * (T(1) + std::abs(pcb[0])) || t2 > pcb[1] + T(1e-12) * (T(1) + std::abs(pcb[1])))
                            refuse<T>("a tangent edge that cannot be continued", {ei.edge});
                    }
                    // the offsets meet along it
                    for (int i = 0; i <= 8; ++i)
                    {
                        const T t = t1 + (t2 - t1) * T(i) / T(8);
                        const auto qa = pa->value(t), qb = pb->value(t);
                        if (norm(fa.offset->value(qa[0], qa[1]) - fb.offset->value(qb[0], qb[1])) > tol)
                            refuse<T>("the offsets part along a tangent edge", {ei.edge});
                    }
                    Edge<T> ne;
                    ne.curve = std::make_shared<CurveOnSurface<T, 3>>(pa, fa.offset);
                    ne.u1 = t1;
                    ne.u2 = t2;
                    ne.v1 = new_vertex.at(e.v1.index);
                    ne.v2 = new_vertex.at(e.v2.index);
                    ne.tol = std::max(tol, e.tol);
                    new_edge[ei.edge.index] = {m.add(std::move(ne)), false};
                    new_pc[ei.uses[0].coedge] = pa;
                    new_pc[ei.uses[1].coedge] = pb;
                    continue;
                }
                // sharp: the branch of the intersection of the offsets through the new vertices
                if (fa.face == fb.face)
                    refuse<T>("a sharp edge of a face with itself", {ei.edge});
                const auto res = intersect(*fa.offset, fa.box, *fb.offset, fb.box, io);
                const bool closed = e.v1 == e.v2;
                const auto pm = e.curve->value((e.u1 + e.u2) / T(2));
                const IntersectionCurve<T> *best = nullptr;
                T best_gap = std::numeric_limits<T>::infinity();
                T start{}, end{};      // the new edge on the branch (end past b on a loop run across its ends)
                bool reversed = false; // the branch runs from the new v2 to the new v1
                const T reach = T(10) * std::max(tol, io.tol);
                for (const auto &c : res.curves)
                {
                    const T span = c.b - c.a;
                    const auto [s1, d1] = extrema_curve_point(*c.curve, a.p, T(1e-12) * (T(1) + std::abs(span)));
                    const auto [s2, d2] = extrema_curve_point(*c.curve, b.p, T(1e-12) * (T(1) + std::abs(span)));
                    if (d1 > reach || d2 > reach)
                        continue;
                    // a loop: a closed branch, or one whose ends meet (the box of a closed surface, a
                    // NURBS cylinder, cuts the loop at its seam)
                    const bool loop = c.closed || norm(c.curve->value(c.a) - c.curve->value(c.b)) <= reach;
                    if (closed && !loop)
                        continue;
                    // on a loop, a parameter in [a, b) (a vertex at its ends is at a), a point anywhere
                    const T eps = T(1e-12) * (T(1) + std::abs(c.a) + std::abs(c.b));
                    auto wrap = [&](T t) { return loop && t >= c.b - eps ? t - span : t; };
                    auto at = [&](T t) {
                        if (loop)
                        {
                            t = c.a + std::fmod(t - c.a, span);
                            if (t < c.a)
                                t += span;
                        }
                        return c.curve->value(t);
                    };
                    // the arcs between the vertices: the whole loop for a closed edge, on a loop either
                    // way round it; the one whose middle is nearest that of the old edge
                    std::vector<std::tuple<T, T, bool>> arcs;
                    const T w1 = wrap(s1), w2 = wrap(s2);
                    if (closed)
                        arcs.emplace_back(w1, w1 + span, false);
                    else if (!loop)
                        arcs.emplace_back(std::min(s1, s2), std::max(s1, s2), s1 > s2);
                    else
                    {
                        arcs.emplace_back(w1, w2 > w1 ? w2 : w2 + span, false);
                        arcs.emplace_back(w2, w1 > w2 ? w1 : w1 + span, true);
                    }
                    for (const auto &[t1, t2, rev] : arcs)
                    {
                        const T gap = norm(at((t1 + t2) / T(2)) - pm);
                        if (gap < best_gap)
                            best_gap = gap, best = &c, start = t1, end = t2, reversed = rev;
                    }
                }
                if (!best)
                    refuse<T>("no branch of the intersection of the offsets joins the vertices of an edge", {ei.edge});
                Edge<T> ne;
                ne.tol = std::max({tol, e.tol, best->deviation});
                std::shared_ptr<Curve<T, 2>> pca = best->pc1, pcb = best->pc2;
                ne.curve = best->curve;
                ne.u1 = start;
                ne.u2 = end;
                if (end > best->b + T(1e-12) * (T(1) + std::abs(best->b)))
                {
                    // across the ends of a loop: the loop run from the start of the edge
                    if (const auto *ec = dynamic_cast<const ElementaryCurve<T> *>(best->curve.get()))
                    {
                        ne.curve = intersect_detail::restricted(*ec, start, end);
                        const auto *ea = dynamic_cast<const ElementarySurface<T> *>(fa.offset.get());
                        const auto *eb = dynamic_cast<const ElementarySurface<T> *>(fb.offset.get());
                        if (!ea || !eb)
                            refuse<T>("a sharp edge on a loop of a surface that is not a primitive", {ei.edge});
                        const auto ka = intersect_detail::primitive_pcurve(*ne.curve, start, end, *ea, fa.box, io.tol);
                        const auto kb = intersect_detail::primitive_pcurve(*ne.curve, start, end, *eb, fb.box, io.tol);
                        pca = ka.pcurve;
                        pcb = kb.pcurve;
                        ne.tol = std::max({ne.tol, ka.deviation, kb.deviation});
                    }
                    else
                    {
                        ne.curve = from_parameter(best->curve, best->a, best->b, start);
                        pca = from_parameter(best->pc1, best->a, best->b, start);
                        pcb = from_parameter(best->pc2, best->a, best->b, start);
                    }
                }
                if (closed)
                {
                    // the loop runs the way of the edge?
                    const T tm = (ne.u1 + ne.u2) / T(2);
                    const auto dn = ne.curve->value(tm, 1);
                    const auto [s, dist] = extrema_curve_point(*e.curve, ne.curve->value(tm), T(1e-12));
                    reversed = dn * e.curve->value(s, 1) < T(0);
                    ne.v1 = ne.v2 = new_vertex.at(e.v1.index);
                }
                else
                {
                    if (!(ne.u2 - ne.u1 > T(0)))
                        refuse<T>("a face vanishes (an edge reduced to a point)", {ei.edge});
                    // the new edge runs the way of the old one, or the face turned back
                    const auto dn = ne.curve->value((ne.u1 + ne.u2) / T(2), 1);
                    const T dir = (reversed ? T(-1) : T(1)) * (dn * e.curve->value((e.u1 + e.u2) / T(2), 1));
                    if (!(dir > T(0)))
                        refuse<T>("a face vanishes (an edge turned back)", {ei.edge});
                    ne.v1 = new_vertex.at((reversed ? e.v2 : e.v1).index);
                    ne.v2 = new_vertex.at((reversed ? e.v1 : e.v2).index);
                }
                new_edge[ei.edge.index] = {m.add(std::move(ne)), reversed};
                new_pc[ei.uses[0].coedge] = pca;
                new_pc[ei.uses[1].coedge] = pcb;
            }
            std::set<CoEdgeRef> dropped; // pole edges left out
            for (const auto &[eid, use] : degenerate)
            {
                const auto &[k, ref] = use;
                if (new_edge.contains(eid.index))
                {
                    new_pc[ref] = shifted(faces[k], pcurve_of(ref));
                    continue;
                }
                const auto &e = m.edge(eid);
                const auto pc = shifted(faces[k], pcurve_of(ref));
                const auto &nv = moved.at(e.v1.index);
                bool point = true;
                for (int i = 0; i <= 8 && point; ++i)
                {
                    const auto q = pc->value(e.u1 + (e.u2 - e.u1) * T(i) / T(8));
                    point = norm(faces[k].offset->value(q[0], q[1]) - nv.p) <= tol;
                }
                if (!point)
                {
                    // a sharp vertex leaves the pole of the offset (a sphere cut through its poles):
                    // the face no longer reaches it, its pole edge is left out of its wire
                    if (!nv.sharp)
                        refuse<T>("a pole whose offset is not a point (the apex of a cone)", {eid});
                    dropped.insert(ref);
                    continue;
                }
                Edge<T> ne;
                ne.u1 = e.u1;
                ne.u2 = e.u2;
                ne.v1 = ne.v2 = new_vertex.at(e.v1.index);
                ne.tol = std::max(tol, e.tol);
                ne.degenerate = true;
                new_edge[eid.index] = {m.add(std::move(ne)), false};
                new_pc[ref] = pc;
            }

            std::vector<ShellId> new_shells;
            for (auto sh : shells)
            {
                Shell ns;
                ns.closed = true;
                for (const auto &fu : m.shell(sh).faces)
                {
                    const auto &f = faces[face_index.at(fu.face.index)];
                    const auto &face = m.face(fu.face);
                    Face<T> nf{f.offset, {}, std::max(tol, face.tol), face.natural_bounds && !f.extended && f.dv == T(0)};
                    for (auto w : face.wires)
                    {
                        Wire<T> nw;
                        nw.closed = true;
                        const auto &wire = m.wire(w);
                        bool left_out = false;
                        for (std::uint32_t i = 0; i < wire.coedges.size(); ++i)
                        {
                            if (dropped.contains(CoEdgeRef{w, i}))
                            {
                                left_out = true;
                                continue;
                            }
                            const auto &ce = wire.coedges[i];
                            const auto &[ne, rev] = new_edge.at(ce.edge.index);
                            nw.coedges.push_back(CoEdge<T>{ne, rev ? reverse(ce.orient) : ce.orient, new_pc.at(CoEdgeRef{w, i})});
                        }
                        if (left_out)
                        {
                            // the neighbours of a pole edge left out meet in (u,v): the new vertex is a
                            // regular point of the offset
                            const auto [u1, u2, v1, v2] = f.box;
                            const T gap = T(1e-7) * std::max({T(1), u2 - u1, v2 - v1});
                            for (std::size_t i = 0; i < nw.coedges.size(); ++i)
                            {
                                const auto &c0 = nw.coedges[i], &c1 = nw.coedges[(i + 1) % nw.coedges.size()];
                                if (norm(coedge_uv(m, c0, m.edge(c0.edge).u2) - coedge_uv(m, c1, m.edge(c1.edge).u1)) > gap)
                                    refuse<T>("a pole edge left out where the pcurves do not meet", {fu.face});
                            }
                        }
                        nf.wires.push_back(m.add(std::move(nw)));
                    }
                    const auto nid = m.add(std::move(nf));
                    if (f.d != T(0))
                    {
                        // the offset regular on the new face (its continuation was not checked)
                        const auto nb = classify_detail::face_uv_box(m, nid);
                        if (!offset_surface(f.surface, f.d, tol / T(10), std::array<T, 4>{nb[0], nb[1], nb[2] - f.dv, nb[3] - f.dv}))
                            refuse<T>("the offset folds (beyond a radius of curvature)", {fu.face});
                    }
                    ns.faces.push_back(FaceUse{nid, fu.orient});
                }
                new_shells.push_back(m.add(std::move(ns)));
            }
            Solid nsol;
            nsol.outer = new_shells.front();
            nsol.voids.assign(new_shells.begin() + 1, new_shells.end());
            const auto out = m.add(std::move(nsol));
            const auto report = check(m, ShapeId{out});
            if (!report.ok())
            {
                const auto &first = report.entries.front();
                refuse<T>(std::string("the offset solid is not valid: ") + to_string(first.issue) + " " + first.detail, {first.shape});
            }
            // the offset intersecting itself: two faces sharing no vertex that meet (the
            // section of each pair, its debris erased)
            const auto nf = explore<FaceId>(m, ShapeId{out});
            std::vector<std::set<std::uint32_t>> fv(nf.size());
            for (std::size_t i = 0; i < nf.size(); ++i)
                for (auto v : explore<VertexId>(m, ShapeId{nf[i]}))
                    fv[i].insert(v.index);
            SectionOptions<T> sopt;
            sopt.tol = tol;
            sopt.intersect = io;
            // the boxes of the faces once (`section` filters its pairs the same way)
            std::vector<BoundingBox<T>> boxes;
            for (auto f : nf)
                boxes.push_back(section_detail::face_box(m, f, tol));
            const auto built = detail::arena_size(m);
            for (std::size_t i = 0; i < nf.size(); ++i)
                for (std::size_t j = i + 1; j < nf.size(); ++j)
                {
                    if (!section_detail::boxes_overlap(boxes[i], boxes[j]) || std::ranges::any_of(fv[i], [&](std::uint32_t v) { return fv[j].contains(v); }))
                        continue;
                    const auto r = section(m, ShapeId{nf[i]}, ShapeId{nf[j]}, sopt);
                    const bool meet = r && (!r->edges.empty() || !r->points.empty() || !r->coincident.empty() || !r->boundary_edges.empty());
                    detail::rollback(m, built);
                    if (meet)
                        refuse<T>("the offset intersects itself", {nf[i], nf[j]});
                }
            return out;
        }

        template <std::floating_point T, typename F>
        auto guarded(Model<T> &m, const char *what, F &&f) -> BuildResult<SolidId>
        {
            const auto size = detail::arena_size(m);
            try
            {
                return f();
            }
            catch (const BuildFailure &e)
            {
                detail::rollback(m, size);
                auto err = detail::rolled_back(size, e.error());
                if (err.message.starts_with("offset_solid: "))
                    err.message = std::string(what) + ": " + err.message.substr(14);
                return std::unexpected(std::move(err));
            }
            catch (const std::exception &e)
            {
                detail::rollback(m, size);
                return detail::fail(BuildErrc::NonGeneric, std::string(what) + ": " + e.what());
            }
        }
    }

    /**
     * @brief The solid `s` offset by `d` along its outward normal (d > 0 grows it, d < 0 shrinks
     * it), with sharp joints: each face on its offset surface, each sharp edge the intersection
     * of the offsets of its faces, continued past the faces as needed (`thick.h`). The solid is
     * not modified; on failure the model is left as it was.
     */
    template <std::floating_point T>
    [[nodiscard]] auto offset_solid(Model<T> &m, SolidId s, std::type_identity_t<T> d, const OffsetOptions<T> &o = {}) -> BuildResult<SolidId>
    {
        if (!m.alive(s))
            return detail::fail(BuildErrc::InvalidId, "offset_solid: dead or invalid solid", {s});
        if (!detail::valid_tolerance(o.tol) || !std::isfinite(d))
            return detail::fail(BuildErrc::InvalidBounds, "offset_solid: a tolerance not positive or a distance not finite", {s});
        if (std::abs(d) <= o.tol)
            return detail::fail(BuildErrc::InvalidBounds, "offset_solid: a distance within the tolerance", {s});
        return thick_detail::guarded(m, "offset_solid", [&] { return BuildResult<SolidId>{thick_detail::offset_faces(m, s, [d](FaceId) { return d; }, o)}; });
    }

    /**
     * @brief The solid `s` hollowed out to the thickness `t` and opened by the faces `removed`:
     * `s` minus its offset by -t, the removed faces offset by 0 (`cut`; they lie on the faces of
     * `s`, and the boolean of coincident faces opens the shell there). A negative t grows outward:
     * the offset by |t| minus `s`. Without a removed face, a closed hollow solid. On failure the
     * model is left as it was.
     */
    template <std::floating_point T>
    [[nodiscard]] auto make_thick(Model<T> &m, SolidId s, std::span<const FaceId> removed, std::type_identity_t<T> t, const OffsetOptions<T> &o = {})
        -> BuildResult<SolidId>
    {
        if (!m.alive(s))
            return detail::fail(BuildErrc::InvalidId, "make_thick: dead or invalid solid", {s});
        if (!detail::valid_tolerance(o.tol) || !std::isfinite(t))
            return detail::fail(BuildErrc::InvalidBounds, "make_thick: a tolerance not positive or a thickness not finite", {s});
        if (std::abs(t) <= o.tol)
            return detail::fail(BuildErrc::InvalidBounds, "make_thick: a thickness within the tolerance", {s});
        const auto faces = explore<FaceId>(m, ShapeId{s});
        for (auto f : removed)
            if (!std::ranges::contains(faces, f))
                return detail::fail(BuildErrc::InvalidId, "make_thick: a removed face is not a face of the solid", {f});
        std::set<std::uint32_t> gone;
        for (auto f : removed)
            gone.insert(f.index);
        return thick_detail::guarded(m, "make_thick", [&]() -> BuildResult<SolidId> {
            const T d = -t;
            const auto before = detail::arena_size(m);
            const auto other = thick_detail::offset_faces(m, s, [&](FaceId f) { return gone.contains(f.index) ? T(0) : d; }, o);
            BooleanOptions<T> bo;
            bo.tol = o.tol;
            auto r = t > T(0) ? cut(m, ShapeId{s}, ShapeId{other}, bo) : cut(m, ShapeId{other}, ShapeId{s}, bo);
            if (!r)
            {
                auto err = r.error();
                err.message = "make_thick: " + err.message;
                throw BuildFailure(std::move(err));
            }
            if (r->solids.size() != 1)
                thick_detail::refuse<T>("the thick solid is not one solid", {r->shape});
            // the offset was a step: erased
            boolean_detail::erase_shape(m, ShapeId{other});
            boolean_detail::erase_unused(m, before[1], before[0]);
            return r->solids.front();
        });
    }
}
