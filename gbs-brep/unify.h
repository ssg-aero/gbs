#pragma once
/**
 * @file unify.h
 * @brief Faces of one domain unified (stage 4, `docs/sources/design/brep_palier4.md`, § 0 ter):
 * adjacent faces of a set lying on one surface (the same object, or equal primitives) become one
 * face, whatever the seam of the parametrization; then the pieces of an edge split at a vertex
 * that no longer bounds anything are joined back.
 *
 * - Two faces are of one domain when the parameters of one map onto those of the other by an
 *   affine map (the identity for the same surface object; a rigid map of the plane between two
 *   coplanar planes, a shift of the angle between two coaxial cylinders), checked on samples, and
 *   when their outward normals (orientation in the shell included) agree.
 * - The pcurves of the second face are mapped into the parameters of the first (exact for the
 *   B-splines: an affine map commutes with the basis), shifted by whole periods so that the edges
 *   they share lie at the same place; the shared edges are dropped and the remaining co-edges
 *   chained into loops by their vertices and their parameters.
 * - A periodic parameter gets a window (box) centred on the merged face, so that a face across the
 *   seam of the surface needs no seam: the face lives on a copy of the primitive with that box
 *   (`ElementarySurface::with_box`). A region that goes round the period keeps its seam.
 * - A plane in B-spline form is put on the Plane primitive to merge with a coplanar face
 *   (`to_plane`); a face of a sphere around a pole is put on the sphere turned away from it
 *   (`reframe_pole`); the pieces of edges and of seams left are joined (`join_edges`, `join_seams`).
 */
#include <gbs-brep/builders.h>
#include <gbs-brep/closure.h>
#include <gbs-brep/explore.h>
#include <gbs-brep/model.h>
#include <gbs/bscurve.h>
#include <gbs/primitives.h>
#include <gbs/transform.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numbers>
#include <memory>
#include <optional>
#include <set>
#include <typeinfo>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs::brep
{
    /// What `unify_faces` did.
    struct UnifyReport
    {
        std::size_t faces_merged{}; ///< pairs of faces merged into one
        std::size_t edges_joined{}; ///< pairs of edge pieces joined into one edge
        std::size_t faces_planed{}; ///< faces on a plane in B-spline form put on the Plane primitive (`to_plane`)
        std::size_t poles_reframed{}; ///< faces of a sphere holding a pole put on a turned sphere (`reframe_pole`)
        /// {new face, first piece, second piece}, in order; the second piece is invalid for a face
        /// only put on another surface (`to_plane`, `reframe_pole`)
        std::vector<std::array<FaceId, 3>> history;
    };

    namespace unify_detail
    {
        /// uv -> A uv + b
        template <std::floating_point T>
        struct Affine2
        {
            std::array<T, 4> a{T(1), T(0), T(0), T(1)}; // row major
            point<T, 2> b{};

            point<T, 2> operator()(const point<T, 2> &p) const { return {a[0] * p[0] + a[1] * p[1] + b[0], a[2] * p[0] + a[3] * p[1] + b[1]}; }
            [[nodiscard]] T det() const { return a[0] * a[3] - a[1] * a[2]; }
            [[nodiscard]] bool identity() const
            {
                return a == std::array<T, 4>{T(1), T(0), T(0), T(1)} && b == point<T, 2>{};
            }
        };

        /// The pcurve mapped by A (exact for B-splines); nullptr for a curve that cannot be.
        template <std::floating_point T>
        std::shared_ptr<Curve<T, 2>> mapped(const std::shared_ptr<Curve<T, 2>> &c, const Affine2<T> &A)
        {
            if (!c || A.identity())
                return c;
            auto f = [&](const point<T, 2> &p) { return A(p); };
            if (auto r = std::dynamic_pointer_cast<BSCurveRational<T, 2>>(c))
            {
                auto n = std::make_shared<BSCurveRational<T, 2>>(*r);
                gbs::transform(*n, f);
                return n;
            }
            if (auto b = std::dynamic_pointer_cast<BSCurve<T, 2>>(c))
            {
                auto n = std::make_shared<BSCurve<T, 2>>(*b);
                gbs::transform(*n, f);
                return n;
            }
            return nullptr;
        }

        /// A co-edge of a merged face, with its ends in the parameters of the merged face.
        template <std::floating_point T>
        struct Use
        {
            CoEdge<T> ce;
            VertexId vs, ve;     // start and end vertex in the sense of the co-edge
            point<T, 2> s, e;    // start and end parameters
        };

        template <std::floating_point T>
        Use<T> use_of(const Model<T> &m, const CoEdge<T> &ce)
        {
            const auto &ed = m.edge(ce.edge);
            const bool fwd = ce.orient == Orientation::Forward;
            const auto p1 = ce.pcurve->value(ed.u1), p2 = ce.pcurve->value(ed.u2);
            return {ce, fwd ? ed.v1 : ed.v2, fwd ? ed.v2 : ed.v1, fwd ? p1 : p2, fwd ? p2 : p1};
        }

        template <std::floating_point T>
        T dist(const point<T, 2> &p, const point<T, 2> &q)
        {
            return std::hypot(p[0] - q[0], p[1] - q[1]);
        }

        /// Points of the boundary of a face in its parameters (ends and middles of its co-edges).
        template <std::floating_point T>
        std::vector<point<T, 2>> boundary_samples(const Model<T> &m, const Face<T> &f, int n = 4)
        {
            std::vector<point<T, 2>> r;
            for (auto w : f.wires)
                for (const auto &ce : m.wire(w).coedges)
                {
                    if (!ce.pcurve)
                        continue;
                    const auto &ed = m.edge(ce.edge);
                    for (int k = 0; k < n; ++k)
                        r.push_back(ce.pcurve->value(ed.u1 + (ed.u2 - ed.u1) * T(k) / T(n)));
                }
            return r;
        }

        /// The periods of the parameters of a surface on its box (0: not periodic).
        template <std::floating_point T>
        std::array<T, 2> periods(const Surface<T, 3> &s, T tol)
        {
            const auto c = surface_closure(s, tol);
            const auto b = s.bounds();
            return {c.closed_u ? b[1] - b[0] : T(0), c.closed_v ? b[3] - b[2] : T(0)};
        }

        /// The map of the parameters of fb onto those of fa when both lie on one surface.
        template <std::floating_point T>
        std::optional<Affine2<T>> parameter_map(const Model<T> &m, const Face<T> &fa, const Face<T> &fb, T tol)
        {
            if (fa.surface == fb.surface)
                return Affine2<T>{};
            const auto *ea = dynamic_cast<const ElementarySurface<T> *>(fa.surface.get());
            const auto *eb = dynamic_cast<const ElementarySurface<T> *>(fb.surface.get());
            if (!ea || !eb || typeid(*ea) != typeid(*eb))
                return std::nullopt;
            const auto pa = periods(*ea, tol);
            const auto uvb = boundary_samples(m, fb);
            if (uvb.size() < 3)
                return std::nullopt;
            // a first map from the derivatives at the first sample: dSa(q0) L = dSb(uv0)
            const auto uv0 = uvb.front();
            const auto q0 = ea->parameters(fb.surface->value(uv0[0], uv0[1]));
            Affine2<T> guess;
            {
                const auto au = ea->value(q0[0], q0[1], 1, 0), av = ea->value(q0[0], q0[1], 0, 1);
                const auto bu = fb.surface->value(uv0[0], uv0[1], 1, 0), bv = fb.surface->value(uv0[0], uv0[1], 0, 1);
                const T g11 = au * au, g12 = au * av, g22 = av * av, g = g11 * g22 - g12 * g12;
                if (!(std::abs(g) > T(0)))
                    return std::nullopt;
                auto col = [&](const point<T, 3> &d) -> std::array<T, 2> {
                    const T x = au * d, y = av * d;
                    return {(g22 * x - g12 * y) / g, (g11 * y - g12 * x) / g};
                };
                const auto cu = col(bu), cv = col(bv);
                guess.a = {cu[0], cv[0], cu[1], cv[1]};
                guess.b = q0 - point<T, 2>{guess.a[0] * uv0[0] + guess.a[1] * uv0[1], guess.a[2] * uv0[0] + guess.a[3] * uv0[1]};
            }
            std::vector<point<T, 2>> uva;
            for (const auto &uv : uvb)
            {
                const auto p = fb.surface->value(uv[0], uv[1]);
                auto q = ea->parameters(p);
                if (norm(ea->value(q[0], q[1]) - p) > T(10) * tol)
                    return std::nullopt; // (not on the surface of fa)
                // periodic parameters unwrapped next to the first map (a face going round the
                // period has samples at both ends of it)
                const auto g = guess(uv);
                for (int i = 0; i < 2; ++i)
                    if (pa[i] > T(0))
                        q[i] += pa[i] * std::round((g[i] - q[i]) / pa[i]);
                uva.push_back(q);
            }
            // least squares: uva = A uvb + b, one coordinate at a time (normal equations 3 x 3)
            Affine2<T> A;
            std::array<std::array<T, 3>, 3> N{};
            std::array<std::array<T, 3>, 2> rhs{};
            for (std::size_t k = 0; k < uvb.size(); ++k)
            {
                const std::array<T, 3> x{uvb[k][0], uvb[k][1], T(1)};
                for (int i = 0; i < 3; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                        N[i][j] += x[i] * x[j];
                    rhs[0][i] += x[i] * uva[k][0];
                    rhs[1][i] += x[i] * uva[k][1];
                }
            }
            auto solve = [](std::array<std::array<T, 3>, 3> M, std::array<T, 3> y) -> std::optional<std::array<T, 3>> {
                for (int c = 0; c < 3; ++c)
                {
                    int piv = c;
                    for (int r = c + 1; r < 3; ++r)
                        if (std::abs(M[r][c]) > std::abs(M[piv][c]))
                            piv = r;
                    if (!(std::abs(M[piv][c]) > T(1e-300)))
                        return std::nullopt;
                    std::swap(M[c], M[piv]);
                    std::swap(y[c], y[piv]);
                    for (int r = 0; r < 3; ++r)
                        if (r != c)
                        {
                            const T f = M[r][c] / M[c][c];
                            for (int k = 0; k < 3; ++k)
                                M[r][k] -= f * M[c][k];
                            y[r] -= f * y[c];
                        }
                }
                return std::array<T, 3>{y[0] / M[0][0], y[1] / M[1][1], y[2] / M[2][2]};
            };
            const auto r0 = solve(N, rhs[0]), r1 = solve(N, rhs[1]);
            if (!r0 || !r1)
                return std::nullopt;
            A.a = {(*r0)[0], (*r0)[1], (*r1)[0], (*r1)[1]};
            A.b = {(*r0)[2], (*r1)[2]};
            if (!(std::abs(A.det()) > T(1e-12)))
                return std::nullopt;
            // checked on the samples and between them, within tol: the edges of b keep their
            // tolerance on the surface of a (faces a few tolerances apart stay two faces)
            for (std::size_t k = 0; k < uvb.size(); ++k)
            {
                const auto mid = (uvb[k] + uvb[(k + 1) % uvb.size()]) * T(0.5);
                for (const auto &uv : {uvb[k], mid})
                {
                    const auto q = A(uv);
                    if (norm(ea->value(q[0], q[1]) - fb.surface->value(uv[0], uv[1])) > tol)
                        return std::nullopt;
                }
            }
            return A;
        }

        /// Signed area of a loop of uses in the parameter plane (its pcurves sampled).
        template <std::floating_point T>
        T loop_area(const Model<T> &m, const std::vector<Use<T>> &loop)
        {
            std::vector<point<T, 2>> p;
            for (const auto &u : loop)
            {
                const auto &ed = m.edge(u.ce.edge);
                for (int k = 0; k < 8; ++k)
                {
                    const T s = T(k) / T(8);
                    const T t = u.ce.orient == Orientation::Forward ? ed.u1 + (ed.u2 - ed.u1) * s : ed.u2 - (ed.u2 - ed.u1) * s;
                    p.push_back(u.ce.pcurve->value(t));
                }
            }
            T a{};
            for (std::size_t i = 0, j = p.size() - 1; i < p.size(); j = i++)
                a += p[j][0] * p[i][1] - p[i][0] * p[j][1];
            return a / T(2);
        }

        /// The box of the parameters of a face (its boundary sampled).
        template <std::floating_point T>
        std::array<T, 4> face_box(const Model<T> &m, const Face<T> &f)
        {
            std::array<T, 4> box{std::numeric_limits<T>::max(), std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max(),
                                 std::numeric_limits<T>::lowest()};
            for (const auto &p : boundary_samples(m, f, 16))
                box = {std::min(box[0], p[0]), std::max(box[1], p[0]), std::min(box[2], p[1]), std::max(box[3], p[1])};
            return box;
        }

        /// A surface whose map is affine on the box [u1, u2] x [v1, v2] of a face, S(u, v) =
        /// S(u1, v1) + (u - u1) X + (v - v1) Y (a plane in B-spline form, as other systems write
        /// planes, or the planar part of a B-spline): the plane, and the map of the parameters of
        /// the surface onto those of the plane; nothing for any other surface.
        template <std::floating_point T>
        std::optional<std::pair<Frame3<T>, Affine2<T>>> affine_plane(const Surface<T, 3> &s, const std::array<T, 4> &box, T tol)
        {
            if (dynamic_cast<const ElementarySurface<T> *>(&s))
                return std::nullopt;
            const auto [u1, u2, v1, v2] = box;
            if (!(u2 > u1) || !(v2 > v1))
                return std::nullopt;
            const auto o = s.value(u1, v1);
            const auto X = (s.value(u2, v1) - o) / (u2 - u1), Y = (s.value(u1, v2) - o) / (v2 - v1);
            const auto z = cross(X, Y);
            if (!(norm(z) > T(1e-12) * norm(X) * norm(Y)))
                return std::nullopt;
            constexpr int n = 6;
            for (int i = 0; i <= n; ++i)
                for (int j = 0; j <= n; ++j)
                {
                    const T u = u1 + (u2 - u1) * T(i) / T(n), v = v1 + (v2 - v1) * T(j) / T(n);
                    if (norm(s.value(u, v) - (o + (u - u1) * X + (v - v1) * Y)) > tol)
                        return std::nullopt;
                }
            const Frame3<T> f{ax2<T, 3>{o, z, X}};
            // plane parameters: ((p - o) . x, (p - o) . y), p = o + (u - u1) X + (v - v1) Y
            Affine2<T> A;
            A.a = {X * f.x, Y * f.x, X * f.y, Y * f.y};
            A.b = {-(A.a[0] * u1 + A.a[1] * v1), -(A.a[2] * u1 + A.a[3] * v1)};
            return std::pair{f, A};
        }

        /// Whether two frames carry one plane (normals parallel, either way, origins on it).
        template <std::floating_point T>
        bool coplanar(const Frame3<T> &f, const Frame3<T> &g, T tol)
        {
            return norm(cross(f.z, g.z)) <= T(1e-9) && std::abs((g.o - f.o) * f.z) <= tol;
        }

        /// The plane of a face: its own, or that of its surface in B-spline form.
        template <std::floating_point T>
        std::optional<Frame3<T>> plane_of(const Model<T> &m, const Face<T> &f, T tol)
        {
            if (const auto *p = dynamic_cast<const Plane<T> *>(f.surface.get()))
                return p->frame();
            if (auto r = affine_plane(*f.surface, face_box(m, f), tol))
                return r->first;
            return std::nullopt;
        }
    }

    /**
     * @brief Puts a face on a plane in B-spline form (on the box of the face) onto the Plane primitive it is (its pcurves
     * mapped, exactly): the face can then merge with faces of that plane beyond the box of its
     * surface. Returns the new face, which replaces the old one in the shells; nothing when the
     * surface is not such a plane.
     */
    template <std::floating_point T>
    auto to_plane(Model<T> &m, const FaceUse &fu, T tol) -> std::optional<FaceUse>
    {
        using namespace unify_detail;
        const Face<T> f = m.face(fu.face);
        if (!f.surface)
            return std::nullopt;
        const auto ap = affine_plane(*f.surface, face_box(m, f), tol);
        if (!ap)
            return std::nullopt;
        const auto &[frame, A] = *ap;
        std::vector<Wire<T>> wires;
        std::array<T, 4> box{std::numeric_limits<T>::max(), std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max(),
                             std::numeric_limits<T>::lowest()};
        for (auto w : f.wires)
        {
            Wire<T> nw = m.wire(w);
            for (auto &ce : nw.coedges)
            {
                ce.pcurve = mapped(ce.pcurve, A);
                if (!ce.pcurve)
                    return std::nullopt;
                const auto &ed = m.edge(ce.edge);
                for (int k = 0; k <= 16; ++k)
                {
                    const auto p = ce.pcurve->value(ed.u1 + (ed.u2 - ed.u1) * T(k) / T(16));
                    box = {std::min(box[0], p[0]), std::max(box[1], p[0]), std::min(box[2], p[1]), std::max(box[3], p[1])};
                }
            }
            wires.push_back(std::move(nw));
        }
        // the box of the face, with a margin for the bulge of curved pcurves between the samples
        const T margin = T(0.05) * std::max(box[1] - box[0], box[3] - box[2]);
        if (!(margin > T(0)))
            return std::nullopt;
        box = {box[0] - margin, box[1] + margin, box[2] - margin, box[3] + margin};
        Face<T> nf{std::make_shared<Plane<T>>(frame, box), {}, f.tol, false};
        for (auto &w : wires)
            nf.wires.push_back(m.add(std::move(w)));
        const auto name = std::string(m.name(fu.face));
        const auto id = m.add(std::move(nf));
        if (!name.empty())
            m.setName(id, name);
        for (auto sid : m.template ids<ShellId>())
            for (auto &u : m.shell(sid).faces)
                if (u.face == fu.face)
                    u.face = id;
        for (auto w : f.wires)
            m.erase(w);
        m.erase(fu.face);
        return FaceUse{id, fu.orient};
    }

    /**
     * @brief Merges the faces a and b (adjacent, of one domain, facing the same way) into one face,
     * which replaces a in the shells (b is removed from them); returns its use, or nothing when the
     * faces are not of one domain, share no edge at the same place, or would go round a period.
     * The faces a and b are erased; the edges they shared stay in the model, unused by them.
     */
    template <std::floating_point T>
    auto merge_faces(Model<T> &m, const FaceUse &a, const FaceUse &b, T tol) -> std::optional<FaceUse>
    {
        using namespace unify_detail;
        if (a.face == b.face)
            return std::nullopt;
        const Face<T> fa = m.face(a.face), fb = m.face(b.face);
        if (!fa.surface || !fb.surface)
            return std::nullopt;
        auto A = parameter_map(m, fa, fb, tol);
        if (!A)
            return std::nullopt;
        // outward normals agree (the orientation in the shell included)
        {
            const auto samples = boundary_samples(m, fb, 1);
            if (samples.empty())
                return std::nullopt;
            const auto uv = samples.front();
            const auto q = (*A)(uv);
            const auto nb = cross(fb.surface->value(uv[0], uv[1], 1, 0), fb.surface->value(uv[0], uv[1], 0, 1));
            const auto na = cross(fa.surface->value(q[0], q[1], 1, 0), fa.surface->value(q[0], q[1], 0, 1));
            const T s = na * nb;
            if (!(std::abs(s) > T(0)) || ((s > T(0)) != (a.orient == b.orient)))
                return std::nullopt;
        }
        const auto per = periods(*fa.surface, tol);
        const auto box = fa.surface->bounds();
        const T eps_uv = T(1e-6) * std::max({box[1] - box[0], box[3] - box[2], T(1)});

        std::vector<Use<T>> ua, ub;
        for (auto w : fa.wires)
            for (const auto &ce : m.wire(w).coedges)
            {
                if (!ce.pcurve)
                    return std::nullopt;
                ua.push_back(use_of(m, ce));
            }
        std::vector<CoEdge<T>> cb;
        for (auto w : fb.wires)
            for (const auto &ce : m.wire(w).coedges)
            {
                if (!ce.pcurve)
                    return std::nullopt;
                cb.push_back(ce);
            }
        // a shared edge aligns the periods: b shifted by whole periods so that it lies at its place in a
        std::map<std::uint32_t, int> in_a, in_b;
        for (const auto &u : ua)
            ++in_a[u.ce.edge.index];
        for (const auto &ce : cb)
            ++in_b[ce.edge.index];
        std::optional<std::uint32_t> anchor;
        for (const auto &[e, n] : in_a)
            if (n == 1 && in_b.contains(e) && in_b[e] == 1)
            {
                anchor = e;
                break;
            }
        if (!anchor)
            return std::nullopt;
        {
            const auto &ea = *std::ranges::find_if(ua, [&](const auto &u) { return u.ce.edge.index == *anchor; });
            const auto &eb = *std::ranges::find_if(cb, [&](const auto &c) { return c.edge.index == *anchor; });
            const auto &ed = m.edge(ea.ce.edge);
            const T tm = (ed.u1 + ed.u2) / T(2);
            const auto d = ea.ce.pcurve->value(tm) - (*A)(eb.pcurve->value(tm));
            for (int i = 0; i < 2; ++i)
            {
                const T k = per[i] > T(0) ? std::round(d[i] / per[i]) : T(0);
                if (std::abs(d[i] - k * per[i]) > T(100) * eps_uv)
                    return std::nullopt; // (the edge is not at the same place in both faces)
                A->b[i] += k * per[i];
            }
        }
        const bool flip = A->det() < T(0); // (b's parameters run the other way: its co-edges too)
        for (auto ce : cb)
        {
            ce.pcurve = mapped(ce.pcurve, *A);
            if (!ce.pcurve)
                return std::nullopt;
            if (flip)
                ce.orient = reverse(ce.orient);
            ub.push_back(use_of(m, ce));
        }
        // the shared edges at the same place in both faces are dropped
        std::vector<Use<T>> all;
        std::vector<char> dropped_b(ub.size(), 0);
        for (const auto &u : ua)
        {
            bool drop = false;
            for (std::size_t k = 0; k < ub.size() && !drop; ++k)
                if (!dropped_b[k] && ub[k].ce.edge == u.ce.edge && in_a[u.ce.edge.index] == 1 && in_b[u.ce.edge.index] == 1 &&
                    dist(ub[k].s, u.e) <= T(100) * eps_uv && dist(ub[k].e, u.s) <= T(100) * eps_uv)
                {
                    dropped_b[k] = 1;
                    drop = true;
                }
            if (!drop)
                all.push_back(u);
        }
        for (std::size_t k = 0; k < ub.size(); ++k)
            if (!dropped_b[k])
                all.push_back(ub[k]);
        // loops: each co-edge followed by the one starting where it ends (vertex and parameters)
        std::vector<std::vector<Use<T>>> loops;
        std::vector<char> used(all.size(), 0);
        for (std::size_t first = 0; first < all.size(); ++first)
        {
            if (used[first])
                continue;
            std::vector<Use<T>> loop;
            std::size_t cur = first;
            for (std::size_t guard = 0; guard <= all.size(); ++guard)
            {
                used[cur] = 1;
                loop.push_back(all[cur]);
                if (all[cur].ve == all[first].vs && dist(all[cur].e, all[first].s) <= T(100) * eps_uv)
                    break;
                std::optional<std::size_t> next;
                T best = T(100) * eps_uv;
                for (std::size_t k = 0; k < all.size(); ++k)
                    if (!used[k] && all[k].vs == all[cur].ve)
                        if (const T d = dist(all[k].s, all[cur].e); d <= best)
                        {
                            best = d;
                            next = k;
                        }
                if (!next)
                    return std::nullopt; // (the boundary does not chain: the region goes round a period)
                cur = *next;
            }
            if (!(loop.back().ve == loop.front().vs))
                return std::nullopt;
            loops.push_back(std::move(loop));
        }
        // one outer loop (positive area), the holes after it
        std::vector<T> areas;
        for (const auto &l : loops)
            areas.push_back(loop_area(m, l));
        if (std::ranges::count_if(areas, [](T x) { return x > T(0); }) != 1)
            return std::nullopt;
        const auto outer = static_cast<std::size_t>(std::ranges::max_element(areas) - areas.begin());
        // the window of the periodic parameters, centred on the face
        std::array<T, 4> lo_hi{std::numeric_limits<T>::max(), std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max(),
                               std::numeric_limits<T>::lowest()};
        for (const auto &l : loops)
            for (const auto &u : l)
            {
                const auto &ed = m.edge(u.ce.edge);
                for (int k = 0; k <= 8; ++k)
                {
                    const auto p = u.ce.pcurve->value(ed.u1 + (ed.u2 - ed.u1) * T(k) / T(8));
                    lo_hi = {std::min(lo_hi[0], p[0]), std::max(lo_hi[1], p[0]), std::min(lo_hi[2], p[1]), std::max(lo_hi[3], p[1])};
                }
            }
        std::array<T, 4> nb = box;
        for (int i = 0; i < 2; ++i)
        {
            const T lo = lo_hi[2 * i], hi = lo_hi[2 * i + 1];
            const bool outside = lo < box[2 * i] - eps_uv || hi > box[2 * i + 1] + eps_uv;
            if (per[i] > T(0))
            {
                // a face across the seam gets a window centred on it; one that goes round the
                // period starts its window on itself: at a vertex where it pinches (the double
                // point of a curve), or along an edge that becomes its seam
                if (!outside)
                    continue;
                if (!(hi - lo <= per[i] + T(100) * eps_uv))
                    return std::nullopt;
                const T mid = hi - lo < per[i] - T(100) * eps_uv ? (lo + hi) / T(2) : lo + per[i] / T(2);
                nb[2 * i] = mid - per[i] / T(2), nb[2 * i + 1] = mid + per[i] / T(2);
            }
            else if (fa.surface != fb.surface) // (pieces of one surface: the box held them both)
                nb[2 * i] = std::min(box[2 * i], lo), nb[2 * i + 1] = std::max(box[2 * i + 1], hi);
        }
        auto surface = fa.surface;
        if (nb != box)
        {
            const auto *e = dynamic_cast<const ElementarySurface<T> *>(fa.surface.get());
            if (!e)
                return std::nullopt; // (only a primitive can move its window)
            surface = e->with_box(nb);
        }

        // ---- write
        auto wire_of = [&](const std::vector<Use<T>> &l) {
            Wire<T> w;
            w.closed = true;
            for (const auto &u : l)
                w.coedges.push_back(u.ce);
            return m.add(std::move(w));
        };
        Face<T> nf{surface, {wire_of(loops[outer])}, std::max(fa.tol, fb.tol), false};
        for (std::size_t l = 0; l < loops.size(); ++l)
            if (l != outer)
                nf.wires.push_back(wire_of(loops[l]));
        const auto name = std::string(m.name(a.face));
        const auto id = m.add(std::move(nf));
        if (!name.empty())
            m.setName(id, name);
        for (auto sid : m.template ids<ShellId>())
        {
            auto &uses = m.shell(sid).faces;
            for (auto &u : uses)
                if (u.face == a.face)
                    u.face = id;
            std::erase_if(uses, [&](const FaceUse &u) { return u.face == b.face; });
        }
        for (auto w : fa.wires)
            m.erase(w);
        for (auto w : fb.wires)
            m.erase(w);
        m.erase(a.face);
        m.erase(b.face);
        return FaceUse{id, a.orient};
    }

    /**
     * @brief Puts a face of a sphere that holds a pole (a degenerate edge, and the seam to it) onto
     * the same sphere turned so that the face lies around its equator, away from the poles and
     * the seam: the face then has neither seam nor degenerate edge. The face must lie in an open
     * hemisphere (its boundary and the pole), around c, the mean direction of its boundary: the
     * axis of the new sphere is normal to c, its seam opposite c. The pcurves are projected anew
     * (`make_face_use`). Returns the new face, which replaces the old one in the shells; nothing
     * when the face is not such a face.
     */
    template <std::floating_point T>
    auto reframe_pole(Model<T> &m, const FaceUse &fu, T tol) -> std::optional<FaceUse>
    {
        const Face<T> f = m.face(fu.face);
        const auto *sph = dynamic_cast<const SphericalSurface<T> *>(f.surface.get());
        if (!sph)
            return std::nullopt;
        std::map<std::uint32_t, int> uses;
        std::size_t poles{};
        for (auto w : f.wires)
            for (const auto &ce : m.wire(w).coedges)
            {
                ++uses[ce.edge.index];
                poles += m.edge(ce.edge).degenerate;
            }
        if (poles != 1)
            return std::nullopt;
        // the loops of the face without the seam and the pole, in the sense of the face use
        const bool rev = fu.orient == Orientation::Reversed;
        std::vector<std::vector<CoEdge<T>>> loops;
        for (auto w : f.wires)
        {
            std::vector<CoEdge<T>> l;
            for (const auto &ce : m.wire(w).coedges)
                if (!m.edge(ce.edge).degenerate && uses[ce.edge.index] == 1)
                    l.push_back(CoEdge<T>{ce.edge, rev ? reverse(ce.orient) : ce.orient, nullptr});
            if (l.empty())
                return std::nullopt;
            if (rev)
                std::ranges::reverse(l);
            loops.push_back(std::move(l));
        }
        // the hemisphere: around the mean direction of the boundary, which holds it and the pole
        const auto &fr = sph->frame();
        const T R = sph->radius();
        std::vector<std::vector<point<T, 3>>> pts;
        point<T, 3> c{};
        for (const auto &l : loops)
        {
            auto &ps = pts.emplace_back();
            for (const auto &ce : l)
            {
                const auto &ed = m.edge(ce.edge);
                for (int k = 0; k < 8; ++k)
                {
                    const T s = T(k) / T(8);
                    const T t = ce.orient == Orientation::Forward ? ed.u1 + (ed.u2 - ed.u1) * s : ed.u2 - (ed.u2 - ed.u1) * s;
                    const auto d = (ed.curve->value(t) - fr.o) / R;
                    ps.push_back(d);
                    c = c + d;
                }
            }
        }
        if (!(norm(c) > T(1e-9)))
            return std::nullopt;
        c = c / norm(c);
        point<T, 3> pole{};
        for (auto w : f.wires)
            for (const auto &ce : m.wire(w).coedges)
                if (m.edge(ce.edge).degenerate)
                    pole = (m.vertex(m.edge(ce.edge).v1).pnt - fr.o) / R;
        constexpr T margin = T(0.05); // (cos of the angle kept from the great circle)
        if (!(pole * c > margin) || !std::ranges::all_of(pts, [&](const auto &ps) { return std::ranges::all_of(ps, [&](const auto &d) { return d * c > margin; }); }))
            return std::nullopt;
        // the outer loop: counter-clockwise about the normal of the use (the loops are in its
        // sense), in the plane tangent at c (central projection: one to one on the hemisphere)
        auto zc = fr.z - (fr.z * c) * c;
        if (!(norm(zc) > T(1e-6)))
            zc = norm(cross(c, fr.x)) > T(0.5) ? cross(c, fr.x) : cross(c, fr.y);
        zc = zc / norm(zc);
        const auto yc = cross(zc, c);
        std::vector<T> areas;
        for (const auto &ps : pts)
        {
            T a{};
            for (std::size_t i = 0, j = ps.size() - 1; i < ps.size(); j = i++)
            {
                const auto p = ps[j] / (ps[j] * c), q = ps[i] / (ps[i] * c);
                a += (p * yc) * (q * zc) - (q * yc) * (p * zc);
            }
            areas.push_back(rev ? -a / T(2) : a / T(2));
        }
        if (std::ranges::count_if(areas, [](T a) { return a > T(0); }) != 1)
            return std::nullopt;
        const auto outer = static_cast<std::size_t>(std::ranges::max_element(areas) - areas.begin());
        // ---- write: the face rebuilt on the turned sphere, its window centred on c
        constexpr T pi = std::numbers::pi_v<T>;
        auto surface = std::make_shared<SphericalSurface<T>>(R, ax2<T, 3>{fr.o, zc, c}, std::array<T, 4>{-pi, pi, -pi / T(2), pi / T(2)});
        auto wire_of = [&](std::vector<CoEdge<T>> l) {
            Wire<T> w;
            w.closed = true;
            w.coedges = std::move(l);
            return m.add(std::move(w));
        };
        const auto ow = wire_of(loops[outer]);
        std::vector<WireId> holes;
        for (std::size_t k = 0; k < loops.size(); ++k)
            if (k != outer)
                holes.push_back(wire_of(loops[k]));
        MakeFaceOptions<T> opts;
        opts.tol = std::max(f.tol, tol);
        auto r = make_face_use(m, std::static_pointer_cast<Surface<T, 3>>(surface), ow, holes, opts);
        if (!r)
        {
            m.erase(ow);
            for (auto w : holes)
                m.erase(w);
            return std::nullopt;
        }
        const auto name = std::string(m.name(fu.face));
        if (!name.empty())
            m.setName(r->face, name);
        // (the wires were given in the sense of the use: the new use keeps the senses in the shell)
        const FaceUse nu = *r;
        for (auto sid : m.template ids<ShellId>())
            for (auto &u : m.shell(sid).faces)
                if (u.face == fu.face)
                    u = nu;
        for (auto w : f.wires)
            m.erase(w);
        m.erase(fu.face);
        return nu;
    }

    /**
     * @brief Joins back the pieces of an edge: two edges on the same curve, end to end at a vertex
     * that no other edge of `faces` uses, consecutive with the same pcurve in every wire of
     * `faces` that uses them, become one edge (the inverse of `split_edge`). `faces` must hold
     * every face using these edges (a whole shell or solid). Returns the number of joins.
     */
    template <std::floating_point T>
    std::size_t join_edges(Model<T> &m, std::span<const FaceUse> faces)
    {
        std::size_t joined{};
        for (bool changed = true; changed;)
        {
            changed = false;
            std::map<std::uint32_t, std::set<std::uint32_t>> edges_at; // vertex -> edges
            std::set<std::uint32_t> wires;
            for (const auto &fu : faces)
                if (m.alive(fu.face))
                    for (auto w : m.face(fu.face).wires)
                    {
                        wires.insert(w.index);
                        for (const auto &ce : m.wire(w).coedges)
                        {
                            const auto &ed = m.edge(ce.edge);
                            edges_at[ed.v1.index].insert(ce.edge.index);
                            edges_at[ed.v2.index].insert(ce.edge.index);
                        }
                    }
            for (const auto &[v, es] : edges_at)
            {
                if (es.size() != 2)
                    continue;
                EdgeId e1{*es.begin()}, e2{*std::next(es.begin())};
                auto ed1 = m.edge(e1), ed2 = m.edge(e2);
                if (ed1.degenerate || ed2.degenerate || !ed1.curve || ed1.curve != ed2.curve)
                    continue;
                // end to end on the curve at v: e1 = [u1, u], e2 = [u, u2]
                if (ed1.v1.index == v && ed2.v2.index == v)
                {
                    std::swap(e1, e2);
                    std::swap(ed1, ed2);
                }
                if (!(ed1.v2.index == v && ed2.v1.index == v) || ed1.v1 == ed1.v2 || ed2.v1 == ed2.v2 ||
                    !(std::abs(ed1.u2 - ed2.u1) <= T(1e-12) * (T(1) + std::abs(ed1.u2))))
                    continue;
                // consecutive with one pcurve in every wire using them
                bool ok = true;
                for (auto wi : wires)
                {
                    const auto &cos = m.wire(WireId{wi}).coedges;
                    for (std::size_t i = 0; i < cos.size() && ok; ++i)
                    {
                        if (cos[i].edge != e1)
                            continue;
                        const auto &prev = cos[(i + cos.size() - 1) % cos.size()], &next = cos[(i + 1) % cos.size()];
                        const auto &other = cos[i].orient == Orientation::Forward ? next : prev;
                        ok = other.edge == e2 && other.orient == cos[i].orient && other.pcurve == cos[i].pcurve;
                    }
                    ok = ok && std::ranges::count_if(cos, [&](const auto &c) { return c.edge == e1; }) ==
                                   std::ranges::count_if(cos, [&](const auto &c) { return c.edge == e2; });
                }
                if (!ok)
                    continue;
                // ---- write: e1 covers both, e2 leaves the wires
                auto &x = m.edge(e1);
                x.u2 = ed2.u2;
                x.v2 = ed2.v2;
                x.tol = std::max(ed1.tol, ed2.tol);
                // (the ends hold the edge: a vertex tolerance is not below its edges')
                for (auto v : {x.v1, x.v2})
                    m.vertex(v).tol = std::max(m.vertex(v).tol, x.tol);
                for (auto wi : wires)
                    std::erase_if(m.wire(WireId{wi}).coedges, [&](const auto &c) { return c.edge == e2; });
                ++joined;
                changed = true;
                break;
            }
        }
        return joined;
    }

    /**
     * @brief Joins the pieces of a seam: two edges that are both the seam of one face of
     * `faces` on a primitive (each run twice by its boundary, on the iso-lines c and c + period
     * of one parameter), end to end at a vertex no other edge of `faces` uses, become one edge on
     * the iso-curve of the primitive, its pcurves the two iso-lines (exact). Pieces of a seam come
     * from merged faces whose edges lie on different curves, which `join_edges` keeps apart.
     * Returns the number of joins.
     */
    template <std::floating_point T>
    std::size_t join_seams(Model<T> &m, std::span<const FaceUse> faces)
    {
        std::size_t joined{};
        for (bool changed = true; changed;)
        {
            changed = false;
            std::map<std::uint32_t, std::set<std::uint32_t>> edges_at, faces_of;
            for (const auto &fu : faces)
                for (auto w : m.face(fu.face).wires)
                    for (const auto &ce : m.wire(w).coedges)
                    {
                        const auto &ed = m.edge(ce.edge);
                        edges_at[ed.v1.index].insert(ce.edge.index);
                        edges_at[ed.v2.index].insert(ce.edge.index);
                        faces_of[ce.edge.index].insert(fu.face.index);
                    }
            for (const auto &fu : faces)
            {
                const Face<T> f = m.face(fu.face);
                const auto *srf = dynamic_cast<const ElementarySurface<T> *>(f.surface.get());
                if (!srf)
                    continue;
                // the seam edges of the face: run twice, by this face only, along an iso-line
                struct Seam
                {
                    int iso;       // the constant parameter
                    T c1, c2;      // its values on the two runs
                    T t1, t2;      // the other parameter at the ends of the edge (its u1, u2)
                };
                std::map<std::uint32_t, std::vector<const CoEdge<T> *>> runs;
                for (auto w : f.wires)
                    for (const auto &ce : m.wire(w).coedges)
                        runs[ce.edge.index].push_back(&ce);
                std::map<std::uint32_t, Seam> seams;
                for (const auto &[e, rs] : runs)
                {
                    const auto &ed = m.edge(EdgeId{e});
                    if (rs.size() != 2 || ed.degenerate || faces_of[e].size() != 1 || !rs[0]->pcurve || !rs[1]->pcurve)
                        continue;
                    const T tm = (ed.u1 + ed.u2) / T(2);
                    const auto a1 = rs[0]->pcurve->value(ed.u1), a2 = rs[0]->pcurve->value(ed.u2), am = rs[0]->pcurve->value(tm);
                    const auto b1 = rs[1]->pcurve->value(ed.u1), bm = rs[1]->pcurve->value(tm);
                    const T eps = T(1e-9) * (T(1) + std::abs(a1[0]) + std::abs(a1[1]));
                    for (int i = 0; i < 2; ++i)
                        if (std::abs(a1[i] - a2[i]) <= eps && std::abs(a1[i] - am[i]) <= eps && std::abs(b1[i] - bm[i]) <= eps &&
                            std::abs(a1[1 - i] - b1[1 - i]) <= eps)
                            seams[e] = Seam{i, a1[i], b1[i], a1[1 - i], a2[1 - i]};
                }
                for (const auto &[e1, s1] : seams)
                {
                    for (const auto &[e2, s2] : seams)
                    {
                        if (e2 <= e1 || s1.iso != s2.iso || std::abs(std::min(s1.c1, s1.c2) - std::min(s2.c1, s2.c2)) > T(1e-9) * (T(1) + std::abs(s1.c1)))
                            continue;
                        const auto &d1 = m.edge(EdgeId{e1}), &d2 = m.edge(EdgeId{e2});
                        // the vertex between them, used by these two edges only
                        std::optional<VertexId> v;
                        for (auto a : {d1.v1, d1.v2})
                            for (auto b : {d2.v1, d2.v2})
                                if (a == b && edges_at[a.index].size() == 2)
                                    v = a;
                        if (!v)
                            continue;
                        // the other parameter on the joined edge: from the far end of e1 to the far end of e2
                        const T far1 = d1.v1 == *v ? s1.t2 : s1.t1, far2 = d2.v1 == *v ? s2.t2 : s2.t1;
                        const T mid1 = d1.v1 == *v ? s1.t1 : s1.t2, mid2 = d2.v1 == *v ? s2.t1 : s2.t2;
                        if (std::abs(mid1 - mid2) > T(1e-9) * (T(1) + std::abs(mid1)) || !((far1 - mid1) * (far2 - mid2) < T(0)))
                            continue; // (not end to end along the iso-line)
                        const VertexId vfar1 = d1.v1 == *v ? d1.v2 : d1.v1, vfar2 = d2.v1 == *v ? d2.v2 : d2.v1;
                        const bool up = far1 < far2; // the new edge runs from the lower end
                        const T lo = up ? far1 : far2, hi = up ? far2 : far1;
                        const VertexId vlo = up ? vfar1 : vfar2, vhi = up ? vfar2 : vfar1;
                        const T cmin = std::min(s1.c1, s1.c2);
                        auto curve = s1.iso == 0 ? srf->iso_u(cmin, lo, hi) : srf->iso_v(cmin, lo, hi);
                        if (!curve)
                            continue;
                        const auto at = [&](T t) { return s1.iso == 0 ? f.surface->value(cmin, t) : f.surface->value(t, cmin); };
                        const T tm = (lo + hi) / T(2);
                        if (norm(curve->value(lo) - at(lo)) > f.tol || norm(curve->value(tm) - at(tm)) > f.tol ||
                            norm(curve->value(hi) - at(hi)) > f.tol)
                            continue; // (an iso-curve not parametrized by the other parameter)
                        const auto e = m.add(Edge<T>{std::static_pointer_cast<Curve<T, 3>>(curve), lo, hi, vlo, vhi, std::max(d1.tol, d2.tol)});
                        for (auto v : {vlo, vhi})
                            m.vertex(v).tol = std::max(m.vertex(v).tol, m.edge(e).tol);
                        auto line = [&](T c) {
                            const point<T, 2> p = s1.iso == 0 ? point<T, 2>{c, lo} : point<T, 2>{lo, c}, q = s1.iso == 0 ? point<T, 2>{c, hi} : point<T, 2>{hi, c};
                            return std::static_pointer_cast<Curve<T, 2>>(
                                std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{p, q}, std::vector<T>{lo, lo, hi, hi}, 1));
                        };
                        // ---- write: in each run, the two co-edges become one (they follow each other)
                        for (auto w : f.wires)
                        {
                            auto &cos = m.wire(w).coedges;
                            for (std::size_t i = 0; i < cos.size(); ++i)
                            {
                                if (cos[i].edge.index != e1)
                                    continue;
                                // the run of e1 along increasing t keeps the constant of its pcurve
                                const auto &ce = cos[i];
                                const T c = ce.pcurve->value(d1.u1)[s1.iso];
                                const bool fwd1 = ce.orient == Orientation::Forward;
                                // sense of the run in t: e1 from its start to its end in its own sense
                                const T from = fwd1 ? (s1.t1) : (s1.t2), to = fwd1 ? s1.t2 : s1.t1;
                                const bool increasing = to > from;
                                cos[i] = CoEdge<T>{e, increasing ? Orientation::Forward : Orientation::Reversed, line(c)};
                            }
                            std::erase_if(cos, [&](const auto &c) { return c.edge.index == e2; });
                        }
                        ++joined;
                        changed = true;
                        break;
                    }
                    if (changed)
                        break;
                }
                if (changed)
                    break;
            }
        }
        return joined;
    }

    /**
     * @brief Unifies the faces of one domain in `faces` (`merge_faces`, pairs that share an edge,
     * until none merges; a plane in B-spline form put on its plane first when its neighbour is
     * coplanar, `to_plane`), puts the faces of spheres around a pole on a turned sphere
     * (`reframe_pole`), then joins the pieces of edges and of seams left (`join_edges`,
     * `join_seams`). `faces` is updated: the new faces replace the old ones.
     */
    template <std::floating_point T>
    UnifyReport unify_faces(Model<T> &m, std::vector<FaceUse> &faces, T tol)
    {
        UnifyReport r;
        std::set<std::pair<std::uint32_t, std::uint32_t>> failed; // (pairs of faces tried in vain)
        for (bool changed = true; changed;)
        {
            changed = false;
            std::map<std::uint32_t, std::vector<std::size_t>> by_edge;
            for (std::size_t k = 0; k < faces.size(); ++k)
                for (auto e : explore<EdgeId>(m, ShapeId{faces[k].face}))
                    by_edge[e.index].push_back(k);
            for (const auto &[e, ks] : by_edge)
            {
                if (ks.size() != 2 || ks[0] == ks[1] || failed.contains({faces[ks[0]].face.index, faces[ks[1]].face.index}))
                    continue;
                auto merged = merge_faces(m, faces[ks[0]], faces[ks[1]], tol);
                if (!merged)
                {
                    // a plane in B-spline form next to a face of the same plane: both put on the
                    // primitive, then merged
                    const auto pa = unify_detail::plane_of(m, m.face(faces[ks[0]].face), tol),
                               pb = unify_detail::plane_of(m, m.face(faces[ks[1]].face), tol);
                    if (!pa || !pb || !unify_detail::coplanar(*pa, *pb, tol) ||
                        (dynamic_cast<const Plane<T> *>(m.face(faces[ks[0]].face).surface.get()) &&
                         dynamic_cast<const Plane<T> *>(m.face(faces[ks[1]].face).surface.get())))
                    {
                        failed.insert({faces[ks[0]].face.index, faces[ks[1]].face.index});
                        continue;
                    }
                    for (auto k : ks)
                        if (auto p = to_plane(m, faces[k], tol))
                        {
                            r.history.push_back({p->face, faces[k].face, FaceId{}});
                            ++r.faces_planed;
                            faces[k] = *p;
                        }
                    merged = merge_faces(m, faces[ks[0]], faces[ks[1]], tol);
                    changed = true; // (the faces changed: the edges are explored again)
                    if (!merged)
                        break;
                }
                {
                    r.history.push_back({merged->face, faces[ks[0]].face, faces[ks[1]].face});
                    faces[ks[0]] = *merged;
                    faces.erase(faces.begin() + static_cast<std::ptrdiff_t>(ks[1]));
                    ++r.faces_merged;
                    changed = true;
                    break;
                }
            }
        }
        for (auto &fu : faces)
            if (auto p = reframe_pole(m, fu, tol))
            {
                r.history.push_back({p->face, fu.face, FaceId{}});
                ++r.poles_reframed;
                fu = *p;
            }
        r.edges_joined = join_edges(m, faces) + join_seams(m, faces);
        return r;
    }
}
