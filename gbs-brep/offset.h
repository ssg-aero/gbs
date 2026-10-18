#pragma once
/**
 * @file offset.h
 * @brief Offset of a face: the face at the distance d along the normal of its surface.
 *
 * The offset surface keeps the parametrization of the surface, so that the pcurves of the face
 * are reused as they are (a cone shifts its v by -d sin a: its pcurves are shifted alike):
 * - primitives: exact offset primitive (plane moved, cylinder, sphere and cone of another radius,
 *   torus of another minor radius), checked against S + d n at a few points;
 * - other surfaces: an offset surface S + d n with derivatives by finite differences.
 * The edges are rebuilt as the pcurves on the offset surface (`CurveOnSurface`), the vertices at
 * their ends; the face of the offset keeps the topology and the tolerance of the face, except that
 * a pole whose offset is not a point (the apex of a cone) becomes a closed edge.
 */
#include <gbs-brep/builders.h>
#include <gbs-brep/classify.h>
#include <gbs/curveonsurface.h>
#include <gbs/primitives.h>

#include <cmath>
#include <map>
#include <memory>
#include <optional>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs::brep
{
    /// S(u, v) + d n(u, v), first derivatives by central finite differences.
    template <std::floating_point T>
    class OffsetSurface : public Surface<T, 3>
    {
        std::shared_ptr<Surface<T, 3>> s_;
        T d_;

        point<T, 3> at(T u, T v) const
        {
            const auto n = cross(s_->value(u, v, 1, 0), s_->value(u, v, 0, 1));
            return s_->value(u, v) + (d_ / norm(n)) * n;
        }

    public:
        OffsetSurface(std::shared_ptr<Surface<T, 3>> s, T d) : s_{std::move(s)}, d_{d} {}
        [[nodiscard]] auto basis() const -> const std::shared_ptr<Surface<T, 3>> & { return s_; }
        [[nodiscard]] T distance() const { return d_; }
        auto bounds() const -> std::array<T, 4> override { return s_->bounds(); }
        auto has_second_derivatives() const -> bool override { return false; }
        auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, 3> override
        {
            if (du == 0 && dv == 0)
                return at(u, v);
            const auto b = bounds();
            const T hu = T(1e-6) * (b[1] - b[0]), hv = T(1e-6) * (b[3] - b[2]);
            if (du == 1 && dv == 0)
            {
                const T a = std::max(b[0], u - hu), c = std::min(b[1], u + hu);
                return (at(c, v) - at(a, v)) / (c - a);
            }
            if (du == 0 && dv == 1)
            {
                const T a = std::max(b[2], v - hv), c = std::min(b[3], v + hv);
                return (at(u, c) - at(u, a)) / (c - a);
            }
            throw std::runtime_error("OffsetSurface: only first derivatives");
        }
    };

    namespace offset_detail
    {
        /// Exact offset of a primitive and the shift of its v, when S' (u, v + dv) = S(u, v) + d n.
        template <std::floating_point T>
        auto primitive_offset(const ElementarySurface<T> &s, T d, T tol) -> std::optional<std::pair<std::shared_ptr<Surface<T, 3>>, T>>
        {
            const auto &f = s.frame();
            const auto box = s.bounds();
            std::vector<std::pair<std::shared_ptr<ElementarySurface<T>>, T>> cands;
            if (dynamic_cast<const Plane<T> *>(&s))
                for (T sg : {T(1), T(-1)})
                {
                    auto g = f;
                    g.o = f.o + sg * d * f.z;
                    cands.push_back({std::make_shared<Plane<T>>(g, box), T(0)});
                }
            else if (const auto *c = dynamic_cast<const CylindricalSurface<T> *>(&s))
                for (T sg : {T(1), T(-1)})
                {
                    if (c->radius() + sg * d > T(0))
                        cands.push_back({std::make_shared<CylindricalSurface<T>>(c->radius() + sg * d, f, box), T(0)});
                }
            else if (const auto *sp = dynamic_cast<const SphericalSurface<T> *>(&s))
                for (T sg : {T(1), T(-1)})
                {
                    if (sp->radius() + sg * d > T(0))
                        cands.push_back({std::make_shared<SphericalSurface<T>>(sp->radius() + sg * d, f, box), T(0)});
                }
            else if (const auto *to = dynamic_cast<const ToroidalSurface<T> *>(&s))
                for (T sg : {T(1), T(-1)})
                {
                    if (to->minor_radius() + sg * d > T(0))
                        cands.push_back({std::make_shared<ToroidalSurface<T>>(to->major_radius(), to->minor_radius() + sg * d, f, box), T(0)});
                }
            else if (const auto *co = dynamic_cast<const ConicalSurface<T> *>(&s))
            {
                const T a = co->semi_angle();
                for (T sg : {T(1), T(-1)})
                {
                    const T dv = -sg * d * std::sin(a), R = co->radius() + sg * d / std::cos(a);
                    if (R >= T(0))
                        cands.push_back({std::make_shared<ConicalSurface<T>>(R, a, f, std::array<T, 4>{box[0], box[1], box[2] + dv, box[3] + dv}), dv});
                }
            }
            // the candidate where S' (u, v + dv) = S(u, v) + d n
            for (const auto &[c, dv] : cands)
            {
                bool ok = true;
                for (T fu : {T(0.3), T(0.7)})
                    for (T fv : {T(0.4), T(0.6)})
                    {
                        const T u = box[0] + fu * (box[1] - box[0]), v = box[2] + fv * (box[3] - box[2]);
                        const auto n = cross(s.value(u, v, 1, 0), s.value(u, v, 0, 1));
                        const auto target = s.value(u, v) + (d / norm(n)) * n;
                        ok = ok && norm(c->value(u, v + dv) - target) <= tol;
                    }
                if (ok)
                    return std::pair<std::shared_ptr<Surface<T, 3>>, T>{c, dv};
            }
            return std::nullopt;
        }

        template <std::floating_point T>
        auto shifted(const std::shared_ptr<Curve<T, 2>> &pc, T dv) -> std::shared_ptr<Curve<T, 2>>
        {
            if (dv == T(0))
                return pc;
            if (const auto *bs = dynamic_cast<const BSCurve<T, 2> *>(pc.get()))
            {
                auto poles = bs->poles();
                for (auto &q : poles)
                    q[1] += dv;
                return std::make_shared<BSCurve<T, 2>>(poles, bs->knotsFlats(), bs->degree());
            }
            return nullptr;
        }
    }

    /**
     * @brief The face at the distance `d` of the face `f`, along the normal of its surface
     * (Su x Sv; a negative d goes the other way). The face `f` is not modified.
     */
    template <std::floating_point T>
    [[nodiscard]] auto offset_face(Model<T> &m, FaceId f, std::type_identity_t<T> d, std::type_identity_t<T> tol = brep_default_tolerance<T>)
        -> BuildResult<FaceId>
    {
        if (!m.alive(f))
            return detail::fail(BuildErrc::InvalidId, "offset_face: dead or invalid face", {f});
        const auto face = m.face(f);
        std::shared_ptr<Surface<T, 3>> srf;
        T dv{};
        if (const auto *es = dynamic_cast<const ElementarySurface<T> *>(face.surface.get()))
            if (auto po = offset_detail::primitive_offset(*es, d, tol))
                std::tie(srf, dv) = *po;
        if (!srf)
            srf = std::make_shared<OffsetSurface<T>>(face.surface, d);

        std::map<std::uint32_t, VertexId> vertices;
        std::map<std::uint32_t, EdgeId> edges;
        auto vertex = [&](VertexId v, const point<T, 3> &p) {
            if (auto it = vertices.find(v.index); it != vertices.end())
                return it->second;
            const auto id = m.add(Vertex<T>{p, std::max(tol, m.vertex(v).tol)});
            vertices.emplace(v.index, id);
            return id;
        };
        Face<T> nf{srf, {}, face.tol, face.natural_bounds};
        for (auto w : face.wires)
        {
            Wire<T> nw;
            nw.closed = m.wire(w).closed;
            for (const auto &ce : m.wire(w).coedges)
            {
                auto pc = offset_detail::shifted(ce.pcurve, dv);
                if (!pc)
                    return detail::fail(BuildErrc::PCurveApproximation, "offset_face: a pcurve cannot be shifted", {ce.edge});
                const auto &e = m.edge(ce.edge);
                auto it = edges.find(ce.edge.index);
                if (it == edges.end())
                {
                    Edge<T> ne;
                    ne.u1 = e.u1;
                    ne.u2 = e.u2;
                    ne.tol = std::max(tol, e.tol);
                    // the end points on the offset surface: where the pcurve of this co-edge starts and ends
                    const auto q1 = pc->value(e.u1), q2 = pc->value(e.u2);
                    const auto p1 = srf->value(q1[0], q1[1]);
                    ne.v1 = vertex(e.v1, p1);
                    ne.v2 = vertex(e.v2, srf->value(q2[0], q2[1]));
                    // a pole stays a pole only if its offset is a point: the apex of a cone becomes a circle
                    ne.degenerate = e.degenerate;
                    if (e.degenerate)
                        for (int k = 1; k <= 8; ++k)
                        {
                            const auto q = pc->value(e.u1 + (e.u2 - e.u1) * T(k) / T(8));
                            ne.degenerate = ne.degenerate && norm(srf->value(q[0], q[1]) - p1) <= ne.tol;
                        }
                    if (!ne.degenerate)
                        ne.curve = std::make_shared<CurveOnSurface<T, 3>>(pc, srf);
                    it = edges.emplace(ce.edge.index, m.add(std::move(ne))).first;
                }
                nw.coedges.push_back(CoEdge<T>{it->second, ce.orient, pc});
            }
            nf.wires.push_back(m.add(std::move(nw)));
        }
        return m.add(std::move(nf));
    }
}
