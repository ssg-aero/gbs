#pragma once
/**
 * @file offset.h
 * @brief Offset of a face: the face at the distance d along the normal of its surface.
 *
 * The offset surface keeps the parametrization of the surface, so that the pcurves of the face
 * are reused as they are (a cone shifts its v by -d sin a: its pcurves are shifted alike):
 * - primitives: exact offset primitive (plane moved, cylinder, sphere and cone of another radius,
 *   torus of another minor radius), checked against S + d n at a few points;
 * - an offset surface: the offset of its basis at the sum of the distances;
 * - other surfaces: an offset surface S + d n with exact derivatives (`gbs/offsets.h`).
 * The edges are rebuilt as the pcurves on the offset surface (`CurveOnSurface`), the vertices at
 * their ends (a vertex of the face gives one per point it moves to); the face of the offset keeps
 * the topology and the tolerance of the face, except that a pole whose offset is not a point (the
 * apex of a cone) becomes an edge: a circle, or an arc between two vertices for a sector. An
 * offset that folds on the face (beyond a radius of curvature: `offset_margin`, or no primitive of
 * the same kind for a primitive) is refused.
 */
#include <gbs-brep/builders.h>
#include <gbs-brep/classify.h>
#include <gbs/bssanalysis.h>
#include <gbs/curveonsurface.h>
#include <gbs/offsets.h>
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
    using gbs::OffsetSurface;

    namespace offset_detail
    {
        template <std::floating_point T>
        auto shifted(const std::shared_ptr<Curve<T, 2>> &pc, T dv) -> std::shared_ptr<Curve<T, 2>>
        {
            return shifted_pcurve(pc, point<T, 2>{T(0), dv});
        }
    }

    /**
     * @brief The face at the distance `d` of the face `f`, along the normal of its surface
     * (Su x Sv; a negative d goes the other way). The face `f` is not modified; on failure the
     * model is left as it was. An offset that folds on the face (d beyond a radius of curvature
     * of its surface, measured on the box of its pcurves) gives `InvalidBounds`.
     */
    template <std::floating_point T>
    [[nodiscard]] auto offset_face(Model<T> &m, FaceId f, std::type_identity_t<T> d, std::type_identity_t<T> tol = brep_default_tolerance<T>)
        -> BuildResult<FaceId>
    {
        if (!m.alive(f))
            return detail::fail(BuildErrc::InvalidId, "offset_face: dead or invalid face", {f});
        const auto face = m.face(f);
        // exact for a primitive or an offset surface; refused where the offset folds on the face
        const auto off = offset_surface(face.surface, d, tol, classify_detail::face_uv_box(m, f));
        if (!off)
            return detail::fail(BuildErrc::InvalidBounds, "offset_face: the offset folds (beyond a radius of curvature)", {f});
        const auto &srf = off->surface;
        const T dv = off->dv;

        const auto size = detail::arena_size(m);
        // the offset vertices of each vertex of the face: one per point it moves to (the two ends of
        // the arc an apex of a sector of a cone becomes)
        std::map<std::uint32_t, std::vector<VertexId>> vertices;
        std::map<std::uint32_t, EdgeId> edges;
        auto vertex = [&](VertexId v, const point<T, 3> &p) {
            const T t = std::max(tol, m.vertex(v).tol);
            auto &done = vertices[v.index];
            for (auto id : done)
                if (norm(m.vertex(id).pnt - p) <= t)
                    return id;
            done.push_back(m.add(Vertex<T>{p, t}));
            return done.back();
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
                {
                    detail::rollback(m, size);
                    return detail::fail(BuildErrc::PCurveApproximation, "offset_face: a pcurve cannot be shifted", {ce.edge});
                }
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
