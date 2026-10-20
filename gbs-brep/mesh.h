#pragma once
/**
 * @file mesh.h
 * @brief Triangle meshes of the faces, shells and solids of the native BREP (stage 3,
 * `docs/sources/design/brep_palier3_mesh.md`).
 *
 * This part: the options, the discretization of the edges and the boundaries of the faces in
 * their parameter space.
 *
 * - Each edge is discretized once, on the parameter of its curve, until the chordal deviation,
 *   the turn of the tangent and the length of its segments are within the options; every face
 *   using the edge gets the same points, hence the same vertices of the mesh (watertight).
 * - The boundary of a face is a list of loops of points of (u, v), the outer one first, each
 *   point with its vertex of the mesh: the pcurves of the co-edges at the parameters of their
 *   edges, in the order and the sense of the wire. A seam appears twice (at u1 and at u2) with
 *   the same vertices; all the points of a degenerate edge (a pole, an apex) are one vertex.
 */
#include <gbs-brep/builders.h>
#include <gbs-brep/closure.h>
#include <gbs-brep/explore.h>
#include <gbs-brep/model.h>
#include <gbs/bscurve.h>
#include <gbs/intersect/curve_pieces.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <vector>

#ifdef GBS_USE_MODULES
    import vecop;
#endif

namespace gbs::brep
{
    /// Precision of a mesh.
    template <std::floating_point T>
    struct MeshOptions
    {
        T deflection{};   ///< largest distance between the mesh and the geometry (0: 1e-3 of the diagonal of the shape)
        T angle{T(0.5)};  ///< largest turn of the tangent along a segment, of the normal across a triangle (radians)
        T max_edge{};     ///< longest segment or side of a triangle (0: none)
    };

    namespace mesh_detail
    {
        /// Angle between two vectors (0 when one is null).
        template <std::floating_point T>
        T angle_between(const point<T, 3> &a, const point<T, 3> &b)
        {
            const T na = norm(a), nb = norm(b);
            if (!(na > T(0)) || !(nb > T(0)))
                return T(0);
            return std::acos(std::clamp((a * b) / (na * nb), T(-1), T(1)));
        }

        /**
         * A bound of the turn of the tangent along a piece: the sum of the angles between the
         * consecutive legs of its control polygon (the Bézier poles; the samples of a sampled
         * piece).
         */
        template <std::floating_point T>
        T turn(const CurvePiece<T, 3> &pc)
        {
            std::vector<point<T, 3>> legs;
            for (std::size_t k = 0; k + 1 < pc.poles.size(); ++k)
            {
                const auto a = intersect_detail::cartesian<T, 3>(pc.poles[k]), b = intersect_detail::cartesian<T, 3>(pc.poles[k + 1]);
                if (norm(b - a) > T(0))
                    legs.push_back(b - a);
            }
            T t{};
            for (std::size_t k = 0; k + 1 < legs.size(); ++k)
                t += angle_between(legs[k], legs[k + 1]);
            return t;
        }
    }

    /**
     * @brief Parameters of the points of the curve c on [a, b]: the ends of its pieces (its Bézier
     * segments, `CurvePieces`), each piece halved while it may be farther than `deflection` from
     * its chord (its poles bound it: convex hull), its tangent may turn by more than `angle` (the
     * turn of its control polygon), or its chord is longer than `max_edge`. A curve which is not a
     * B-spline is sampled instead. Throws `std::invalid_argument` for a deflection or an angle
     * that is not positive.
     */
    template <std::floating_point T>
    [[nodiscard]] auto curve_parameters(const Curve<T, 3> &c, T a, T b, const MeshOptions<T> &o) -> std::vector<T>
    {
        if (!(o.deflection > T(0)) || !(o.angle > T(0)))
            throw std::invalid_argument("mesh: the deflection and the angle must be positive");
        CurvePieces<T, 3> pieces{c};
        const T min_step = T(1e-9) * (b - a);
        auto split = [&](const CurvePiece<T, 3> &pc) {
            if (!(pc.b - pc.a > min_step))
                return false;
            return pc.flatness > o.deflection || mesh_detail::turn(pc) > o.angle || (o.max_edge > T(0) && pc.chord() > o.max_edge);
        };
        std::vector<T> r{a};
        for (auto &first : pieces.pieces(a, b))
        {
            // depth first: the parameters come out in order
            std::vector<CurvePiece<T, 3>> stack{std::move(first)};
            while (!stack.empty())
            {
                auto pc = std::move(stack.back());
                stack.pop_back();
                if (split(pc))
                {
                    auto [l, rr] = pieces.split(pc);
                    stack.push_back(std::move(rr));
                    stack.push_back(std::move(l));
                }
                else
                    r.push_back(pc.b);
            }
        }
        return r;
    }

    /// A point of a boundary of a face: its parameters, its vertex in the mesh.
    template <std::floating_point T>
    struct BoundaryPoint
    {
        point<T, 2> uv;
        std::uint32_t vertex;
    };

    /**
     * @brief The vertices of a mesh along the edges of a model: each model vertex and each edge
     * discretized once, shared by the faces that use them.
     */
    template <std::floating_point T>
    class MeshBuilder
    {
    public:
        /// The points of an edge: parameters on its curve (or on its pcurve when degenerate), vertices of the mesh.
        struct EdgePoints
        {
            std::vector<T> t;
            std::vector<std::uint32_t> vertices;
        };

        /// The options as given: a deflection and an angle that are not positive throw (see `resolved`).
        MeshBuilder(const Model<T> &m, const MeshOptions<T> &o) : m_{m}, o_{o}
        {
            if (!(o.deflection > T(0)) || !(o.angle > T(0)))
                throw std::invalid_argument("mesh: the deflection and the angle must be positive");
        }
        /// The options resolved for the shape s: a deflection of 0 is 1e-3 of its diagonal.
        MeshBuilder(const Model<T> &m, const ShapeId &s, const MeshOptions<T> &o);

        [[nodiscard]] auto options() const -> const MeshOptions<T> & { return o_; }
        [[nodiscard]] auto vertices() const -> const std::vector<point<T, 3>> & { return vertices_; }

        /// The vertex of the mesh at a vertex of the model.
        std::uint32_t vertex(VertexId v)
        {
            if (v.index >= of_vertex_.size())
                of_vertex_.resize(v.index + 1, none);
            auto &x = of_vertex_[v.index];
            if (x == none)
            {
                x = static_cast<std::uint32_t>(vertices_.size());
                vertices_.push_back(m_.vertex(v).pnt);
            }
            return x;
        }

        /// The points of a non-degenerate edge.
        auto edge(EdgeId e) -> const EdgePoints &
        {
            if (auto *known = find(e))
                return *known;
            const auto &ed = m_.edge(e);
            EdgePoints p;
            p.t = curve_parameters(*ed.curve, ed.u1, ed.u2, o_);
            p.vertices.reserve(p.t.size());
            p.vertices.push_back(vertex(ed.v1));
            for (std::size_t k = 1; k + 1 < p.t.size(); ++k)
            {
                p.vertices.push_back(static_cast<std::uint32_t>(vertices_.size()));
                vertices_.push_back(ed.curve->value(p.t[k]));
            }
            p.vertices.push_back(vertex(ed.v2));
            return store(e, std::move(p));
        }

        /**
         * The points of a degenerate edge (a pole, an apex): on its pcurve, as many as the angle
         * asks for a full turn around the point, in proportion of the part of the face it spans
         * (its whole side: a full turn); all of them the vertex of the edge.
         */
        auto degenerate_edge(EdgeId e, const Curve<T, 2> &pc, const std::array<T, 4> &box) -> const EdgePoints &
        {
            if (auto *known = find(e))
                return *known;
            const auto &ed = m_.edge(e);
            const auto a = pc.value(ed.u1), b = pc.value(ed.u2);
            const T du = std::abs(b[0] - a[0]) / std::max(box[1] - box[0], std::numeric_limits<T>::min());
            const T dv = std::abs(b[1] - a[1]) / std::max(box[3] - box[2], std::numeric_limits<T>::min());
            const T turn = std::max(du, dv) * T(2) * std::numbers::pi_v<T>;
            const auto n = static_cast<std::size_t>(std::max(T(2), std::ceil(turn / o_.angle)));
            EdgePoints p;
            const auto v = vertex(ed.v1);
            for (std::size_t k = 0; k <= n; ++k)
            {
                p.t.push_back(ed.u1 + (ed.u2 - ed.u1) * T(k) / T(n));
                p.vertices.push_back(v);
            }
            return store(e, std::move(p));
        }

        /**
         * @brief The boundary of the face f in its parameter space: one loop per wire (the outer one
         * first), the points of each co-edge in the order and the sense of the wire, the last
         * point of a co-edge being the first of the next.
         */
        auto boundary(FaceId f) -> BuildResult<std::vector<std::vector<BoundaryPoint<T>>>>
        {
            const auto &face = m_.face(f);
            std::optional<std::array<T, 4>> box; // (for the degenerate edges only)
            std::vector<std::vector<BoundaryPoint<T>>> loops;
            for (auto w : face.wires)
            {
                std::vector<BoundaryPoint<T>> loop;
                for (const auto &ce : m_.wire(w).coedges)
                {
                    if (!ce.pcurve)
                        return detail::fail(BuildErrc::MissingPCurve, "mesh: a co-edge of the face has no pcurve", {f, ce.edge});
                    const auto &ed = m_.edge(ce.edge);
                    if (!ed.same_parameter)
                        return detail::fail(BuildErrc::NotSameParameter, "mesh: an edge is not parametrized like its pcurves", {f, ce.edge});
                    if (ed.degenerate && !box)
                        box = face_parameter_box(f);
                    const auto &pts = ed.degenerate ? degenerate_edge(ce.edge, *ce.pcurve, *box) : edge(ce.edge);
                    const std::size_t n = pts.t.size();
                    for (std::size_t i = 0; i + 1 < n; ++i)
                    {
                        const std::size_t k = ce.orient == Orientation::Forward ? i : n - 1 - i;
                        loop.push_back({ce.pcurve->value(pts.t[k]), pts.vertices[k]});
                    }
                }
                loops.push_back(std::move(loop));
            }
            return loops;
        }

    private:
        const Model<T> &m_;
        MeshOptions<T> o_;
        static constexpr std::uint32_t none = std::numeric_limits<std::uint32_t>::max();
        std::vector<point<T, 3>> vertices_;
        std::vector<std::uint32_t> of_vertex_;           ///< by model vertex index
        std::vector<std::optional<EdgePoints>> edges_;   ///< by model edge index

        const EdgePoints *find(EdgeId e) const { return e.index < edges_.size() && edges_[e.index] ? &*edges_[e.index] : nullptr; }
        const EdgePoints &store(EdgeId e, EdgePoints p)
        {
            if (e.index >= edges_.size())
                edges_.resize(e.index + 1);
            return *(edges_[e.index] = std::move(p));
        }

        /// The box of the pcurves of the face (its parameter space in use).
        std::array<T, 4> face_parameter_box(FaceId f) const
        {
            std::array<T, 4> b{std::numeric_limits<T>::max(), std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max(),
                               std::numeric_limits<T>::lowest()};
            for (auto w : m_.face(f).wires)
                for (const auto &ce : m_.wire(w).coedges)
                {
                    if (!ce.pcurve)
                        continue;
                    const auto &ed = m_.edge(ce.edge);
                    for (int k = 0; k <= 8; ++k)
                    {
                        const auto uv = ce.pcurve->value(ed.u1 + (ed.u2 - ed.u1) * T(k) / T(8));
                        b = {std::min(b[0], uv[0]), std::max(b[1], uv[0]), std::min(b[2], uv[1]), std::max(b[3], uv[1])};
                    }
                }
            return b;
        }
    };

    /// The options with their defaults resolved for the shape s: the deflection 1e-3 of its diagonal.
    template <std::floating_point T>
    [[nodiscard]] auto resolved(const Model<T> &m, const ShapeId &s, MeshOptions<T> o) -> MeshOptions<T>
    {
        if (!(o.deflection > T(0)))
        {
            const auto b = bounding_box(m, s);
            o.deflection = T(1e-3) * std::max(norm(b.max - b.min), std::numeric_limits<T>::min());
        }
        return o;
    }

    template <std::floating_point T>
    MeshBuilder<T>::MeshBuilder(const Model<T> &m, const ShapeId &s, const MeshOptions<T> &o) : MeshBuilder{m, resolved(m, s, o)}
    {
    }
}
