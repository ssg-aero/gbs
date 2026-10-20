#pragma once
/**
 * @file mesh.h
 * @brief Triangle meshes of the faces, shells and solids of the native BREP (stage 3,
 * `docs/sources/design/brep_palier3_mesh.md`).
 *
 * This part: the options, the discretization of the edges, the boundaries of the faces in their
 * parameter space and the triangulation of each face.
 *
 * - Each edge is discretized once, on the parameter of its curve, until the chordal deviation,
 *   the turn of the tangent and the length of its segments are within the options; every face
 *   using the edge gets the same points, hence the same vertices of the mesh (watertight).
 * - The boundary of a face is a list of loops of points of (u, v), the outer one first, each
 *   point with its vertex of the mesh: the pcurves of the co-edges at the parameters of their
 *   edges, in the order and the sense of the wire. A seam appears twice (at u1 and at u2) with
 *   the same vertices; all the points of a degenerate edge (a pole, an apex) are one vertex.
 * - A face is triangulated in its parameter space scaled by the mean lengths of S_u and S_v (about
 *   isotropic in 3D): its boundary loops are the constrained edges of a Delaunay triangulation
 *   (`gbs-mesh/cdt.h`), refined inside by circumcentres until each triangle is within the
 *   deflection, the angle and the length of the options. The boundary is never refined: its
 *   points are those of the edges, shared with the other faces.
 */
#include <gbs-brep/builders.h>
#include <gbs-brep/closure.h>
#include <gbs-brep/explore.h>
#include <gbs-brep/model.h>
#include <gbs/bscurve.h>
#include <gbs/intersect/curve_pieces.h>
#include <gbs-mesh/cdt.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <unordered_set>
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

        /// The distance of q to the triangle (a, b, c), which may be degenerate (a pole).
        template <std::floating_point T>
        T distance_to_triangle(const point<T, 3> &q, const point<T, 3> &a, const point<T, 3> &b, const point<T, 3> &c)
        {
            const auto n = cross(b - a, c - a);
            const T nn = norm(n);
            if (nn > T(0))
            {
                // inside the prism of the triangle: the distance to its plane
                const bool in = (cross(b - a, q - a) * n) >= T(0) && (cross(c - b, q - b) * n) >= T(0) && (cross(a - c, q - c) * n) >= T(0);
                if (in)
                    return std::abs((q - a) * n) / nn;
            }
            return std::min({intersect_detail::distance_to_segment(q, a, b), intersect_detail::distance_to_segment(q, b, c),
                             intersect_detail::distance_to_segment(q, c, a)});
        }

        /// The mean lengths of S_u and S_v on the box {u1, u2, v1, v2} (sampled; 1 for one that vanishes).
        template <std::floating_point T>
        std::array<T, 2> parameter_scales(const Surface<T, 3> &s, const std::array<T, 4> &box)
        {
            constexpr int n = 8;
            T su{}, sv{};
            for (int i = 0; i <= n; ++i)
                for (int j = 0; j <= n; ++j)
                {
                    const T u = box[0] + (box[1] - box[0]) * T(i) / T(n), v = box[2] + (box[3] - box[2]) * T(j) / T(n);
                    su += norm(s.value(u, v, 1, 0));
                    sv += norm(s.value(u, v, 0, 1));
                }
            su /= T((n + 1) * (n + 1));
            sv /= T((n + 1) * (n + 1));
            if (!(su > T(0)) || !std::isfinite(su))
                su = sv > T(0) && std::isfinite(sv) ? sv : T(1);
            if (!(sv > T(0)) || !std::isfinite(sv))
                sv = su;
            return {su, sv};
        }

        /// The centre of the circle through a, b and c (not finite when they are aligned).
        template <std::floating_point T>
        point<T, 2> circumcentre(const point<T, 2> &a, const point<T, 2> &b, const point<T, 2> &c)
        {
            const T bx = b[0] - a[0], by = b[1] - a[1], cx = c[0] - a[0], cy = c[1] - a[1];
            const T d = T(2) * (bx * cy - by * cx);
            const T b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
            return {a[0] + (cy * b2 - by * c2) / d, a[1] + (bx * c2 - cx * b2) / d};
        }

        /// The smallest angle of the triangle (a, b, c) (0 when degenerate).
        template <std::floating_point T>
        T smallest_angle(const point<T, 2> &a, const point<T, 2> &b, const point<T, 2> &c)
        {
            auto at = [](const point<T, 2> &o, const point<T, 2> &p, const point<T, 2> &q) {
                const T x1 = p[0] - o[0], y1 = p[1] - o[1], x2 = q[0] - o[0], y2 = q[1] - o[1];
                return std::atan2(std::abs(x1 * y2 - y1 * x2), x1 * x2 + y1 * y2);
            };
            return std::min({at(a, b, c), at(b, c, a), at(c, a, b)});
        }

        struct TriangleHash
        {
            std::size_t operator()(const std::array<std::uint32_t, 3> &t) const noexcept
            {
                std::size_t h = t[0];
                h = h * 0x9E3779B97F4A7C15ull + t[1];
                return h * 0x9E3779B97F4A7C15ull + t[2];
            }
        };

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

    /// The triangles of a face.
    template <std::floating_point T>
    struct FaceMesh
    {
        std::vector<point<T, 2>> uv;                         ///< its points, in the parameter space of the face
        std::vector<std::uint32_t> vertices;                 ///< the vertex of the mesh at each point
        std::vector<std::array<std::uint32_t, 3>> triangles; ///< points, counter-clockwise in (u, v): along S_u x S_v
        bool complete{true};                                 ///< false: the refinement stopped at the limit of points
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

        /// The largest number of points of a face: beyond it, the face is left as it is (`FaceMesh::complete`).
        static constexpr std::size_t max_face_points = 200'000;

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

        /**
         * @brief The triangles of the face f: its boundary (`boundary`) constrained in a Delaunay
         * triangulation of its parameter space scaled by the mean lengths of S_u and S_v, refined
         * inside until each triangle is within the options:
         *
         * - deflection: the surface at its centroid and at the middles of its inner sides (those
         *   that are not on the boundary) is within the deflection of the triangle in 3D;
         * - angle: the normals of the surface at its corners turn by less than the angle;
         * - length: its inner sides are shorter than `max_edge`.
         *
         * A triangle out of these is given a point at its circumcentre in the scaled space, or at
         * its centroid when the circumcentre is outside the face or near its boundary (Ruppert's
         * encroachment). Triangles whose smallest angle is under 20° are given their circumcentre
         * too, when it is well inside. Triangles with two corners at the same vertex (along a pole)
         * are left out: they have no area. Interior points are new vertices of the mesh.
         */
        auto face(FaceId f) -> BuildResult<FaceMesh<T>>
        {
            auto loops = boundary(f);
            if (!loops)
                return std::unexpected(loops.error());
            const auto &srf = *m_.face(f).surface;
            std::array<T, 4> box{std::numeric_limits<T>::max(), std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max(),
                                 std::numeric_limits<T>::lowest()};
            for (const auto &loop : *loops)
                for (const auto &p : loop)
                    box = {std::min(box[0], p.uv[0]), std::max(box[1], p.uv[0]), std::min(box[2], p.uv[1]), std::max(box[3], p.uv[1])};
            if (!(box[0] < box[1]) || !(box[2] < box[3]))
                return detail::fail(BuildErrc::DegenerateWire, "mesh: the boundary of the face encloses no area", {f});
            const auto [a, b] = mesh_detail::parameter_scales(srf, box);
            auto scaled = [&](const point<T, 2> &uv) { return point<T, 2>{a * uv[0], b * uv[1]}; };
            auto unscaled = [&](const point<T, 2> &xy) { return point<T, 2>{xy[0] / a, xy[1] / b}; };
            using CDT = ConstrainedDelaunay<T>;
            CDT cdt{{a * box[0], a * box[1], b * box[2], b * box[3]}};

            FaceMesh<T> r;
            std::vector<point<T, 3>> normal; // of the surface at each point
            auto add = [&](const point<T, 2> &uv, std::uint32_t vertex) {
                const auto k = cdt.insert(scaled(uv));
                if (k == r.uv.size())
                {
                    r.uv.push_back(uv);
                    r.vertices.push_back(vertex);
                    // (none where S_u x S_v vanishes: a pole, an apex)
                    const auto n = cross(srf.value(uv[0], uv[1], 1, 0), srf.value(uv[0], uv[1], 0, 1));
                    normal.push_back(norm(n) > T(1e-6) * a * b ? n : point<T, 3>{});
                }
                return k;
            };
            try
            {
                std::vector<std::vector<std::uint32_t>> ids;
                for (const auto &loop : *loops)
                {
                    auto &l = ids.emplace_back();
                    for (const auto &p : loop)
                        l.push_back(add(p.uv, p.vertex));
                }
                for (const auto &l : ids)
                    for (std::size_t k = 0; k < l.size(); ++k)
                        cdt.constrain(l[k], l[(k + 1) % l.size()]);
            }
            catch (const std::runtime_error &e)
            {
                return detail::fail(BuildErrc::NonGeneric, std::string{"mesh: the boundary of the face does not triangulate: "} + e.what(), {f});
            }

            // the refinement
            const T tiny = T(1e-9) * std::max(a * (box[1] - box[0]), b * (box[3] - box[2]));
            const T quality = T(20) * std::numbers::pi_v<T> / T(180);
            std::unordered_set<std::array<std::uint32_t, 3>, mesh_detail::TriangleHash> good;
            enum class Need { nothing, geometry, quality };
            auto need = [&](const std::array<std::uint32_t, 3> &t) {
                // along a pole: no area, left out of the mesh (its neighbours cover the surface)
                if (r.vertices[t[0]] == r.vertices[t[1]] || r.vertices[t[1]] == r.vertices[t[2]] || r.vertices[t[2]] == r.vertices[t[0]])
                    return Need::nothing;
                const auto &p0 = vertices_[r.vertices[t[0]]], &p1 = vertices_[r.vertices[t[1]]], &p2 = vertices_[r.vertices[t[2]]];
                const auto centroid = (r.uv[t[0]] + r.uv[t[1]] + r.uv[t[2]]) * (T(1) / T(3));
                // folded against the surface (near a pole, where (u, v) and 3D disagree most)
                if (!(cross(p1 - p0, p2 - p0) * cross(srf.value(centroid[0], centroid[1], 1, 0), srf.value(centroid[0], centroid[1], 0, 1)) > T(0)))
                    return Need::geometry;
                // smaller than the deflection: within it
                if (std::max({norm(p1 - p0), norm(p2 - p1), norm(p0 - p2)}) < o_.deflection)
                    return Need::nothing;
                for (int i = 0; i < 3; ++i)
                {
                    const auto u = t[(i + 1) % 3], w = t[(i + 2) % 3];
                    if (cdt.constrained(u, w))
                        continue;
                    if (o_.max_edge > T(0) && norm(vertices_[r.vertices[u]] - vertices_[r.vertices[w]]) > o_.max_edge)
                        return Need::geometry;
                    const auto mid = (r.uv[u] + r.uv[w]) * T(0.5);
                    if (mesh_detail::distance_to_triangle(srf.value(mid[0], mid[1]), p0, p1, p2) > o_.deflection)
                        return Need::geometry;
                }
                if (mesh_detail::distance_to_triangle(srf.value(centroid[0], centroid[1]), p0, p1, p2) > o_.deflection)
                    return Need::geometry;
                for (int i = 0; i < 3; ++i)
                    if (mesh_detail::angle_between(normal[t[i]], normal[t[(i + 1) % 3]]) > o_.angle)
                        return Need::geometry;
                if (mesh_detail::smallest_angle(scaled(r.uv[t[0]]), scaled(r.uv[t[1]]), scaled(r.uv[t[2]])) < quality)
                    return Need::quality;
                return Need::nothing;
            };
            for (;;)
            {
                struct Candidate
                {
                    std::array<std::uint32_t, 3> of; // (its triangle: the point is skipped once it is gone)
                    point<T, 2> xy;
                    bool forced; // (a centroid: inside, even near the boundary)
                };
                std::vector<Candidate> todo;
                bool unmet = false;
                for (const auto &t : cdt.domain())
                {
                    if (good.contains(t))
                        continue;
                    const auto a0 = cdt.point_at(t[0]), a1 = cdt.point_at(t[1]), a2 = cdt.point_at(t[2]);
                    const T longest = std::sqrt(std::max({dist2(a0, a1), dist2(a1, a2), dist2(a2, a0)}));
                    const auto n = longest > tiny ? need(t) : Need::nothing;
                    if (n == Need::nothing)
                    {
                        good.insert(t);
                        continue;
                    }
                    const auto c = mesh_detail::circumcentre(a0, a1, a2);
                    if (cdt.place(c) == CDT::Place::inside)
                        todo.push_back({t, c, false});
                    else if (n == Need::geometry)
                        todo.push_back({t, (a0 + a1 + a2) * (T(1) / T(3)), true});
                    unmet = unmet || n == Need::geometry;
                }
                std::size_t inserted{};
                for (const auto &[of, xy, forced] : todo)
                {
                    if (cdt.size() >= max_face_points)
                        break;
                    if (!cdt.has_triangle(of[0], of[1], of[2]))
                        continue;
                    const auto where = cdt.place(xy);
                    if (where == CDT::Place::outside || (where == CDT::Place::encroaching && !forced))
                        continue;
                    const auto before = r.uv.size();
                    const auto uv = unscaled(xy);
                    add(uv, static_cast<std::uint32_t>(vertices_.size()));
                    if (r.uv.size() > before)
                    {
                        vertices_.push_back(srf.value(uv[0], uv[1]));
                        ++inserted;
                    }
                }
                if (inserted == 0)
                {
                    r.complete = !unmet;
                    break;
                }
            }

            cdt.remove_outside();
            for (const auto &t : cdt.triangles())
                if (r.vertices[t[0]] != r.vertices[t[1]] && r.vertices[t[1]] != r.vertices[t[2]] && r.vertices[t[2]] != r.vertices[t[0]])
                    r.triangles.push_back(t);
            return r;
        }

    private:
        const Model<T> &m_;
        MeshOptions<T> o_;
        static constexpr std::uint32_t none = std::numeric_limits<std::uint32_t>::max();
        std::vector<point<T, 3>> vertices_;
        std::vector<std::uint32_t> of_vertex_;           ///< by model vertex index
        std::vector<std::optional<EdgePoints>> edges_;   ///< by model edge index

        static T dist2(const point<T, 2> &p, const point<T, 2> &q)
        {
            const T dx = p[0] - q[0], dy = p[1] - q[1];
            return dx * dx + dy * dy;
        }

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
