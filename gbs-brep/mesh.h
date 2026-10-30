#pragma once
/**
 * @file mesh.h
 * @brief Triangle meshes of the faces, shells and solids of the native BREP (stage 3,
 * `docs/sources/design/brep_palier3_mesh.md`).
 *
 * The options, the discretization of the edges, the boundaries of the faces in their parameter
 * space, the triangulation of each face, and the mesh of a shape.
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
 * - The mesh of a shape gathers its faces: the vertices along the edges are shared, the
 *   triangles of a face used `Reversed` in its shell are turned over, so that the mesh of a
 *   solid faces outwards; each corner has the normal of the surface there, oriented likewise.
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
#include <tuple>
#include <unordered_map>
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
        std::size_t max_points{200'000}; ///< most points of a face: beyond, it is left as it is (`FaceMesh::complete`)
        T min_angle{T(20) * std::numbers::pi_v<T> / T(180)}; ///< triangles under it refined, in (u, v) scaled (radians; 0: none)
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

        /**
         * The scales of u and v for the triangulation: the mean length of S_u (S_v) over the
         * spacing the options ask for along u (v). That spacing is the one the chord of a curve of
         * the mean normal curvature along u needs to be within the deflection (sqrt(8 deflection /
         * k)) and to turn by less than the angle (angle / k), at most the size of the face: a
         * cylinder gets triangles long along its generatrices, a sphere isotropic ones. The
         * curvature along u is the turn of the normals between samples over the length between
         * them, both summed on a grid of 9 x 9 (first derivatives only). Without `deflection` (0), the mean lengths only (as isotropic in 3D).
         * Returns {scale of u, scale of v, mean |S_u|, mean |S_v|}.
         */
        template <std::floating_point T>
        std::array<T, 4> parameter_scales(const Surface<T, 3> &s, const std::array<T, 4> &box, T deflection, T angle)
        {
            constexpr int n = 8;
            std::array<std::array<point<T, 3>, n + 1>, n + 1> p, nrm;
            T su{}, sv{};
            point<T, 3> lo{std::numeric_limits<T>::max(), std::numeric_limits<T>::max(), std::numeric_limits<T>::max()}, hi{-lo[0], -lo[1], -lo[2]};
            for (int i = 0; i <= n; ++i)
                for (int j = 0; j <= n; ++j)
                {
                    const T u = box[0] + (box[1] - box[0]) * T(i) / T(n), v = box[2] + (box[3] - box[2]) * T(j) / T(n);
                    const auto du = s.value(u, v, 1, 0), dv = s.value(u, v, 0, 1);
                    su += norm(du);
                    sv += norm(dv);
                    p[i][j] = s.value(u, v);
                    nrm[i][j] = cross(du, dv);
                    for (int k = 0; k < 3; ++k)
                    {
                        lo[k] = std::min(lo[k], p[i][j][k]);
                        hi[k] = std::max(hi[k], p[i][j][k]);
                    }
                }
            su /= T((n + 1) * (n + 1));
            sv /= T((n + 1) * (n + 1));
            if (!(su > T(0)) || !std::isfinite(su))
                su = sv > T(0) && std::isfinite(sv) ? sv : T(1);
            if (!(sv > T(0)) || !std::isfinite(sv))
                sv = su;
            if (!(deflection > T(0)))
                return {su, sv, su, sv};
            // mean normal curvatures along u and v: the total turn of the normals over the total
            // length of the chords (a pole, where the chords vanish and the normals turn, weighs
            // nothing)
            T turn_u{}, len_u{}, turn_v{}, len_v{};
            for (int i = 0; i <= n; ++i)
                for (int j = 0; j <= n; ++j)
                {
                    auto add = [&](const point<T, 3> &n2, const point<T, 3> &p2, T &turn, T &len) {
                        len += norm(p2 - p[i][j]);
                        turn += angle_between(nrm[i][j], n2); // (0 where a normal vanishes)
                    };
                    if (i < n)
                        add(nrm[i + 1][j], p[i + 1][j], turn_u, len_u);
                    if (j < n)
                        add(nrm[i][j + 1], p[i][j + 1], turn_v, len_v);
                }
            const T ku = len_u > T(0) ? turn_u / len_u : T(0), kv = len_v > T(0) ? turn_v / len_v : T(0);
            const T size = std::max(norm(hi - lo), deflection);
            auto spacing = [&](T k) {
                if (!(k > T(0)))
                    return size;
                return std::min({std::sqrt(T(8) * deflection / k), angle / k, size});
            };
            return {su / spacing(ku), sv / spacing(kv), su, sv};
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

        /**
         * The parameter that does not move a point at a pole or an apex (S_u x S_v null there,
         * `pole`): the one whose derivative vanishes, against the mean lengths su and sv of S_u
         * and S_v on the face (0: u, 1: v); -1 elsewhere, or when both vanish.
         */
        template <std::floating_point T>
        std::int8_t free_parameter(const point<T, 3> &du, const point<T, 3> &dv, T su, T sv, bool pole)
        {
            if (!pole)
                return -1;
            const T ru = norm(du) / su, rv = norm(dv) / sv;
            if (!(std::max(ru, rv) > T(1e-6)))
                return -1;
            return ru < rv ? 0 : 1;
        }

        /**
         * The point of (u, v) with the barycentric weights w on the corners uv, where the
         * triangle is in 3D: a corner at a pole (`free[i]` the parameter that does not move it,
         * -1 elsewhere) is one point whatever that parameter, which it takes from the other
         * corners. The side from a pole to a corner follows the iso line of that corner through
         * the pole (a generatrix of a cone, a meridian of a sphere); the centroid of a triangle at
         * a pole lies between the iso lines of its other two corners.
         */
        template <std::floating_point T, std::size_t N>
        point<T, 2> combination(const std::array<point<T, 2>, N> &uv, const std::array<std::int8_t, N> &free, const std::array<T, N> &w)
        {
            point<T, 2> r{};
            for (std::size_t k = 0; k < 2; ++k)
            {
                T s{}, sw{}, all{};
                for (std::size_t i = 0; i < N; ++i)
                {
                    all += w[i] * uv[i][k];
                    if (free[i] != static_cast<std::int8_t>(k))
                    {
                        s += w[i] * uv[i][k];
                        sw += w[i];
                    }
                }
                r[k] = sw > T(0) ? s / sw : all;
            }
            return r;
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
     * @brief Parameters of the points of the curve c on [a, b]: its pieces (its Bézier segments,
     * `CurvePieces`), each halved while it may be farther than `deflection` from its chord (its
     * poles bound it: convex hull), its tangent may turn by more than `angle` (the turn of its
     * control polygon), or its chord is longer than `max_edge`; then consecutive pieces merged
     * while the same bounds on their poles together allow it, so that a curve of many knots (an
     * approximated intersection) gets no more points than its shape asks for. A curve which is not
     * a B-spline is sampled instead (its pieces not merged). Throws `std::invalid_argument` for a
     * deflection or an angle that is not positive.
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
        std::vector<CurvePiece<T, 3>> leaves;
        for (auto &first : pieces.pieces(a, b))
        {
            // depth first: the pieces come out in order
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
                    leaves.push_back(std::move(pc));
            }
        }
        // the pieces [i, j) together: their poles within the deflection of the chord, the turn of
        // their control polygons (junctions included) within the angle
        auto together = [&](std::size_t i, std::size_t j) {
            const auto &p = leaves[i].first, &q = leaves[j - 1].last;
            if (o.max_edge > T(0) && norm(q - p) > o.max_edge)
                return false;
            std::vector<point<T, 3>> legs;
            point<T, 3> prev = p;
            for (std::size_t k = i; k < j; ++k)
            {
                if (!leaves[k].bezier)
                    return false;
                for (const auto &h : leaves[k].poles)
                {
                    const auto x = intersect_detail::cartesian<T, 3>(h);
                    if (intersect_detail::distance_to_segment(x, p, q) > o.deflection)
                        return false;
                    if (norm(x - prev) > T(0))
                        legs.push_back(x - prev);
                    prev = x;
                }
            }
            T t{};
            for (std::size_t k = 0; k + 1 < legs.size(); ++k)
                if ((t += mesh_detail::angle_between(legs[k], legs[k + 1])) > o.angle)
                    return false;
            return true;
        };
        std::vector<T> r{a};
        // the run from i: its length doubled while together, then bisected (O(k log k) tests of
        // poles for a run of k pieces; any run kept passes `together`)
        for (std::size_t i = 0; i < leaves.size();)
        {
            const std::size_t most = leaves.size() - i; // (the longest run from i)
            std::size_t good = 1, len = 2;
            while (len <= most && together(i, i + len))
            {
                good = len;
                len *= 2;
            }
            for (std::size_t bad = std::min(len, most + 1); bad - good > 1;)
            {
                const std::size_t mid = good + (bad - good) / 2;
                (together(i, i + mid) ? good : bad) = mid;
            }
            r.push_back(leaves[i + good - 1].b);
            i += good;
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
        std::vector<std::int8_t> free;                       ///< at a pole or an apex, the parameter that does not move the point (0: u, 1: v); -1 elsewhere
        std::vector<std::array<std::uint32_t, 3>> triangles; ///< points, counter-clockwise in (u, v): along S_u x S_v
        bool complete{true};                                 ///< false: the refinement stopped at `MeshOptions::max_points`
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

        /**
         * The points of a non-degenerate edge: `curve_parameters` of its curve, then each segment
         * halved while, on a face of the edge, the surface at the quarters or the middle of the
         * chord of its pcurve is farther than the deflection from its chord in 3D (beyond the gap
         * between the pcurve and the edge, measured at the same parameters). The boundary of a face is
         * never refined by its triangulation: the pcurve chords must already follow the surface.
         */
        auto edge(EdgeId e) -> const EdgePoints &
        {
            if (auto *known = find(e))
                return *known;
            const auto &ed = m_.edge(e);
            EdgePoints p;
            p.t = on_faces(e, curve_parameters(*ed.curve, ed.u1, ed.u2, o_));
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
         *
         * A corner at a pole or an apex is one 3D point whatever its free parameter (u at the pole
         * of a sphere, at the apex of a cone): a triangle there is measured where it is in 3D, that
         * corner taking the free parameter of the other two (`mesh_detail::combination`), and its
         * centroid, when given, is moved between their iso lines.
         * A boundary that does not triangulate (crossing itself in (u, v)) gives `NonGeneric`.
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
            // (the triangulation throws on a boundary crossing itself, or a walk that does not end)
            const auto first = vertices_.size();
            try
            {
                return triangulate(*loops, srf, box);
            }
            catch (const std::exception &e)
            {
                vertices_.resize(first); // (the interior points of the face)
                return detail::fail(BuildErrc::NonGeneric, std::string{"mesh: the face does not triangulate: "} + e.what(), {f});
            }
        }

    private:
        /// `face` once its boundary is known: throws when the triangulation fails.
        FaceMesh<T> triangulate(const std::vector<std::vector<BoundaryPoint<T>>> &loops, const Surface<T, 3> &srf, const std::array<T, 4> &box)
        {
            const auto [a, b, su, sv] = mesh_detail::parameter_scales(srf, box, o_.deflection, o_.angle);
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
                    const auto du = srf.value(uv[0], uv[1], 1, 0), dv = srf.value(uv[0], uv[1], 0, 1);
                    const auto n = cross(du, dv);
                    const bool pole = !(norm(n) > T(1e-6) * su * sv);
                    normal.push_back(pole ? point<T, 3>{} : n);
                    r.free.push_back(mesh_detail::free_parameter(du, dv, su, sv, pole));
                }
                return k;
            };
            auto corners = [&](const std::array<std::uint32_t, 3> &t) { return std::array{r.uv[t[0]], r.uv[t[1]], r.uv[t[2]]}; };
            auto frees = [&](const std::array<std::uint32_t, 3> &t) { return std::array{r.free[t[0]], r.free[t[1]], r.free[t[2]]}; };
            for (const auto &loop : loops)
            {
                std::vector<std::uint32_t> l;
                for (const auto &p : loop)
                    l.push_back(add(p.uv, p.vertex));
                for (std::size_t k = 0; k < l.size(); ++k)
                    cdt.constrain(l[k], l[(k + 1) % l.size()]);
            }

            // the refinement
            const T tiny = T(1e-9) * std::max(a * (box[1] - box[0]), b * (box[3] - box[2]));
            // triangles found within the options (their corners from the smallest: one key for
            // each triangle); the surface at the middles of the inner sides (shared by two
            // triangles), a side from a pole along the iso line of its other end
            std::unordered_set<std::array<std::uint32_t, 3>, mesh_detail::TriangleHash> good;
            std::unordered_map<std::uint64_t, point<T, 3>> middles;
            auto middle = [&](std::uint32_t u, std::uint32_t w) -> const point<T, 3> & {
                const auto key = (std::uint64_t{std::min(u, w)} << 32) | std::max(u, w);
                auto it = middles.find(key);
                if (it == middles.end())
                {
                    const auto mid = mesh_detail::combination<T, 2>({r.uv[u], r.uv[w]}, {r.free[u], r.free[w]}, {T(0.5), T(0.5)});
                    it = middles.emplace(key, srf.value(mid[0], mid[1])).first;
                }
                return it->second;
            };
            // (folded: a defect of geometry given its centroid)
            enum class Need { nothing, folded, geometry, quality };
            auto need = [&](const std::array<std::uint32_t, 3> &t) {
                // along a pole: no area, left out of the mesh (its neighbours cover the surface)
                if (r.vertices[t[0]] == r.vertices[t[1]] || r.vertices[t[1]] == r.vertices[t[2]] || r.vertices[t[2]] == r.vertices[t[0]])
                    return Need::nothing;
                const auto &p0 = vertices_[r.vertices[t[0]]], &p1 = vertices_[r.vertices[t[1]]], &p2 = vertices_[r.vertices[t[2]]];
                constexpr T third = T(1) / T(3);
                const auto centroid = mesh_detail::combination<T, 3>(corners(t), frees(t), {third, third, third});
                // folded against the surface (near a pole, where (u, v) and 3D disagree most: a
                // corner at a pole is one 3D point, the triangle folded when the iso lines of its
                // other two corners are in the wrong order)
                if (!(cross(p1 - p0, p2 - p0) * cross(srf.value(centroid[0], centroid[1], 1, 0), srf.value(centroid[0], centroid[1], 0, 1)) > T(0)))
                    return Need::folded;
                for (int i = 0; i < 3; ++i)
                {
                    const auto u = t[(i + 1) % 3], w = t[(i + 2) % 3];
                    if (cdt.constrained(u, w))
                        continue;
                    if (o_.max_edge > T(0) && norm(vertices_[r.vertices[u]] - vertices_[r.vertices[w]]) > o_.max_edge)
                        return Need::geometry;
                    if (mesh_detail::distance_to_triangle(middle(u, w), p0, p1, p2) > o_.deflection)
                        return Need::geometry;
                }
                if (mesh_detail::distance_to_triangle(srf.value(centroid[0], centroid[1]), p0, p1, p2) > o_.deflection)
                    return Need::geometry;
                for (int i = 0; i < 3; ++i)
                    if (mesh_detail::angle_between(normal[t[i]], normal[t[(i + 1) % 3]]) > o_.angle)
                        return Need::geometry;
                if (o_.min_angle > T(0) && mesh_detail::smallest_angle(scaled(r.uv[t[0]]), scaled(r.uv[t[1]]), scaled(r.uv[t[2]])) < o_.min_angle)
                    return Need::quality;
                return Need::nothing;
            };
            // the point given to a triangle out of the options whose circumcentre is not inside:
            // its centroid, or, at a pole, the point of its median from the pole whose free
            // parameter stays between those of the other two corners (the centroid when the pole
            // corner is between them). The fan around the pole keeps the order of its iso lines: a
            // centroid pulled towards a pole corner beyond them (the end of a seam) would fold the
            // triangle it makes with the nearer one. (A folded triangle keeps its centroid: no
            // point inserted in it unfolds it, the cavity of the point does.)
            auto inner = [&](const std::array<std::uint32_t, 3> &t, const std::array<point<T, 2>, 3> &p) {
                for (std::size_t i = 0; i < 3; ++i)
                {
                    const auto k = r.free[t[i]];
                    const auto j = (i + 1) % 3, l = (i + 2) % 3;
                    if (k < 0 || r.free[t[j]] >= 0 || r.free[t[l]] >= 0)
                        continue;
                    const auto mid = (p[j] + p[l]) * T(0.5);
                    const auto x = static_cast<std::size_t>(k);
                    const T d = std::abs(p[i][x] - mid[x]), half = std::abs(p[j][x] - p[l][x]) / T(2);
                    return mid + (p[i] - mid) * std::min(T(1) / T(3), d > T(0) ? half / (T(2) * d) : T(1));
                }
                return (p[0] + p[1] + p[2]) * (T(1) / T(3));
            };
            for (;;)
            {
                struct Candidate
                {
                    std::array<std::uint32_t, 3> of; // (its triangle: the point is skipped once it is gone)
                    point<T, 2> xy;
                    bool forced; // (inside the triangle, even near the boundary)
                };
                std::vector<Candidate> todo;
                bool unmet = false;
                for (auto t : cdt.domain())
                {
                    std::ranges::rotate(t, std::ranges::min_element(t));
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
                    if (cdt.place(c, t[0]) == CDT::Place::inside)
                        todo.push_back({t, c, false});
                    else if (n == Need::folded)
                        todo.push_back({t, (a0 + a1 + a2) * (T(1) / T(3)), true});
                    else if (n == Need::geometry)
                        todo.push_back({t, inner(t, {a0, a1, a2}), true});
                    unmet = unmet || n != Need::quality;
                }
                std::size_t inserted{};
                for (const auto &[of, xy, forced] : todo)
                {
                    if (cdt.size() >= o_.max_points)
                        break;
                    if (!cdt.has_triangle(of[0], of[1], of[2]))
                        continue;
                    const auto where = cdt.place(xy, of[0]);
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

        /// The surfaces and pcurves of the co-edges of each edge (all the faces of the model), built once.
        std::vector<std::vector<std::pair<const Surface<T, 3> *, const Curve<T, 2> *>>> uses_;
        bool indexed_{};

        std::vector<T> on_faces(EdgeId e, std::vector<T> t)
        {
            if (!indexed_)
            {
                uses_.resize(m_.template capacity<ShapeType::Edge>());
                for (auto f : m_.template ids<FaceId>())
                    for (auto w : m_.face(f).wires)
                        for (const auto &ce : m_.wire(w).coedges)
                            if (ce.pcurve && ce.edge.index < uses_.size())
                                uses_[ce.edge.index].push_back({m_.face(f).surface.get(), ce.pcurve.get()});
                indexed_ = true;
            }
            if (e.index >= uses_.size() || uses_[e.index].empty())
                return t;
            const auto &c = *m_.edge(e).curve;
            // (the pcurves follow the surface within their approximation only: that gap, at the
            // ends and the middle of the segment, is not the chord's)
            // (the chord of the pcurve sampled at its quarters, its middle; the gap measured at
            // the ends and at the same parameters of the edge)
            constexpr std::array<T, 3> at{T(0.25), T(0.5), T(0.75)};
            auto off = [&](T a, T b) {
                const auto pa = c.value(a), pb = c.value(b);
                std::array<point<T, 3>, 3> pt;
                for (std::size_t i = 0; i < 3; ++i)
                    pt[i] = c.value(a + (b - a) * at[i]);
                for (const auto &[srf, pc] : uses_[e.index])
                {
                    const auto ua = pc->value(a), ub = pc->value(b);
                    T gap = std::max(norm(srf->value(ua[0], ua[1]) - pa), norm(srf->value(ub[0], ub[1]) - pb));
                    for (std::size_t i = 0; i < 3; ++i)
                    {
                        const auto ui = pc->value(a + (b - a) * at[i]);
                        gap = std::max(gap, norm(srf->value(ui[0], ui[1]) - pt[i]));
                    }
                    for (T s : at)
                    {
                        const auto uv = ua + (ub - ua) * s;
                        if (intersect_detail::distance_to_segment(srf->value(uv[0], uv[1]), pa, pb) > o_.deflection + gap)
                            return true;
                    }
                }
                return false;
            };
            constexpr int max_depth = 16;
            std::vector<T> r{t.front()};
            for (std::size_t k = 0; k + 1 < t.size(); ++k)
            {
                std::vector<std::tuple<T, T, int>> stack{{t[k], t[k + 1], 0}};
                while (!stack.empty())
                {
                    const auto [a, b, depth] = stack.back();
                    stack.pop_back();
                    if (depth < max_depth && off(a, b))
                    {
                        const T mid = (a + b) / T(2);
                        stack.push_back({mid, b, depth + 1});
                        stack.push_back({a, mid, depth + 1});
                    }
                    else
                        r.push_back(b);
                }
            }
            return r;
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

    /// A triangle mesh of a shape.
    template <std::floating_point T>
    struct Mesh
    {
        /// The triangles of a face: triangles[first, first + count).
        struct FaceRange
        {
            FaceId face;
            std::uint32_t first, count;
        };
        std::vector<point<T, 3>> vertices;
        std::vector<std::array<std::uint32_t, 3>> triangles; ///< counter-clockwise seen from outside (along the oriented normal)
        std::vector<FaceRange> faces;                        ///< in the order of the triangles
        std::vector<point<T, 2>> uv;                         ///< per corner (3 per triangle): (u, v) on its face
        std::vector<point<T, 3>> normals;                    ///< per corner: unit normal of the surface, oriented
        std::vector<FaceId> incomplete;                      ///< faces whose refinement stopped at `MeshOptions::max_points`
        std::vector<FaceId> folded;                          ///< faces with a triangle against the normal of the surface at a corner (a defect)
    };

    /**
     * @brief The mesh of the shape s: each of its faces triangulated (`MeshBuilder::face`), the
     * edges shared. A face is oriented by its first use in a shell of s (`Reversed`: its triangles
     * and normals turned over), `Forward` when it is in no shell. The normal at a corner is that of
     * the surface; where it has none (a pole, an apex), that of the surface just inside the
     * triangle as it is in 3D (`mesh_detail::combination`), or of the triangle; a triangle
     * against the normal of the surface at a corner keeps it and its face is reported `folded`
     * (a defect of the mesh or of the orientation of the face in its shell). The vertices
     * are those used by the triangles, in their order of creation. The first face that cannot be
     * meshed gives its error; options out of their range give InvalidTolerance (deflection,
     * max_edge) or InvalidOptions (angles).
     */
    template <std::floating_point T>
    [[nodiscard]] auto mesh(const Model<T> &m, const ShapeId &s, const MeshOptions<T> &o = {}) -> BuildResult<Mesh<T>>
    {
        auto finite_from = [](T x, T low) { return std::isfinite(x) && x >= low; };
        if (!finite_from(o.deflection, T(0)) || !finite_from(o.max_edge, T(0)))
            return detail::fail(BuildErrc::InvalidTolerance, "mesh: the deflection and max_edge must be finite and >= 0 (0: default)");
        if (!(o.angle > T(0)) || !std::isfinite(o.angle) || !finite_from(o.min_angle, T(0)))
            return detail::fail(BuildErrc::InvalidOptions, "mesh: the angle must be finite and > 0, min_angle finite and >= 0");
        MeshBuilder<T> mb{m, s, o};
        // the faces, oriented by their shells
        std::vector<FaceUse> uses;
        std::vector<std::uint8_t> seen(m.template capacity<ShapeType::Face>(), 0);
        for (auto sh : explore<ShellId>(m, s))
            for (const auto &fu : m.shell(sh).faces)
                if (!seen[fu.face.index])
                {
                    seen[fu.face.index] = 1;
                    uses.push_back(fu);
                }
        for (auto f : explore<FaceId>(m, s))
            if (!seen[f.index])
            {
                seen[f.index] = 1;
                uses.push_back({f, Orientation::Forward});
            }

        Mesh<T> r;
        std::vector<std::array<std::uint32_t, 3>> tris; // in the vertices of the builder
        for (const auto &fu : uses)
        {
            auto fm = mb.face(fu.face);
            if (!fm)
                return std::unexpected(fm.error());
            const auto &srf = *m.face(fu.face).surface;
            const bool reversed = fu.orient == Orientation::Reversed;
            r.faces.push_back({fu.face, static_cast<std::uint32_t>(tris.size()), static_cast<std::uint32_t>(fm->triangles.size())});
            if (!fm->complete)
                r.incomplete.push_back(fu.face);
            bool folded = false;
            for (auto t : fm->triangles)
            {
                if (reversed)
                    std::swap(t[1], t[2]);
                const auto &p0 = mb.vertices()[fm->vertices[t[0]]], &p1 = mb.vertices()[fm->vertices[t[1]]], &p2 = mb.vertices()[fm->vertices[t[2]]];
                const auto nt = cross(p1 - p0, p2 - p0);
                constexpr T third = T(1) / T(3);
                const auto centroid = mesh_detail::combination<T, 3>({fm->uv[t[0]], fm->uv[t[1]], fm->uv[t[2]]},
                                                                     {fm->free[t[0]], fm->free[t[1]], fm->free[t[2]]}, {third, third, third});
                // S_u x S_v at uv, oriented; none where it vanishes against max(|S_u|, |S_v|)^2
                auto normal_at = [&](const point<T, 2> &uv) -> std::optional<point<T, 3>> {
                    const auto su = srf.value(uv[0], uv[1], 1, 0), sv = srf.value(uv[0], uv[1], 0, 1);
                    const auto n = cross(su, sv);
                    const T scale = std::max(norm(su), norm(sv));
                    if (!(norm(n) > T(1e-6) * scale * scale))
                        return std::nullopt;
                    return reversed ? n * T(-1) : n;
                };
                for (auto k : t)
                {
                    const auto &uv = fm->uv[k];
                    // at a pole or an apex: the surface just inside the triangle (its limit on
                    // that side, on the iso line of its centroid), else the triangle's
                    auto n = normal_at(uv);
                    if (!n)
                        n = normal_at(mesh_detail::combination<T, 2>({uv, centroid}, {fm->free[k], -1}, {T(1) - T(1e-3), T(1e-3)}));
                    if (!n)
                        n = nt;
                    else if (!(*n * nt > T(0)))
                        folded = true;
                    const T l = norm(*n);
                    r.uv.push_back(uv);
                    r.normals.push_back(l > T(0) ? *n / l : *n);
                }
                tris.push_back({fm->vertices[t[0]], fm->vertices[t[1]], fm->vertices[t[2]]});
            }
            if (folded)
                r.folded.push_back(fu.face);
        }
        // the vertices used, in their order
        constexpr auto unused = std::numeric_limits<std::uint32_t>::max();
        std::vector<std::uint32_t> index(mb.vertices().size(), unused);
        for (const auto &t : tris)
            for (auto v : t)
                index[v] = 0;
        for (std::size_t v = 0; v < index.size(); ++v)
            if (index[v] != unused)
            {
                index[v] = static_cast<std::uint32_t>(r.vertices.size());
                r.vertices.push_back(mb.vertices()[v]);
            }
        r.triangles.reserve(tris.size());
        for (const auto &t : tris)
            r.triangles.push_back({index[t[0]], index[t[1]], index[t[2]]});
        return r;
    }

    /// The area of a mesh.
    template <std::floating_point T>
    [[nodiscard]] T area(const Mesh<T> &mesh)
    {
        T a{};
        for (const auto &t : mesh.triangles)
            a += norm(cross(mesh.vertices[t[1]] - mesh.vertices[t[0]], mesh.vertices[t[2]] - mesh.vertices[t[0]]));
        return a / T(2);
    }

    /// The volume enclosed by a closed mesh (positive when it faces outwards).
    template <std::floating_point T>
    [[nodiscard]] T volume(const Mesh<T> &mesh)
    {
        T v{};
        for (const auto &t : mesh.triangles)
            v += mesh.vertices[t[0]] * cross(mesh.vertices[t[1]], mesh.vertices[t[2]]);
        return v / T(6);
    }

    /**
     * The sides of the triangles that are not used once in each sense: none for a closed,
     * consistently oriented mesh; the boundary of an open one.
     */
    [[nodiscard]] inline std::size_t open_sides(const std::vector<std::array<std::uint32_t, 3>> &triangles)
    {
        std::unordered_map<std::uint64_t, int> sides; // directed: +1 each use, the reverse -1
        auto key = [](std::uint32_t a, std::uint32_t b) { return (std::uint64_t{a} << 32) | b; };
        for (const auto &t : triangles)
            for (int k = 0; k < 3; ++k)
            {
                const auto a = t[k], b = t[(k + 1) % 3];
                if (a < b)
                    ++sides[key(a, b)];
                else
                    --sides[key(b, a)];
            }
        std::size_t n{};
        for (const auto &[k, c] : sides)
            n += static_cast<std::size_t>(std::abs(c));
        return n;
    }
    template <std::floating_point T>
    [[nodiscard]] std::size_t open_sides(const Mesh<T> &mesh)
    {
        return open_sides(mesh.triangles);
    }

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
