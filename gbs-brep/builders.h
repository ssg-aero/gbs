#pragma once

/**
 * @file builders.h
 * @brief Native BREP core, stage 1: builders of vertices, edges and wires.
 *
 * Design: docs/sources/design/brep_core.md, section 6.0; architecture note
 * docs/sources/design/brep_pr03_builders_wire.md.
 *
 * Every builder returns a `BuildResult<Id>` (`std::expected<Id, BuildError>`):
 * an expected failure (curve off a vertex, edges that do not chain…) is a
 * value, not an exception. A builder that fails leaves the model unchanged.
 * `unwrap()` turns an error into a `BRepError` exception for callers that
 * prefer exceptions.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <expected>
#include <memory>
#include <numeric>
#include <span>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <gbs-brep/model.h>
#include <gbs-brep/closure.h>
#include <gbs/bscbuild.h>
#include <gbs/curveonsurface.h>

namespace gbs::brep
{
    // =========================================================================
    // Result type
    // =========================================================================

    /// Reason of a builder failure.
    enum class BuildErrc : std::uint8_t
    {
        InvalidId,            ///< dead or invalid identifier in the input
        InvalidTolerance,     ///< tolerance not finite or not > 0
        InvalidPoint,         ///< point with a non finite coordinate
        NullCurve,            ///< no curve given
        UnboundedCurve,       ///< infinite parameter range (e.g. a Line) without explicit bounds
        InvalidBounds,        ///< u1 >= u2, or bounds outside the curve's range
        DegenerateEdge,       ///< the edge would be shorter than the tolerance
        VertexOffCurve,       ///< a given vertex is farther from the curve end than allowed
        EmptyInput,           ///< no edge given
        DuplicateEdge,        ///< the same edge given twice
        DegenerateEdgeInWire, ///< degenerate edges belong to face wires, built by face builders
        Branching,            ///< more than two edges meet at a vertex
        Disconnected,         ///< the edges do not form a single chain
        NullSurface,          ///< no surface given
        UnboundedSurface,     ///< infinite or empty parametric rectangle
        DegenerateSurface,    ///< the whole surface collapses (all four sides degenerate)
    };

    [[nodiscard]] inline constexpr const char *to_string(BuildErrc c) noexcept
    {
        switch (c)
        {
        case BuildErrc::InvalidId: return "invalid id";
        case BuildErrc::InvalidTolerance: return "invalid tolerance";
        case BuildErrc::InvalidPoint: return "invalid point";
        case BuildErrc::NullCurve: return "null curve";
        case BuildErrc::UnboundedCurve: return "unbounded curve";
        case BuildErrc::InvalidBounds: return "invalid bounds";
        case BuildErrc::DegenerateEdge: return "degenerate edge";
        case BuildErrc::VertexOffCurve: return "vertex off curve";
        case BuildErrc::EmptyInput: return "empty input";
        case BuildErrc::DuplicateEdge: return "duplicate edge";
        case BuildErrc::DegenerateEdgeInWire: return "degenerate edge in wire";
        case BuildErrc::Branching: return "branching";
        case BuildErrc::Disconnected: return "disconnected";
        case BuildErrc::NullSurface: return "null surface";
        case BuildErrc::UnboundedSurface: return "unbounded surface";
        case BuildErrc::DegenerateSurface: return "degenerate surface";
        }
        std::unreachable();
    }

    /// Why a builder failed, with the entities involved when there are some.
    struct BuildError
    {
        BuildErrc code;
        std::string message;
        std::vector<ShapeId> shapes;
    };

    template <typename V>
    using BuildResult = std::expected<V, BuildError>;

    /// Value of a successful result; throws `BRepError` with the error message otherwise.
    template <typename V>
    [[nodiscard]] auto unwrap(BuildResult<V> r) -> V
    {
        if (!r)
            throw BRepError(std::string(to_string(r.error().code)) + ": " + r.error().message);
        return std::move(*r);
    }

    namespace detail
    {
        inline auto fail(BuildErrc c, std::string msg, std::vector<ShapeId> shapes = {}) -> std::unexpected<BuildError>
        {
            return std::unexpected(BuildError{c, std::move(msg), std::move(shapes)});
        }

        template <std::floating_point T>
        [[nodiscard]] bool valid_tolerance(T tol) noexcept
        {
            return std::isfinite(tol) && tol > T(0);
        }

        template <std::floating_point T, std::size_t dim>
        [[nodiscard]] bool finite(const point<T, dim> &p) noexcept
        {
            return std::ranges::all_of(p, [](T x) { return std::isfinite(x); });
        }

        /// Parameters beyond this magnitude denote an unbounded curve (Line::bounds() returns the type limits).
        template <std::floating_point T>
        inline constexpr T unbounded_parameter = std::numeric_limits<T>::max() / T(16);

        /// Endpoint of an edge's curve at its v1 (first = true) or v2 side.
        template <std::floating_point T>
        [[nodiscard]] auto edge_end_point(const Edge<T> &e, bool first) -> point<T, 3>
        {
            return e.curve->value(first ? e.u1 : e.u2);
        }

        /// Raises the tolerance of the vertices of `e` so that each contains its curve end and is >= e.tol.
        template <std::floating_point T>
        void fit_vertex_tolerances(Model<T> &m, const Edge<T> &e)
        {
            if (e.degenerate || !e.curve)
            {
                auto &v = m.vertex(e.v1);
                v.tol = std::max(v.tol, e.tol);
                return;
            }
            for (bool first : {true, false})
            {
                auto &v = m.vertex(first ? e.v1 : e.v2);
                v.tol = std::max({v.tol, e.tol, distance(edge_end_point(e, first), v.pnt)});
            }
        }
    } // namespace detail

    // =========================================================================
    // Vertex
    // =========================================================================

    template <std::floating_point T>
    [[nodiscard]] auto make_vertex(Model<T> &m, const std::type_identity_t<point<T, 3>> &p, std::type_identity_t<T> tol = brep_default_tolerance<T>) -> BuildResult<VertexId>
    {
        if (!detail::valid_tolerance(tol))
            return detail::fail(BuildErrc::InvalidTolerance, "vertex tolerance must be finite and > 0");
        if (!detail::finite(p))
            return detail::fail(BuildErrc::InvalidPoint, "vertex point must be finite");
        return m.add(Vertex<T>{p, tol});
    }

    // =========================================================================
    // Edges
    // =========================================================================

    /**
     * @brief Edge on `curve` between `u1` and `u2`, with new vertices at both
     * ends (a single one if the ends are within `tol`: closed edge).
     */
    template <std::floating_point T>
    [[nodiscard]] auto make_edge(Model<T> &m, std::type_identity_t<std::shared_ptr<Curve<T, 3>>> curve, std::type_identity_t<T> u1, std::type_identity_t<T> u2, std::type_identity_t<T> tol = brep_default_tolerance<T>) -> BuildResult<EdgeId>
    {
        if (!curve)
            return detail::fail(BuildErrc::NullCurve, "make_edge needs a curve");
        if (!detail::valid_tolerance(tol))
            return detail::fail(BuildErrc::InvalidTolerance, "edge tolerance must be finite and > 0");
        if (!std::isfinite(u1) || !std::isfinite(u2) ||
            std::abs(u1) >= detail::unbounded_parameter<T> || std::abs(u2) >= detail::unbounded_parameter<T>)
            return detail::fail(BuildErrc::UnboundedCurve, "edge bounds must be finite; give explicit bounds for an unbounded curve");
        auto [c1, c2] = curve->bounds();
        if (!(u1 < u2) || u1 < c1 - knot_eps<T> || u2 > c2 + knot_eps<T>)
            return detail::fail(BuildErrc::InvalidBounds, "edge bounds must satisfy curve.u1 <= u1 < u2 <= curve.u2");

        const auto p1 = curve->value(u1);
        const auto p2 = curve->value(u2);
        if (!detail::finite(p1) || !detail::finite(p2))
            return detail::fail(BuildErrc::InvalidPoint, "curve evaluates to a non finite point at the edge bounds");

        // An edge whose curve stays within tol of its start is degenerate: reject.
        bool collapsed = true;
        for (T s : {T(0.25), T(0.5), T(0.75), T(1)})
            if (distance(curve->value(u1 + s * (u2 - u1)), p1) > tol)
                collapsed = false;
        if (collapsed)
            return detail::fail(BuildErrc::DegenerateEdge, "edge curve stays within tolerance of its start point");

        const T d = distance(p1, p2);
        const auto v1 = m.add(Vertex<T>{p1, tol}); // a closed edge (d <= tol) gets a single vertex
        const auto v2 = d <= tol ? v1 : m.add(Vertex<T>{p2, tol});
        return m.add(Edge<T>{std::move(curve), u1, u2, v1, v2, tol});
    }

    /// Edge on the whole parameter range of a bounded `curve`.
    template <std::floating_point T>
    [[nodiscard]] auto make_edge(Model<T> &m, std::type_identity_t<std::shared_ptr<Curve<T, 3>>> curve, std::type_identity_t<T> tol = brep_default_tolerance<T>) -> BuildResult<EdgeId>
    {
        if (!curve)
            return detail::fail(BuildErrc::NullCurve, "make_edge needs a curve");
        auto [u1, u2] = curve->bounds();
        return make_edge(m, std::move(curve), u1, u2, tol);
    }

    /**
     * @brief Edge on `curve` between `u1` and `u2` ending on existing vertices.
     * Each vertex must lie within `max(tol, vertex.tol)` of the curve end; its
     * tolerance is then raised to contain the curve end and to stay >= `tol`.
     * `v1 == v2` builds a closed edge.
     */
    template <std::floating_point T>
    [[nodiscard]] auto make_edge(Model<T> &m, std::type_identity_t<std::shared_ptr<Curve<T, 3>>> curve, std::type_identity_t<T> u1, std::type_identity_t<T> u2, VertexId v1, VertexId v2,
                                 std::type_identity_t<T> tol = brep_default_tolerance<T>) -> BuildResult<EdgeId>
    {
        if (!m.alive(v1) || !m.alive(v2))
            return detail::fail(BuildErrc::InvalidId, "make_edge: dead or invalid vertex", {v1, v2});
        if (!curve)
            return detail::fail(BuildErrc::NullCurve, "make_edge needs a curve");
        if (!detail::valid_tolerance(tol))
            return detail::fail(BuildErrc::InvalidTolerance, "edge tolerance must be finite and > 0");
        if (!std::isfinite(u1) || !std::isfinite(u2) ||
            std::abs(u1) >= detail::unbounded_parameter<T> || std::abs(u2) >= detail::unbounded_parameter<T>)
            return detail::fail(BuildErrc::UnboundedCurve, "edge bounds must be finite");
        auto [c1, c2] = curve->bounds();
        if (!(u1 < u2) || u1 < c1 - knot_eps<T> || u2 > c2 + knot_eps<T>)
            return detail::fail(BuildErrc::InvalidBounds, "edge bounds must satisfy curve.u1 <= u1 < u2 <= curve.u2");

        const auto p1 = curve->value(u1);
        const auto p2 = curve->value(u2);
        const T d1 = distance(p1, m.vertex(v1).pnt);
        const T d2 = distance(p2, m.vertex(v2).pnt);
        if (d1 > std::max(tol, m.vertex(v1).tol))
            return detail::fail(BuildErrc::VertexOffCurve, "curve start is " + std::to_string(d1) + " away from v1", {v1});
        if (d2 > std::max(tol, m.vertex(v2).tol))
            return detail::fail(BuildErrc::VertexOffCurve, "curve end is " + std::to_string(d2) + " away from v2", {v2});

        Edge<T> e{std::move(curve), u1, u2, v1, v2, tol};
        detail::fit_vertex_tolerances(m, e);
        return m.add(std::move(e));
    }

    /// Straight edge between two points (degree-1 B-spline parametrized by arc length).
    template <std::floating_point T>
    [[nodiscard]] auto make_edge(Model<T> &m, const std::type_identity_t<point<T, 3>> &p1, const std::type_identity_t<point<T, 3>> &p2, std::type_identity_t<T> tol = brep_default_tolerance<T>) -> BuildResult<EdgeId>
    {
        if (!detail::finite(p1) || !detail::finite(p2))
            return detail::fail(BuildErrc::InvalidPoint, "segment ends must be finite");
        if (!detail::valid_tolerance(tol))
            return detail::fail(BuildErrc::InvalidTolerance, "edge tolerance must be finite and > 0");
        if (distance(p1, p2) <= tol)
            return detail::fail(BuildErrc::DegenerateEdge, "segment shorter than the tolerance");
        auto crv = std::make_shared<BSCurve<T, 3>>(build_segment(p1, p2));
        auto [u1, u2] = crv->bounds();
        return make_edge(m, std::move(crv), u1, u2, tol);
    }

    /// Straight edge between two existing vertices.
    template <std::floating_point T>
    [[nodiscard]] auto make_edge(Model<T> &m, VertexId v1, VertexId v2, std::type_identity_t<T> tol = brep_default_tolerance<T>) -> BuildResult<EdgeId>
    {
        if (!m.alive(v1) || !m.alive(v2))
            return detail::fail(BuildErrc::InvalidId, "make_edge: dead or invalid vertex", {v1, v2});
        const auto p1 = m.vertex(v1).pnt;
        const auto p2 = m.vertex(v2).pnt;
        if (v1 == v2 || distance(p1, p2) <= std::max(tol, m.vertex(v1).tol + m.vertex(v2).tol))
            return detail::fail(BuildErrc::DegenerateEdge, "segment between confused vertices", {v1, v2});
        auto crv = std::make_shared<BSCurve<T, 3>>(build_segment(p1, p2));
        auto [u1, u2] = crv->bounds();
        return make_edge(m, std::move(crv), u1, u2, v1, v2, tol);
    }

    /// Degenerate edge reduced to the vertex `v` (pole of a sphere, apex of a cone), parameter range [u1, u2].
    template <std::floating_point T>
    [[nodiscard]] auto make_degenerate_edge(Model<T> &m, VertexId v, std::type_identity_t<T> u1 = T(0), std::type_identity_t<T> u2 = T(1)) -> BuildResult<EdgeId>
    {
        if (!m.alive(v))
            return detail::fail(BuildErrc::InvalidId, "make_degenerate_edge: dead or invalid vertex", {v});
        if (!std::isfinite(u1) || !std::isfinite(u2) || !(u1 < u2))
            return detail::fail(BuildErrc::InvalidBounds, "degenerate edge range must satisfy u1 < u2");
        return m.add(Edge<T>{nullptr, u1, u2, v, v, m.vertex(v).tol, true, true});
    }

    // =========================================================================
    // Wire
    // =========================================================================

    /**
     * @brief Chains `edges`, given in any order and any sense, into a wire.
     *
     * 1. Vertices of the edges closer than `max(tol, tol_a + tol_b)` are merged
     *    (the one with the smallest index survives, its point is not moved,
     *    its tolerance grows to contain the merged curve ends). Every edge of
     *    the model pointing to a merged vertex is redirected to the survivor;
     *    merged vertices are erased (tombstones).
     * 2. The edges must then form one simple chain: no vertex shared by more
     *    than two edges, one connected component. The wire is closed iff the
     *    chain is a cycle.
     * 3. The first edge of the input keeps the Forward sense; for a cycle the
     *    wire starts with it.
     *
     * Edges are never modified geometrically and no pcurve is set (free wire).
     * On failure the model is unchanged.
     */
    template <std::floating_point T>
    [[nodiscard]] auto make_wire(Model<T> &m, std::span<const EdgeId> edges, std::type_identity_t<T> tol = brep_default_tolerance<T>) -> BuildResult<WireId>
    {
        if (edges.empty())
            return detail::fail(BuildErrc::EmptyInput, "make_wire needs at least one edge");
        if (!detail::valid_tolerance(tol))
            return detail::fail(BuildErrc::InvalidTolerance, "wire tolerance must be finite and > 0");
        for (std::size_t i{}; i < edges.size(); ++i)
        {
            if (!m.alive(edges[i]))
                return detail::fail(BuildErrc::InvalidId, "make_wire: dead or invalid edge", {edges[i]});
            if (m.edge(edges[i]).degenerate)
                return detail::fail(BuildErrc::DegenerateEdgeInWire, "degenerate edges are inserted by face builders", {edges[i]});
            for (std::size_t j{}; j < i; ++j)
                if (edges[j] == edges[i])
                    return detail::fail(BuildErrc::DuplicateEdge, "edge given twice", {edges[i]});
        }

        // ---- 1. vertex classes (union-find on the distinct end vertices) ----
        std::vector<VertexId> verts;
        for (auto eid : edges)
            for (auto v : {m.edge(eid).v1, m.edge(eid).v2})
                if (std::ranges::find(verts, v) == verts.end())
                    verts.push_back(v);
        std::ranges::sort(verts);

        std::vector<std::size_t> parent(verts.size());
        std::iota(parent.begin(), parent.end(), std::size_t{0});
        auto find = [&](std::size_t i) {
            while (parent[i] != i)
                i = parent[i] = parent[parent[i]];
            return i;
        };
        for (std::size_t i{}; i < verts.size(); ++i)
            for (std::size_t j{i + 1}; j < verts.size(); ++j)
            {
                const auto &a = m.vertex(verts[i]);
                const auto &b = m.vertex(verts[j]);
                if (distance(a.pnt, b.pnt) <= std::max(tol, a.tol + b.tol))
                {
                    auto ri = find(i), rj = find(j);
                    if (ri != rj)
                        parent[std::max(ri, rj)] = std::min(ri, rj); // smallest index survives (verts is sorted)
                }
            }
        auto rep = [&](VertexId v) {
            auto i = static_cast<std::size_t>(std::ranges::lower_bound(verts, v) - verts.begin());
            return verts[find(i)];
        };

        // ---- 2. chain on the representatives (no mutation yet) -------------
        std::unordered_map<VertexId, std::vector<std::size_t>> incident;
        for (std::size_t i{}; i < edges.size(); ++i)
        {
            const auto &e = m.edge(edges[i]);
            incident[rep(e.v1)].push_back(i);
            incident[rep(e.v2)].push_back(i); // a closed edge is incident twice to its vertex
        }
        std::vector<VertexId> ends;
        for (const auto &[v, inc] : incident)
        {
            if (inc.size() > 2)
                return detail::fail(BuildErrc::Branching, "more than two edges meet at a vertex", {v});
            if (inc.size() == 1)
                ends.push_back(v);
        }
        if (ends.size() != 0 && ends.size() != 2)
            return detail::fail(BuildErrc::Disconnected, "edges do not form a single chain");
        const bool closed = ends.empty();

        std::vector<CoEdge<T>> chain;
        chain.reserve(edges.size());
        std::vector<std::uint8_t> used(edges.size(), 0);
        auto step = [&](std::size_t i, VertexId from) {
            const auto &e = m.edge(edges[i]);
            used[i] = 1;
            const bool forward = rep(e.v1) == from;
            chain.push_back(CoEdge<T>{edges[i], forward ? Orientation::Forward : Orientation::Reversed, nullptr});
            return forward ? rep(e.v2) : rep(e.v1);
        };
        VertexId cur;
        if (closed)
            cur = step(0, rep(m.edge(edges[0]).v1)); // the first edge leads, Forward
        else
            cur = std::ranges::min(ends);
        while (chain.size() < edges.size())
        {
            auto it = std::ranges::find_if(incident[cur], [&](std::size_t i) { return !used[i]; });
            if (it == incident[cur].end())
                break;
            cur = step(*it, cur);
        }
        if (chain.size() != edges.size())
            return detail::fail(BuildErrc::Disconnected, "edges do not form a single chain");

        // the first input edge keeps its sense: reverse an open chain if needed
        auto first = std::ranges::find_if(chain, [&](const CoEdge<T> &ce) { return ce.edge == edges[0]; });
        if (first->orient == Orientation::Reversed)
        {
            std::ranges::reverse(chain);
            for (auto &ce : chain)
                ce.orient = reverse(ce.orient);
        }

        // ---- 3. apply the vertex merge to the model -----------------------
        std::vector<std::pair<VertexId, VertexId>> merged; // absorbed -> survivor
        for (auto v : verts)
            if (auto r = rep(v); r != v)
                merged.emplace_back(v, r);
        if (!merged.empty())
        {
            auto redirect = [&](VertexId v) {
                for (const auto &[a, s] : merged)
                    if (a == v)
                        return s;
                return v;
            };
            for (auto eid : m.template ids<EdgeId>())
            {
                auto &e = m.edge(eid);
                const auto n1 = redirect(e.v1), n2 = redirect(e.v2);
                if (n1 != e.v1 || n2 != e.v2)
                {
                    e.v1 = n1;
                    e.v2 = n2;
                    detail::fit_vertex_tolerances(m, e);
                }
            }
            for (const auto &[a, s] : merged)
            {
                auto &sv = m.vertex(s);
                sv.tol = std::max(sv.tol, m.vertex(a).tol); // keep the absorbed ball's contract
                m.erase(a);
            }
        }

        return m.add(Wire<T>{std::move(chain), closed});
    }

    template <std::floating_point T>
    [[nodiscard]] auto make_wire(Model<T> &m, const std::vector<EdgeId> &edges, std::type_identity_t<T> tol = brep_default_tolerance<T>) -> BuildResult<WireId>
    {
        return make_wire(m, std::span<const EdgeId>{edges}, tol);
    }

    // =========================================================================
    // Face with natural bounds
    // =========================================================================

    /**
     * @brief Face covering the whole parametric rectangle of `srf`.
     *
     * The outer wire has four co-edges, counter-clockwise in (u,v):
     * v = v1 (u increasing), u = u2 (v increasing), v = v2 (u decreasing),
     * u = u1 (v decreasing). Each pcurve is a degree-1 B-spline parametrized
     * like its edge; each 3D curve is the exact `CurveOnSurface` of the side's
     * pcurve. `surface_closure(srf, tol)` decides the topology:
     * - a degenerate side becomes a degenerate edge (pole, apex);
     * - a closed direction becomes a seam: one edge used twice, Forward with
     *   the pcurve on one side and Reversed with the pcurve on the other;
     * - corners closer than `tol` share a vertex.
     * The whole surface collapsing to a point or a curve is rejected.
     */
    template <std::floating_point T>
    [[nodiscard]] auto make_face(Model<T> &m, std::type_identity_t<std::shared_ptr<Surface<T, 3>>> srf,
                                 std::type_identity_t<T> tol = brep_default_tolerance<T>) -> BuildResult<FaceId>
    {
        if (!srf)
            return detail::fail(BuildErrc::NullSurface, "make_face needs a surface");
        if (!detail::valid_tolerance(tol))
            return detail::fail(BuildErrc::InvalidTolerance, "face tolerance must be finite and > 0");
        const auto [u1, u2, v1, v2] = srf->bounds();
        for (T x : {u1, u2, v1, v2})
            if (!std::isfinite(x) || std::abs(x) >= detail::unbounded_parameter<T>)
                return detail::fail(BuildErrc::UnboundedSurface, "surface bounds must be finite");
        if (!(u1 < u2) || !(v1 < v2))
            return detail::fail(BuildErrc::UnboundedSurface, "surface parametric rectangle is empty");

        const auto c = surface_closure(*srf, tol);
        if ((c.degenerate_u1 && c.degenerate_u2) && (c.degenerate_v1 && c.degenerate_v2))
            return detail::fail(BuildErrc::DegenerateSurface, "surface collapses to a point");

        // ---- corners, counter-clockwise: 0 (u1,v1) 1 (u2,v1) 2 (u2,v2) 3 (u1,v2)
        const std::array<point<T, 2>, 4> uv{{{u1, v1}, {u2, v1}, {u2, v2}, {u1, v2}}};
        std::array<point<T, 3>, 4> p;
        for (std::size_t i{}; i < 4; ++i)
        {
            p[i] = srf->value(uv[i][0], uv[i][1]);
            if (!detail::finite(p[i]))
                return detail::fail(BuildErrc::InvalidPoint, "surface evaluates to a non finite corner");
        }
        std::array<std::size_t, 4> cls{0, 1, 2, 3};
        for (std::size_t j{1}; j < 4; ++j)
            for (std::size_t i{}; i < j; ++i)
                if (cls[i] == i && distance(p[i], p[j]) <= tol)
                {
                    cls[j] = i;
                    break;
                }

        // ---- sides: (start corner, end corner) in the edge's increasing-parameter sense
        struct Side
        {
            std::size_t c_start, c_end;
            bool degenerate;
            Orientation orient; // co-edge sense in the counter-clockwise wire
            T t1, t2;           // edge parameter range (u for v-isos, v for u-isos)
        };
        const std::array<Side, 4> sides{{
            {0, 1, c.degenerate_v1, Orientation::Forward, u1, u2},  // v = v1
            {1, 2, c.degenerate_u2, Orientation::Forward, v1, v2},  // u = u2
            {3, 2, c.degenerate_v2, Orientation::Reversed, u1, u2}, // v = v2
            {0, 3, c.degenerate_u1, Orientation::Reversed, v1, v2}, // u = u1
        }};
        auto pcurve = [&](const Side &s) -> std::shared_ptr<Curve<T, 2>> {
            return std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{uv[s.c_start], uv[s.c_end]},
                                                   std::vector<T>{s.t1, s.t1, s.t2, s.t2}, 1);
        };

        // ---- write (nothing can fail from here on)
        std::array<VertexId, 4> vtx;
        for (std::size_t i{}; i < 4; ++i)
            if (cls[i] == i)
            {
                T vt = tol;
                for (std::size_t j{}; j < 4; ++j)
                    if (cls[j] == i)
                        vt = std::max(vt, distance(p[i], p[j]));
                vtx[i] = m.add(Vertex<T>{p[i], vt});
            }
        for (std::size_t i{}; i < 4; ++i)
            vtx[i] = vtx[cls[i]];

        std::array<std::shared_ptr<Curve<T, 2>>, 4> pc;
        for (std::size_t k{}; k < 4; ++k)
            pc[k] = pcurve(sides[k]);
        std::array<EdgeId, 4> edge;
        auto make_side_edge = [&](std::size_t k) {
            const auto &s = sides[k];
            const auto va = vtx[s.c_start], vb = vtx[s.c_end];
            if (s.degenerate)
                return m.add(Edge<T>{nullptr, s.t1, s.t2, va, va, m.vertex(va).tol, true, true});
            return m.add(Edge<T>{std::make_shared<CurveOnSurface<T, 3>>(pc[k], srf), s.t1, s.t2, va, vb, tol});
        };
        edge[0] = make_side_edge(0);
        edge[3] = make_side_edge(3);
        edge[1] = c.closed_u ? edge[3] : make_side_edge(1); // seam u = u1 == u = u2
        edge[2] = c.closed_v ? edge[0] : make_side_edge(2); // seam v = v1 == v = v2

        Wire<T> w;
        for (std::size_t k{}; k < 4; ++k)
            w.coedges.push_back(CoEdge<T>{edge[k], sides[k].orient, pc[k]});
        w.closed = true;
        const auto wid = m.add(std::move(w));
        return m.add(Face<T>{std::move(srf), {wid}, tol, true});
    }

} // namespace gbs::brep
