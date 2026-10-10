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
#include <optional>
#include <span>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <gbs-brep/model.h>
#include <gbs-brep/closure.h>
#include <gbs-brep/explore.h>
#include <gbs-brep/pcurve.h>
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
        WireNotClosed,        ///< a face boundary must be a closed wire
        WireInUse,            ///< the wire already bounds a face (or is given twice)
        EdgeOffSurface,       ///< an edge is farther than pcurve_tol from the surface
        CrossesSeam,          ///< an edge or a wire goes across the seam of a closed surface
        PCurveApproximation,  ///< the pcurve could not reach pcurve_tol
        DegenerateWire,       ///< a boundary encloses no area in the parametric space
        HoleOutsideOuter,     ///< a hole is not inside the outer boundary
        ShellNotClosed,       ///< a solid needs closed shells
        ShellNotOrientable,   ///< a solid needs consistently oriented shells
        ZeroVolume,           ///< the shell encloses no volume
        InvalidFile,          ///< a file that cannot be read or parsed, or contains no shape
        UnsupportedEntity,    ///< an entity of the file that the reader does not support (strict reading)
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
        case BuildErrc::WireNotClosed: return "wire not closed";
        case BuildErrc::WireInUse: return "wire in use";
        case BuildErrc::EdgeOffSurface: return "edge off surface";
        case BuildErrc::CrossesSeam: return "crosses seam";
        case BuildErrc::PCurveApproximation: return "pcurve approximation";
        case BuildErrc::DegenerateWire: return "degenerate wire";
        case BuildErrc::HoleOutsideOuter: return "hole outside outer";
        case BuildErrc::ShellNotClosed: return "shell not closed";
        case BuildErrc::ShellNotOrientable: return "shell not orientable";
        case BuildErrc::ZeroVolume: return "zero volume";
        case BuildErrc::InvalidFile: return "invalid file";
        case BuildErrc::UnsupportedEntity: return "unsupported entity";
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

        /// Applies vertex merges (absorbed -> survivor): every edge of the model is redirected
        /// to the survivors, survivors keep the absorbed balls, absorbed vertices are erased.
        template <std::floating_point T>
        void merge_vertices(Model<T> &m, const std::vector<std::pair<VertexId, VertexId>> &merged)
        {
            if (merged.empty())
                return;
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
                    fit_vertex_tolerances(m, e);
                }
            }
            for (const auto &[a, s] : merged)
            {
                auto &sv = m.vertex(s);
                sv.tol = std::max(sv.tol, m.vertex(a).tol); // keep the absorbed ball's contract
                m.erase(a);
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
        detail::merge_vertices(m, merged);

        return m.add(Wire<T>{std::move(chain), closed});
    }

    template <std::floating_point T>
    [[nodiscard]] auto make_wire(Model<T> &m, const std::vector<EdgeId> &edges, std::type_identity_t<T> tol = brep_default_tolerance<T>) -> BuildResult<WireId>
    {
        return make_wire(m, std::span<const EdgeId>{edges}, tol);
    }

    /// A co-edge of an ordered wire: an edge and the sense it is traversed in.
    struct OrientedEdge
    {
        EdgeId edge;
        Orientation orient{Orientation::Forward};
    };

    /**
     * @brief Wire of co-edges given in order and sense (topology read from a
     * file, e.g. a STEP edge_loop), without reordering.
     *
     * - Each co-edge must end where the next one starts: same vertex, or
     *   vertices closer than `max(tol, tol_a + tol_b)`, which are then merged
     *   as in `make_wire`. The wire is closed iff the last co-edge ends where
     *   the first one starts.
     * - An edge may be used twice, in opposite senses: the seam of a closed
     *   surface (cylinder: bottom circle, seam, top circle reversed, seam
     *   reversed). A third use, or two uses in the same sense, is rejected.
     * - Degenerate edges are rejected: face builders insert them at poles.
     *
     * No pcurve is set (free wire). On failure the model is unchanged.
     */
    template <std::floating_point T>
    [[nodiscard]] auto make_wire_ordered(Model<T> &m, std::span<const OrientedEdge> coedges, std::type_identity_t<T> tol = brep_default_tolerance<T>) -> BuildResult<WireId>
    {
        if (coedges.empty())
            return detail::fail(BuildErrc::EmptyInput, "make_wire_ordered needs at least one co-edge");
        if (!detail::valid_tolerance(tol))
            return detail::fail(BuildErrc::InvalidTolerance, "wire tolerance must be finite and > 0");
        for (std::size_t i{}; i < coedges.size(); ++i)
        {
            const auto &ce = coedges[i];
            if (!m.alive(ce.edge))
                return detail::fail(BuildErrc::InvalidId, "make_wire_ordered: dead or invalid edge", {ce.edge});
            if (m.edge(ce.edge).degenerate)
                return detail::fail(BuildErrc::DegenerateEdgeInWire, "degenerate edges are inserted by face builders", {ce.edge});
            std::size_t uses = 0;
            for (std::size_t j{}; j < i; ++j)
                if (coedges[j].edge == ce.edge)
                {
                    ++uses;
                    if (coedges[j].orient == ce.orient)
                        return detail::fail(BuildErrc::DuplicateEdge, "edge used twice in the same sense", {ce.edge});
                }
            if (uses > 1)
                return detail::fail(BuildErrc::DuplicateEdge, "edge used more than twice", {ce.edge});
        }

        auto start = [&](const OrientedEdge &ce) { const auto &e = m.edge(ce.edge); return ce.orient == Orientation::Forward ? e.v1 : e.v2; };
        auto end = [&](const OrientedEdge &ce) { const auto &e = m.edge(ce.edge); return ce.orient == Orientation::Forward ? e.v2 : e.v1; };
        auto close = [&](VertexId a, VertexId b) {
            const auto &va = m.vertex(a), &vb = m.vertex(b);
            return a == b || distance(va.pnt, vb.pnt) <= std::max(tol, va.tol + vb.tol);
        };

        // ---- junctions: vertex pairs to merge (union-find, smallest id survives)
        std::unordered_map<VertexId, VertexId> parent;
        auto find = [&](VertexId v) {
            while (parent.contains(v) && parent[v] != v)
                v = parent[v];
            return v;
        };
        auto unite = [&](VertexId a, VertexId b) {
            a = find(a), b = find(b);
            if (a != b)
                parent[std::max(a, b)] = std::min(a, b);
        };
        for (std::size_t i{}; i + 1 < coedges.size(); ++i)
        {
            const auto a = end(coedges[i]), b = start(coedges[i + 1]);
            if (!close(a, b))
                return detail::fail(BuildErrc::Disconnected, "co-edge " + std::to_string(i) + " does not end where the next one starts",
                                    {coedges[i].edge, coedges[i + 1].edge});
            unite(a, b);
        }
        const bool closed = close(end(coedges.back()), start(coedges.front()));
        if (closed)
            unite(end(coedges.back()), start(coedges.front()));

        // ---- write
        std::vector<std::pair<VertexId, VertexId>> merged;
        for (const auto &[v, p] : parent)
            if (auto r = find(v); r != v)
                merged.emplace_back(v, r);
        std::ranges::sort(merged);
        detail::merge_vertices(m, merged);

        Wire<T> w;
        w.closed = closed;
        for (const auto &ce : coedges)
            w.coedges.push_back(CoEdge<T>{ce.edge, ce.orient, nullptr});
        return m.add(std::move(w));
    }

    template <std::floating_point T>
    [[nodiscard]] auto make_wire_ordered(Model<T> &m, const std::vector<OrientedEdge> &coedges, std::type_identity_t<T> tol = brep_default_tolerance<T>) -> BuildResult<WireId>
    {
        return make_wire_ordered(m, std::span<const OrientedEdge>{coedges}, tol);
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

    // =========================================================================
    // Face bounded by wires on a surface
    // =========================================================================

    namespace detail
    {
        template <std::floating_point T>
        auto check_surface_bounds(const Surface<T, 3> &srf) -> std::optional<BuildError>
        {
            const auto [u1, u2, v1, v2] = srf.bounds();
            for (T x : {u1, u2, v1, v2})
                if (!std::isfinite(x) || std::abs(x) >= unbounded_parameter<T>)
                    return BuildError{BuildErrc::UnboundedSurface, "surface bounds must be finite", {}};
            if (!(u1 < u2) || !(v1 < v2))
                return BuildError{BuildErrc::UnboundedSurface, "surface parametric rectangle is empty", {}};
            return std::nullopt;
        }

        /// Start / end (u,v) of a co-edge, its sense applied.
        template <std::floating_point T>
        auto coedge_uv_end(const Model<T> &m, const CoEdge<T> &ce, bool end) -> point<T, 2>
        {
            const auto &e = m.edge(ce.edge);
            return coedge_uv(m, ce, end ? e.u2 : e.u1);
        }
    } // namespace detail

    namespace detail
    {
        template <std::floating_point T>
        struct BuiltFace
        {
            FaceId face;
            bool outer_reversed; ///< the outer wire was turned to be counter-clockwise in (u,v)
        };

        /// Core of the face-from-wires builders; see `make_face` and `make_face_use`.
        template <std::floating_point T>
        auto make_face_from_wires(Model<T> &m, std::shared_ptr<Surface<T, 3>> srf, std::span<const WireId> wires,
                                  const MakeFaceOptions<T> &opts) -> BuildResult<BuiltFace<T>>
        {
            if (!srf)
                return fail(BuildErrc::NullSurface, "make_face needs a surface");
            if (!valid_tolerance(opts.tol) || !valid_tolerance(opts.pcurve_tol))
                return fail(BuildErrc::InvalidTolerance, "face and pcurve tolerances must be finite and > 0");
            if (auto err = check_surface_bounds(*srf))
                return std::unexpected(*err);

            for (std::size_t i{}; i < wires.size(); ++i)
            {
                const auto wid = wires[i];
                if (!m.alive(wid))
                    return fail(BuildErrc::InvalidId, "make_face: dead or invalid wire", {wid});
                for (std::size_t j{}; j < i; ++j)
                    if (wires[j] == wid)
                        return fail(BuildErrc::WireInUse, "wire given twice", {wid});
                if (!m.wire(wid).closed || !is_closed(m, wid))
                    return fail(BuildErrc::WireNotClosed, "face boundaries must be closed wires", {wid});
                for (const auto &ce : m.wire(wid).coedges)
                {
                    if (ce.pcurve)
                        return fail(BuildErrc::WireInUse, "wire already carries pcurves (bounds a face)", {wid});
                    if (m.edge(ce.edge).degenerate)
                        return fail(BuildErrc::DegenerateEdgeInWire, "degenerate edge in a face boundary wire", {ce.edge});
                }
            }
            for (auto fid : m.template ids<FaceId>())
                for (auto w : m.face(fid).wires)
                    if (std::ranges::find(wires, w) != wires.end())
                        return fail(BuildErrc::WireInUse, "wire already bounds a face", {w, fid});

            const auto cl = surface_closure(*srf, opts.tol);
            const auto [su1, su2, sv1, sv2] = srf->bounds();
            const std::array<T, 2> lo{su1, sv1}, hi{su2, sv2};
            const std::array<T, 2> period{su2 - su1, sv2 - sv1};
            const std::array<bool, 2> closed{cl.closed_u, cl.closed_v};

            // degenerate edges inserted at poles; erased if the build fails
            std::vector<EdgeId> inserted;
            auto rollback = [&](BuildError err) -> BuildResult<BuiltFace<T>> {
                for (auto e : inserted)
                    m.erase(e);
                return std::unexpected(std::move(err));
            };

            // ---- 1. pcurves, computed on copies of the co-edges (the wires are not touched yet)
            std::vector<std::vector<CoEdge<T>>> cw;
            std::vector<std::vector<std::uint8_t>> exact;
            std::unordered_map<EdgeId, T> deviation;
            auto fit_coedge = [&](CoEdge<T> &ce, const std::optional<point<T, 2>> &hint, std::uint8_t &is_exact) -> std::optional<BuildError> {
                const auto &e = m.edge(ce.edge);
                if (auto pc = extract_pcurve(e, srf))
                {
                    ce.pcurve = std::move(pc);
                    deviation.try_emplace(ce.edge, T(0));
                    is_exact = 1;
                    return std::nullopt;
                }
                auto r = project_pcurve(e, *srf, cl, hint, opts);
                if (!r)
                {
                    switch (r.error())
                    {
                    case PCurveErrc::OffSurface:
                        return BuildError{BuildErrc::EdgeOffSurface, "edge is not on the surface within pcurve_tol", {ce.edge}};
                    case PCurveErrc::CrossesSeam:
                        return BuildError{BuildErrc::CrossesSeam, "edge goes across the seam; split it at the seam", {ce.edge}};
                    default:
                        return BuildError{BuildErrc::PCurveApproximation, "pcurve did not reach pcurve_tol", {ce.edge}};
                    }
                }
                ce.pcurve = r->pcurve;
                auto &d = deviation[ce.edge];
                d = std::max(d, r->deviation);
                is_exact = 0;
                return std::nullopt;
            };
            for (auto wid : wires)
            {
                cw.push_back(m.wire(wid).coedges);
                exact.emplace_back(cw.back().size(), 0);
                for (std::size_t i{}; i < cw.back().size(); ++i)
                    if (auto err = fit_coedge(cw.back()[i], std::nullopt, exact.back()[i]))
                        return std::unexpected(*err);
            }

            // ---- 2. side of the co-edges lying on a seam
            // A co-edge entirely on the seam of a closed direction k (u = lo or u = hi) can be on either
            // side. Its side follows the neighbour co-edge it touches, unless that junction is a pole
            // (k is free there). Otherwise it is opposite to the other use of the same edge (both
            // sides of a seam), or lo by default.
            const auto seam_eps = [&](std::size_t k) { return T(1e-7) * period[k]; };
            auto on_seam = [&](const CoEdge<T> &ce, std::size_t k) {
                const auto &e = m.edge(ce.edge);
                for (int j = 0; j <= 4; ++j)
                {
                    const T x = ce.pcurve->value(e.u1 + (e.u2 - e.u1) * T(j) / T(4))[k];
                    if (std::abs(x - lo[k]) > seam_eps(k) && std::abs(x - hi[k]) > seam_eps(k))
                        return false;
                }
                return true;
            };
            std::vector<std::vector<std::array<std::uint8_t, 2>>> ambiguous(cw.size());
            for (std::size_t wi{}; wi < cw.size(); ++wi)
                for (std::size_t i{}; i < cw[wi].size(); ++i)
                {
                    auto &a = ambiguous[wi].emplace_back(std::array<std::uint8_t, 2>{0, 0});
                    for (std::size_t k = 0; k < 2; ++k)
                        a[k] = closed[k] && !exact[wi][i] && on_seam(cw[wi][i], k);
                }
            auto put_on_side = [&](std::size_t wi, std::size_t i, std::size_t k, T side) -> std::optional<BuildError> {
                auto &ce = cw[wi][i];
                auto hint = coedge_uv_end(m, ce, false);
                if (std::abs(hint[k] - side) > seam_eps(k)) // on the other side: project again with the side as hint
                {
                    hint[k] = side;
                    if (auto err = fit_coedge(ce, hint, exact[wi][i]))
                        return err;
                }
                ambiguous[wi][i][k] = 0;
                return std::nullopt;
            };
            auto nearest_side = [&](std::size_t k, T x) { return std::abs(x - lo[k]) <= std::abs(x - hi[k]) ? lo[k] : hi[k]; };
            for (bool progress = true; progress;)
            {
                progress = false;
                for (std::size_t wi{}; wi < cw.size(); ++wi)
                {
                    const auto n = cw[wi].size();
                    for (std::size_t i{}; i < n; ++i)
                        for (std::size_t k = 0; k < 2; ++k)
                        {
                            if (!ambiguous[wi][i][k])
                                continue;
                            const std::size_t ip = (i + n - 1) % n, in = (i + 1) % n;
                            std::optional<T> side;
                            if (!ambiguous[wi][ip][k])
                                if (const auto a = coedge_uv_end(m, cw[wi][ip], true); !free_coordinate(*srf, a, k))
                                    side = nearest_side(k, a[k]);
                            if (!side && !ambiguous[wi][in][k])
                                if (const auto b = coedge_uv_end(m, cw[wi][in], false); !free_coordinate(*srf, b, k))
                                    side = nearest_side(k, b[k]);
                            if (!side)
                                continue;
                            if (auto err = put_on_side(wi, i, k, *side))
                                return rollback(*err);
                            progress = true;
                        }
                }
            }
            for (std::size_t wi{}; wi < cw.size(); ++wi)
                for (std::size_t i{}; i < cw[wi].size(); ++i)
                    for (std::size_t k = 0; k < 2; ++k)
                    {
                        if (!ambiguous[wi][i][k])
                            continue;
                        T side = lo[k];
                        for (std::size_t wj{}; wj < cw.size(); ++wj)
                            for (std::size_t j{}; j < cw[wj].size(); ++j)
                                if ((wj != wi || j != i) && cw[wj][j].edge == cw[wi][i].edge && !ambiguous[wj][j][k])
                                    side = nearest_side(k, coedge_uv_end(m, cw[wj][j], false)[k]) == lo[k] ? hi[k] : lo[k];
                        if (auto err = put_on_side(wi, i, k, side))
                            return rollback(*err);
                    }

            // ---- 3. degenerate co-edges where two co-edges meet at a pole with different (u,v)
            for (std::size_t wi{}; wi < cw.size(); ++wi)
            {
                auto &coedges = cw[wi];
                for (std::size_t i{}; i < coedges.size(); ++i)
                {
                    const auto &ce = coedges[i];
                    const auto &next = coedges[(i + 1) % coedges.size()];
                    if (m.edge(ce.edge).degenerate)
                        continue;
                    const auto a = coedge_uv_end(m, ce, true), b = coedge_uv_end(m, next, false);
                    std::optional<std::size_t> moving;
                    bool other_equal = true;
                    for (std::size_t k = 0; k < 2; ++k)
                        if (std::abs(a[k] - b[k]) > T(1e-7) * period[k])
                        {
                            if (moving)
                                other_equal = false;
                            moving = k;
                        }
                    if (!moving || !other_equal)
                        continue;
                    const std::size_t k = *moving;
                    if (!free_coordinate(*srf, a, k) || !free_coordinate(*srf, b, k))
                        continue; // not a pole: left to the continuity check
                    const bool forward = b[k] > a[k];
                    const T t1 = std::min(a[k], b[k]), t2 = std::max(a[k], b[k]);
                    auto pole = make_degenerate_edge(m, coedge_end(m, ce), t1, t2);
                    if (!pole)
                        return rollback(pole.error());
                    inserted.push_back(*pole);
                    const auto &p1 = forward ? a : b, &p2 = forward ? b : a;
                    auto pc = std::make_shared<BSCurve<T, 2>>(points_vector<T, 2>{p1, p2}, std::vector<T>{t1, t1, t2, t2}, 1);
                    coedges.insert(coedges.begin() + static_cast<std::ptrdiff_t>(i) + 1,
                                   CoEdge<T>{*pole, forward ? Orientation::Forward : Orientation::Reversed, std::move(pc)});
                    ++i;
                }
            }

            // ---- 4. continuity, orientation, holes
            bool outer_reversed = false;
            for (std::size_t wi{}; wi < cw.size(); ++wi)
            {
                auto &coedges = cw[wi];
                for (std::size_t k = 0; k < 2; ++k)
                    if (closed[k])
                        for (std::size_t i{}; i < coedges.size(); ++i)
                        {
                            const auto a = coedge_uv_end(m, coedges[i], true);
                            const auto b = coedge_uv_end(m, coedges[(i + 1) % coedges.size()], false);
                            if (std::abs(a[k] - b[k]) > period[k] / 2)
                                return rollback(BuildError{BuildErrc::CrossesSeam, "wire goes around the seam of a closed surface", {wires[wi]}});
                        }

                const T area = uv_signed_area(m, std::span<const CoEdge<T>>{coedges});
                if (!(std::abs(area) > T(1e-12) * period[0] * period[1]))
                    return rollback(BuildError{BuildErrc::DegenerateWire, "boundary encloses no area in the parametric space", {wires[wi]}});
                if ((area > T(0)) != (wi == 0)) // outer counter-clockwise, holes clockwise
                {
                    std::ranges::reverse(coedges);
                    for (auto &ce : coedges)
                        ce.orient = reverse(ce.orient);
                    if (wi == 0)
                        outer_reversed = true;
                }
            }

            const auto outer_poly = uv_polygon(m, std::span<const CoEdge<T>>{cw.front()});
            for (std::size_t wi{1}; wi < wires.size(); ++wi)
            {
                const auto hole_poly = uv_polygon(m, std::span<const CoEdge<T>>{cw[wi]});
                if (!std::ranges::all_of(hole_poly, [&](const point<T, 2> &p) { return inside(outer_poly, p); }))
                    return rollback(BuildError{BuildErrc::HoleOutsideOuter, "hole not inside the outer boundary", {wires[wi]}});
            }

            // ---- write (nothing can fail from here on)
            for (std::size_t wi{}; wi < wires.size(); ++wi)
                m.wire(wires[wi]).coedges = std::move(cw[wi]);
            for (const auto &[eid, dev] : deviation)
            {
                auto &e = m.edge(eid);
                e.tol = std::max({e.tol, dev, opts.tol});
                fit_vertex_tolerances(m, e);
            }
            const auto f = m.add(Face<T>{std::move(srf), std::vector<WireId>(wires.begin(), wires.end()), opts.tol, false});
            return BuiltFace<T>{f, outer_reversed};
        }
    } // namespace detail

    /**
     * @brief Face on `srf` bounded by the closed free wire `outer` and the
     * closed free wires `holes`.
     *
     * For every co-edge a pcurve is computed: extracted exactly when the edge
     * is a `CurveOnSurface` on this very surface, projected and interpolated
     * otherwise (`project_pcurve`, accuracy `opts.pcurve_tol`). On a closed
     * surface the samples are made continuous across the seam; an edge or a
     * wire going across the seam is rejected (to be split first). A co-edge
     * lying on the seam takes the side of the co-edge it touches; when both
     * uses of a seam edge are in the face (wire from `make_wire_ordered`),
     * they are put on opposite sides. Where two consecutive co-edges meet at a
     * pole (sphere, cone apex) with different (u,v), a degenerate co-edge is
     * inserted. The outer wire is turned counter-clockwise in (u,v), the holes
     * clockwise, by reversing the order and the senses of their co-edges; each
     * hole must lie inside the outer boundary. The wires are then attached to
     * the face (their co-edges receive the pcurves); edge tolerances are raised
     * to the pcurve deviation and to `opts.tol`, vertex tolerances follow.
     * On failure the model is unchanged.
     */
    template <std::floating_point T>
    [[nodiscard]] auto make_face(Model<T> &m, std::type_identity_t<std::shared_ptr<Surface<T, 3>>> srf, WireId outer,
                                 std::span<const WireId> holes, std::type_identity_t<MakeFaceOptions<T>> opts = {}) -> BuildResult<FaceId>
    {
        std::vector<WireId> wires{outer};
        wires.insert(wires.end(), holes.begin(), holes.end());
        auto r = detail::make_face_from_wires(m, std::move(srf), std::span<const WireId>{wires}, opts);
        if (!r)
            return std::unexpected(r.error());
        return r->face;
    }

    /**
     * @brief Face bounded by wires given in the sense of the face, as a face
     * use that keeps those senses (STEP: the outer loop is counter-clockwise
     * seen from the face normal, which opposes the surface normal when
     * `advanced_face.same_sense` is false).
     *
     * Builds the face as `make_face`, which turns the outer wire
     * counter-clockwise in (u,v) (that is, around the surface normal). The use
     * is Forward if the wire already was, Reversed otherwise: composed with the
     * use, every co-edge keeps the sense it was given, so the edges of
     * neighbouring faces keep opposite senses in a shell.
     */
    template <std::floating_point T>
    [[nodiscard]] auto make_face_use(Model<T> &m, std::type_identity_t<std::shared_ptr<Surface<T, 3>>> srf, WireId outer,
                                     std::span<const WireId> holes, std::type_identity_t<MakeFaceOptions<T>> opts = {}) -> BuildResult<FaceUse>
    {
        std::vector<WireId> wires{outer};
        wires.insert(wires.end(), holes.begin(), holes.end());
        auto r = detail::make_face_from_wires(m, std::move(srf), std::span<const WireId>{wires}, opts);
        if (!r)
            return std::unexpected(r.error());
        return FaceUse{r->face, r->outer_reversed ? Orientation::Reversed : Orientation::Forward};
    }

    template <std::floating_point T>
    [[nodiscard]] auto make_face_use(Model<T> &m, std::type_identity_t<std::shared_ptr<Surface<T, 3>>> srf, WireId outer,
                                     const std::vector<WireId> &holes = {}, std::type_identity_t<MakeFaceOptions<T>> opts = {}) -> BuildResult<FaceUse>
    {
        return make_face_use(m, std::move(srf), outer, std::span<const WireId>{holes}, opts);
    }

    /// Face on `srf` bounded by the closed free wire `outer`, without holes.
    template <std::floating_point T>
    [[nodiscard]] auto make_face(Model<T> &m, std::type_identity_t<std::shared_ptr<Surface<T, 3>>> srf, WireId outer,
                                 std::type_identity_t<MakeFaceOptions<T>> opts = {}) -> BuildResult<FaceId>
    {
        return make_face(m, std::move(srf), outer, std::span<const WireId>{}, opts);
    }

    template <std::floating_point T>
    [[nodiscard]] auto make_face(Model<T> &m, std::type_identity_t<std::shared_ptr<Surface<T, 3>>> srf, WireId outer,
                                 const std::vector<WireId> &holes, std::type_identity_t<MakeFaceOptions<T>> opts = {}) -> BuildResult<FaceId>
    {
        return make_face(m, std::move(srf), outer, std::span<const WireId>{holes}, opts);
    }

    /**
     * @brief Face on `srf` bounded by closed loops of 2D curves given in its
     * parametric space (the stage 2 path: trimming curves from intersections).
     * Each 2D curve becomes an edge whose 3D curve is the exact
     * `CurveOnSurface`, the loops become wires (`make_wire`, vertex merge
     * within `opts.tol`), then `make_face(srf, outer, holes)` extracts the
     * pcurves exactly. On failure every entity created is erased.
     */
    template <std::floating_point T>
    [[nodiscard]] auto make_face(Model<T> &m, std::type_identity_t<std::shared_ptr<Surface<T, 3>>> srf,
                                 const std::vector<std::shared_ptr<Curve<T, 2>>> &outer,
                                 const std::vector<std::vector<std::shared_ptr<Curve<T, 2>>>> &holes = {},
                                 std::type_identity_t<MakeFaceOptions<T>> opts = {}) -> BuildResult<FaceId>
    {
        if (!srf)
            return detail::fail(BuildErrc::NullSurface, "make_face needs a surface");
        const auto n_vtx = static_cast<std::uint32_t>(m.template capacity<ShapeType::Vertex>());
        const auto n_edg = static_cast<std::uint32_t>(m.template capacity<ShapeType::Edge>());
        const auto n_wir = static_cast<std::uint32_t>(m.template capacity<ShapeType::Wire>());
        auto rollback = [&](BuildError err) -> BuildResult<FaceId> {
            for (auto i = n_wir; i < m.template capacity<ShapeType::Wire>(); ++i)
                m.erase(WireId{i});
            for (auto i = n_edg; i < m.template capacity<ShapeType::Edge>(); ++i)
                m.erase(EdgeId{i});
            for (auto i = n_vtx; i < m.template capacity<ShapeType::Vertex>(); ++i)
                m.erase(VertexId{i});
            return std::unexpected(std::move(err));
        };

        auto make_loop = [&](const std::vector<std::shared_ptr<Curve<T, 2>>> &loop) -> BuildResult<WireId> {
            std::vector<EdgeId> edges;
            for (const auto &pc : loop)
            {
                if (!pc)
                    return detail::fail(BuildErrc::NullCurve, "null 2D curve in a loop");
                auto e = make_edge(m, std::make_shared<CurveOnSurface<T, 3>>(pc, srf), opts.tol);
                if (!e)
                    return std::unexpected(e.error());
                edges.push_back(*e);
            }
            return make_wire(m, edges, opts.tol);
        };

        auto wo = make_loop(outer);
        if (!wo)
            return rollback(wo.error());
        std::vector<WireId> wh;
        for (const auto &loop : holes)
        {
            auto w = make_loop(loop);
            if (!w)
                return rollback(w.error());
            wh.push_back(*w);
        }
        auto f = make_face(m, srf, *wo, std::span<const WireId>{wh}, opts);
        if (!f)
            return rollback(f.error());
        return f;
    }

    // =========================================================================
    // Solid and compound
    // =========================================================================

    /**
     * @brief Solid bounded by the closed shell `outer` and the closed shells
     * `voids` (cavities).
     *
     * Every shell must be closed and consistently oriented. The outer shell is
     * turned so that its normals point outwards (positive signed volume), each
     * cavity so that its normals point into the cavity (negative signed
     * volume), by reversing all the face uses of a shell when needed. The
     * inclusion of the cavities in the outer shell is not checked at stage 1.
     * On failure the model is unchanged.
     */
    template <std::floating_point T>
    [[nodiscard]] auto make_solid(Model<T> &m, ShellId outer, std::span<const ShellId> voids, std::size_t n_volume = 64) -> BuildResult<SolidId>
    {
        std::vector<ShellId> shells{outer};
        shells.insert(shells.end(), voids.begin(), voids.end());
        std::vector<T> volume;
        for (std::size_t i{}; i < shells.size(); ++i)
        {
            const auto sid = shells[i];
            if (!m.alive(sid))
                return detail::fail(BuildErrc::InvalidId, "make_solid: dead or invalid shell", {sid});
            for (std::size_t j{}; j < i; ++j)
                if (shells[j] == sid)
                    return detail::fail(BuildErrc::InvalidId, "make_solid: shell given twice", {sid});
            for (const auto &fu : m.shell(sid).faces)
                if (!m.alive(fu.face))
                    return detail::fail(BuildErrc::InvalidId, "make_solid: dead face in shell", {sid, fu.face});
            if (!is_closed(m, sid))
                return detail::fail(BuildErrc::ShellNotClosed, "make_solid needs closed shells", {sid});
            if (!is_orientable(m, sid))
                return detail::fail(BuildErrc::ShellNotOrientable, "make_solid needs consistently oriented shells", {sid});
            const T v = signed_volume(m, sid, n_volume);
            if (!(std::abs(v) > T(0)) || !std::isfinite(v))
                return detail::fail(BuildErrc::ZeroVolume, "shell encloses no volume", {sid});
            volume.push_back(v);
        }
        // ---- write
        for (std::size_t i{}; i < shells.size(); ++i)
        {
            const bool want_positive = i == 0;
            if ((volume[i] > T(0)) != want_positive)
                for (auto &fu : m.shell(shells[i]).faces)
                    fu.orient = reverse(fu.orient);
            m.shell(shells[i]).closed = true;
        }
        return m.add(Solid{outer, std::vector<ShellId>(voids.begin(), voids.end())});
    }

    template <std::floating_point T>
    [[nodiscard]] auto make_solid(Model<T> &m, ShellId outer, const std::vector<ShellId> &voids = {}, std::size_t n_volume = 64) -> BuildResult<SolidId>
    {
        return make_solid(m, outer, std::span<const ShellId>{voids}, n_volume);
    }

    /// Compound of any live shapes (a compound may contain compounds).
    template <std::floating_point T>
    [[nodiscard]] auto make_compound(Model<T> &m, std::span<const ShapeId> shapes) -> BuildResult<CompoundId>
    {
        for (const auto &s : shapes)
            if (!m.alive(s))
                return detail::fail(BuildErrc::InvalidId, "make_compound: dead or invalid shape", {s});
        return m.add(Compound{std::vector<ShapeId>(shapes.begin(), shapes.end())});
    }

    template <std::floating_point T>
    [[nodiscard]] auto make_compound(Model<T> &m, const std::vector<ShapeId> &shapes) -> BuildResult<CompoundId>
    {
        return make_compound(m, std::span<const ShapeId>{shapes});
    }

} // namespace gbs::brep
