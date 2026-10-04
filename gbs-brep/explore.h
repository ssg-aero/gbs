#pragma once

/**
 * @file explore.h
 * @brief Native BREP core, stage 1: topological explorer, adjacency index and
 * base queries (closure, manifoldness, orientation consistency, bounding box).
 *
 * Design: docs/sources/design/brep_core.md, section 5.
 */

#include <algorithm>
#include <cstdint>
#include <span>
#include <unordered_map>
#include <utility>
#include <vector>

#include <gbs-brep/model.h>

namespace gbs::brep
{
    // =========================================================================
    // Explorer: sub-entities of a shape, by type
    // =========================================================================

    namespace detail
    {
        template <std::floating_point T>
        struct Visitor
        {
            const Model<T> &m;
            std::array<std::vector<std::uint8_t>, 7> seen; // one visited table per ShapeType

            explicit Visitor(const Model<T> &model) : m{model}
            {
                seen[0].assign(m.template capacity<ShapeType::Vertex>(), 0);
                seen[1].assign(m.template capacity<ShapeType::Edge>(), 0);
                seen[2].assign(m.template capacity<ShapeType::Wire>(), 0);
                seen[3].assign(m.template capacity<ShapeType::Face>(), 0);
                seen[4].assign(m.template capacity<ShapeType::Shell>(), 0);
                seen[5].assign(m.template capacity<ShapeType::Solid>(), 0);
                seen[6].assign(m.template capacity<ShapeType::Compound>(), 0);
            }

            template <ShapeType K>
            bool first_visit(Id<K> h)
            {
                auto &tbl = seen[static_cast<std::size_t>(K)];
                if (!h.valid() || h.index >= tbl.size() || tbl[h.index])
                    return false;
                tbl[h.index] = 1;
                return true;
            }

            // Depth-first descent; `out` receives every first-visited id of kind Sub.
            template <typename Sub>
            void visit(const ShapeId &id, std::vector<Sub> &out)
            {
                std::visit([&](auto h) { visit_id(h, out); }, id);
            }

            template <ShapeType K, typename Sub>
            void visit_id(Id<K> h, std::vector<Sub> &out)
            {
                if (!m.alive(h) || !first_visit(h))
                    return;
                if constexpr (K == Sub::type)
                    out.push_back(h);

                if constexpr (K == ShapeType::Compound)
                {
                    for (const auto &s : m.compound(h).shapes)
                        visit(s, out);
                }
                else if constexpr (K == ShapeType::Solid)
                {
                    const auto &so = m.solid(h);
                    visit_id(so.outer, out);
                    for (auto v : so.voids)
                        visit_id(v, out);
                }
                else if constexpr (K == ShapeType::Shell)
                {
                    for (const auto &fu : m.shell(h).faces)
                        visit_id(fu.face, out);
                }
                else if constexpr (K == ShapeType::Face)
                {
                    for (auto w : m.face(h).wires)
                        visit_id(w, out);
                }
                else if constexpr (K == ShapeType::Wire)
                {
                    for (const auto &ce : m.wire(h).coedges)
                        visit_id(ce.edge, out);
                }
                else if constexpr (K == ShapeType::Edge)
                {
                    const auto &e = m.edge(h);
                    visit_id(e.v1, out);
                    visit_id(e.v2, out);
                }
            }
        };
    } // namespace detail

    /**
     * @brief Sub-entities of kind `Sub` of `shape`, without duplicates, in
     * depth-first discovery order. A seam edge appears once. The shape itself
     * is included when it is of kind `Sub`.
     *
     * Example: `explore<EdgeId>(m, FaceId{3})`, `explore<FaceId>(m, solid)`.
     */
    template <typename Sub, std::floating_point T>
    [[nodiscard]] auto explore(const Model<T> &m, const ShapeId &shape) -> std::vector<Sub>
    {
        std::vector<Sub> out;
        detail::Visitor<T> v{m};
        v.visit(shape, out);
        return out;
    }

    // =========================================================================
    // Adjacency index (ancestors), built on demand
    // =========================================================================

    /// Position of a co-edge: the wire and the index in `Wire::coedges`.
    struct CoEdgeRef
    {
        WireId wire;
        std::uint32_t index{npos};
        friend constexpr auto operator<=>(const CoEdgeRef &, const CoEdgeRef &) = default;
    };

    /**
     * @brief Upward adjacency of the entities under a root shape, like
     * `TopExp::MapShapesAndAncestors`. Built once, read many times; it is
     * invalidated by any mutation of the model.
     */
    template <std::floating_point T>
    class TopologyIndex
    {
        std::vector<std::vector<EdgeId>> m_edges_of_vertex;
        std::vector<std::vector<CoEdgeRef>> m_coedges_of_edge;
        std::vector<std::vector<FaceId>> m_faces_of_edge;
        std::vector<FaceId> m_face_of_wire;
        std::vector<std::vector<ShellId>> m_shells_of_face;
        std::vector<EdgeId> m_edges;
        std::vector<FaceId> m_faces;

        template <typename V>
        static auto at(const std::vector<V> &tbl, std::uint32_t i) -> std::span<const typename V::value_type>
        {
            if (i >= tbl.size())
                return {};
            return tbl[i];
        }

        static void push_unique(auto &vec, auto id)
        {
            if (std::find(vec.begin(), vec.end(), id) == vec.end())
                vec.push_back(id);
        }

    public:
        TopologyIndex(const Model<T> &m, const ShapeId &root)
        {
            m_edges_of_vertex.resize(m.template capacity<ShapeType::Vertex>());
            m_coedges_of_edge.resize(m.template capacity<ShapeType::Edge>());
            m_faces_of_edge.resize(m.template capacity<ShapeType::Edge>());
            m_face_of_wire.resize(m.template capacity<ShapeType::Wire>());
            m_shells_of_face.resize(m.template capacity<ShapeType::Face>());

            m_edges = explore<EdgeId>(m, root);
            m_faces = explore<FaceId>(m, root);

            for (auto eid : m_edges)
            {
                const auto &e = m.edge(eid);
                if (e.v1.valid()) push_unique(m_edges_of_vertex[e.v1.index], eid);
                if (e.v2.valid()) push_unique(m_edges_of_vertex[e.v2.index], eid);
            }
            for (auto fid : m_faces)
            {
                const auto &f = m.face(fid);
                for (auto wid : f.wires)
                {
                    if (!m.alive(wid))
                        continue;
                    m_face_of_wire[wid.index] = fid;
                    const auto &w = m.wire(wid);
                    for (std::uint32_t i{}; i < w.coedges.size(); ++i)
                    {
                        auto eid = w.coedges[i].edge;
                        if (!m.alive(eid))
                            continue;
                        m_coedges_of_edge[eid.index].push_back(CoEdgeRef{wid, i});
                        push_unique(m_faces_of_edge[eid.index], fid);
                    }
                }
            }
            for (auto sid : explore<ShellId>(m, root))
                for (const auto &fu : m.shell(sid).faces)
                    if (m.alive(fu.face))
                        push_unique(m_shells_of_face[fu.face.index], sid);
        }

        /// Edges incident to a vertex.
        [[nodiscard]] auto edges_of(VertexId v) const -> std::span<const EdgeId> { return at(m_edges_of_vertex, v.index); }
        /// Faces using an edge: 0 (free edge of a wire), 1 (boundary), 2 (manifold), more (non-manifold).
        [[nodiscard]] auto faces_of(EdgeId e) const -> std::span<const FaceId> { return at(m_faces_of_edge, e.index); }
        /// Co-edges using an edge (two in the same wire for a seam).
        [[nodiscard]] auto coedges_of(EdgeId e) const -> std::span<const CoEdgeRef> { return at(m_coedges_of_edge, e.index); }
        /// Face owning a wire (invalid for a free wire).
        [[nodiscard]] auto face_of(WireId w) const -> FaceId { return w.index < m_face_of_wire.size() ? m_face_of_wire[w.index] : FaceId{}; }
        /// Shells using a face.
        [[nodiscard]] auto shells_of(FaceId f) const -> std::span<const ShellId> { return at(m_shells_of_face, f.index); }

        [[nodiscard]] const std::vector<EdgeId> &edges() const noexcept { return m_edges; }
        [[nodiscard]] const std::vector<FaceId> &faces() const noexcept { return m_faces; }
    };

    // =========================================================================
    // Wire queries
    // =========================================================================

    /// True if consecutive co-edges share a vertex (orientation applied) and the last one closes on the first.
    template <std::floating_point T>
    [[nodiscard]] bool is_closed(const Model<T> &m, WireId wid)
    {
        const auto &w = m.wire(wid);
        if (w.coedges.empty())
            return false;
        for (std::size_t i{}; i < w.coedges.size(); ++i)
        {
            const auto &ce = w.coedges[i];
            const auto &next = w.coedges[(i + 1) % w.coedges.size()];
            if (!m.alive(ce.edge) || !m.alive(next.edge))
                return false;
            if (coedge_end(m, ce) != coedge_start(m, next))
                return false;
        }
        return true;
    }

    /// True if consecutive co-edges share a vertex, the chain being open or closed.
    template <std::floating_point T>
    [[nodiscard]] bool is_chained(const Model<T> &m, WireId wid)
    {
        const auto &w = m.wire(wid);
        if (w.coedges.empty())
            return false;
        for (std::size_t i{}; i + 1 < w.coedges.size(); ++i)
        {
            if (!m.alive(w.coedges[i].edge) || !m.alive(w.coedges[i + 1].edge))
                return false;
            if (coedge_end(m, w.coedges[i]) != coedge_start(m, w.coedges[i + 1]))
                return false;
        }
        return m.alive(w.coedges.back().edge);
    }

    /// 3D point of the i-th co-edge of a wire at parameter t of its edge, orientation applied
    /// (t grows along the co-edge: t = u1 is the co-edge start, t = u2 its end).
    template <std::floating_point T>
    [[nodiscard]] auto coedge_point(const Model<T> &m, WireId wid, std::size_t i, T t) -> point<T, 3>
    {
        const auto &ce = m.wire(wid).coedges.at(i);
        const auto &e = m.edge(ce.edge);
        const T u = ce.orient == Orientation::Forward ? t : e.u1 + e.u2 - t;
        return edge_point(m, ce.edge, u);
    }

    // =========================================================================
    // Shell queries
    // =========================================================================

    /// One usage of an edge by a shell: the face and the effective sense
    /// (co-edge sense composed with the face-use sense).
    struct EdgeUse
    {
        FaceId face;
        CoEdgeRef coedge;
        Orientation orient;
    };

    /**
     * @brief Usages of every non-degenerate edge of a shell, keyed by edge.
     * Building block of `is_closed`, `is_manifold`, `is_orientable`,
     * `free_edges` and of the sewing.
     */
    template <std::floating_point T>
    [[nodiscard]] auto shell_edge_uses(const Model<T> &m, ShellId sid) -> std::unordered_map<EdgeId, std::vector<EdgeUse>>
    {
        std::unordered_map<EdgeId, std::vector<EdgeUse>> uses;
        for (const auto &fu : m.shell(sid).faces)
        {
            if (!m.alive(fu.face))
                continue;
            for (auto wid : m.face(fu.face).wires)
            {
                if (!m.alive(wid))
                    continue;
                const auto &w = m.wire(wid);
                for (std::uint32_t i{}; i < w.coedges.size(); ++i)
                {
                    const auto &ce = w.coedges[i];
                    if (!m.alive(ce.edge) || m.edge(ce.edge).degenerate)
                        continue;
                    uses[ce.edge].push_back(EdgeUse{fu.face, CoEdgeRef{wid, i}, compose(ce.orient, fu.orient)});
                }
            }
        }
        return uses;
    }

    /// Every edge is used at most twice.
    template <std::floating_point T>
    [[nodiscard]] bool is_manifold(const Model<T> &m, ShellId sid)
    {
        for (const auto &[eid, u] : shell_edge_uses(m, sid))
            if (u.size() > 2)
                return false;
        return true;
    }

    /// Every (non-degenerate) edge is used exactly twice.
    template <std::floating_point T>
    [[nodiscard]] bool is_closed(const Model<T> &m, ShellId sid)
    {
        auto uses = shell_edge_uses(m, sid);
        if (uses.empty())
            return false;
        for (const auto &[eid, u] : uses)
            if (u.size() != 2)
                return false;
        return true;
    }

    /**
     * @brief True if the shell is consistently oriented: the two usages of
     * every twice-used edge have opposite effective senses. Edges used once
     * are ignored, edges used more than twice make the shell non orientable.
     * A shell with one flipped face fails here; the sewing orientation pass
     * repairs that case, a Möbius strip cannot be repaired.
     */
    template <std::floating_point T>
    [[nodiscard]] bool is_orientable(const Model<T> &m, ShellId sid)
    {
        for (const auto &[eid, u] : shell_edge_uses(m, sid))
        {
            if (u.size() > 2)
                return false;
            if (u.size() == 2 && u[0].orient == u[1].orient)
                return false;
        }
        return true;
    }

    /// Edges used exactly once by the shell, in increasing index order.
    template <std::floating_point T>
    [[nodiscard]] auto free_edges(const Model<T> &m, ShellId sid) -> std::vector<EdgeId>
    {
        std::vector<EdgeId> res;
        for (const auto &[eid, u] : shell_edge_uses(m, sid))
            if (u.size() == 1)
                res.push_back(eid);
        std::sort(res.begin(), res.end());
        return res;
    }

    /// Edges used more than twice by the shell, in increasing index order.
    template <std::floating_point T>
    [[nodiscard]] auto non_manifold_edges(const Model<T> &m, ShellId sid) -> std::vector<EdgeId>
    {
        std::vector<EdgeId> res;
        for (const auto &[eid, u] : shell_edge_uses(m, sid))
            if (u.size() > 2)
                res.push_back(eid);
        std::sort(res.begin(), res.end());
        return res;
    }

    // =========================================================================
    // Bounding box
    // =========================================================================

    template <std::floating_point T>
    struct BoundingBox
    {
        point<T, 3> min{std::numeric_limits<T>::max(), std::numeric_limits<T>::max(), std::numeric_limits<T>::max()};
        point<T, 3> max{std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest()};

        [[nodiscard]] bool empty() const noexcept { return min[0] > max[0]; }

        void add(const point<T, 3> &p, T tol = T(0))
        {
            for (std::size_t i{}; i < 3; ++i)
            {
                min[i] = std::min(min[i], p[i] - tol);
                max[i] = std::max(max[i], p[i] + tol);
            }
        }

        void add(const BoundingBox &b)
        {
            if (b.empty())
                return;
            add(b.min);
            add(b.max);
        }

        /// Grows the box by `d` in every direction.
        void inflate(T d)
        {
            for (std::size_t i{}; i < 3; ++i)
            {
                min[i] -= d;
                max[i] += d;
            }
        }

        [[nodiscard]] bool contains(const point<T, 3> &p) const noexcept
        {
            for (std::size_t i{}; i < 3; ++i)
                if (p[i] < min[i] || p[i] > max[i])
                    return false;
            return true;
        }

        [[nodiscard]] bool intersects(const BoundingBox &b) const noexcept
        {
            if (empty() || b.empty())
                return false;
            for (std::size_t i{}; i < 3; ++i)
                if (b.max[i] < min[i] || b.min[i] > max[i])
                    return false;
            return true;
        }

        [[nodiscard]] auto diagonal() const noexcept -> T { return empty() ? T(0) : norm(max - min); }
    };

    /**
     * @brief Bounding box of a shape: vertices inflated by their tolerance,
     * `n_samples` points on each non-degenerate edge inflated by the edge
     * tolerance, and for a face with natural bounds an `n_samples × n_samples`
     * grid of its surface (a face bounded by its wires only is covered by
     * its edges; a bulging trimmed face may exceed the box, which stays a
     * cheap conservative estimate for the sewing filter).
     */
    template <std::floating_point T>
    [[nodiscard]] auto bounding_box(const Model<T> &m, const ShapeId &shape, std::size_t n_samples = 10) -> BoundingBox<T>
    {
        BoundingBox<T> box;
        n_samples = std::max<std::size_t>(n_samples, 2);

        for (auto vid : explore<VertexId>(m, shape))
        {
            const auto &v = m.vertex(vid);
            box.add(v.pnt, v.tol);
        }
        for (auto eid : explore<EdgeId>(m, shape))
        {
            const auto &e = m.edge(eid);
            if (e.degenerate || !e.curve)
                continue;
            for (std::size_t i{}; i < n_samples; ++i)
            {
                const T u = e.u1 + (e.u2 - e.u1) * T(i) / T(n_samples - 1);
                box.add(e.curve->value(u), e.tol);
            }
        }
        for (auto fid : explore<FaceId>(m, shape))
        {
            const auto &f = m.face(fid);
            if (!f.natural_bounds || !f.surface)
                continue;
            auto [u1, u2, v1, v2] = f.surface->bounds();
            for (std::size_t i{}; i < n_samples; ++i)
                for (std::size_t j{}; j < n_samples; ++j)
                {
                    const T u = u1 + (u2 - u1) * T(i) / T(n_samples - 1);
                    const T v = v1 + (v2 - v1) * T(j) / T(n_samples - 1);
                    box.add(f.surface->value(u, v), f.tol);
                }
        }
        return box;
    }

} // namespace gbs::brep
