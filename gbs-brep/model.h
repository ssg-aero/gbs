#pragma once

/**
 * @file model.h
 * @brief Native BREP core, stage 1: data model.
 *
 * Design: docs/sources/design/brep_core.md, sections 2 to 4.
 *
 * The topology is stored in an arena, `gbs::brep::Model<T>`: one contiguous
 * table per entity type, entities designated by typed handles (`VertexId`,
 * `EdgeId`, ...). Relations between entities are handles, never pointers.
 * Geometry (curves, surfaces) stays outside the arena as `std::shared_ptr`,
 * shared with the rest of gbs.
 *
 * Orientation is carried by usage: a `CoEdge` is an edge seen from a wire
 * (edge + sense + pcurve), a `FaceUse` is a face seen from a shell.
 */

#include <array>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <gbs/gbslib.h>
#include <gbs/gbsconstants.h>
#include <gbs/exceptions.h>
#include <gbs/bscurve.h>
#include <gbs/bssurf.h>

namespace gbs::brep
{
    // =========================================================================
    // Handles
    // =========================================================================

    /// Sentinel index of an invalid handle.
    inline constexpr std::uint32_t npos = std::numeric_limits<std::uint32_t>::max();

    /// Entity kinds, in descending order of the topological hierarchy.
    enum class ShapeType : std::uint8_t
    {
        Vertex,
        Edge,
        Wire,
        Face,
        Shell,
        Solid,
        Compound
    };

    /**
     * @brief Typed handle to an entity of a `Model`: an index into the table
     * of entities of kind `K`. Cheap to copy, compare and hash. Carries no
     * ownership and no lifetime: it is valid as long as the model that issued
     * it has not erased the entity nor been compacted.
     */
    template <ShapeType K>
    struct Handle
    {
        static constexpr ShapeType type = K;
        std::uint32_t index{npos};

        [[nodiscard]] constexpr bool valid() const noexcept { return index != npos; }
        friend constexpr auto operator<=>(const Handle &, const Handle &) = default;
    };

    using VertexId   = Handle<ShapeType::Vertex>;
    using EdgeId     = Handle<ShapeType::Edge>;
    using WireId     = Handle<ShapeType::Wire>;
    using FaceId     = Handle<ShapeType::Face>;
    using ShellId    = Handle<ShapeType::Shell>;
    using SolidId    = Handle<ShapeType::Solid>;
    using CompoundId = Handle<ShapeType::Compound>;

    /// Handle of any kind: what a `Compound` holds and what explorers accept.
    using ShapeId = std::variant<VertexId, EdgeId, WireId, FaceId, ShellId, SolidId, CompoundId>;

    [[nodiscard]] inline constexpr ShapeType shape_type(const ShapeId &id) noexcept
    {
        return std::visit([](auto h) { return decltype(h)::type; }, id);
    }

    [[nodiscard]] inline constexpr std::uint32_t shape_index(const ShapeId &id) noexcept
    {
        return std::visit([](auto h) { return h.index; }, id);
    }

    [[nodiscard]] inline constexpr bool valid(const ShapeId &id) noexcept
    {
        return shape_index(id) != npos;
    }

    // =========================================================================
    // Orientation
    // =========================================================================

    /// Sense of a usage (co-edge in a wire, face in a shell).
    enum class Orientation : std::uint8_t
    {
        Forward,
        Reversed
    };

    [[nodiscard]] inline constexpr Orientation reverse(Orientation o) noexcept
    {
        return o == Orientation::Forward ? Orientation::Reversed : Orientation::Forward;
    }

    /// Composition of two senses: Reversed ∘ Reversed == Forward.
    [[nodiscard]] inline constexpr Orientation compose(Orientation a, Orientation b) noexcept
    {
        return a == b ? Orientation::Forward : Orientation::Reversed;
    }

    // =========================================================================
    // Entities (plain data; behaviour lives in free functions)
    // =========================================================================

    template <std::floating_point T>
    struct Vertex
    {
        point<T, 3> pnt{};
        T tol{brep_default_tolerance<T>}; ///< radius of the confusion ball, > 0
    };

    template <std::floating_point T>
    struct Edge
    {
        std::shared_ptr<Curve<T, 3>> curve; ///< nullptr iff degenerate
        T u1{}, u2{};                       ///< bounds on curve, u1 < u2; the edge runs from v1 to v2 as u grows
        VertexId v1, v2;                    ///< v1 == v2 allowed (closed or degenerate edge)
        T tol{brep_default_tolerance<T>};   ///< radius of the tube around curve
        bool degenerate{false};             ///< edge reduced to a point (pole, apex)
        bool same_parameter{true};          ///< pcurves are parametrized like curve, up to tol
    };

    template <std::floating_point T>
    struct CoEdge
    {
        EdgeId edge;
        Orientation orient{Orientation::Forward}; ///< Forward: traversed from v1 to v2
        std::shared_ptr<Curve<T, 2>> pcurve;      ///< in the (u,v) space of the owning face; nullptr for a free wire
    };

    template <std::floating_point T>
    struct Wire
    {
        std::vector<CoEdge<T>> coedges; ///< ordered, chained end to end (according to each orientation)
        bool closed{false};
    };

    template <std::floating_point T>
    struct Face
    {
        std::shared_ptr<Surface<T, 3>> surface;
        std::vector<WireId> wires;        ///< wires[0] = outer boundary, then holes; all closed
        T tol{brep_default_tolerance<T>};
        bool natural_bounds{false};       ///< wires[0] is the parametric rectangle of surface
    };

    struct FaceUse
    {
        FaceId face;
        Orientation orient{Orientation::Forward}; ///< Forward: the surface normal points outside the shell
    };

    struct Shell
    {
        std::vector<FaceUse> faces;
        bool closed{false}; ///< every edge used exactly twice, in opposite senses
    };

    struct Solid
    {
        ShellId outer;
        std::vector<ShellId> voids; ///< cavities (STEP brep_with_voids); stored only at stage 1
    };

    struct Compound
    {
        std::vector<ShapeId> shapes;
    };

    // =========================================================================
    // Id remapping (result of compact() and append())
    // =========================================================================

    /**
     * @brief Old index -> new index for each entity table. An entry equal to
     * `npos` means the entity was dead (compact) or not copied (append).
     */
    struct IdRemap
    {
        /// One table per ShapeType, indexed by static_cast<std::size_t>(type).
        std::array<std::vector<std::uint32_t>, 7> tables;

        [[nodiscard]] auto table(ShapeType t) -> std::vector<std::uint32_t> & { return tables[static_cast<std::size_t>(t)]; }
        [[nodiscard]] auto table(ShapeType t) const -> const std::vector<std::uint32_t> & { return tables[static_cast<std::size_t>(t)]; }

        template <ShapeType K>
        [[nodiscard]] Handle<K> map(Handle<K> old) const
        {
            const auto &tbl = table(K);
            if (!old.valid() || old.index >= tbl.size())
                return Handle<K>{};
            return Handle<K>{tbl[old.index]};
        }

        [[nodiscard]] ShapeId map(const ShapeId &old) const
        {
            return std::visit([this](auto h) -> ShapeId { return map(h); }, old);
        }
    };

    // =========================================================================
    // Model
    // =========================================================================

    /**
     * @brief Arena owning every topological entity of a BREP model.
     *
     * Entities are addressed by typed handles. `erase()` only marks an entity
     * dead (tombstone) and never invalidates other handles; `compact()`
     * removes the dead entries, renumbers and returns the `IdRemap` to apply
     * to handles held outside the model. `erase()` does not cascade: the
     * caller is responsible for the consistency of the entities that still
     * reference the erased one (builders do it).
     */
    template <std::floating_point T>
    class Model
    {
        template <typename E>
        struct Table
        {
            std::vector<E> items;
            std::vector<std::uint8_t> alive;
            std::size_t n_alive{0};

            auto add(E e) -> std::uint32_t
            {
                items.push_back(std::move(e));
                alive.push_back(1);
                ++n_alive;
                return static_cast<std::uint32_t>(items.size() - 1);
            }
            [[nodiscard]] bool is_alive(std::uint32_t i) const noexcept { return i < items.size() && alive[i]; }
            void erase(std::uint32_t i)
            {
                if (is_alive(i))
                {
                    alive[i] = 0;
                    --n_alive;
                }
            }
            /// new index for each old index (npos for dead), then drops dead items
            auto compact() -> std::vector<std::uint32_t>
            {
                std::vector<std::uint32_t> remap(items.size(), npos);
                std::vector<E> kept;
                kept.reserve(n_alive);
                for (std::size_t i{}; i < items.size(); ++i)
                {
                    if (alive[i])
                    {
                        remap[i] = static_cast<std::uint32_t>(kept.size());
                        kept.push_back(std::move(items[i]));
                    }
                }
                items = std::move(kept);
                alive.assign(items.size(), 1);
                return remap;
            }
        };

        Table<Vertex<T>> m_vertices;
        Table<Edge<T>> m_edges;
        Table<Wire<T>> m_wires;
        Table<Face<T>> m_faces;
        Table<Shell> m_shells;
        Table<Solid> m_solids;
        Table<Compound> m_compounds;

        template <ShapeType K>
        auto table() -> auto &
        {
            if constexpr (K == ShapeType::Vertex) return m_vertices;
            else if constexpr (K == ShapeType::Edge) return m_edges;
            else if constexpr (K == ShapeType::Wire) return m_wires;
            else if constexpr (K == ShapeType::Face) return m_faces;
            else if constexpr (K == ShapeType::Shell) return m_shells;
            else if constexpr (K == ShapeType::Solid) return m_solids;
            else return m_compounds;
        }

        template <ShapeType K>
        auto table() const -> const auto &
        {
            return const_cast<Model *>(this)->template table<K>();
        }

        static const char *type_name(ShapeType t) noexcept
        {
            switch (t)
            {
            case ShapeType::Vertex: return "vertex";
            case ShapeType::Edge: return "edge";
            case ShapeType::Wire: return "wire";
            case ShapeType::Face: return "face";
            case ShapeType::Shell: return "shell";
            case ShapeType::Solid: return "solid";
            default: return "compound";
            }
        }

        template <ShapeType K>
        void check(Handle<K> h) const
        {
            if (!table<K>().is_alive(h.index))
            {
                throw BRepError(std::string("invalid ") + type_name(K) + " handle " +
                                (h.valid() ? std::to_string(h.index) : std::string("npos")));
            }
        }

        // Rewrites every handle stored in the entities through `remap`.
        void apply_remap(const IdRemap &remap)
        {
            for (auto &e : m_edges.items)
            {
                e.v1 = remap.map(e.v1);
                e.v2 = remap.map(e.v2);
            }
            for (auto &w : m_wires.items)
                for (auto &ce : w.coedges)
                    ce.edge = remap.map(ce.edge);
            for (auto &f : m_faces.items)
                for (auto &w : f.wires)
                    w = remap.map(w);
            for (auto &sh : m_shells.items)
                for (auto &fu : sh.faces)
                    fu.face = remap.map(fu.face);
            for (auto &so : m_solids.items)
            {
                so.outer = remap.map(so.outer);
                for (auto &v : so.voids)
                    v = remap.map(v);
            }
            for (auto &c : m_compounds.items)
                for (auto &s : c.shapes)
                    s = remap.map(s);
        }

    public:
        Model() = default;

        // ---- access (throw BRepError on an invalid or dead handle) ---------

        [[nodiscard]] auto vertex(VertexId h) const -> const Vertex<T> & { check(h); return m_vertices.items[h.index]; }
        [[nodiscard]] auto vertex(VertexId h) -> Vertex<T> & { check(h); return m_vertices.items[h.index]; }
        [[nodiscard]] auto edge(EdgeId h) const -> const Edge<T> & { check(h); return m_edges.items[h.index]; }
        [[nodiscard]] auto edge(EdgeId h) -> Edge<T> & { check(h); return m_edges.items[h.index]; }
        [[nodiscard]] auto wire(WireId h) const -> const Wire<T> & { check(h); return m_wires.items[h.index]; }
        [[nodiscard]] auto wire(WireId h) -> Wire<T> & { check(h); return m_wires.items[h.index]; }
        [[nodiscard]] auto face(FaceId h) const -> const Face<T> & { check(h); return m_faces.items[h.index]; }
        [[nodiscard]] auto face(FaceId h) -> Face<T> & { check(h); return m_faces.items[h.index]; }
        [[nodiscard]] auto shell(ShellId h) const -> const Shell & { check(h); return m_shells.items[h.index]; }
        [[nodiscard]] auto shell(ShellId h) -> Shell & { check(h); return m_shells.items[h.index]; }
        [[nodiscard]] auto solid(SolidId h) const -> const Solid & { check(h); return m_solids.items[h.index]; }
        [[nodiscard]] auto solid(SolidId h) -> Solid & { check(h); return m_solids.items[h.index]; }
        [[nodiscard]] auto compound(CompoundId h) const -> const Compound & { check(h); return m_compounds.items[h.index]; }
        [[nodiscard]] auto compound(CompoundId h) -> Compound & { check(h); return m_compounds.items[h.index]; }

        // ---- creation -------------------------------------------------------

        auto add(Vertex<T> v) -> VertexId { return VertexId{m_vertices.add(std::move(v))}; }
        auto add(Edge<T> e) -> EdgeId { return EdgeId{m_edges.add(std::move(e))}; }
        auto add(Wire<T> w) -> WireId { return WireId{m_wires.add(std::move(w))}; }
        auto add(Face<T> f) -> FaceId { return FaceId{m_faces.add(std::move(f))}; }
        auto add(Shell s) -> ShellId { return ShellId{m_shells.add(std::move(s))}; }
        auto add(Solid s) -> SolidId { return SolidId{m_solids.add(std::move(s))}; }
        auto add(Compound c) -> CompoundId { return CompoundId{m_compounds.add(std::move(c))}; }

        // ---- liveness, counts, iteration -----------------------------------

        template <ShapeType K>
        [[nodiscard]] bool alive(Handle<K> h) const noexcept { return table<K>().is_alive(h.index); }

        [[nodiscard]] bool alive(const ShapeId &id) const noexcept
        {
            return std::visit([this](auto h) { return alive(h); }, id);
        }

        /// Number of live entities of kind K.
        template <ShapeType K>
        [[nodiscard]] std::size_t count() const noexcept { return table<K>().n_alive; }

        [[nodiscard]] std::size_t count(ShapeType t) const noexcept
        {
            switch (t)
            {
            case ShapeType::Vertex: return m_vertices.n_alive;
            case ShapeType::Edge: return m_edges.n_alive;
            case ShapeType::Wire: return m_wires.n_alive;
            case ShapeType::Face: return m_faces.n_alive;
            case ShapeType::Shell: return m_shells.n_alive;
            case ShapeType::Solid: return m_solids.n_alive;
            default: return m_compounds.n_alive;
            }
        }

        /// Size of the table of kind K, dead entries included (upper bound of a valid index).
        template <ShapeType K>
        [[nodiscard]] std::size_t capacity() const noexcept { return table<K>().items.size(); }

        /// Handles of the live entities of kind K, in increasing index order.
        template <typename Id>
        [[nodiscard]] auto ids() const -> std::vector<Id>
        {
            const auto &tbl = table<Id::type>();
            std::vector<Id> res;
            res.reserve(tbl.n_alive);
            for (std::uint32_t i{}; i < tbl.items.size(); ++i)
                if (tbl.alive[i])
                    res.push_back(Id{i});
            return res;
        }

        [[nodiscard]] bool empty() const noexcept
        {
            return m_vertices.n_alive + m_edges.n_alive + m_wires.n_alive + m_faces.n_alive +
                       m_shells.n_alive + m_solids.n_alive + m_compounds.n_alive == 0;
        }

        // ---- removal --------------------------------------------------------

        /// Marks the entity dead. Other handles stay valid; no cascade.
        template <ShapeType K>
        void erase(Handle<K> h) { table<K>().erase(h.index); }

        void erase(const ShapeId &id)
        {
            std::visit([this](auto h) { erase(h); }, id);
        }

        /**
         * @brief Removes the dead entries and renumbers every table. Handles
         * stored inside the model are rewritten; handles held outside must be
         * translated with the returned `IdRemap` (dead ones map to invalid).
         */
        auto compact() -> IdRemap
        {
            IdRemap remap;
            remap.table(ShapeType::Vertex) = m_vertices.compact();
            remap.table(ShapeType::Edge) = m_edges.compact();
            remap.table(ShapeType::Wire) = m_wires.compact();
            remap.table(ShapeType::Face) = m_faces.compact();
            remap.table(ShapeType::Shell) = m_shells.compact();
            remap.table(ShapeType::Solid) = m_solids.compact();
            remap.table(ShapeType::Compound) = m_compounds.compact();
            apply_remap(remap);
            return remap;
        }

        /**
         * @brief Copies the live entities of `other` into this model. Geometry
         * is shared (same `shared_ptr`), not duplicated. Returns the map from
         * the handles of `other` to the handles in this model.
         */
        auto append(const Model &other) -> IdRemap
        {
            IdRemap remap;
            auto copy_table = [](const auto &src, auto &dst) {
                std::vector<std::uint32_t> r(src.items.size(), npos);
                for (std::uint32_t i{}; i < src.items.size(); ++i)
                    if (src.alive[i])
                        r[i] = dst.add(src.items[i]);
                return r;
            };
            const auto first_edge = static_cast<std::uint32_t>(m_edges.items.size());
            const auto first_wire = static_cast<std::uint32_t>(m_wires.items.size());
            const auto first_face = static_cast<std::uint32_t>(m_faces.items.size());
            const auto first_shell = static_cast<std::uint32_t>(m_shells.items.size());
            const auto first_solid = static_cast<std::uint32_t>(m_solids.items.size());
            const auto first_compound = static_cast<std::uint32_t>(m_compounds.items.size());

            remap.table(ShapeType::Vertex) = copy_table(other.m_vertices, m_vertices);
            remap.table(ShapeType::Edge) = copy_table(other.m_edges, m_edges);
            remap.table(ShapeType::Wire) = copy_table(other.m_wires, m_wires);
            remap.table(ShapeType::Face) = copy_table(other.m_faces, m_faces);
            remap.table(ShapeType::Shell) = copy_table(other.m_shells, m_shells);
            remap.table(ShapeType::Solid) = copy_table(other.m_solids, m_solids);
            remap.table(ShapeType::Compound) = copy_table(other.m_compounds, m_compounds);

            // Rewrite the handles of the copied entities only.
            for (auto i = first_edge; i < m_edges.items.size(); ++i)
            {
                auto &e = m_edges.items[i];
                e.v1 = remap.map(e.v1);
                e.v2 = remap.map(e.v2);
            }
            for (auto i = first_wire; i < m_wires.items.size(); ++i)
                for (auto &ce : m_wires.items[i].coedges)
                    ce.edge = remap.map(ce.edge);
            for (auto i = first_face; i < m_faces.items.size(); ++i)
                for (auto &w : m_faces.items[i].wires)
                    w = remap.map(w);
            for (auto i = first_shell; i < m_shells.items.size(); ++i)
                for (auto &fu : m_shells.items[i].faces)
                    fu.face = remap.map(fu.face);
            for (auto i = first_solid; i < m_solids.items.size(); ++i)
            {
                auto &so = m_solids.items[i];
                so.outer = remap.map(so.outer);
                for (auto &v : so.voids)
                    v = remap.map(v);
            }
            for (auto i = first_compound; i < m_compounds.items.size(); ++i)
                for (auto &s : m_compounds.items[i].shapes)
                    s = remap.map(s);
            return remap;
        }
    };

    // =========================================================================
    // Small helpers on usages
    // =========================================================================

    /// Vertex where the co-edge starts, its orientation taken into account.
    template <std::floating_point T>
    [[nodiscard]] auto coedge_start(const Model<T> &m, const CoEdge<T> &ce) -> VertexId
    {
        const auto &e = m.edge(ce.edge);
        return ce.orient == Orientation::Forward ? e.v1 : e.v2;
    }

    /// Vertex where the co-edge ends, its orientation taken into account.
    template <std::floating_point T>
    [[nodiscard]] auto coedge_end(const Model<T> &m, const CoEdge<T> &ce) -> VertexId
    {
        const auto &e = m.edge(ce.edge);
        return ce.orient == Orientation::Forward ? e.v2 : e.v1;
    }

    /// Point of the edge's 3D curve at parameter u (throws on a degenerate edge).
    template <std::floating_point T>
    [[nodiscard]] auto edge_point(const Model<T> &m, EdgeId h, T u) -> point<T, 3>
    {
        const auto &e = m.edge(h);
        if (e.degenerate || !e.curve)
            return m.vertex(e.v1).pnt;
        return e.curve->value(u);
    }

    /// Point of a co-edge's pcurve at parameter t of the edge, orientation applied.
    /// A Reversed co-edge is traversed from u2 down to u1 as t grows from u1 to u2.
    template <std::floating_point T>
    [[nodiscard]] auto coedge_uv(const Model<T> &m, const CoEdge<T> &ce, T t) -> point<T, 2>
    {
        if (!ce.pcurve)
            throw BRepError("co-edge without pcurve");
        const auto &e = m.edge(ce.edge);
        const T u = ce.orient == Orientation::Forward ? t : e.u1 + e.u2 - t;
        return ce.pcurve->value(u);
    }

} // namespace gbs::brep

// Hash support so handles can key unordered containers.
template <gbs::brep::ShapeType K>
struct std::hash<gbs::brep::Handle<K>>
{
    std::size_t operator()(const gbs::brep::Handle<K> &h) const noexcept
    {
        return std::hash<std::uint32_t>{}(h.index) ^ (static_cast<std::size_t>(K) << 29);
    }
};
