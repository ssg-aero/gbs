#pragma once

/**
 * @file check.h
 * @brief Native BREP core, stage 1: minimal validity check ("BRepCheck-lite").
 *
 * Design: docs/sources/design/brep_core.md, sections 2.6 and 5.4;
 * architecture note docs/sources/design/brep_pr06_check.md.
 *
 * `check(model, shape)` visits every entity under `shape` once and reports
 * the violated invariants of the model, without stopping at the first one.
 * It never modifies the model and never throws on an inconsistent model
 * (dead references are reported, not followed).
 *
 * Out of scope (on purpose): self-intersection of a wire in (u,v),
 * intersection between faces, validity of the surfaces themselves.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <utility>
#include <string>
#include <vector>

#include <gbs-brep/model.h>
#include <gbs-brep/explore.h>
#include <gbs-brep/closure.h>

namespace gbs::brep
{
    /// A violated invariant.
    enum class Issue : std::uint8_t
    {
        DeadReference,          ///< an entity refers to a dead or invalid id
        InvalidTolerance,       ///< tolerance not finite or not > 0
        InvalidPoint,           ///< vertex point not finite
        EdgeWithoutCurve,       ///< non-degenerate edge with no 3D curve
        InvalidDegenerateEdge,  ///< degenerate edge with a curve or with two different vertices
        EdgeBoundsInverted,     ///< u1 >= u2
        VertexOffCurve,         ///< curve end outside its vertex tolerance ball
        ToleranceOrderViolated, ///< vertex.tol < edge.tol or edge.tol < face.tol
        WireNotChained,         ///< consecutive co-edges do not share a vertex
        WireNotClosed,          ///< closed flag set (or face boundary) but the chain does not close
        FaceWithoutSurface,     ///< no surface
        FaceWithoutBoundary,    ///< no wire
        CoEdgeWithoutPCurve,    ///< co-edge of a face wire without pcurve
        PCurveRangeMismatch,    ///< pcurve parameter range does not cover the edge range
        SameParameterViolated,  ///< S(pcurve(t)) farther than edge.tol from C(t)
        OuterWireNotDirect,     ///< outer boundary not counter-clockwise in (u,v)
        InnerWireNotIndirect,   ///< hole not clockwise in (u,v)
        InnerWireOutsideOuter,  ///< hole not inside the outer boundary
        InnerWiresOverlap,      ///< two holes overlap
        NonManifoldEdge,        ///< edge used by more than two co-edges of a shell
        ShellNotClosed,         ///< closed flag set, or outer/void shell of a solid, but some edge is used once
        ShellNotOrientable,     ///< the two uses of an edge have the same effective sense
        SolidNotOutward,        ///< the outer shell of a solid has a negative signed volume
        VoidNotInward,          ///< a cavity shell of a solid has a positive signed volume
    };

    [[nodiscard]] inline constexpr const char *to_string(Issue i) noexcept
    {
        switch (i)
        {
        case Issue::DeadReference: return "dead reference";
        case Issue::InvalidTolerance: return "invalid tolerance";
        case Issue::InvalidPoint: return "invalid point";
        case Issue::EdgeWithoutCurve: return "edge without curve";
        case Issue::InvalidDegenerateEdge: return "invalid degenerate edge";
        case Issue::EdgeBoundsInverted: return "edge bounds inverted";
        case Issue::VertexOffCurve: return "vertex off curve";
        case Issue::ToleranceOrderViolated: return "tolerance order violated";
        case Issue::WireNotChained: return "wire not chained";
        case Issue::WireNotClosed: return "wire not closed";
        case Issue::FaceWithoutSurface: return "face without surface";
        case Issue::FaceWithoutBoundary: return "face without boundary";
        case Issue::CoEdgeWithoutPCurve: return "co-edge without pcurve";
        case Issue::PCurveRangeMismatch: return "pcurve range mismatch";
        case Issue::SameParameterViolated: return "same parameter violated";
        case Issue::OuterWireNotDirect: return "outer wire not direct";
        case Issue::InnerWireNotIndirect: return "inner wire not indirect";
        case Issue::InnerWireOutsideOuter: return "inner wire outside outer";
        case Issue::InnerWiresOverlap: return "inner wires overlap";
        case Issue::NonManifoldEdge: return "non-manifold edge";
        case Issue::ShellNotClosed: return "shell not closed";
        case Issue::ShellNotOrientable: return "shell not orientable";
        case Issue::SolidNotOutward: return "solid not outward";
        case Issue::VoidNotInward: return "void not inward";
        }
        std::unreachable();
    }

    /// One finding: the entity at fault, the issue, a human-readable detail.
    struct CheckEntry
    {
        ShapeId shape;
        Issue issue;
        std::string detail;
    };

    struct CheckReport
    {
        std::vector<CheckEntry> entries;

        [[nodiscard]] bool ok() const noexcept { return entries.empty(); }
        [[nodiscard]] std::size_t count(Issue i) const
        {
            return static_cast<std::size_t>(std::ranges::count(entries, i, &CheckEntry::issue));
        }
        [[nodiscard]] bool has(Issue i, const ShapeId &s) const
        {
            return std::ranges::any_of(entries, [&](const CheckEntry &e) { return e.issue == i && e.shape == s; });
        }
    };

    struct CheckOptions
    {
        std::size_t n_samples = 9;    ///< parameters per co-edge for SameParameter
        std::size_t n_uv_polygon = 16; ///< points per co-edge for the (u,v) orientation and inclusion tests
        bool geometry = true;         ///< false: topological checks only (no curve / surface evaluation)
    };

    namespace detail
    {
        template <std::floating_point T>
        struct Checker
        {
            const Model<T> &m;
            CheckOptions opts;
            CheckReport report;

            void add(const ShapeId &s, Issue i, std::string detail = {})
            {
                report.entries.push_back(CheckEntry{s, i, std::move(detail)});
            }

            static bool bad_tol(T t) { return !std::isfinite(t) || !(t > T(0)); }
            static bool within(T d, T tol) { return d <= tol * (T(1) + T(1e-9)) + std::numeric_limits<T>::min(); }

            void vertex(VertexId id)
            {
                const auto &v = m.vertex(id);
                if (bad_tol(v.tol))
                    add(id, Issue::InvalidTolerance);
                if (!std::ranges::all_of(v.pnt, [](T x) { return std::isfinite(x); }))
                    add(id, Issue::InvalidPoint);
            }

            void edge(EdgeId id)
            {
                const auto &e = m.edge(id);
                if (bad_tol(e.tol))
                    add(id, Issue::InvalidTolerance);
                if (!(e.u1 < e.u2))
                    add(id, Issue::EdgeBoundsInverted);
                if (!m.alive(e.v1) || !m.alive(e.v2))
                {
                    add(id, Issue::DeadReference, "vertex");
                    return;
                }
                for (auto v : {e.v1, e.v2})
                    if (m.vertex(v).tol < e.tol)
                        add(id, Issue::ToleranceOrderViolated, "vertex tolerance below edge tolerance");
                if (e.degenerate)
                {
                    if (e.curve || e.v1 != e.v2)
                        add(id, Issue::InvalidDegenerateEdge);
                    return;
                }
                if (!e.curve)
                {
                    add(id, Issue::EdgeWithoutCurve);
                    return;
                }
                if (!opts.geometry || !(e.u1 < e.u2))
                    return;
                for (bool first : {true, false})
                {
                    const auto &v = m.vertex(first ? e.v1 : e.v2);
                    const T d = distance(e.curve->value(first ? e.u1 : e.u2), v.pnt);
                    if (!within(d, v.tol))
                        add(id, Issue::VertexOffCurve, (first ? "v1 at " : "v2 at ") + std::to_string(d));
                }
            }

            /// Returns false when the wire refers to dead edges (no further check possible).
            bool wire(WireId id, bool must_close)
            {
                const auto &w = m.wire(id);
                for (const auto &ce : w.coedges)
                    if (!m.alive(ce.edge))
                    {
                        add(id, Issue::DeadReference, "edge");
                        return false;
                    }
                if (w.coedges.empty() || !is_chained(m, id))
                    add(id, Issue::WireNotChained);
                else if ((w.closed || must_close) && !is_closed(m, id))
                    add(id, Issue::WireNotClosed);
                return true;
            }

            void face(FaceId id)
            {
                const auto &f = m.face(id);
                if (bad_tol(f.tol))
                    add(id, Issue::InvalidTolerance);
                if (!f.surface)
                    add(id, Issue::FaceWithoutSurface);
                if (f.wires.empty())
                    add(id, Issue::FaceWithoutBoundary);

                bool all_pcurves = !f.wires.empty();
                for (auto wid : f.wires)
                {
                    if (!m.alive(wid))
                    {
                        add(id, Issue::DeadReference, "wire");
                        all_pcurves = false;
                        continue;
                    }
                    if (!wire(wid, true))
                    {
                        all_pcurves = false;
                        continue;
                    }
                    for (const auto &ce : m.wire(wid).coedges)
                    {
                        const auto &e = m.edge(ce.edge);
                        if (e.tol < f.tol)
                            add(ce.edge, Issue::ToleranceOrderViolated, "edge tolerance below face tolerance");
                        if (!ce.pcurve)
                        {
                            add(wid, Issue::CoEdgeWithoutPCurve);
                            all_pcurves = false;
                            continue;
                        }
                        const auto [p1, p2] = ce.pcurve->bounds();
                        if (p1 > e.u1 + knot_eps<T> || p2 < e.u2 - knot_eps<T>)
                        {
                            add(wid, Issue::PCurveRangeMismatch);
                            all_pcurves = false;
                            continue;
                        }
                        if (opts.geometry && f.surface && (e.u1 < e.u2))
                            same_parameter(wid, ce, *f.surface);
                    }
                }
                if (opts.geometry && all_pcurves && f.surface)
                    orientation(id);
            }

            void same_parameter(WireId wid, const CoEdge<T> &ce, const Surface<T, 3> &srf)
            {
                const auto &e = m.edge(ce.edge);
                const std::size_t n = std::max<std::size_t>(opts.n_samples, 2);
                T worst{};
                for (std::size_t i{}; i < n; ++i)
                {
                    const T t = e.u1 + (e.u2 - e.u1) * T(i) / T(n - 1);
                    const auto uv = ce.pcurve->value(t);
                    const auto p3 = e.degenerate || !e.curve ? m.vertex(e.v1).pnt : e.curve->value(t);
                    worst = std::max(worst, distance(srf.value(uv[0], uv[1]), p3));
                }
                if (!within(worst, e.tol))
                    add(wid, Issue::SameParameterViolated, "edge " + std::to_string(ce.edge.index) + " deviates " + std::to_string(worst));
            }

            void orientation(FaceId id)
            {
                const auto &f = m.face(id);
                std::vector<std::vector<point<T, 2>>> polys;
                for (auto wid : f.wires)
                    polys.push_back(uv_polygon(m, std::span<const CoEdge<T>>{m.wire(wid).coedges}, opts.n_uv_polygon));
                if (!(signed_area(polys[0]) > T(0)))
                    add(f.wires[0], Issue::OuterWireNotDirect);
                for (std::size_t i{1}; i < polys.size(); ++i)
                {
                    if (!(signed_area(polys[i]) < T(0)))
                        add(f.wires[i], Issue::InnerWireNotIndirect);
                    if (!std::ranges::all_of(polys[i], [&](const point<T, 2> &p) { return inside(polys[0], p); }))
                        add(f.wires[i], Issue::InnerWireOutsideOuter);
                    for (std::size_t j{1}; j < i; ++j)
                        if (std::ranges::any_of(polys[i], [&](const point<T, 2> &p) { return inside(polys[j], p); }) ||
                            std::ranges::any_of(polys[j], [&](const point<T, 2> &p) { return inside(polys[i], p); }))
                            add(f.wires[i], Issue::InnerWiresOverlap, "with wire " + std::to_string(f.wires[j].index));
                }
            }

            void shell(ShellId id, bool must_close)
            {
                const auto &sh = m.shell(id);
                for (const auto &fu : sh.faces)
                    if (!m.alive(fu.face))
                    {
                        add(id, Issue::DeadReference, "face");
                        return;
                    }
                for (const auto &[eid, uses] : shell_edge_uses(m, id))
                {
                    if (uses.size() > 2)
                        add(eid, Issue::NonManifoldEdge, std::to_string(uses.size()) + " uses");
                    else if (uses.size() == 2 && uses[0].orient == uses[1].orient)
                        add(id, Issue::ShellNotOrientable, "edge " + std::to_string(eid.index));
                }
                if ((sh.closed || must_close) && !is_closed(m, id))
                    add(id, Issue::ShellNotClosed);
            }

            /// Orientation of a solid's shells, only when they are valid enough to integrate.
            void solid(SolidId id)
            {
                const auto &so = m.solid(id);
                auto usable = [&](ShellId s) {
                    if (!m.alive(s) || !is_closed(m, s) || !is_orientable(m, s))
                        return false;
                    for (const auto &fu : m.shell(s).faces)
                    {
                        if (!m.alive(fu.face) || !m.face(fu.face).surface)
                            return false;
                        for (auto w : m.face(fu.face).wires)
                            if (!m.alive(w) || std::ranges::any_of(m.wire(w).coedges, [&](const CoEdge<T> &ce) { return !m.alive(ce.edge) || !ce.pcurve; }))
                                return false;
                    }
                    return true;
                };
                if (usable(so.outer) && !(signed_volume(m, so.outer) > T(0)))
                    add(id, Issue::SolidNotOutward);
                for (auto v : so.voids)
                    if (usable(v) && !(signed_volume(m, v) < T(0)))
                        add(v, Issue::VoidNotInward);
            }
        };
    } // namespace detail

    /**
     * @brief Checks every entity under `shape` (each once) and reports all the
     * violated invariants. A shell is required to be closed when its `closed`
     * flag is set or when it bounds a solid; a wire when its `closed` flag is
     * set or when it bounds a face.
     */
    template <std::floating_point T>
    [[nodiscard]] auto check(const Model<T> &m, const ShapeId &shape, CheckOptions opts = {}) -> CheckReport
    {
        detail::Checker<T> c{m, opts, {}};
        if (!m.alive(shape))
        {
            c.add(shape, Issue::DeadReference, "root");
            return c.report;
        }
        // shells bounding a solid, and wires bounding a face, must close
        std::vector<ShellId> solid_shells;
        for (auto sid : explore<SolidId>(m, shape))
        {
            const auto &so = m.solid(sid);
            if (!m.alive(so.outer))
                c.add(sid, Issue::DeadReference, "outer shell");
            else
                solid_shells.push_back(so.outer);
            for (auto v : so.voids)
                if (!m.alive(v))
                    c.add(sid, Issue::DeadReference, "void shell");
                else
                    solid_shells.push_back(v);
        }
        for (auto cid : explore<CompoundId>(m, shape))
            for (const auto &s : m.compound(cid).shapes)
                if (!m.alive(s))
                    c.add(cid, Issue::DeadReference, "compound member");

        for (auto vid : explore<VertexId>(m, shape))
            c.vertex(vid);
        for (auto eid : explore<EdgeId>(m, shape))
            c.edge(eid);
        const auto faces = explore<FaceId>(m, shape);
        std::vector<WireId> face_wires;
        for (auto fid : faces)
        {
            c.face(fid);
            for (auto w : m.face(fid).wires)
                face_wires.push_back(w);
        }
        for (auto wid : explore<WireId>(m, shape))
            if (std::ranges::find(face_wires, wid) == face_wires.end())
                c.wire(wid, false); // free wire
        for (auto sid : explore<ShellId>(m, shape))
            c.shell(sid, std::ranges::find(solid_shells, sid) != solid_shells.end());
        if (opts.geometry)
            for (auto sid : explore<SolidId>(m, shape))
                c.solid(sid);
        return c.report;
    }

} // namespace gbs::brep
