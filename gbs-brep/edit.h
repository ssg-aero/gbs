#pragma once

/**
 * @file edit.h
 * @brief Native BREP core: local edits of an existing model (first operation
 * of stage 2, needed by the STEP reader to insert seams).
 *
 * Architecture note: docs/sources/design/step_pr07_seam.md.
 */

#include <gbs-brep/builders.h>

#include <utility>
#include <vector>

namespace gbs::brep
{
    /// Result of split_edge.
    struct SplitEdge
    {
        EdgeId first;   ///< the original edge, now from v1 to the new vertex
        EdgeId second;  ///< new edge from the new vertex to the former v2
        VertexId vertex; ///< new vertex at the split point
    };

    /**
     * @brief Splits the edge `e` at its parameter `u` (u1 < u < u2).
     *
     * The edge keeps its id and becomes [u1, u]; a new edge covers [u, u2];
     * both share the 3D curve and meet at a new vertex (tolerance of the edge).
     * Every wire of the model using `e` is updated: a Forward co-edge becomes
     * (first, second), a Reversed one (second, first), both keeping the
     * co-edge's pcurve, which is parametrized like the edge. Name and external
     * id are copied to the new edge. A closed edge becomes two edges between
     * its vertex and the new one.
     * On failure the model is unchanged.
     */
    template <std::floating_point T>
    [[nodiscard]] auto split_edge(Model<T> &m, EdgeId e, std::type_identity_t<T> u) -> BuildResult<SplitEdge>
    {
        if (!m.alive(e))
            return detail::fail(BuildErrc::InvalidId, "split_edge: dead or invalid edge", {e});
        const auto ed = m.edge(e);
        if (ed.degenerate || !ed.curve)
            return detail::fail(BuildErrc::DegenerateEdge, "split_edge: degenerate edge", {e});
        const T eps = T(1e-9) * (ed.u2 - ed.u1);
        if (!std::isfinite(u) || !(u > ed.u1 + eps) || !(u < ed.u2 - eps))
            return detail::fail(BuildErrc::InvalidBounds, "split_edge: parameter not strictly inside the edge", {e});
        const auto p = ed.curve->value(u);
        if (!detail::finite(p))
            return detail::fail(BuildErrc::InvalidPoint, "split_edge: curve evaluates to a non finite point", {e});

        // ---- write
        const auto v = m.add(Vertex<T>{p, ed.tol});
        auto &first = m.edge(e);
        first.u2 = u;
        first.v2 = v;
        const auto second = m.add(Edge<T>{ed.curve, u, ed.u2, v, ed.v2, ed.tol});
        if (auto n = m.name(e); !n.empty())
            m.setName(second, std::string(n));
        m.setExternalId(second, m.externalId(e));

        for (auto wid : m.template ids<WireId>())
        {
            auto &cos = m.wire(wid).coedges;
            for (std::size_t i = 0; i < cos.size(); ++i)
            {
                if (cos[i].edge != e)
                    continue;
                const auto ce = cos[i];
                const CoEdge<T> half{second, ce.orient, ce.pcurve};
                // Forward: e then second ; Reversed: second then e
                cos.insert(cos.begin() + static_cast<std::ptrdiff_t>(ce.orient == Orientation::Forward ? i + 1 : i), half);
                ++i;
            }
        }
        return SplitEdge{e, second, v};
    }
}
