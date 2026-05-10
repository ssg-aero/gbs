#pragma once

#include "halfEdgeMeshData.h"
#include "halfEdgeMeshGetters.h"
#include "halfEdgeMeshEditors.h"
#include "baseGeom.h"
#include "baseIntersection.h"

#include <list>
#include <memory>

namespace gbs
{

/**
 * @brief Recovers the constrained edge (v_p, v_q) in a 2D triangulation by
 *        flipping the diagonals of triangles crossed by the segment p-q.
 *
 * Implements the Anglada / Sloan edge-recovery step of a Constrained
 * Delaunay Triangulation. Assumes the input triangulation already contains
 * v_p and v_q as vertices, and that the segment p-q does not pass through
 * any other vertex of the triangulation (no Steiner-point handling here).
 *
 * @tparam T Floating-point type.
 * @param faces_lst Triangulation face list (will be mutated in place).
 * @param v_p Origin vertex of the constrained edge.
 * @param v_q End vertex of the constrained edge.
 * @param tol Tolerance forwarded to seg_seg_strict_intersection / orient_2d.
 * @return Shared pointer to the half-edge from v_p to v_q after recovery,
 *         or nullptr on failure (boundary hit, degenerate input, no progress).
 */
    template <std::floating_point T>
    auto recoverEdge(
        auto& faces_lst,
        const std::shared_ptr<HalfEdgeVertex<T, 2>>& v_p,
        const std::shared_ptr<HalfEdgeVertex<T, 2>>& v_q,
        T tol = T{1e-10}) -> std::shared_ptr<HalfEdge<T, 2>>
    {
        if (!v_p || !v_q || v_p == v_q)
            return nullptr;

        // Edge already present? Cheap fast path.
        if (auto e = findHalfEdge(v_p, v_q))
            return e;
        if (auto e = findHalfEdge(v_q, v_p))
            return e->opposite ? e->opposite : nullptr;

        const auto& p = v_p->coords;
        const auto& q = v_q->coords;

        // Step 1: collect the chain of half-edges of `faces_lst` whose interior
        // is strictly crossed by segment p-q. We start from a face attached to
        // v_p whose opposite-to-v_p edge crosses p-q, then walk through
        // successive faces by following the opposite half-edge.
        std::list<std::shared_ptr<HalfEdge<T, 2>>> crossing;

        for (const auto& f : getFacesAttachedToVertex(v_p))
        {
            // The triangle's three edges; pick the one not incident to v_p.
            for (const auto& e : getTriangleEdges(f))
            {
                if (e->vertex == v_p || e->previous->vertex == v_p)
                    continue;
                const auto& a = e->previous->vertex->coords;
                const auto& b = e->vertex->coords;
                if (seg_seg_strict_intersection(p, q, a, b, tol))
                {
                    crossing.push_back(e);
                    break;
                }
            }
            if (!crossing.empty())
                break;
        }

        if (crossing.empty())
            return nullptr;

        // Walk the chain across opposite half-edges until we reach v_q.
        // The cap is generous enough for any reasonable polygon while still
        // failing fast when the navigation has broken.
        constexpr int max_chain_steps = 4096;
        for (int step = 0; step < max_chain_steps; ++step)
        {
            auto e = crossing.back();
            if (!e->opposite)
                return nullptr; // chain hit the mesh boundary

            auto e_in = e->opposite;        // same edge, viewed from the next face
            auto apex = e_in->next->vertex; // vertex of next face not on the shared edge

            if (apex == v_q)
                break; // segment p-q ends at this face's apex — chain complete

            auto e_left  = e_in->next;     // from v_a (= e->vertex) to apex
            auto e_right = e_in->previous; // from apex to v_b (= e->previous->vertex)

            const auto& la = e_left->previous->vertex->coords;
            const auto& lb = e_left->vertex->coords;
            const auto& ra = e_right->previous->vertex->coords;
            const auto& rb = e_right->vertex->coords;

            if (seg_seg_strict_intersection(p, q, la, lb, tol))
                crossing.push_back(e_left);
            else if (seg_seg_strict_intersection(p, q, ra, rb, tol))
                crossing.push_back(e_right);
            else
                return nullptr; // p-q passes through a vertex (unsupported here)
        }

        // Step 2: Anglada — flip each crossing diagonal whose quadrilateral is
        // convex; non-convex quads are skipped and revisited on the next pass.
        // Each successful flip either yields the target edge p-q, removes the
        // diagonal from the crossing list (if it no longer crosses p-q), or
        // replaces the diagonal in-place (the recycled half-edge is now the new
        // diagonal between the two former apexes).
        std::shared_ptr<HalfEdge<T, 2>> recovered{};

        const std::size_t hard_cap = crossing.size() * crossing.size() + 8;
        for (std::size_t pass = 0; !crossing.empty() && pass < hard_cap; ++pass)
        {
            bool any_flipped = false;
            auto it = crossing.begin();
            while (it != crossing.end())
            {
                auto e = *it;
                if (!e->opposite)
                {
                    it = crossing.erase(it);
                    continue;
                }
                auto f1 = e->face;
                auto f2 = e->opposite->face;

                // Quad corners: shared edge from v_a to v_b, apex of f1, apex of f2.
                auto v_a     = e->previous->vertex;
                auto v_b     = e->vertex;
                auto v_apex1 = e->next->vertex;
                auto v_apex2 = e->opposite->next->vertex;

                // Convex iff the alternative splitting (apex1 - apex2) lies inside
                // the quad: both half-quads must be CCW.
                const auto& A = v_a->coords;
                const auto& B = v_b->coords;
                const auto& C = v_apex1->coords;
                const auto& D = v_apex2->coords;
                if (orient_2d(C, B, D) <= tol || orient_2d(D, A, C) <= tol)
                {
                    ++it;
                    continue; // non-convex or degenerate, defer
                }

                flip(f1, f2);

                // After flip, e is recycled into the new diagonal between the apexes.
                auto na = e->previous->vertex;
                auto nb = e->vertex;

                if ((na == v_p && nb == v_q) || (na == v_q && nb == v_p))
                {
                    recovered = (na == v_p) ? e : e->opposite;
                    it = crossing.erase(it);
                }
                else if (seg_seg_strict_intersection(
                             p, q, na->coords, nb->coords, tol))
                {
                    ++it; // still crossing — keep, will retry
                }
                else
                {
                    it = crossing.erase(it);
                }
                any_flipped = true;
            }
            if (!any_flipped)
                return nullptr; // stuck — every remaining quad is non-convex
        }

        if (!recovered)
            recovered = findHalfEdge(v_p, v_q);

        return recovered;
    }

} // namespace gbs
