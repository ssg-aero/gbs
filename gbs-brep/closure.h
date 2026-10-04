#pragma once

/**
 * @file closure.h
 * @brief Native BREP core, stage 1: geometric queries used by the face
 * builders: closure / degeneracy of a surface's boundary, signed area of a
 * wire in a face's parametric space.
 *
 * Design: docs/sources/design/brep_core.md, sections 3.3 and 6.1;
 * architecture note docs/sources/design/brep_pr04_face_natural.md.
 */

#include <algorithm>
#include <array>
#include <cstddef>

#include <gbs-brep/model.h>

namespace gbs::brep
{
    /**
     * @brief How the parametric rectangle of a surface closes on itself.
     *
     * `closed_u`: S(u1, v) == S(u2, v) for every v (seam along u = u1).
     * `degenerate_u1`: the iso u = u1 is a single point (pole, apex); same for
     * the three other sides.
     */
    struct SurfaceClosure
    {
        bool closed_u{false}, closed_v{false};
        bool degenerate_u1{false}, degenerate_u2{false};
        bool degenerate_v1{false}, degenerate_v2{false};

        [[nodiscard]] constexpr bool degenerate() const noexcept
        {
            return degenerate_u1 || degenerate_u2 || degenerate_v1 || degenerate_v2;
        }
    };

    /**
     * @brief Detects closure and degenerate sides numerically, by sampling
     * the four boundary isos at `n_samples` points and comparing within `tol`.
     * Works for any `Surface<T,3>` (NURBS, revolution, offset…) without
     * relying on the surface type; a side detected as degenerate is never
     * reported as a seam.
     */
    template <std::floating_point T>
    [[nodiscard]] auto surface_closure(const Surface<T, 3> &srf, T tol, std::size_t n_samples = 9) -> SurfaceClosure
    {
        n_samples = std::max<std::size_t>(n_samples, 3);
        const auto [u1, u2, v1, v2] = srf.bounds();
        auto param = [n_samples](T a, T b, std::size_t i) { return a + (b - a) * T(i) / T(n_samples - 1); };

        T d_closed_u{}, d_closed_v{}, d_u1{}, d_u2{}, d_v1{}, d_v2{};
        const auto p_u1v1 = srf.value(u1, v1), p_u2v1 = srf.value(u2, v1), p_u1v2 = srf.value(u1, v2);
        for (std::size_t i{}; i < n_samples; ++i)
        {
            const T v = param(v1, v2, i);
            const T u = param(u1, u2, i);
            const auto a = srf.value(u1, v), b = srf.value(u2, v);
            const auto c = srf.value(u, v1), d = srf.value(u, v2);
            d_closed_u = std::max(d_closed_u, distance(a, b));
            d_closed_v = std::max(d_closed_v, distance(c, d));
            d_u1 = std::max(d_u1, distance(a, p_u1v1));
            d_u2 = std::max(d_u2, distance(b, p_u2v1));
            d_v1 = std::max(d_v1, distance(c, p_u1v1));
            d_v2 = std::max(d_v2, distance(d, p_u1v2));
        }
        SurfaceClosure c;
        c.degenerate_u1 = d_u1 <= tol;
        c.degenerate_u2 = d_u2 <= tol;
        c.degenerate_v1 = d_v1 <= tol;
        c.degenerate_v2 = d_v2 <= tol;
        c.closed_u = d_closed_u <= tol && !c.degenerate_u1 && !c.degenerate_u2;
        c.closed_v = d_closed_v <= tol && !c.degenerate_v1 && !c.degenerate_v2;
        return c;
    }

    /**
     * @brief Signed area enclosed by a wire in the (u,v) space of its face,
     * by the shoelace formula on `n_per_coedge` points per co-edge, the
     * co-edge senses applied. Positive for a counter-clockwise (outer)
     * boundary, negative for a clockwise one (hole). Every co-edge must carry
     * a pcurve.
     */
    template <std::floating_point T>
    [[nodiscard]] auto uv_signed_area(const Model<T> &m, WireId wid, std::size_t n_per_coedge = 16) -> T
    {
        n_per_coedge = std::max<std::size_t>(n_per_coedge, 2);
        std::vector<point<T, 2>> poly;
        for (const auto &ce : m.wire(wid).coedges)
        {
            const auto &e = m.edge(ce.edge);
            for (std::size_t i{}; i + 1 < n_per_coedge; ++i) // the last point is the next co-edge's first
                poly.push_back(coedge_uv(m, ce, e.u1 + (e.u2 - e.u1) * T(i) / T(n_per_coedge - 1)));
        }
        T a{};
        for (std::size_t i{}; i < poly.size(); ++i)
        {
            const auto &p = poly[i];
            const auto &q = poly[(i + 1) % poly.size()];
            a += p[0] * q[1] - q[0] * p[1];
        }
        return a / T(2);
    }

} // namespace gbs::brep
