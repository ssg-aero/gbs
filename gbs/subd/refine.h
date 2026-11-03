#pragma once
/**
 * @file refine.h
 * @brief One step of Catmull-Clark subdivision of a cage, with the rules of OpenSubdiv's
 * Catmark scheme for infinitely sharp edges and corners.
 *
 * - face point: centroid of the face;
 * - edge point: (v₀ + v₁ + f₀ + f₁) / 4 for a smooth edge, its middle for a sharp one (crease
 *   or boundary);
 * - vertex point: ((n - 2) V + Σ eⱼ / n + Σ fⱼ / n) / n for a smooth vertex or a dart of n
 *   edges (eⱼ its neighbours, fⱼ the face points around it), (P₋ + 6 V + P₊) / 8 on a crease,
 *   V at a corner.
 *
 * The refined cage lists the vertex points first (vertex v stays v), then the face points, then
 * the edge points, in the order of `Topology`; each face of n sides gives n quadrangles, one per
 * corner, in the order of OpenSubdiv: the quadrangle of corner i of an n-gon is (vertex point of
 * i, edge point of side i, face point, edge point of side i - 1), its (u, v) from vertex i; the
 * quadrangle of corner i of a quadrangle starts at its corner of lowest (u, v) in the parent, so
 * that it runs along the parent's directions, the child of corner i at its position i. The halves
 * of a crease are creases, the corners stay corners. Design: docs/sources/design/subd_pr01_cage.md.
 */
#include <gbs/subd/cage.h>

#ifdef GBS_USE_MODULES
    import vecop;
#else
    #include <gbs/vecop.ixx>
#endif

#include <array>
#include <cstddef>
#include <vector>

namespace gbs::subd
{
    /// Child of vertex v, face f, edge e in a cage refined once (`refine`).
    struct ChildIndices
    {
        std::size_t n_vertices{}, n_faces{};
        auto of_vertex(std::size_t v) const -> std::size_t { return v; }
        auto of_face(std::size_t f) const -> std::size_t { return n_vertices + f; }
        auto of_edge(std::size_t e) const -> std::size_t { return n_vertices + n_faces + e; }
    };

    /// Cage refined once by Catmull-Clark; `topo` must be the topology of `cage`.
    template <std::floating_point T>
    auto refine(const Cage<T> &cage, const Topology &topo) -> Cage<T>
    {
        const auto nv = topo.n_vertices(), nf = topo.n_faces(), ne = topo.n_edges();
        const ChildIndices child{nv, nf};
        const auto &P = cage.vertices;
        Cage<T> r;
        r.vertices.resize(nv + nf + ne);

        for (std::size_t f{}; f < nf; ++f)
        {
            point<T, 3> c{};
            const auto h0 = topo.face_begin(f), n = topo.face_size(f);
            for (auto h = h0; h < h0 + n; ++h)
                c += P[topo.origin(h)];
            r.vertices[child.of_face(f)] = c / T(n);
        }
        for (std::size_t e{}; e < ne; ++e)
        {
            const auto [a, b] = topo.edge_vertices(e);
            if (topo.is_sharp(e))
                r.vertices[child.of_edge(e)] = (P[a] + P[b]) / T(2);
            else
            {
                const auto h = topo.edge_half_edge(e);
                r.vertices[child.of_edge(e)] = (P[a] + P[b] + r.vertices[child.of_face(topo.face(h))] + r.vertices[child.of_face(topo.face(topo.twin(h)))]) / T(4);
            }
        }
        for (std::size_t v{}; v < nv; ++v)
        {
            auto &q = r.vertices[child.of_vertex(v)];
            switch (topo.rule(v))
            {
            case VertexRule::Corner:
                q = P[v];
                break;
            case VertexRule::Crease:
            {
                point<T, 3> s{};
                for (const auto e : topo.vertex_edges(v))
                    if (topo.is_sharp(e))
                    {
                        const auto &ev = topo.edge_vertices(e);
                        s += P[ev[0] == v ? ev[1] : ev[0]];
                    }
                q = (T(6) * P[v] + s) / T(8);
                break;
            }
            case VertexRule::Smooth:
            case VertexRule::Dart:
            {
                point<T, 3> se{}, sf{};
                for (const auto e : topo.vertex_edges(v))
                {
                    const auto &ev = topo.edge_vertices(e);
                    se += P[ev[0] == v ? ev[1] : ev[0]];
                }
                for (const auto f : topo.vertex_faces(v))
                    sf += r.vertices[child.of_face(f)];
                const auto n = T(topo.valence(v));
                q = (n - T(2)) / n * P[v] + (se + sf) / (n * n);
                break;
            }
            }
        }

        r.faces.reserve(topo.n_half_edges());
        constexpr std::array<std::size_t, 4> start{0, 3, 2, 1}; // corner of lowest (u, v) of the child of corner i of a quadrangle
        for (std::size_t f{}; f < nf; ++f)
        {
            const auto h0 = topo.face_begin(f), n = topo.face_size(f);
            for (auto h = h0; h < h0 + n; ++h)
            {
                const std::array<std::size_t, 4> q{child.of_vertex(topo.origin(h)), child.of_edge(topo.edge(h)), child.of_face(f), child.of_edge(topo.edge(topo.prev(h)))};
                const auto s = n == 4 ? start[h - h0] : 0;
                r.faces.push_back({q[s], q[(s + 1) % 4], q[(s + 2) % 4], q[(s + 3) % 4]});
            }
        }
        for (std::size_t e{}; e < ne; ++e)
            if (topo.is_crease(e))
            {
                const auto [a, b] = topo.edge_vertices(e);
                r.creases.push_back({child.of_vertex(a), child.of_edge(e)});
                r.creases.push_back({child.of_edge(e), child.of_vertex(b)});
            }
        for (std::size_t v{}; v < nv; ++v)
            if (topo.is_corner_tagged(v))
                r.corners.push_back(child.of_vertex(v));
        return r;
    }

    /**
     * @brief Cage refined `levels` times by Catmull-Clark (OpenSubdiv's Catmark rules).
     *
     * The cage is checked first (`CageError`). A face of n sides becomes n quadrangles at the
     * first level, 4 at each next one.
     */
    template <std::floating_point T>
    auto refine(const Cage<T> &cage, std::size_t levels = 1) -> Cage<T>
    {
        Topology topo{cage};
        if (levels == 0)
            return cage;
        auto r = refine(cage, topo);
        for (std::size_t l = 1; l < levels; ++l)
            r = refine(r, Topology{r});
        return r;
    }
}
