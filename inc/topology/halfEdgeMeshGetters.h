#pragma once

#include "halfEdgeMeshData.h"

namespace gbs
{
    template <typename T, size_t dim>
    auto getFaceEdge(const std::shared_ptr<HalfEdgeFace<T, dim>> &face, const std::shared_ptr<HalfEdgeVertex<T, dim>> &vertex) -> std::shared_ptr<HalfEdge<T, dim>>
    {
        auto edge = face->edge;
        while (edge->vertex != vertex)
        {
            edge = edge->next;
            if (edge == face->edge) // loop completed
            {
                return nullptr;
            }
        }
        return edge;
    }

    template <typename T, size_t dim>
    auto getFaceEdges(const HalfEdgeFace<T, dim> &face)
    {
        std::list<std::shared_ptr<HalfEdge<T, dim>>> edges_lst;
        auto edge = face.edge;
        while (edge)
        {
            edges_lst.push_back(edge);
            edge = edge->next;
            if (edge == face.edge)
            {
                break;
            }
        }
        return edges_lst;
    }

    template <typename T, size_t dim>
    auto getFaceEdges(const std::shared_ptr<HalfEdgeFace<T, dim>> &face)
    {
        return getFaceEdges(*face);
    }

    template <typename T, size_t dim>
    auto getFaceVertices(const HalfEdgeFace<T, dim> &face)
    {
        std::list<std::shared_ptr<HalfEdgeVertex<T, dim>>> vtx_lst;
        auto edge = face.edge;
        while (edge)
        {
            vtx_lst.push_back(edge->vertex);
            edge = edge->next;
            if (edge == face.edge)
            {
                break;
            }
        }
        return vtx_lst;
    }

    template <typename T, size_t dim>
    auto getFaceVertices(const std::shared_ptr<HalfEdgeFace<T, dim>> &face)
    {
        return getFaceVertices(*face);
    }

    template <typename T, size_t dim>
    auto getFaceCoords(const std::shared_ptr<HalfEdgeFace<T, dim>> &face)
    {
        std::list<std::array<T, dim>> coords_lst;
        auto edge = face->edge;
        while (edge)
        {
            coords_lst.push_back(edge->vertex->coords);
            edge = edge->next;
            if (edge == face->edge)
            {
                break;
            }
        }
        return coords_lst;
    }

    /**
     * @brief Get the Common Edge of h_f1 with h_f2, nullptr if none
     *
     * @tparam T
     * @tparam dim
     * @param h_f1
     * @param h_f2
     * @return std::shared_ptr<HalfEdge<T, dim>>
     */
    template <typename T, size_t dim>
    auto getCommonEdge(const std::shared_ptr<HalfEdgeFace<T, dim>> &h_f1, const std::shared_ptr<HalfEdgeFace<T, dim>> &h_f2) -> std::shared_ptr<HalfEdge<T, dim>>
    {
        for (const auto &h_e : getFaceEdges(h_f1))
        {
            if (h_e->opposite && h_e->opposite->face == h_f2)
                return h_e;
        }
        return nullptr;
    }
    /**
     * @brief Get the Common Edges of h_f1 with h_f2, (nullptr, nullptr) if none
     *
     * @tparam T
     * @tparam dim
     * @param h_f1
     * @param h_f2
     * @return std::pair< std::shared_ptr<HalfEdge<T, dim>>, std::shared_ptr<HalfEdge<T, dim>> >
     */
    template <typename T, size_t dim>
    auto getCommonEdges(const std::shared_ptr<HalfEdgeFace<T, dim>> &h_f1, const std::shared_ptr<HalfEdgeFace<T, dim>> &h_f2) -> std::pair<std::shared_ptr<HalfEdge<T, dim>>, std::shared_ptr<HalfEdge<T, dim>>>
    {
        if (auto h_e1 = getCommonEdge(h_f1, h_f2))
        {
            return std::make_pair(h_e1, h_e1->opposite);
        }
        return std::make_pair(nullptr, nullptr);
    }

    template <typename T, size_t dim>
    auto getPreviousFace(const std::shared_ptr<HalfEdge<T, dim>> &edge) -> std::shared_ptr<HalfEdgeFace<T, dim>>
    {
        if (edge->next)
        {
            auto opp = edge->next->opposite;
            if (opp)
            {
                return opp->face;
            }
        }
        return nullptr;
    }

    template <typename T, size_t dim>
    auto getNextFace(const std::shared_ptr<HalfEdge<T, dim>> &edge) -> std::shared_ptr<HalfEdgeFace<T, dim>>
    {
        auto opp = edge->opposite;
        if (opp)
        {
            return opp->face;
        }
        return nullptr;
    }

    // template <typename T, size_t dim>
    // auto getFaces(const std::shared_ptr<HalfEdge<T, dim>> &edge) -> std::list< std::shared_ptr< HalfEdgeFace<T,dim> > >
    // {
    //     std::list< std::shared_ptr< HalfEdgeFace<T,dim> > > face;

    // }

    template <typename T, size_t dim>
    auto getFacesAttachedToVertex(const std::shared_ptr<HalfEdgeVertex<T, dim>> &h_v)
    {
        assert(h_v->edge);

        std::list<std::shared_ptr<HalfEdgeFace<T, dim>>> neighbors;
        auto start = h_v->edge;

        auto current = start;
        do
        {
            neighbors.push_front(current->face);
            if (current->opposite)
            {
                current = current->opposite->previous;
            }
            else
            {
                current = nullptr;
            }

        } while (current && current != start);

        if (current != start && start->next->opposite)
        {
            current = start->next->opposite;
            do
            {
                neighbors.push_back(current->face);
                current = current->next->opposite;
            } while (current && current != start);
        }

        return neighbors;
    }

    template <typename T, size_t dim>
    auto getNeighboringFaces(const std::shared_ptr<HalfEdgeFace<T, dim>> &h_f)
    {
        std::list<std::shared_ptr<HalfEdgeFace<T, dim>>> neighbors;
        auto edges = getFaceEdges(h_f);
        for (const auto &h_e : edges)
        {
            if (h_e->opposite)
            {
                assert(h_e->opposite->face);
                neighbors.push_back((h_e->opposite->face));
            }
        }
        return neighbors;
    }

    template <typename T, size_t dim>
    auto getFacesBoundary(const std::list<std::shared_ptr<HalfEdgeFace<T, dim>>> &h_f_lst)
    {
        std::list<std::shared_ptr<HalfEdge<T, dim>>> boundary;
        auto begin = h_f_lst.begin();
        auto end = h_f_lst.end();
        for (const auto &h_f : h_f_lst)
        {
            const auto &h_e_lst = getFaceEdges(h_f);
            for (const auto &h_e : h_e_lst)
            {
                if (
                    !h_e->opposite // whole mesh boundary
                    ||
                    end == std::find(begin, end, h_e->opposite->face) // opposite face is not within th selection
                )
                {
                    boundary.push_back(h_e);
                }
            }
        }

        return boundary;
    }

    template <typename T, size_t dim>
    auto getOrientedFacesBoundaries(const std::list<std::shared_ptr<HalfEdgeFace<T, dim>>> &h_f_lst)
    {
        auto boundary = getFacesBoundary(h_f_lst);
        return takeClosedLoops(boundary);
    }

    template <typename T, size_t dim>
    auto getOrientedFacesBoundary(const std::list<std::shared_ptr<HalfEdgeFace<T, dim>>> &h_f_lst)
    {
        auto boundary = getFacesBoundary(h_f_lst);
        return takeClosedLoops(boundary).front();
    }

    template <typename T, size_t dim>
    auto getVerticesMapFromFaces(const auto &faces_lst ) -> std::map< std::shared_ptr< HalfEdgeVertex<T,dim> >, size_t >
    {
        std::map< std::shared_ptr< HalfEdgeVertex<T,dim> >, size_t > vertices_map;
        size_t index{};
        for( const auto &f : faces_lst)
        {
            auto vtx_lst = getFaceVertices(f);
            for( const auto &vtx : vtx_lst)
            {
                if(!vertices_map.contains(vtx))
                {
                    vertices_map[vtx] = index;
                    index++;
                }
            }
        }

        return vertices_map;
    }

    template <typename T, size_t dim>
    auto getVerticesVectorFromFaces(const auto &faces_lst ) -> std::vector< std::shared_ptr< HalfEdgeVertex<T,dim> > >
    {
        auto vertices_map = getVerticesMapFromFaces<T,dim>(faces_lst);

        std::vector< std::shared_ptr< HalfEdgeVertex<T,dim> > > vertices(vertices_map.size());
        std::transform(
            std::execution::par,
            vertices_map.begin(), vertices_map.end(),
            vertices.begin(),
            [](const auto &pair){return pair.first;}
        );

        return vertices;
    }

    template <typename T, size_t dim>
    auto getEdgesMap(const auto &faces_lst ) -> std::map< std::shared_ptr< HalfEdgeVertex<T,dim> >, size_t >
    {
        std::map< std::shared_ptr< HalfEdgeEdge<T,dim> >, size_t > edges_map;
        size_t index{};
        for( const auto &f : faces_lst)
        {
            auto hed_lst = getFaceEdges(f);
            for( const auto &hed : hed_lst)
            {
                if(!edges_map.contains(hed))
                {
                    edges_map[hed] = index;
                    index++;
                }
            }
        }

        return edges_map;
    }
}