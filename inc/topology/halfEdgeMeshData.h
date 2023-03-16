#pragma once
#include <memory>
#include <algorithm>
namespace gbs
{
    template <typename T, size_t dim>
    struct HalfEdge;
    template <typename T, size_t dim>
    struct HalfEdgeVertex;
    template <typename T, size_t dim>
    struct HalfEdgeFace;

    template <typename T, size_t dim>
    struct HalfEdgeVertex
    {
        std::array<T, dim> coords;
        std::shared_ptr<HalfEdge<T, dim>> edge;
    };

    template <typename T, size_t dim>
    struct HalfEdgeFace
    {
        std::shared_ptr<HalfEdge<T, dim>> edge;
    };

    template <typename T, size_t dim>
    struct HalfEdge
    {
        std::shared_ptr<HalfEdgeVertex<T, dim>> vertex;
        std::shared_ptr<HalfEdgeFace<T, dim>> face;
        std::shared_ptr<HalfEdge<T, dim>> next;
        std::shared_ptr<HalfEdge<T, dim>> previous;
        std::shared_ptr<HalfEdge<T, dim>> opposite;

        // HalfEdge(const std::array<T, dim> &coords) : vertex{std::make_shared<HalfEdgeVertex<T, dim>>({coords, this})} {}
    };

    /**
     * @brief Builds half edge from vertex and tag, if the vertex is free it is tagged as belonging to the edge
     *
     * @tparam T
     * @tparam dim
     * @param vertex
     * @return std::shared_ptr< HalfEdge<T, dim> >
     */

    template <typename T, size_t dim>
    auto make_shared_h_edge(const std::shared_ptr<HalfEdgeVertex<T, dim>> &vertex, const std::shared_ptr<HalfEdgeFace<T, dim>> &face = nullptr) -> std::shared_ptr<HalfEdge<T, dim>>
    {
        auto hedge = std::make_shared<HalfEdge<T, dim>>(
            HalfEdge<T, dim>{
                .vertex = vertex,
                .face = face});

        if (!vertex->edge)
        {
            vertex->edge = hedge;
        }

        return hedge;
    }

    template <typename T, size_t dim>
    auto make_shared_h_vertex(const std::array<T, dim> &coords)
    {
        return std::make_shared<HalfEdgeVertex<T, dim>>(HalfEdgeVertex<T, dim>{coords, nullptr});
    }

    /**
     * @brief Builds half edge and its vertices from coordinate
     *
     * @tparam T
     * @tparam dim
     * @param coords
     * @return std::shared_ptr< HalfEdge<T, dim> >
     */
    template <typename T, size_t dim>
    auto make_shared_h_edge(const std::array<T, dim> &coords) -> std::shared_ptr<HalfEdge<T, dim>>
    {
        auto vertex = std::make_shared<HalfEdgeVertex<T, dim>>(HalfEdgeVertex<T, dim>{coords, nullptr});
        return make_shared_h_edge<T, dim>(vertex);
    }

    template <typename T, size_t dim, typename _Container>
    auto make_shared_h_vertices(const _Container &coords)
    {
        auto n = coords.size();
        std::vector<std::shared_ptr<HalfEdgeVertex<T, dim>>> vertices_cloud(n);
        std::transform(
            coords.begin(),
            coords.end(),
            vertices_cloud.begin(),
            [](const auto &xyz)
            { return make_shared_h_vertex<T, dim>(xyz); });
        return vertices_cloud;
    }

    template <typename T, size_t dim, typename _Container>
    auto make_shared_h_edges(const _Container &coords)
    {
        auto n = coords.size();
        std::vector<std::shared_ptr<HalfEdge<T, dim>>> vertices_cloud(n);
        std::transform(
            coords.begin(),
            coords.end(),
            vertices_cloud.begin(),
            [](const auto &xyz)
            { return make_shared_h_edge<T, dim>(xyz); });
        return vertices_cloud;
    }

    template <typename T, size_t dim, typename _It>
    void make_loop(const _It &begin, const _It &end, std::shared_ptr<HalfEdgeFace<T, dim>> &p_face)
    {
        auto n = std::distance(begin, end);
        assert(p_face);
        assert(n > 1);

        p_face->edge = (*begin);

        auto it = begin;
        while (it != end)
        {
            (*it)->face = p_face;
            (*it)->next = std::next(it) != end ? *std::next(it) : *begin;
            (*it)->previous = it != begin ? *std::prev(it) : *std::next(begin, n - 1);
            std::advance(it, 1);
        }
    }

    /**
     * @brief Builds a face from a list of half edges. Edges are tagged as belonging to the newly created face.
     *  A nullptr is returned if an half edge already belongs to a face.
     *
     * @tparam T
     * @tparam dim
     * @tparam _It
     * @param begin
     * @param end
     * @return std::shared_ptr< HalfEdgeFace<T, dim> >
     */
    template <typename T, size_t dim, typename _It>
    auto make_shared_h_face(const _It &begin, const _It &end) -> std::shared_ptr<HalfEdgeFace<T, dim>>
    {

        auto n = std::distance(begin, end);
        if (n < 2)
        {
            return nullptr;
        }

        auto p_face = std::make_shared<HalfEdgeFace<T, dim>>();
        make_loop(begin, end, p_face);

        return p_face;
    }

    /**
     * @brief Builds a face from a list of half edges. Edges are tagged as belonging to the newly created face.
     *  A nullptr is returned if an half edge already belongs to a face.
     *
     * @tparam T
     * @tparam dim
     * @param lst
     * @return std::shared_ptr< HalfEdgeFace<T, dim> >
     */
    template <typename T, size_t dim>
    auto make_shared_h_face(const auto &lst) -> std::shared_ptr<HalfEdgeFace<T, dim>>
    {
        return make_shared_h_face<T, dim>(lst.begin(), lst.end());
    }

    template <typename T, size_t dim>
    auto make_opposite(const std::shared_ptr<HalfEdgeVertex<T, dim>> &vertex, const std::shared_ptr<HalfEdge<T, dim>> &edge)
    {
        auto opposite = make_shared_h_edge(edge->previous->vertex);
        opposite->opposite = edge;
        edge->opposite = opposite;
        return opposite;
    }

    template <typename T, size_t dim>
    inline auto make_HalfEdges(std::vector< std::array<T,dim> > &coords)
    {
        
        long long n = coords.size();
        std::vector<std::shared_ptr< HalfEdge<T,dim> > > h_edges(n);
        std::transform(
            coords.begin(), coords.end(),
            h_edges.begin(),
            [](const auto &X){
                return std::make_shared< HalfEdge<T,dim> >(
                HalfEdge<T,dim>{
                    .vertex = make_shared_h_vertex(X)
                }
            );}
        );

        auto nm = n-1;
        h_edges.front()->next = h_edges[1];
        h_edges.front()->previous = h_edges.back();
        for( long long i{1}; i < nm; i++)
        {
            h_edges[i]->previous = h_edges[i-1];
            h_edges[i]->next = h_edges[i+1];
        }
        h_edges.back()->next =h_edges.front();
        h_edges.back()->previous =  h_edges[nm-1];

        return h_edges;
    }

}