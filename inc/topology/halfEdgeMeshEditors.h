#pragma once
#include "halfEdgeMeshData.h"
#include "halfEdgeMeshGetters.h"
#include <list>

namespace gbs
{
    template <typename T, size_t dim>
    auto add_face(
        const std::shared_ptr<HalfEdgeFace<T, dim>> &face,
        const std::shared_ptr<HalfEdge<T, dim>> &edge,
        const std::array<T, dim> &coords) -> std::shared_ptr<HalfEdgeFace<T, dim>>
    {
        if (!edge || edge->opposite || edge->face != face)
        {
            return nullptr;
        }

        auto opposite = make_opposite(edge->previous->vertex, edge);

        auto lst = {make_shared_h_edge(edge->vertex), opposite, make_shared_h_edge(coords)};

        return make_shared_h_face<T, dim>(lst);
    }

    template <typename T, size_t dim>
    auto add_face(
        const std::shared_ptr<HalfEdge<T, dim>> &edge,
        const std::array<T, dim> &coords) -> std::shared_ptr<HalfEdgeFace<T, dim>>
    {
        if (!edge || edge->opposite)
        {
            return nullptr;
        }

        auto opposite = make_opposite(edge->previous->vertex, edge);

        auto lst = {make_shared_h_edge(edge->vertex), opposite, make_shared_h_edge(coords)};

        return make_shared_h_face<T, dim>(lst);
    }

    template <typename T, size_t dim>
    void associate(std::shared_ptr<HalfEdgeVertex<T, dim>> &h_v, std::shared_ptr<HalfEdge<T, dim>> &h_e)
    {
        h_e->vertex = h_v;
        h_v->edge = h_e;
    }

    template <typename T, size_t dim>
    void flip(std::shared_ptr<HalfEdgeFace<T, dim>> &h_f1, std::shared_ptr<HalfEdgeFace<T, dim>> &h_f2)
    {
        assert(getFaceEdges(h_f1).size() == 3);
        assert(getFaceEdges(h_f2).size() == 3);
        auto [h_e1_1, h_e1_2] = getCommonEdges(h_f1, h_f2);
        if (h_e1_1 == nullptr || h_e1_2 == nullptr)
            return;
        auto h_e2_1 = h_e1_1->next;
        auto h_e3_1 = h_e2_1->next;
        auto h_e2_2 = h_e1_2->next;
        auto h_e3_2 = h_e2_2->next;

        auto h_v1 = h_e1_1->vertex;
        auto h_v2 = h_e2_1->vertex;
        auto h_v3 = h_e3_1->vertex;
        auto h_v4 = h_e2_2->vertex;

        associate(h_v4, h_e1_1);
        associate(h_v2, h_e1_2);

        std::list<std::shared_ptr<HalfEdge<T, dim>>> lst1{h_e1_1, h_e3_2, h_e2_1};
        std::list<std::shared_ptr<HalfEdge<T, dim>>> lst2{h_e1_2, h_e3_1, h_e2_2};

        make_loop(lst1.begin(), lst1.end(), h_f1);
        make_loop(lst2.begin(), lst2.end(), h_f2);
    }
    /**
     * @brief Link 2 half edges to as one edge
     *
     * @tparam T
     * @tparam dim
     * @param h_e1
     * @param h_e2
     */
    template <typename T, size_t dim>
    void link_edges(const std::shared_ptr<HalfEdge<T, dim>> &h_e1, const std::shared_ptr<HalfEdge<T, dim>> &h_e2)
    {
        assert(h_e1);
        assert(h_e2);
        h_e1->opposite = h_e2;
        h_e2->opposite = h_e1;
    }

    template <typename T, size_t dim>
    auto add_vertex(const std::shared_ptr<HalfEdgeFace<T, dim>> &h_f, const std::shared_ptr<HalfEdgeVertex<T, dim>> &h_v)
    {
        return add_vertex(getFaceEdges(h_f), h_v);
    }

    template <typename T, size_t dim>
    auto add_vertex(const std::list<std::shared_ptr<HalfEdge<T, dim>>> &h_e_lst, const std::shared_ptr<HalfEdgeVertex<T, dim>> &h_v)
    {
        std::list<std::shared_ptr<HalfEdgeFace<T, dim>>> h_f_lst;

        std::shared_ptr<HalfEdge<T, dim>> h_e_prev{};
        h_v->edge = nullptr; // first new half edge takes ownership
        for (const auto &h_e : h_e_lst)
        {
            assert(h_e->previous);
            auto h_e1 = make_shared_h_edge(h_v);
            auto h_e2 = make_shared_h_edge(h_e->previous->vertex);
            if (h_e_prev)
            {
                link_edges(h_e2, h_e_prev);
            }
            h_e_prev = h_e1;
            auto lst = {h_e, h_e1, h_e2};
            h_f_lst.push_back(
                make_shared_h_face<T, dim>(lst));
            assert(is_ccw(h_f_lst.back()));
        }
        link_edges(h_f_lst.back()->edge->next, h_f_lst.front()->edge->previous);
        return h_f_lst;
    }

    /**
     * @brief Remove faces attached to vertex vtx
     *
     * @tparam T
     * @tparam dim
     * @param faces_lst Face container
     * @param vtx       The vertex
     * @return size_t   The number of removed faces
     */
    template <typename T, size_t dim>
    auto remove_faces(auto &faces_lst, const std::shared_ptr<HalfEdgeVertex<T, dim>> &vtx) -> size_t
    {
        size_t count{};
        auto it_hf = faces_lst.begin();

        while (it_hf != faces_lst.end())
        {
            auto h_e_lst = getFaceEdges(*it_hf);
            auto it = std::find_if(h_e_lst.begin(), h_e_lst.end(),
                                   [&vtx](const auto &he)
                                   { return vtx == he->vertex; });
            if (h_e_lst.end() != it)
            {
                for (auto &h_e : h_e_lst) // no more opposite on adjacent faces
                {
                    if (h_e->opposite)
                    {
                        h_e->opposite->opposite = nullptr;
                    }
                }
                it_hf = faces_lst.erase(it_hf);
                count++;
            }
            else
            {
                std::advance(it_hf, 1);
            }
        }
        return count;
    }

    template <typename T, size_t dim>
    auto takeClosedLoops(std::list<std::shared_ptr<HalfEdge<T, dim>>> &boundary)
    {
        std::list<std::list<std::shared_ptr<HalfEdge<T, dim>>>> boundaries_oriented;

        while (boundary.size())
        {
            std::list<std::shared_ptr<HalfEdge<T, dim>>> boundary_oriented;

            auto previous = boundary.end();

            boundary_oriented.push_front(boundary.front());
            boundary.erase(boundary.begin());

            while (previous != boundary.begin())
            {
                auto tail = boundary_oriented.front()->previous->vertex;
                auto it = std::find_if(
                    boundary.begin(), boundary.end(),
                    [tail](const auto &e)
                    {
                        return e->vertex == tail;
                    });
                if (it != boundary.end())
                {
                    boundary_oriented.push_front(*it);
                    boundary.erase(it);
                }
                else
                {
                    break;
                }
            }

            boundaries_oriented.push_back(boundary_oriented);
        }

        return boundaries_oriented;
    }

    template <typename T, typename _Container1, typename _Container2>
    auto takeExternalFaces(_Container1 &faces_lst, const _Container2 &boundary)
    {
        std::list<std::shared_ptr<HalfEdgeFace<T, 2>>> external_faces;
        auto it = faces_lst.begin();
        while (it != faces_lst.end())
        {
            it = std::find_if(
                faces_lst.begin(), faces_lst.end(),
                [&boundary](const auto &hf)
                {
                    auto coords = getFaceCoords(hf);
                    std::array<T, 2> G{};
                    for (const auto &xy : coords)
                    {
                        G = G + xy;
                    }
                    G = G / T(coords.size());
                    return !is_inside(G, boundary);

                    // for(const auto &xy : coords)
                    // {
                    //     if(!is_inside(xy,boundary))
                    //     {
                    //         return true;
                    //     }
                    // }
                    // return false;
                });
            if (it != faces_lst.end())
            {
                for (auto &h_e : getFaceEdges(*it))
                {
                    if (h_e->opposite)
                    {
                        h_e->opposite->opposite = nullptr;
                    }
                }
                external_faces.push_back(*it);
                it = faces_lst.erase(it);
            }
        }
        return external_faces;
    }

    template <typename T, typename _Container1, typename _Container2>
    auto takeInternalFaces(_Container1 &faces_lst, const _Container2 &boundary)
    {
        std::list<std::shared_ptr<HalfEdgeFace<T, 2>>> external_faces;
        auto it = faces_lst.begin();
        while (it != faces_lst.end())
        {
            it = std::find_if(
                faces_lst.begin(), faces_lst.end(),
                [&boundary](const auto &hf)
                {
                    auto coords = getFaceCoords(hf);
                    std::array<T, 2> G{};
                    for (const auto &xy : coords)
                    {
                        G = G + xy;
                    }
                    G = G / T(coords.size());
                    return is_inside(G, boundary);
                });
            if (it != faces_lst.end())
            {
                for (auto &h_e : getFaceEdges(*it))
                {
                    if (h_e->opposite)
                    {
                        h_e->opposite->opposite = nullptr;
                    }
                }
                external_faces.push_back(*it);
                it = faces_lst.erase(it);
            }
        }
        return external_faces;
    }

    void reverseBoundary(auto &boundary)
    {
        std::reverse(boundary.begin(), boundary.end());
        std::transform(
            boundary.begin(), boundary.end(),
            boundary.begin(),
            [](const auto &he)
            {
                std::swap(he->previous, he->next);
                return he;
            });
    }

    template <typename _Container, typename _It>
    auto eraseFace(_Container &h_f_lst, const _It &it)
    {

        auto face_edges = getFaceEdges(*it);
        for (auto &he : face_edges)
        {
            if (he->opposite)
            {
                he->opposite->opposite = nullptr;
            }
        }
        return h_f_lst.erase(it);
    }
}