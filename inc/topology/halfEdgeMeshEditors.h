#pragma once
#include "halfEdgeMeshData.h"
#include "halfEdgeMeshGetters.h"
#include <list>

namespace gbs
{
/**
 * @brief Adds a new face to an existing half-edge face by connecting it to a given half-edge.
 *
 * This function takes a shared pointer to an existing HalfEdgeFace, a shared pointer to a HalfEdge,
 * and a new coordinate. It creates a new face by connecting the existing face to the new coordinate
 * using the given half-edge.
 *
 * @tparam T The floating-point type used for coordinates.
 * @tparam dim The dimensionality of the space.
 * @param face A shared pointer to the existing HalfEdgeFace.
 * @param edge A shared pointer to the HalfEdge that will be used to connect the new face.
 * @param coords The new coordinate, represented as std::array<T, dim>, that will be part of the new face.
 * @return A shared pointer to the new HalfEdgeFace created by connecting the existing face to the new coordinate,
 *         or nullptr if the edge is invalid or already has an opposite edge or a different face.
 */
    template <std::floating_point T, size_t dim>
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
/**
 * @brief Adds a new face to an existing half-edge by connecting it to a given half-edge.
 *
 * This function takes a shared pointer to a HalfEdge and a new coordinate. It creates a new face
 * by connecting the existing half-edge to the new coordinate using the given half-edge.
 *
 * @tparam T The floating-point type used for coordinates.
 * @tparam dim The dimensionality of the space.
 * @param edge A shared pointer to the HalfEdge that will be used to connect the new face.
 * @param coords The new coordinate, represented as std::array<T, dim>, that will be part of the new face.
 * @return A shared pointer to the new HalfEdgeFace created by connecting the existing half-edge to the new coordinate,
 *         or nullptr if the edge is invalid or already has an opposite edge.
 */
    template <std::floating_point T, size_t dim>
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
/**
 * @brief Associates a HalfEdgeVertex with a HalfEdge by setting their respective pointers.
 *
 * This function takes shared pointers to a HalfEdgeVertex and a HalfEdge, and sets the vertex
 * pointer of the HalfEdge and the edge pointer of the HalfEdgeVertex to point to each other.
 *
 * @tparam T The floating-point type used for coordinates.
 * @tparam dim The dimensionality of the space.
 * @param h_v A shared pointer to the HalfEdgeVertex that will be associated with the HalfEdge.
 * @param h_e A shared pointer to the HalfEdge that will be associated with the HalfEdgeVertex.
 */
    template <std::floating_point T, size_t dim>
    void associate(std::shared_ptr<HalfEdgeVertex<T, dim>> &h_v, std::shared_ptr<HalfEdge<T, dim>> &h_e)
    {
        h_e->vertex = h_v;
        h_v->edge = h_e;
    }
/**
 * @brief Flips the connection between two triangular HalfEdgeFaces that share a common edge.
 *
 * This function takes shared pointers to two HalfEdgeFaces and performs a flip operation if they share
 * a common edge. A flip operation swaps the connection between two triangles that share an edge, effectively
 * changing their connectivity in the mesh. This function assumes that both input faces are triangles (have 3 edges).
 *
 * @tparam T The floating-point type used for coordinates.
 * @tparam dim The dimensionality of the space.
 * @param h_f1 A shared pointer to the first HalfEdgeFace that will be involved in the flip operation.
 * @param h_f2 A shared pointer to the second HalfEdgeFace that will be involved in the flip operation.
 */
    template <std::floating_point T, size_t dim>
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
 * @brief Links two HalfEdge objects by setting their opposite pointers to each other.
 *
 * This function takes shared pointers to two HalfEdge objects and sets their opposite pointers to point
 * to each other. This establishes a bidirectional relationship between the two HalfEdge objects.
 *
 * @tparam T The floating-point type used for coordinates.
 * @tparam dim The dimensionality of the space.
 * @param h_e1 A shared pointer to the first HalfEdge that will be linked to the second HalfEdge.
 * @param h_e2 A shared pointer to the second HalfEdge that will be linked to the first HalfEdge.
 */
    template <std::floating_point T, size_t dim>
    void link_edges(const std::shared_ptr<HalfEdge<T, dim>> &h_e1, const std::shared_ptr<HalfEdge<T, dim>> &h_e2)
    {
        assert(h_e1);
        assert(h_e2);
        h_e1->opposite = h_e2;
        h_e2->opposite = h_e1;
    }
/**
 * @brief Adds a vertex to a list of HalfEdge objects by creating new HalfEdge objects and connecting them.
 *
 * This function takes a list of shared pointers to HalfEdge objects and a shared pointer to a HalfEdgeVertex.
 * It creates new HalfEdge objects that connect the existing vertices of the face to the new vertex
 * and vice versa. It then returns a list of shared pointers to the new HalfEdgeFace objects created.
 *
 * @tparam T The floating-point type used for coordinates.
 * @tparam dim The dimensionality of the space.
 * @param h_e_lst A list of shared pointers to HalfEdge objects to which the new vertex will be added.
 * @param h_v A shared pointer to the HalfEdgeVertex that will be added.
 * @return A list of shared pointers to HalfEdgeFace objects containing the new vertex.
 */
    template <std::floating_point T, size_t dim>
    auto add_vertex(const std::list<std::shared_ptr<HalfEdge<T, dim>>> &h_e_lst, const std::shared_ptr<HalfEdgeVertex<T, dim>> &h_v)
    {
        // std::list<std::shared_ptr<HalfEdgeFace<T, dim>>> h_f_lst;

        // std::shared_ptr<HalfEdge<T, dim>> h_e_prev{};
        // h_v->edge = nullptr; // first new half edge takes ownership
        // for (const auto &h_e : h_e_lst)
        // {
        //     assert(h_e->previous);
        //     auto h_e1 = make_shared_h_edge(h_v);
        //     auto h_e2 = make_shared_h_edge(h_e->previous->vertex);
        //     if (h_e_prev)
        //     {
        //         link_edges(h_e2, h_e_prev);
        //     }
        //     h_e_prev = h_e1;
        //     auto lst = {h_e, h_e1, h_e2};
        //     h_f_lst.push_back(
        //         make_shared_h_face<T, dim>(lst));
        //     assert(is_ccw(h_f_lst.back()));
        // }
        // link_edges(h_f_lst.back()->edge->next, h_f_lst.front()->edge->previous);
        // return h_f_lst;
        std::list<std::shared_ptr<HalfEdgeFace<T, dim>>> h_f_lst;

        std::shared_ptr<HalfEdge<T, dim>> h_e_prev{};
        h_v->edge = nullptr;

        std::ranges::for_each(h_e_lst, [&](const auto &h_e) {
            assert(h_e->previous);
            auto h_e1 = make_shared_h_edge(h_v);
            auto h_e2 = make_shared_h_edge(h_e->previous->vertex);

            if (h_e_prev)
            {
                link_edges(h_e2, h_e_prev);
            }

            h_e_prev = h_e1;
            auto lst = {h_e, h_e1, h_e2};

            auto new_face = make_shared_h_face<T, dim>(lst);
            assert(is_ccw(new_face));
            h_f_lst.push_back(new_face);
        });

        link_edges(h_f_lst.back()->edge->next, h_f_lst.front()->edge->previous);

        return h_f_lst;
    }
/**
 * @brief Adds a vertex to a HalfEdgeFace by creating new HalfEdge objects and connecting them.
 *
 * This function takes a shared pointer to a HalfEdgeFace and a shared pointer to a HalfEdgeVertex.
 * It creates new HalfEdge objects that connect the existing vertices of the face to the new vertex
 * and vice versa.
 *
 * @tparam T The floating-point type used for coordinates.
 * @tparam dim The dimensionality of the space.
 * @param h_f A shared pointer to the HalfEdgeFace to which the new vertex will be added.
 * @param h_v A shared pointer to the HalfEdgeVertex that will be added to the face.
 * @return A shared pointer to a HalfEdgeFace that contains the new vertex.
 */
    template <std::floating_point T, size_t dim>
    auto add_vertex(const std::shared_ptr<HalfEdgeFace<T, dim>> &h_f, const std::shared_ptr<HalfEdgeVertex<T, dim>> &h_v)
    {
        return add_vertex(getFaceEdges(h_f), h_v);
    }
/**
 * @brief Removes faces from the list of faces that have a specified vertex and updates the count of removed faces.
 * 
 * @tparam T A floating-point type.
 * @tparam dim The dimension of the half-edge data structure.
 * @param faces_lst A list of shared pointers to HalfEdgeFace objects.
 * @param vtx A shared pointer to a HalfEdgeVertex object that is to be checked for its presence in the faces.
 * @return size_t The count of faces removed from the list.
 */
    template <std::floating_point T, size_t dim>
    auto remove_faces(auto &faces_lst, const std::shared_ptr<HalfEdgeVertex<T, dim>> &vtx) -> size_t
    {
        size_t count{};

        faces_lst.erase(
            std::remove_if(faces_lst.begin(), faces_lst.end(),
                [&vtx, &count](const auto &face)
                {
                    auto h_e_lst = getFaceEdges(face);
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
                        count++;
                        return true;
                    }
                    return false;
                }),
        faces_lst.end());

        return count;
    }
/**
 * @brief Extracts and separates closed loops from a list of half-edges.
 * 
 * @tparam T A floating-point type.
 * @tparam dim The dimension of the half-edge data structure.
 * @param boundary A list of shared pointers to HalfEdge objects representing a sequence of half-edges.
 * @return std::list<std::list<std::shared_ptr<HalfEdge<T, dim>>>> A list of lists of shared pointers to HalfEdge objects,
 *         where each inner list represents a closed loop of half-edges.
 */
    template <std::floating_point T, size_t dim>
    auto takeClosedLoops(std::list<std::shared_ptr<HalfEdge<T, dim>>> &boundary)
    {
        std::list<std::list<std::shared_ptr<HalfEdge<T, dim>>>> boundaries_oriented;

        while (!boundary.empty())
        {
            std::list<std::shared_ptr<HalfEdge<T, dim>>> boundary_oriented;

            auto tail_vertex = [&]()
            {
                return boundary_oriented.front()->previous->vertex;
            };

            boundary_oriented.push_front(std::move(boundary.front()));
            boundary.pop_front();

            while (true)
            {
                auto it = std::ranges::find_if(boundary, [&](const auto &e) { return e->vertex == tail_vertex(); });

                if (it != boundary.end())
                {
                    boundary_oriented.push_front(std::move(*it));
                    boundary.erase(it);
                }
                else
                {
                    break;
                }
            }

            boundaries_oriented.push_back(std::move(boundary_oriented));
        }

        return boundaries_oriented;
    }

/**
 * @brief Erases a face from the container of faces and removes the opposites of its half-edges.
 * 
 * @tparam Container The container type for the list of faces.
 * @tparam Iterator The iterator type pointing to the face to be erased.
 * @param h_f_lst The container of faces.
 * @param it The iterator pointing to the face to be erased.
 * @return An iterator pointing to the element following the erased face.
 */
    template <typename Container, std::forward_iterator Iterator>
    auto eraseFace(Container &h_f_lst, const Iterator &it)
    {
        // Get the half-edges of the face to be erased
        auto face_edges = getFaceEdges(*it);

        // Remove the opposites of the half-edges
        for (auto &he : face_edges)
        {
            if (he->opposite)
            {
                he->opposite->opposite = nullptr;
            }
        }

        return h_f_lst.erase(it);
    }
/**
 * @brief Removes and returns external faces from a list of faces based on a given boundary.
 * 
 * @tparam T Floating point type used for coordinates.
 * @param faces_lst List of faces to search and remove external faces from.
 * @param boundary Container holding the boundary information.
 * @return std::list<std::shared_ptr<HalfEdgeFace<T, 2>>> List of external faces.
 */
    template <std::floating_point T, typename _Container1, typename _Container2>
    auto takeExternalFaces(_Container1 &faces_lst, const _Container2 &boundary)
    {
        std::list<std::shared_ptr<HalfEdgeFace<T, 2>>> external_faces;
        auto it = faces_lst.begin();
        
        // Iterate through faces_lst to find external faces
        while (it != faces_lst.end())
        {
            // Search for the next external face in faces_lst
            it = std::ranges::find_if(faces_lst, [&boundary](const auto &hf)
            {
                return !is_centroid_inside_boundary(hf, boundary);
            });

            if (it != faces_lst.end())
            {
                external_faces.push_back(*it);
                it = eraseFace(faces_lst,it);
            }
        }
        
        return external_faces;
    }

/**
 * @brief Removes and returns internal faces from a list of faces based on a given boundary.
 * 
 * @tparam T Floating point type used for coordinates.
 * @param faces_lst List of faces to search and remove internal faces from.
 * @param boundary Container holding the boundary information.
 * @return std::list<std::shared_ptr<HalfEdgeFace<T, 2>>> List of internal faces.
 */
    template <std::floating_point T, typename _Container1, typename _Container2>
    auto takeInternalFaces(_Container1 &faces_lst, const _Container2 &boundary)
    {
        std::list<std::shared_ptr<HalfEdgeFace<T, 2>>> internal_faces;
        auto it = faces_lst.begin();

        // Iterate through faces_lst to find internal faces
        while (it != faces_lst.end())
        {
            // Search for the next internal face in faces_lst
            it = std::ranges::find_if(faces_lst, [&boundary](const auto &hf)
            {
                return is_centroid_inside_boundary(hf, boundary);
            });

            if (it != faces_lst.end())
            {
                internal_faces.push_back(*it);
                it = eraseFace(faces_lst,it);
            }
        }
        
        return internal_faces;
    }

/**
 * @brief Reverses the order of a given boundary, along with swapping the previous and next half-edges.
 * 
 * @tparam Container The container type for the boundary.
 * @param boundary The boundary container to be reversed.
 */
    template <typename Container>
    void reverseBoundary(Container &boundary)
    {
        // Reverse the order of the boundary
        std::ranges::reverse(boundary);

        // Swap the previous and next half-edges for each half-edge in the boundary
        std::ranges::transform(
            boundary,
            boundary.begin(),
            [](const auto &he)
            {
                std::swap(he->previous, he->next);
                return he;
            });
    }

}