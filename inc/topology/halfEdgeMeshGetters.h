#pragma once

#include <list>
#include <map>
#include "halfEdgeMeshData.h"
#include "baseGeom.h"
namespace gbs
{
/**
 * @brief Finds the half-edge of a face originating from a given vertex.
 * 
 * @tparam T Floating point type used for coordinates.
 * @tparam dim Dimension of the half-edge data structure.
 * @param face The face to search for the half-edge.
 * @param vertex The vertex from which the desired half-edge originates.
 * @return std::shared_ptr<HalfEdge<T, dim>> The half-edge originating from the given vertex, or nullptr if not found.
 */
    template <std::floating_point T, size_t dim>
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

/**
 * @brief Retrieves the three edges of a triangle face.
 * 
 * @tparam T The type of the coordinates (must be a floating-point type).
 * @tparam dim The dimension of the space (must be 2 or 3).
 * @param face The face whose edges to retrieve.
 * @return An array of three HalfEdge objects representing the edges of the triangle.
 */
    template <std::floating_point T, size_t dim>
    auto getTriangleEdges(const HalfEdgeFace<T, dim> &face)
    {
        assert(face.edge->next->next->next == face.edge);
        return std::array<std::shared_ptr<HalfEdge<T, dim>>,3> {
            face.edge->next,
            face.edge,
            face.edge->previous
        };
    }

/**
 * @brief Retrieves the three edges of a triangle face given a shared pointer to the face.
 * 
 * @tparam T The type of the coordinates (must be a floating-point type).
 * @tparam dim The dimension of the space (must be 2 or 3).
 * @param face A shared pointer to the face whose edges to retrieve.
 * @return An array of three HalfEdge objects representing the edges of the triangle.
 */
    template <std::floating_point T, size_t dim>
    auto getTriangleEdges(const std::shared_ptr<HalfEdgeFace<T, dim>> &face)
    {
        return getTriangleEdges(*face);
    }

/**
 * @brief Retrieves a list of edges belonging to a given face.
 * 
 * @tparam T Floating point type used for coordinates.
 * @tparam dim Dimension of the half-edge data structure.
 * @param face The face for which the edges are to be retrieved.
 * @return std::list<std::shared_ptr<HalfEdge<T, dim>>> A list of edges belonging to the face.
 */
    template <std::floating_point T, size_t dim>
    auto getFaceEdges(const HalfEdgeFace<T, dim> &face)
    {
        std::list<std::shared_ptr<HalfEdge<T, dim>>> edges_lst;
        auto edge = face.edge;
        assert(edge);
        while (edge)
        {
            edges_lst.push_back(edge);
            edge = edge->next;
            assert(edge);
            if (edge == face.edge)
            {
                break;
            }
        }
        return edges_lst;
    }

/**
 * @brief Retrieves a list of edges belonging to a given face.
 * 
 * @tparam T Floating point type used for coordinates.
 * @tparam dim Dimension of the half-edge data structure.
 * @param face A shared pointer to the face for which the edges are to be retrieved.
 * @return auto A list of edges belonging to the face.
 */
    template <std::floating_point T, size_t dim>
    auto getFaceEdges(const std::shared_ptr<HalfEdgeFace<T, dim>> &face)
    {
        return getFaceEdges(*face);
    }

/**
 * @brief Retrieves a list of vertices belonging to a given face.
 * 
 * @tparam T Floating point type used for coordinates.
 * @tparam dim Dimension of the half-edge data structure.
 * @param face The face for which the vertices are to be retrieved.
 * @return auto A list of vertices belonging to the face.
 */
    template <std::floating_point T, size_t dim>
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

/**
 * @brief Retrieves a list of vertices belonging to a given face.
 * 
 * @tparam T Floating point type used for coordinates.
 * @tparam dim Dimension of the half-edge data structure.
 * @param face A shared pointer to the face for which the vertices are to be retrieved.
 * @return auto A list of vertices belonging to the face.
 */
    template <std::floating_point T, size_t dim>
    auto getFaceVertices(const std::shared_ptr<HalfEdgeFace<T, dim>> &face)
    {
        return getFaceVertices(*face);
    }

/**
 * @brief Retrieves the three vertices of a triangle face.
 * 
 * @tparam T The type of the coordinates (must be a floating-point type).
 * @tparam dim The dimension of the space (must be 2 or 3).
 * @param face The face whose edges to retrieve.
 * @return An array of three HalfEdgeVertices objects representing the vertices of the triangle.
 */
    template <std::floating_point T, size_t dim>
    auto getTriangleVertices(const HalfEdgeFace<T, dim> &face)
    {
        assert(face.edge->next->next->next == face.edge);
        return std::array<std::shared_ptr<HalfEdgeVertex<T, dim>>,3> {
            face.edge->next->vertex,
            face.edge->vertex,
            face.edge->previous->vertex
        };
    }

/**
 * @brief Retrieves the three vertices of a triangle face given a shared pointer to the face.
 * 
 * @tparam T The type of the coordinates (must be a floating-point type).
 * @tparam dim The dimension of the space (must be 2 or 3).
 * @param face A shared pointer to the face whose edges to retrieve.
 * @return An array of three HalfEdgeVertices objects representing the vertices of the triangle.
 */
    template <std::floating_point T, size_t dim>
    auto getTriangleVertices(const std::shared_ptr<HalfEdgeFace<T, dim>> &face)
    {
        return getTriangleVertices(*face);
    }

/**
 * @brief Retrieves a list of coordinates of the vertices belonging to a given face.
 * 
 * @tparam T Floating point type used for coordinates.
 * @tparam dim Dimension of the half-edge data structure.
 * @param face A shared pointer to the face for which the coordinates are to be retrieved.
 * @return auto A list of coordinates of the vertices belonging to the face.
 */
    template <std::floating_point T, size_t dim>
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
 * @brief Retrieves the three points of a triangle face.
 * 
 * @tparam T The type of the coordinates (must be a floating-point type).
 * @tparam dim The dimension of the space (must be 2 or 3).
 * @param face The face whose edges to retrieve.
 * @return An array of three std::array<T,3> objects representing the coordinates of the triangle.
 */
    template <std::floating_point T, size_t dim>
    auto getTriangleCoords(const HalfEdgeFace<T, dim> &face)
    {
        assert(face.edge->next->next->next == face.edge);
        return std::array<HalfEdgeVertex<T, dim>,3> {
            face.edge->next->vertex->coords,
            face.edge->vertex->coords,
            face.edge->previous->vertex->coords
        };
    }

/**
 * @brief Retrieves the three points of a triangle face given a shared pointer to the face.
 * 
 * @tparam T The type of the coordinates (must be a floating-point type).
 * @tparam dim The dimension of the space (must be 2 or 3).
 * @param face A shared pointer to the face whose edges to retrieve.
 * @return An array of three std::array<T,3> objects representing the coordinates of the triangle.
 */
    template <std::floating_point T, size_t dim>
    auto getTriangleCoords(const std::shared_ptr<HalfEdgeFace<T, dim>> &face)
    {
        return getTriangleCoords(*face);
    }
/**
 * @brief Finds the common edge between two faces, if it exists.
 * 
 * @tparam T Floating point type used for coordinates.
 * @tparam dim Dimension of the half-edge data structure.
 * @param h_f1 First face.
 * @param h_f2 Second face.
 * @return auto A shared pointer to the common edge, or nullptr if none exists.
 */
    template <std::floating_point T, size_t dim>
    auto getCommonEdge(const std::shared_ptr<HalfEdgeFace<T, dim>> &h_f1, const std::shared_ptr<HalfEdgeFace<T, dim>> &h_f2) -> std::shared_ptr<HalfEdge<T, dim>>
    {
        auto face_edges = getFaceEdges(h_f1);
            auto common_edge_iter = std::ranges::find_if(face_edges, [&](const auto &h_e) {
                return h_e->opposite && h_e->opposite->face == h_f2;
            });

        return common_edge_iter != face_edges.end() ? *common_edge_iter : nullptr;

    }

/**
 * @brief Finds the common edges between two faces, if they exist.
 * 
 * @tparam T Floating point type used for coordinates.
 * @tparam dim Dimension of the half-edge data structure.
 * @param h_f1 First face.
 * @param h_f2 Second face.
 * @return auto A pair of shared pointers to the common edges, or nullptrs if none exist.
 */
    template <std::floating_point T, size_t dim>
    auto getCommonEdges(const std::shared_ptr<HalfEdgeFace<T, dim>> &h_f1, const std::shared_ptr<HalfEdgeFace<T, dim>> &h_f2) -> std::pair<std::shared_ptr<HalfEdge<T, dim>>, std::shared_ptr<HalfEdge<T, dim>>>
    {
        if (auto h_e1 = getCommonEdge(h_f1, h_f2))
        {
            return std::make_pair(h_e1, h_e1->opposite);
        }
        return std::make_pair(nullptr, nullptr);
    }

/**
 * @brief Gets the previous face, if it exists, from a given edge.
 * 
 * @tparam T Floating point type used for coordinates.
 * @tparam dim Dimension of the half-edge data structure.
 * @param edge The edge from which the previous face is to be found.
 * @return auto A shared pointer to the previous face, or nullptr if none exists.
 */
    template <std::floating_point T, size_t dim>
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

/**
 * @brief Gets the next face, if it exists, from a given edge.
 * 
 * @tparam T Floating point type used for coordinates.
 * @tparam dim Dimension of the half-edge data structure.
 * @param edge The edge from which the next face is to be found.
 * @return auto A shared pointer to the next face, or nullptr if none exists.
 */
    template <std::floating_point T, size_t dim>
    auto getNextFace(const std::shared_ptr<HalfEdge<T, dim>> &edge) -> std::shared_ptr<HalfEdgeFace<T, dim>>
    {
        auto opp = edge->opposite;
        if (opp)
        {
            return opp->face;
        }
        return nullptr;
    }

/**
 * @brief Gets the list of faces attached to a given vertex.
 * 
 * @tparam T Floating point type used for coordinates.
 * @tparam dim Dimension of the half-edge data structure.
 * @param h_v Shared pointer to the vertex.
 * @return auto List of shared pointers to the attached faces.
 */
    template <std::floating_point T, size_t dim>
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

/**
 * @brief Gets the list of neighboring faces of a given face.
 * 
 * @tparam T Floating point type used for coordinates.
 * @tparam dim Dimension of the half-edge data structure.
 * @param h_f Shared pointer to the face.
 * @return auto List of shared pointers to the neighboring faces.
 */
    template <std::floating_point T, size_t dim>
    auto getNeighboringFaces(const std::shared_ptr<HalfEdgeFace<T, dim>> &h_f)
    {
        std::list<std::shared_ptr<HalfEdgeFace<T, dim>>> neighbors;
        auto edges = getFaceEdges(h_f);

        // Transform the list of edges to a list of neighboring faces
        std::transform(edges.begin(), edges.end(), std::back_inserter(neighbors), [](const auto &h_e) {
            return h_e->opposite ? h_e->opposite->face : nullptr;
        });

        // Remove null pointers from the list
        neighbors.remove(nullptr);

        return neighbors;
    }

/**
 * @brief Gets the list of boundary edges of a given list of faces.
 * 
 * @tparam T Floating point type used for coordinates.
 * @tparam dim Dimension of the half-edge data structure.
 * @param h_f_lst List of shared pointers to the faces.
 * @return auto List of shared pointers to the boundary edges.
 */
    template <std::floating_point T, size_t dim>
    auto getFacesBoundary(const std::list<std::shared_ptr<HalfEdgeFace<T, dim>>> &h_f_lst)
    {
        std::list<std::shared_ptr<HalfEdge<T, dim>>> boundary;

        for (const auto &h_f : h_f_lst)
        {
            const auto &h_e_lst = getFaceEdges(h_f);
            for (const auto &h_e : h_e_lst)
            {
                
                if (
                    !h_e->opposite || 
                    std::ranges::none_of(
                        h_f_lst, 
                        [&](const auto &face){ return face == h_e->opposite->face; })
                    )
                {
                    boundary.push_back(h_e);
                }
            }
        }

        return boundary;
    }

/**
 * @brief Retrieves the oriented boundaries of a list of faces.
 * 
 * @tparam T Floating point type used for coordinates.
 * @tparam dim Dimension of the half-edge data structure.
 * @param h_f_lst List of shared pointers to faces.
 * @return auto List of lists containing the oriented boundaries.
 */
    template <std::floating_point T, size_t dim>
    auto getOrientedFacesBoundaries(const std::list<std::shared_ptr<HalfEdgeFace<T, dim>>> &h_f_lst)
    {
        auto boundary = getFacesBoundary(h_f_lst);
        return takeClosedLoops(boundary);
    }

/**
 * @brief Retrieves the oriented boundary of a list of faces.
 * 
 * @tparam T Floating point type used for coordinates.
 * @tparam dim Dimension of the half-edge data structure.
 * @param h_f_lst List of shared pointers to faces.
 * @return auto List containing the oriented boundary.
 */
    template <std::floating_point T, size_t dim>
    auto getOrientedFacesBoundary(const std::list<std::shared_ptr<HalfEdgeFace<T, dim>>> &h_f_lst)
    {
        auto boundary = getFacesBoundary(h_f_lst);
        return takeClosedLoops(boundary).front();
    }

/**
 * @brief Creates a map of unique vertices and their indices from a list of faces.
 * 
 * @tparam T Floating point type used for coordinates.
 * @tparam dim Dimension of the half-edge data structure.
 * @param faces_lst List of shared pointers to faces.
 * @return auto Map of shared pointers to HalfEdgeVertex and their indices.
 */
    template <std::floating_point T, size_t dim>
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

/**
 * @brief Creates a vector of unique vertices from a list of faces.
 * 
 * @tparam T Floating point type.
 * @tparam dim Dimension of the vertices.
 * @param faces_lst List of half-edge faces.
 * @return std::vector<std::shared_ptr<HalfEdgeVertex<T, dim>>> A vector containing unique vertices.
 */
    template <std::floating_point T, size_t dim>
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

/**
 * @brief Creates a map of unique half-edges and their indices from a list of faces.
 * 
 * @tparam T Floating point type.
 * @tparam dim Dimension of the vertices.
 * @param faces_lst List of half-edge faces.
 * @return std::map<std::shared_ptr<HalfEdge<T, dim>>, size_t> A map containing unique half-edges and their indices.
 */
    template <std::floating_point T, size_t dim>
    auto getEdgesMap(const auto &faces_lst) -> std::map<std::shared_ptr<HalfEdge<T, dim>>, size_t>
    {
        std::map<std::shared_ptr<HalfEdge<T, dim>>, size_t> edges_map;
        size_t index{};
        for (const auto &f : faces_lst)
        {
            auto hed_lst = getFaceEdges(f);
            for (const auto &hed : hed_lst)
            {
                if (!edges_map.contains(hed))
                {
                    edges_map[hed] = index;
                    index++;
                }
            }
        }

        return edges_map;
    }

/**
 * @brief Computes a point along the half-edge at a given position.
 * 
 * @tparam T Floating point type.
 * @tparam dim Dimension of the vertices.
 * @param he Shared pointer to the half-edge.
 * @param pos Position of the point along the half-edge (default is 0.5, i.e., midpoint).
 * @return std::array<T, dim> A point along the half-edge.
 */
    template <std::floating_point T, size_t dim>
    auto getEdgePoint(const std::shared_ptr<HalfEdge<T,dim>> &he, T pos = 0.5)
    {
        assert(he && he->vertex && he->previous && he->previous->vertex);
        const auto &p1 = he->previous->vertex->coord;
        const auto &p2 = he->vertex->coord;
        return p1 + pos *(p2-p1);
    }

/**
 * @brief Generates an encompassing mesh around a set of 2D points.
 * 
 * @tparam T Floating point type.
 * @param X_lst Vector of 2D points.
 * @param pc_offset Percentage offset for the encompassing mesh (default is 10%).
 * @return std::list<std::shared_ptr<HalfEdgeFace<T, 2>>> A list of shared pointers to the generated half-edge faces.
 */
    template <std::floating_point T>
    auto getEncompassingMesh(const std::vector<std::array<T, 2>> &X_lst, T pc_offset = 10)
    {
        // Get the minimum and maximum coordinates of the points
        auto [Xmin, Xmax] = getCoordsMinMax(X_lst);

        // Calculate the percentage offsets for the encompassing mesh
        auto pc = pc_offset / 100;
        Xmax = Xmax + pc * (Xmax - Xmin);
        Xmin = Xmin - pc * (Xmax - Xmin);

        // Create half-edges for the encompassing mesh
        auto he1 = make_shared_h_edge<T, 2>({Xmax[0], Xmin[1]});
        auto he2 = make_shared_h_edge<T, 2>({Xmin[0], Xmax[1]});
        auto he3 = make_shared_h_edge<T, 2>(Xmin);

        // Create the first face of the encompassing mesh
        auto lst1 = {he1, he2, he3};
        auto hf1 = make_shared_h_face<T, 2>(lst1);

        // Create the second face and link it to the first face
        auto hf2 = add_face(hf1, he2, Xmax);

        // Return the encompassing mesh as a list of shared pointers to half-edge faces
        return std::list<std::shared_ptr<HalfEdgeFace<T, 2>>>{hf1, hf2};
    }
}