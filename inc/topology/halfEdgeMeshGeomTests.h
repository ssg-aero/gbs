#pragma once
#include "halfEdgeMeshData.h"
#include "baseGeom.h"

namespace gbs
{
/**
 * Check if a point lies inside the circumcircle of a triangle.
 * @tparam T Floating-point type of the coordinates
 * @param pt Coordinates of the point
 * @param t Face of a 2D triangular mesh (must have exactly three vertices)
 * @return Positive value if the point is inside the circumcircle, negative if it's outside, and zero if it's on the circle
 */
    template <std::floating_point T>
    T in_circle(const std::array<T,2> &pt, const HalfEdgeFace<T, 2> &t)
    {
        // Ensure the face has exactly three vertices
        assert(getFaceVertices(t).size() == 3);

        // Extract the three vertices of the triangle
        auto he = t.edge;
        const auto &a = he->vertex->coords;
        he = he->next;
        const auto &b = he->vertex->coords;
        he = he->next;
        const auto &c = he->vertex->coords;

        // Check if the point lies inside the circumcircle of the triangle
        return in_circle(a, b, c, pt);
    }
/**
 * @brief Determines if a point is inside, outside, or on the circumcircle of a triangle.
 * 
 * @tparam T The floating-point type used to represent the point and the triangle.
 * @param pt A 2D array representing the point.
 * @param t A shared pointer to a HalfEdgeFace object representing the triangle.
 * @return The result of calling the helper function in_circle with the same arguments.
 */
    template <std::floating_point T>
    T in_circle(const std::array<T,2> &pt, const std::shared_ptr<HalfEdgeFace<T, 2>> &t)
    {
        assert(t);
        return in_circle(pt,*t);
    }
/**
 * @brief Checks if the vertices of a given HalfEdgeFace are in counter-clockwise order.
 * 
 * @tparam T Floating point type used for coordinates.
 * @param hf A shared pointer to a HalfEdgeFace.
 * @return A positive value if the vertices are in counter-clockwise order,
 *         a negative value if in clockwise order, and zero if the points are collinear.
 */
    template <std::floating_point T>
    T is_ccw(const std::shared_ptr<HalfEdgeFace<T, 2>> &hf)
    {
        auto vertices = getFaceVertices(hf);
        assert(vertices.size() == 3);

        return are_ccw(
            (*std::next(vertices.begin(), 0))->coords,
            (*std::next(vertices.begin(), 1))->coords,
            (*std::next(vertices.begin(), 2))->coords
        );
    }
/**
 * @brief Determines if a point is inside a polygon defined by a list of half-edges.
 * 
 * @tparam T The floating point type used for coordinates.
 * @tparam Container The container type for the list of half-edges.
 * @param xy The point to be tested.
 * @param h_e_lst The list of half-edges defining the polygon.
 * @return true If the point is inside the polygon or on its edges.
 * @return false If the point is outside the polygon.
 */
    template <std::floating_point T, std::ranges::range Container>
    bool is_inside(const std::array<T, 2> &xy, const Container &h_e_lst)
    {
        size_t count{};
        for (const auto &he : h_e_lst)
        {
            const auto &a = he->previous->vertex->coords;
            const auto &b = he->vertex->coords;
            
            // Check if the point is on the polygon's edges
            if (on_seg(a, b, xy))
            {
                return true;
            }
            
            // Check if the segment intersects with the half-open line segment
            if (seg_H_strict_end_intersection(a, b, xy))
            {
                count++;
            }
        }

        // When count is odd, the point is inside the polygon
        return count % 2 == 1;
    }
/**
 * @brief Determines if a half-edge belongs to a locally Delaunay triangulation.
 * 
 * @tparam T The floating point type used for coordinates.
 * @param he The half-edge to be checked.
 * @return A positive value if the half-edge is locally Delaunay.
 * @return Zero if the half-edge is on the boundary.
 * @return A negative value if the half-edge is not locally Delaunay.
 */
    template <std::floating_point T>
    T is_locally_delaunay(const std::shared_ptr<HalfEdge<T, 2>> &he)
    {
        assert(he);

        // If the half-edge is on the boundary, it is considered locally Delaunay
        if (!he->opposite) return 0;

        auto t1 = he->face;
        auto t2 = he->opposite->face;

        assert(t1);
        assert(t2);

        auto vertices1 = getFaceVertices(t1);
        auto vertices2 = getFaceVertices(t2);

        assert(vertices1.size() == 3);
        assert(vertices2.size() == 3);

        auto p1 = he->next->vertex->coords;
        auto p2 = he->opposite->next->vertex->coords;

        // Returns the minimum in_circle value between the two opposite vertices
        return std::min(in_circle(p1, t2), in_circle(p2, t1));
    }
/**
 * @brief Checks if all the faces in a list of faces have a counter-clockwise orientation.
 * 
 * @tparam Container The container type for the list of faces.
 * @param faces_lst The list of faces to be checked.
 * @return true if all the faces have a counter-clockwise orientation.
 * @return false if at least one face does not have a counter-clockwise orientation.
 */
    template <std::ranges::range Container>
    bool are_face_ccw(const Container &faces_lst)
    {
        return std::all_of(faces_lst.begin(), faces_lst.end(),
                        [](const auto &hf) { return is_ccw(hf) > 0; });
    }
/**
 * @brief Checks if the edges in a list of 2D half-edges are in counter-clockwise order.
 * 
 * @tparam HalfEdgeContainer The container type for the list of 2D half-edges.
 * @param edges_lst The list of 2D half-edges.
 * @return true If the edges are in counter-clockwise order.
 * @return false If the edges are not in counter-clockwise order.
 */
    template <typename HalfEdgeContainer>
    bool are_edges_2d_ccw(const HalfEdgeContainer &edges_lst)
    {
        // Calculate the sum of edge cross-products, considering the vertices' coordinates
        auto sum = std::reduce(
            edges_lst.begin(), edges_lst.end(),
            0.,
            [](auto t, const auto &he) {
                auto [x2, y2] = he->vertex->coords;
                auto [x1, y1] = he->previous->vertex->coords;
                return t + (x2 - x1) * (y2 + y1);
            });

        // Check if the sum is negative, indicating that the edges are in counter-clockwise order
        return sum < 0.;
    }
/**
 * @brief Checks if the edges of a face intersect with a given half-edge.
 * 
 * @tparam T Floating point type used for coordinates.
 * @param h_e The half-edge to check for intersection with the face edges.
 * @param h_f The face to check for intersection with the half-edge.
 * @return true If the half-edge intersects any of the face edges.
 * @return false If the half-edge does not intersect any of the face edges.
 */
    template <std::floating_point T>
    bool areFacesEdgesIntersect(const gbs::HalfEdge<T, 2> &h_e, const gbs::HalfEdgeFace<T, 2> &h_f)
    {
        // Ensure the input half-edge has a previous edge
        assert(h_e.previous);

        for (const auto &h_e_ : gbs::getFaceEdges(h_f))
        {
            // Ensure the face edge has a previous edge
            assert(h_e_->previous);

            // Check if the input half-edge intersects with the current face edge
            if (seg_seg_strict_intersection(
                    h_e.previous->vertex->coords,
                    h_e.vertex->coords,
                    h_e_->previous->vertex->coords,
                    h_e_->vertex->coords))
            {
                return true;
            }
        }

        return false;
    }
/**
* @brief Checks if the centroid of a given face is inside the specified boundary.
* 
* @tparam T Floating point type used for coordinates.
* @param face The face to test.
* @param boundary Container holding the boundary information.
* @return bool True if the centroid is inside the boundary, false otherwise.
*/
    template <std::floating_point T, typename _Container>
    bool is_centroid_inside_boundary(const std::shared_ptr<HalfEdgeFace<T, 2>> &face, const _Container &boundary)
    {
        // Calculate the centroid of the face
        auto coords = getFaceCoords(face);
        std::array<T, 2> G{};
        for (const auto &xy : coords)
        {
            G = G + xy;
        }
        G = G / static_cast<T>(coords.size());

        // Check if the centroid is inside the boundary
        return is_inside(G, boundary);
    }
} // namespace gbs