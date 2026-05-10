#pragma once

#include <array>
#include <vector>
#include <algorithm>
#include <span>
#include <unordered_set>

#include <gbs/surfaces>
#include <gbs/bscanalysis.h>
#include <gbs/bssanalysis.h>
#include <gbs/bscbuild.h>

#include "halfEdgeMeshData.h"
#include "halfEdgeMeshGetters.h"
#include "halfEdgeMeshEditors.h"
#include "baseIntersection.h"
#include "halfEdgeMeshGeomTests.h"
#include "halfEdgeMeshSmoothing.h"
#include "edgeRecovery.h"
#include "baseGeom.h"

namespace gbs
{
/**
 * @brief Creates a list of 2D boundary points for a surface mesh.
 *
 * @tparam T A floating-point type used for coordinates and deviation values.
 * @tparam dim The dimension of the surface.
 * @param srf A surface object.
 * @param nu The number of u-direction sampling points (default is 5).
 * @param nv The number of v-direction sampling points (default is 5).
 * @param deviation The maximum allowed deviation for boundary point sampling (default is 0.01).
 * @return A vector of 2D arrays representing the boundary points of the surface mesh.
 */
    template < std::floating_point T, size_t dim>
    auto meshSurfaceBoundary(const Surface<T,dim> &srf, size_t nu = 5, size_t nv = 5, T deviation = 0.01)
    {
        using std::begin;
        using std::end;
        using std::transform;

        auto [u1_, u2_, v1_, v2_] = srf.bounds();

        auto u1 = deviation_based_u_params(srf, u1_, u2_, v1_, nu, deviation);
        auto u2 = deviation_based_u_params(srf, u1_, u2_, v2_, nu, deviation);
        std::reverse(begin(u2), end(u2));
        auto v1 = deviation_based_v_params(srf, v1_, v2_, u1_, nv, deviation);
        std::reverse(begin(v1), end(v1));
        auto v2 = deviation_based_v_params(srf, v1_, v2_, u2_, nv, deviation);

        auto n = u1.size() + u2.size() + v1.size() + v2.size();
        std::vector<std::array<T, 2>> coords(n);

        auto coords_begin = coords.begin();
        coords_begin = transform(begin(u1), end(u1), coords_begin, [v1_](auto u) { return std::array<T, 2>{u, v1_}; });
        coords_begin = transform(++v2.begin(), end(v2), coords_begin, [u2_](auto v) { return std::array<T, 2>{u2_, v}; });
        coords_begin = transform(begin(u2), end(u2), coords_begin, [v2_](auto u) { return std::array<T, 2>{u, v2_}; });
        transform(++v1.begin(), end(v1), coords_begin, [u1_](auto v) { return std::array<T, 2>{u1_, v}; });

        return coords;
    }

/**
 * @brief Inserts a point into a Delaunay triangulation using the Boyer-Watson algorithm.
 *
 * @tparam T A floating-point type used for coordinates and tolerance values.
 * @tparam Container A container type representing the list of HalfEdgeFaces.
 * @param h_f_lst A reference to the container of HalfEdgeFaces forming the initial Delaunay triangulation.
 * @param xy An array containing the coordinates of the point to be inserted.
 * @param tol A floating-point value used as tolerance for the Delaunay condition (default is 1e-10).
 * @return A tuple ( new vertex, deleted HalfEdgeFaces that were violating the Delaunay condition, the newly created HalfEdgeFaces).
 */
    template<std::floating_point T, typename Container>
    auto boyerWatson(Container& h_f_lst, const std::array<T, 2>& xy, T tol = T{1e-10}) {
        using std::begin;
        using std::end;

        auto it{begin(h_f_lst)};
        Container h_f_lst_deleted;

        // Find triangle violation Delaunay condition
        while (it != end(h_f_lst)) {
            it = find_if(it, end(h_f_lst), [xy, tol](const auto& h_f) {
                return in_circle(xy, h_f) > tol;
            });
            if (it != end(h_f_lst)) {
                h_f_lst_deleted.push_back(*it);
                it = h_f_lst.erase(it);
            }
        }

        if (h_f_lst_deleted.empty()) {
            return std::tuple<std::shared_ptr<HalfEdgeVertex<T, 2>>,Container, Container>{}; // if point is outside or confused with an existing point
        }

        assert(are_face_ccw(h_f_lst));

        // Get cavity boundary
        auto h_e_lst = getOrientedFacesBoundary(h_f_lst_deleted);
        assert(are_edges_2d_ccw<T>(h_e_lst));

        // fill cavity
        auto vtx = make_shared_h_vertex(xy);
        auto h_f_lst_new = add_vertex(h_e_lst, vtx);
        assert(are_face_ccw(h_f_lst_new));

        // append new faces
        h_f_lst.insert(end(h_f_lst), begin(h_f_lst_new), end(h_f_lst_new));

        return std::make_tuple(vtx, std::move(h_f_lst_deleted), std::move(h_f_lst_new));
    }

/**
 * @brief Computes the Delaunay triangulation of a 2D point set using the Boyer-Watson algorithm.
 *
 * @tparam T A floating-point type used for coordinates and tolerance.
 * @param coords A container of 2D coordinates representing the input points.
 * @param tol A floating-point value used as tolerance for the Delaunay condition.
 * @return A container of triangulated faces forming the Delaunay triangulation.
 */
    template < std::floating_point T>
    auto delaunay2DBoyerWatson(const auto &coords, T tol)
    {
        auto faces_lst = getEncompassingMesh(coords);
        auto vertices = getVerticesVectorFromFaces<T,2>(faces_lst);
        // insert points
        for(const auto &xy : coords)
        {
            boyerWatson<T>(faces_lst, xy, tol);
        }

        // remove external mesh, i.ei faces attached to initial vertices
        for(const auto &vtx : vertices)
        {
            remove_faces(faces_lst, vtx);
        }

        return faces_lst;
    }
/**
 * @brief Computes the Delaunay triangulation of a 2D point set with inner boundaries using the Boyer-Watson algorithm.
 *
 * @tparam T A floating-point type used for coordinates and tolerance.
 * @param coords_boundary A container of 2D coordinates representing the outer boundary input points.
 * @param coords_inner A container of 2D coordinates representing the inner boundary input points.
 * @param tol A floating-point value used as tolerance for the Delaunay condition.
 * @return A container of triangulated faces forming the Delaunay triangulation with inner boundaries.
 */
    template < std::floating_point T >
    auto delaunay2DBoyerWatson(const auto &coords_boundary,const auto &coords_inner, T tol)
    {
        auto faces_lst = delaunay2DBoyerWatson(coords_boundary, tol);
        auto boundary_inner = make_HalfEdges<T>(coords_inner);
        takeInternalFaces<T>(faces_lst,boundary_inner);
        return faces_lst;
    }
/**
 * @brief Computes the Delaunay triangulation of a 2D point set using the Boyer-Watson algorithm.
 *
 * This function is an overload that accepts a specific container type (std::vector<std::array<T, 2>>) for input points.
 *
 * @tparam T A floating-point type used for coordinates and tolerance.
 * @param coords A std::vector of 2D coordinates (std::array<T, 2>) representing the input points.
 * @param tol A floating-point value used as tolerance for the Delaunay condition (default is 1e-10).
 * @return A container of triangulated faces forming the Delaunay triangulation.
 */
    template < std::floating_point T>
    auto delaunay2DBoyerWatson(const std::vector< std::array<T,2> > &coords, T tol = 1e-10)
    {
        return delaunay2DBoyerWatson<T,std::vector< std::array<T,2> >>(coords, tol);
    }

/**
 * @brief Computes a 2D Constrained Delaunay Triangulation (CDT) of a CCW polygon.
 *
 * The input is a closed simple polygon, given as its vertices in CCW order
 * (closure is implicit — do not repeat the first vertex). The function:
 *   1. builds a Boyer-Watson Delaunay triangulation of the vertices,
 *   2. recovers each polygon edge by edge-flipping (Anglada / Sloan),
 *   3. flood-fills from each constrained half-edge and keeps only the faces
 *      that lie on the interior side (left of every CCW polygon edge).
 *
 * Holes and arbitrary internal constraints are not handled in this first cut.
 * The polygon must be simple (no self-intersection) and no polygon edge may
 * pass through a non-incident vertex.
 *
 * @tparam T Floating-point type.
 * @param coords Polygon vertices in CCW order.
 * @param tol Numerical tolerance.
 * @return Face list of the interior triangles.
 */
    template <std::floating_point T>
    auto delaunay2DConstrained(
        const std::vector<std::array<T, 2>> &outer,
        const std::vector<std::vector<std::array<T, 2>>> &holes,
        T tol = T{1e-10})
    {
        // Each polygon ring is described by a vector of vertices in domain-CCW
        // order: the outer polygon walks CCW around the domain; each hole, as
        // a polygon, is also CCW around itself, but its domain-interior side
        // is OUTSIDE the hole. To get a uniform "interior is always on the
        // left of every constrained edge" invariant — which the flood-fill
        // below relies on — we walk holes in reverse (CW around themselves =
        // CCW around the domain region they bound).
        std::vector<std::vector<std::array<T, 2>>> rings;
        rings.reserve(holes.size() + 1);
        rings.push_back(outer);
        for (const auto &h : holes)
        {
            rings.push_back({h.rbegin(), h.rend()});
        }

        // 1. Boyer-Watson on all ring vertices, tracking which vertex was
        //    created at each insertion so we can name boundary edges later.
        std::vector<std::array<T, 2>> all_coords;
        std::size_t total = 0;
        for (const auto &r : rings) total += r.size();
        all_coords.reserve(total);
        for (const auto &r : rings)
            all_coords.insert(all_coords.end(), r.begin(), r.end());

        auto faces_lst = getEncompassingMesh(all_coords);
        auto super_vertices = getVerticesVectorFromFaces<T, 2>(faces_lst);

        std::vector<std::vector<std::shared_ptr<HalfEdgeVertex<T, 2>>>> ring_vtx(rings.size());
        for (std::size_t r = 0; r < rings.size(); ++r)
        {
            ring_vtx[r].resize(rings[r].size());
            for (std::size_t i = 0; i < rings[r].size(); ++i)
            {
                auto [vtx, deleted, created] = boyerWatson<T>(faces_lst, rings[r][i], tol);
                ring_vtx[r][i] = vtx;
            }
        }

        for (const auto &vtx : super_vertices)
            remove_faces(faces_lst, vtx);
        // remove_faces only severs opposite linkages; vertex->edge pointers can
        // still reference half-edges whose face was deleted. Re-sync them
        // before any navigation (recoverEdge, findHalfEdge) runs.
        repairVertexEdges<T, 2>(faces_lst);

        // 2. Recover constrained edges (outer + holes); collect both directions
        //    of each in a set so the flood-fill can refuse to cross them.
        std::unordered_set<HalfEdge<T, 2> *> constrained;
        for (const auto &vtx : ring_vtx)
        {
            const std::size_t n = vtx.size();
            for (std::size_t i = 0; i < n; ++i)
            {
                auto v_p = vtx[i];
                auto v_q = vtx[(i + 1) % n];
                if (!v_p || !v_q) continue;
                if (auto e = recoverEdge<T>(faces_lst, v_p, v_q, tol))
                {
                    constrained.insert(e.get());
                    if (e->opposite)
                        constrained.insert(e->opposite.get());
                }
            }
        }

        // 3. Flood-fill interior. Seed with each ring half-edge's left face
        //    (CCW orientation ⇒ domain interior on the left of every
        //    constrained edge — true for the outer ring and, after the
        //    reversal above, for each hole as well).
        std::unordered_set<HalfEdgeFace<T, 2> *> interior;
        std::list<std::shared_ptr<HalfEdgeFace<T, 2>>> queue;
        for (const auto &vtx : ring_vtx)
        {
            const std::size_t n = vtx.size();
            for (std::size_t i = 0; i < n; ++i)
            {
                auto e = findHalfEdge(vtx[i], vtx[(i + 1) % n]);
                if (!e || !e->face) continue;
                if (interior.insert(e->face.get()).second)
                    queue.push_back(e->face);
            }
        }
        while (!queue.empty())
        {
            auto f = queue.front();
            queue.pop_front();
            for (const auto &ee : getTriangleEdges(f))
            {
                if (constrained.count(ee.get())) continue;
                if (!ee->opposite) continue;
                auto fn = ee->opposite->face;
                if (fn && interior.insert(fn.get()).second)
                    queue.push_back(fn);
            }
        }

        // 4. Drop exterior faces (includes faces strictly inside each hole).
        faces_lst.remove_if([&interior](const auto &f) {
            return interior.find(f.get()) == interior.end();
        });

        return faces_lst;
    }

    /// Convenience overload: simple polygon, no holes.
    template <std::floating_point T>
    auto delaunay2DConstrained(const std::vector<std::array<T, 2>> &outer, T tol = T{1e-10})
    {
        return delaunay2DConstrained<T>(outer, std::vector<std::vector<std::array<T, 2>>>{}, tol);
    }
/**
 * @brief Refines a 2D Delaunay triangulation of a surface mesh using the Boyer-Watson algorithm.
 *
 * @tparam T A floating-point type used for coordinates, tolerance, and criteria values.
 * @tparam dim A size_t value representing the dimension of the surface.
 * @tparam _Func A callable object type used for distance calculation between mesh and surface.
 * @param srf A surface on which the mesh refinement is performed.
 * @param faces_lst A reference to the container of triangulated faces forming the surface mesh.
 * @param crit_max A floating-point value specifying the maximum allowed distance criterion for refinement.
 * @param max_inner_points A size_t value specifying the maximum number of inner points allowed for refinement (default is 500).
 * @param tol A floating-point value used as tolerance for the Delaunay condition (default is 1e-10).
 * @return A container of triangulated faces forming the refined Delaunay triangulation of the surface mesh.
 */
    template <std::floating_point T, size_t dim, typename _Func>
    auto delaunay2DBoyerWatsonSurfaceMeshRefine(const Surface<T, dim> &srf, auto &faces_lst, T crit_max, size_t max_inner_points = 500, T tol = 1e-10)
    {
        _Func dist_mesh_srf{srf};

        // Store face quality in a map
        std::unordered_map<std::shared_ptr<HalfEdgeFace<T, 2>>, T> face_quality;

        // Function to update face quality map
        auto update_face_quality = [&dist_mesh_srf, &face_quality](const auto &hf) {
            face_quality[hf] = dist_mesh_srf(hf).first;
        };

        // Initialize face quality map
        for (const auto &face : faces_lst)
        {
            update_face_quality(face);
        }

        // Refine mesh to match max criteria
        T crit = std::numeric_limits<T>::max();
        for (size_t i{}; i < max_inner_points && crit > crit_max; i++)
        {
            auto it = std::max_element(
                GBS_PAR_EXEC
                faces_lst.begin(), faces_lst.end(),
                [&face_quality](const auto &hf1, const auto &hf2) {
                    return face_quality[hf1] < face_quality[hf2];
                });

            auto d_G_UV = dist_mesh_srf(*it);
            crit = d_G_UV.first;

            auto [vtx, deleted_faces_list, new_faces_list] = boyerWatson(faces_lst, d_G_UV.second, tol);

            for (const auto &face : deleted_faces_list)
            {
                face_quality.erase(face);
            }

            for (const auto &face : new_faces_list)
            {
                update_face_quality(face);
            }
        }

        return faces_lst;
    }


    template <std::floating_point T, size_t dim, typename _Func>
    auto delaunay2DBoyerWatsonMeshRefine(auto &faces_lst, T crit_max, const _Func &dist_mesh, size_t max_inner_points = 500, T tol = 1e-10)
    {
        // Store face quality in a map
        std::unordered_map<std::shared_ptr<HalfEdgeFace<T, dim>>, T> face_quality;

        // Function to update face quality map
        auto update_face_quality = [&dist_mesh, &face_quality](const auto &hf) {
            face_quality[hf] = dist_mesh(hf).first;
        };

        // Initialize face quality map
        for (const auto &face : faces_lst)
        {
            update_face_quality(face);
        }

        // Refine mesh to match max criteria
        T crit = std::numeric_limits<T>::max();
        for (size_t i{}; i < max_inner_points && crit >= crit_max; i++)
        {
            auto it = std::max_element(
                GBS_PAR_EXEC
                faces_lst.begin(), faces_lst.end(),
                [&face_quality](const auto &hf1, const auto &hf2) {
                    return face_quality[hf1] < face_quality[hf2];
                });

            auto d_G_UV = dist_mesh(*it);
            crit = d_G_UV.first;

            // std::cout << crit << " inserting: " << d_G_UV.second[0] << " " << d_G_UV.second[1] << std::endl;

            auto [vtx, deleted_faces_list, new_faces_list] = boyerWatson(faces_lst, d_G_UV.second, tol);

            for (const auto &face : deleted_faces_list)
            {
                face_quality.erase(face);
            }

            for (const auto &face : new_faces_list)
            {
                update_face_quality(face);
            }

            // Smooth mesh
            laplacian_smoothing(faces_lst, *vtx);
            // Restore Delaunay condition
            auto edges = getEdgesMap<T,dim>(faces_lst);
            for(auto & [ed,i] : edges)
            {
                if(ed->opposite && !is_locally_delaunay(ed))
                {
                    flip(ed->face,ed->opposite->face);
                }
            }
        }

        return faces_lst;
    }

    template <std::floating_point T, size_t dim, typename _Func>
    auto delaunay2DBoyerWatsonMeshRefine(auto &faces_lst, T crit_max, size_t max_inner_points = 500, T tol = 1e-10)
    {
        _Func dist_mesh{};
        return delaunay2DBoyerWatsonMeshRefine<T,dim,_Func>(faces_lst, crit_max, dist_mesh, max_inner_points, tol);
    }

/**
 * @brief Computes the base 2D Delaunay triangulation of a surface using the Boyer-Watson algorithm.
 *
 * @tparam T A floating-point type used for coordinates, tolerance, and deviation values.
 * @tparam dim A size_t value representing the dimension of the surface.
 * @param srf A surface for which the base triangulation is computed.
 * @param nu A size_t value specifying the number of divisions along the U direction (default is 5).
 * @param nv A size_t value specifying the number of divisions along the V direction (default is 5).
 * @param deviation A floating-point value specifying the deviation for mesh surface boundary (default is 0.01).
 * @param tol A floating-point value used as tolerance for the Delaunay condition (default is 1e-10).
 * @return A container of triangulated faces forming the base Delaunay triangulation of the surface.
 */
    template < std::floating_point T, size_t dim>
    auto delaunay2DBoyerWatsonSurfaceBase(const Surface<T,dim> &srf, size_t nu = 5, size_t nv = 5, T deviation = 0.01, T tol = 1e-10)
    {
    
        auto coords = meshSurfaceBoundary(srf, nu ,nv, deviation);
        auto faces_lst = getEncompassingMesh(coords);
        auto vertices = getVerticesVectorFromFaces<T,2>(faces_lst);
        // insert boundary points
        for(const auto &xy : coords)
        {
            boyerWatson<T>(faces_lst, xy, tol);
        }
        // remove external mesh, i.ei faces attached to initial vertices
        // This work because surface coordinate are convex ( rectangle )
        for(const auto &vtx : vertices)
        {
            remove_faces(faces_lst, vtx);
        }
        return faces_lst;
    }
/**
 * @brief Computes the refined 2D Delaunay triangulation of a surface mesh using the Boyer-Watson algorithm.
 *
 * @tparam T A floating-point type used for coordinates, tolerance, deviation, and criteria values.
 * @tparam dim A size_t value representing the dimension of the surface.
 * @tparam _Func A callable object type used for distance calculation between mesh and surface.
 * @param srf A surface for which the refined triangulation is computed.
 * @param crit_max A floating-point value specifying the maximum allowed distance criterion for refinement.
 * @param max_inner_points A size_t value specifying the maximum number of inner points allowed for refinement (default is 500).
 * @param nu A size_t value specifying the number of divisions along the U direction (default is 5).
 * @param nv A size_t value specifying the number of divisions along the V direction (default is 5).
 * @param deviation A floating-point value specifying the deviation for mesh surface boundary (default is 0.01).
 * @param tol A floating-point value used as tolerance for the Delaunay condition (default is 1e-10).
 * @return A container of triangulated faces forming the refined Delaunay triangulation of the surface mesh.
 */
    template < std::floating_point T, size_t dim, typename _Func >
    auto delaunay2DBoyerWatsonSurfaceMesh(const Surface<T,dim> &srf, T crit_max, size_t max_inner_points = 500, size_t nu = 5, size_t nv = 5, T deviation = 0.01, T tol = 1e-10)
    {
        auto faces_lst = delaunay2DBoyerWatsonSurfaceBase(srf, nu, nv, deviation, tol);
        return delaunay2DBoyerWatsonSurfaceMeshRefine<T, dim, _Func>(srf, faces_lst, crit_max, max_inner_points, tol); 
    }
/**
 * @brief Adds an inner boundary to an existing 2D Delaunay triangulation using the Boyer-Watson algorithm.
 *
 * @tparam T A floating-point type used for coordinates and tolerance values.
 * @tparam dim A size_t value representing the dimension of the surface.
 * @param faces_lst A reference to the container of triangulated faces forming the initial Delaunay triangulation.
 * @param coords_inner A vector containing the coordinates of the inner boundary points.
 * @param tol A floating-point value used as tolerance for the Delaunay condition (default is 1e-10).
 * @return A container of triangulated faces forming the Delaunay triangulation with the added inner boundary.
 */
    template < std::floating_point T, size_t dim>
    auto delaunay2DBoyerWatsonAddInnerBound(auto &faces_lst, const std::vector<std::array<T,dim>> &coords_inner, T tol = 1e-10)
    {
        for (const auto &xy : coords_inner)
        {
            boyerWatson<T>(faces_lst, xy, tol);
        }
        auto boundary_inner = make_HalfEdges<T>(coords_inner);
        return takeInternalFaces<T>(faces_lst,boundary_inner, tol);
    }
/**
 * @brief Adds an outer boundary to an existing 2D Delaunay triangulation using the Boyer-Watson algorithm.
 *
 * @tparam T A floating-point type used for coordinates and tolerance values.
 * @tparam dim A size_t value representing the dimension of the surface.
 * @param faces_lst A reference to the container of triangulated faces forming the initial Delaunay triangulation.
 * @param coords_outer A vector containing the coordinates of the outer boundary points.
 * @param tol A floating-point value used as tolerance for the Delaunay condition (default is 1e-10).
 * @return A container of triangulated faces forming the Delaunay triangulation with the added outer boundary.
 */
    template <std::floating_point T, size_t dim>
    auto delaunay2DBoyerWatsonAddOuterBound(auto &faces_lst, const std::vector<std::array<T,dim>> &coords_outer, T tol = 1e-10)
    {
        for (const auto &xy : coords_outer)
        {
            boyerWatson<T>(faces_lst, xy, tol);
        }
        auto boundary_outer = make_HalfEdges<T>(coords_outer);
        return takeExternalFaces<T>(faces_lst,boundary_outer, tol);
    }

// TODO remove
    template < typename T , auto _ExPo = gbs::seq_exec>
    auto base_delaunay2d_mesh(std::vector< std::shared_ptr< HalfEdgeVertex<T,2> > > &vertices_cloud)
    {
        auto n = vertices_cloud.size();

        std::list< std::shared_ptr< HalfEdgeFace<T,2> > > faces_lst;

        for(size_t i{}; i < n; i++)
        {
            auto a = vertices_cloud[i];
            // for(size_t j{i+1}; j < n; j++)
            for(size_t j{}; j < n; j++)
            {
                auto b = vertices_cloud[j];
                for( size_t k{j+1}; k < n ; k++)
                // for( size_t k{}; k < n ; k++)
                {
                    auto c   = vertices_cloud[k];

                    if(orient_2d(a->coords,b->coords,c->coords) > 0.) // avoids also null triangle
                    {
                    
                        auto it_inside_vtx =  std::find_if(
                            _ExPo,
                            vertices_cloud.begin(),
                            vertices_cloud.end(),
                            [a_coords=a->coords,b_coords=b->coords,c_coords=c->coords](const auto &d) { return in_circle(a_coords, b_coords, c_coords, d->coords) > 0.; }
                        );
                        auto is_delauney_triangle {vertices_cloud.end() == it_inside_vtx};

                        if(is_delauney_triangle)
                        {
                            // create potentially new face loop
                            auto he1 = make_shared_h_edge(a);
                            auto he2 = make_shared_h_edge(b);
                            auto he3 = make_shared_h_edge(c);
                            he1->previous = he3;
                            he2->previous = he1;
                            he3->previous = he2;

                            bool cross_tri{false}, is_perm{false};
                            for(const auto &h_face : faces_lst)
                            {
                                auto lst_vtx = {a, b, c};
                                // check if loop is a permutation of previously defined face
                                is_perm = std::is_permutation(lst_vtx.begin(), lst_vtx.end(), getFaceVertices(h_face).begin());
                                if(is_perm)
                                {
                                    break;
                                }
                                // check if edges are intersection face
                                cross_tri = 
                                    areFacesEdgesIntersect(*he1 , *h_face) ||
                                    areFacesEdgesIntersect(*he2 , *h_face) ||
                                    areFacesEdgesIntersect(*he3 , *h_face);
                                if(cross_tri)
                                {
                                    break;
                                }
                            }
                            if(!cross_tri && !is_perm)
                            {
                                auto lst1 = {he1, he2, he3};
                                faces_lst.push_back( make_shared_h_face<T,2>(lst1) );
                            }
                        }
                    
                    }
                }
            }
        }
        return faces_lst;
    }
}