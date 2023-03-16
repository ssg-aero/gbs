#pragma once

#include <array>
#include <vector>
#include <algorithm>
#include <span>

#include <gbs/surfaces>
#include <gbs/bscanalysis.h>

#include "halfEdgeMeshData.h"
#include "halfEdgeMeshGetters.h"
#include "halfEdgeMeshEditors.h"
#include "baseIntersection.h"
#include "halfEdgeMeshGeomTests.h"
#include "baseGeom.h"

namespace gbs
{

    template < typename T, size_t dim>
    auto meshSurfaceBoundary(const std::shared_ptr<Surface<T,dim>> &srf, size_t nu = 5, size_t nv = 5, T deviation = 0.01)
    {
        auto [u1_, u2_, v1_, v2_] = srf->bounds();


        auto isoV1 = CurveOnSurface<T,dim>(
            std::make_shared<BSCurve<T,2>>(build_segment<T,2>({u1_,v1_},{u2_,v1_})),
            srf
        );
        auto isoV2 = CurveOnSurface<T,dim>(
            std::make_shared<BSCurve<T,2>>(build_segment<T,2>({u1_,v2_},{u2_,v2_})),
            srf
        );
        auto isoU1 = CurveOnSurface<T,dim>(
            std::make_shared<BSCurve<T,2>>(build_segment<T,2>({u1_,v1_},{u1_,v2_})),
            srf
        );
        auto isoU2 = CurveOnSurface<T,dim>(
            std::make_shared<BSCurve<T,2>>(build_segment<T,2>({u2_,v1_},{u2_,v2_})),
            srf
        );

        auto u1 = deviation_based_params(isoV1, nu, deviation);
        auto u2 = deviation_based_params(isoV2, nu, deviation);
        std::reverse(u2.begin(), u2.end());
        auto v1 = deviation_based_params(isoU1, nv, deviation);
        std::reverse(v1.begin(), v1.end());
        auto v2 = deviation_based_params(isoU2, nv, deviation);

        auto n = u1.size() + u2.size() + v1.size() + v2.size();
        std::vector<std::array<T, 2>> coords(n);

        auto n_ = u1.size();
        auto coords_begin = coords.begin();
        std::transform(
            u1.begin(), u1.end(),
            coords_begin,
            [v1_](auto u)
            { return std::array<T, 2>{u, v1_}; });

        std::advance(coords_begin, n_);
        n_ = v2.size();
        std::transform(
            std::next(v2.begin()), v2.end(), 
            coords_begin, 
            [u2_](auto v)
            { return std::array<T, 2>{u2_, v}; });

        std::advance(coords_begin, n_);
        n_ = u2.size();
        std::transform(
            u2.begin(), u2.end(),
            coords_begin,
            [v2_](auto u)
            { return std::array<T, 2>{u, v2_}; });
        
        std::advance(coords_begin, n_);
        n_ = v1.size();
        std::transform(
            std::next(v1.begin()), v1.end(), 
            coords_begin, 
            [u1_](auto v)
            { return std::array<T, 2>{u1_, v}; });

        return coords;
    }

    template<typename T, typename _Container>
    void boyerWatson(_Container &h_f_lst, const std::array<T,2> &xy, T tol = 1e-10)
    {
        auto begin = h_f_lst.begin();
        auto it{begin};
        auto end = h_f_lst.end();

        _Container h_f_lst_deleted;
        // Find triangle violation Delaunay condition
        while(it!=end)
        {
            it = std::find_if(
                it, end,
                [xy, tol](const auto &h_f)
                { 
                    return in_circle(xy, h_f) > tol; 
                });
            if(it!=end)
            {
                h_f_lst_deleted.push_back(*it);
                it = h_f_lst.erase( it );
            }
        }
        if(h_f_lst_deleted.size()==0) return; // if point is outside or confused with an existing point
        assert(are_face_ccw(h_f_lst));
        // Get cavity boundary
        auto h_e_lst = getOrientedFacesBoundary(h_f_lst_deleted);
        assert(are_edges_2d_ccw(h_e_lst));
        // fill cavity
        auto h_f_lst_new = add_vertex(h_e_lst, make_shared_h_vertex(xy));
        assert(are_face_ccw(h_f_lst_new));
        // append new faces
        h_f_lst.insert(h_f_lst.end(), h_f_lst_new.begin(), h_f_lst_new.end() );
    }


    // template < typename T> 
    // auto getCoordsMinMax(const std::list< std::shared_ptr<HalfEdgeFace<T, dim>> > &hf_lst)
    // {

    // }

    ///////////////////////


    template < typename T>
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

    template < typename T >
    auto delaunay2DBoyerWatson(const auto &coords_boundary,const auto &coords_inner, T tol)
    {
        auto faces_lst = getEncompassingMesh(coords_boundary);
        auto vertices = getVerticesVectorFromFaces<T,2>(faces_lst);
        // insert points
        for(const auto &xy : coords_boundary)
        {
            boyerWatson<T>(faces_lst, xy, tol);
        }

        // remove external mesh, i.ei faces attached to initial vertices
        for(const auto &vtx : vertices)
        {
            remove_faces(faces_lst, vtx);
        }

        for(const auto &xy : coords_inner)
        {
            boyerWatson<T>(faces_lst, xy, tol);
        }

        return faces_lst;
    }

    template < typename T>
    auto delaunay2DBoyerWatson(const std::vector< std::array<T,2> > &coords, T tol = 1e-10)
    {
        return delaunay2DBoyerWatson<T,std::vector< std::array<T,2> >>(coords, tol);
    }

    template < typename T, size_t dim, typename _Func >
    auto delaunay2DBoyerWatsonSurfaceMesh(const std::shared_ptr<Surface<T,dim>> &srf, T crit_max, size_t max_inner_points = 500, size_t nu = 5, size_t nv = 5, T deviation = 0.01, T tol = 1e-10)
    {
        _Func dist_mesh_srf{*srf};

        auto coords = meshSurfaceBoundary(srf, nu ,nv, deviation);
        auto faces_lst = getEncompassingMesh(coords);
        auto vertices = getVerticesVectorFromFaces<T,2>(faces_lst);
        // insert boundary points
        for(const auto &xy : coords)
        {
            boyerWatson<T>(faces_lst, xy, tol);
        }
        // remove external mesh, i.ei faces attached to initial vertices
        // This work because surface coordiante are convex ( rectangle )
        for(const auto &vtx : vertices)
        {
            remove_faces(faces_lst, vtx);
        }
        // refine mesh to match max criteria
        T crit = 1.;
        for(size_t i{}; i < max_inner_points && crit > crit_max; i++)
        {
            auto it = std::max_element(
                std::execution::par,
                faces_lst.begin(), faces_lst.end(),
                [&dist_mesh_srf, &vertices](const auto &hf1, const auto &hf2)
                {
                    return dist_mesh_srf(hf1).first < dist_mesh_srf(hf2).first;
                }
            );

            auto d_G_UV = dist_mesh_srf(*it);
            crit = d_G_UV.first;
            // std::cout << i << " Max dist: " << crit << " Position: " << std::distance(it, faces_lst.begin())<< " "  << std::endl;

            boyerWatson(faces_lst, d_G_UV.second, tol);
        }

        return faces_lst;
    }

// TODO remove
    template < typename T , auto _ExPo = std::execution::seq>
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
                                is_perm = std::is_permutation(lst_vtx.begin(), lst_vtx.end(), gbs::getFaceVertices(h_face).begin());
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