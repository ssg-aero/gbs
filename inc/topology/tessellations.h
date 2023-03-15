#pragma once

#include <array>
#include <vector>
#include <algorithm>
#include <span>

#include <gbs/surfaces>

#include "halfEdgeMeshData.h"
#include "halfEdgeMeshGetters.h"
#include "halfEdgeMeshEditors.h"
#include "baseIntersection.h"
#include "halfEdgeMeshGeomTests.h"
#include "baseGeom.h"

namespace gbs
{

    template <typename T, typename _Container1, typename _Container2>
    auto take_external_faces(_Container1 &faces_lst, const _Container2 &boundary)
    {
        std::list<std::shared_ptr<HalfEdgeFace<T,2>>> external_faces;
        auto it = faces_lst.begin();
        while(it!=faces_lst.end())
        {
            it = std::find_if(
                faces_lst.begin(), faces_lst.end(),
                [&boundary](const auto &hf)
                {
                    auto coords = getFaceCoords(hf);
                    std::array<T,2> G{};
                    for(const auto &xy : coords)
                    {
                        G = G + xy;
                    }
                    G = G / T(coords.size());
                    return !is_inside(G,boundary);

                    // for(const auto &xy : coords)
                    // {
                    //     if(!is_inside(xy,boundary))
                    //     {
                    //         return true;
                    //     }
                    // }                    
                    // return false;
                }
            );
            if(it!=faces_lst.end())
            {
                for(auto & h_e : getFaceEdges(*it))
                {
                    if(h_e->opposite)
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
    auto take_internal_faces(_Container1 &faces_lst, const _Container2 &boundary)
    {
        std::list<std::shared_ptr<HalfEdgeFace<T,2>>> external_faces;
        auto it = faces_lst.begin();
        while(it!=faces_lst.end())
        {
            it = std::find_if(
                faces_lst.begin(), faces_lst.end(),
                [&boundary](const auto &hf)
                {
                    auto coords = getFaceCoords(hf);
                    std::array<T,2> G{};
                    for(const auto &xy : coords)
                    {
                        G = G + xy;
                    }
                    G = G / T(coords.size());
                    return is_inside(G,boundary);
                }
            );
            if(it!=faces_lst.end())
            {
                for(auto & h_e : getFaceEdges(*it))
                {
                    if(h_e->opposite)
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

    void reverse_boundary(auto &boundary)
    {
        std::reverse(boundary.begin(),boundary.end());
        std::transform(
            boundary.begin(),boundary.end(),
            boundary.begin(),
            [](const auto &he){
                std::swap(he->previous, he->next);
                return he;
            }
        );        
    }

    template<typename _Container,typename _It>
    auto erase_face(_Container &h_f_lst, const _It &it)
    {

        auto face_edges = getFaceEdges(*it);
        for(auto & he : face_edges)
        {
            if(he->opposite)
            {
                he->opposite->opposite = nullptr;
            }
        }
        return h_f_lst.erase( it );
    }

    template<typename T, typename _Container>
    void boyer_watson(_Container &h_f_lst, const std::array<T,2> &xy, T tol = 1e-10)
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

 

    template <typename T>
    auto getEncompassingMesh(const std::vector< std::array<T, 2> > &X_lst, T pc_offset = 10)
    {
        auto [Xmin, Xmax] = getCoordsMinMax(X_lst);

        auto pc = pc_offset / 100;
        Xmax = Xmax + pc*(Xmax-Xmin);
        Xmin = Xmin - pc*(Xmax-Xmin);
        auto he1 = make_shared_h_edge<T,2>({Xmax[0],Xmin[1]});
        auto he2 = make_shared_h_edge<T,2>({Xmin[0],Xmax[1]});
        auto he3 = make_shared_h_edge<T,2>(Xmin);
        auto lst1 = {he1, he2, he3};
        auto hf1 = make_shared_h_face<T,2>(lst1);
        auto hf2 = add_face(hf1,he2,Xmax);

        return std::list< std::shared_ptr<HalfEdgeFace<T, 2>> >{hf1,hf2};
    }

    // template < typename T> 
    // auto getCoordsMinMax(const std::list< std::shared_ptr<HalfEdgeFace<T, dim>> > &hf_lst)
    // {

    // }

//https://www.graphics.rwth-aachen.de/software/openmesh/intro/
/**
 * @brief Get 
 * 
 * @tparam T 
 * @tparam dim 
 * @param vertex 
 * @return auto 
 */
    template <typename T, size_t dim>
    auto getLoop( const std::shared_ptr<HalfEdgeVertex<T, dim>> &vertex)
    {
        std::list< std::shared_ptr< HalfEdge<T, dim> > > vtx_lst;
        auto edge = vertex->edge;

        while (!edge->face)
        {
            /* code */
        }
        

        // while (edge-> && edge != vertex->edge)
        // {
        //     vtx_lst(edge);
            
        // }
        
    }


//     }

//     auto getIncomings()
//     {

//     }

    // template <typename T, size_t dim>
    // auto getEdgeLoop( const std::shared_ptr<HalfEdge<T, dim>> &ed_start) -> std::list< std::shared_ptr<HalfEdge<T, dim>> >
    // {

    //     std::list< std::shared_ptr< HalfEdge<T, dim> > > edges_lst;

    //     auto ed{ed_start};

    //     while(ed)
    //     {
    //         edges_lst.push_back( ed );
    //         ed = ed->next;
    //         if(ed_start == ed)
    //         {
    //             break;
    //         }
    //     }

    //     if(edges_lst.back()->next != edges_lst.front()->previous) // only closed loop
    //     {
    //         edges_lst.clear();
    //     }

    //     return edges_lst;
    // }


    // template <typename T, size_t dim>
    // auto getVertexMainLoop( const std::shared_ptr<HalfEdgeVertex<T, dim>> &vertex) -> std::list< std::shared_ptr<HalfEdge<T, dim>> >
    // {
    //     auto edges_lst = getEdgeLoop(vertex->edge);
    //     // permutation to start with the half edge coming from the vertex
    //     auto front = edges_lst.front();
    //     edges_lst.push_back(  front ); 
    //     edges_lst.pop_front();
    //     return edges_lst;
    // }

    // template <typename T, size_t dim>
    // auto getVertexLoops( const std::shared_ptr<HalfEdgeVertex<T, dim>> &vertex) -> std::list< std::list< std::shared_ptr<HalfEdge<T, dim>> > >
    // {
    //     std::list< std::list< std::shared_ptr< HalfEdge<T, dim> > > > edges_lst;

    //     auto main_loop = getVertexMainLoop(vertex);


    //     auto next_opp = main_loop.front()->opposite;
        
    //     while (next_opp &&  next_opp != main_loop.back() )
    //     {
    //         auto next_loop = getEdgeLoop( next_opp->next );

    //     }
        
        
    //     auto prev_opp = main_loop.back()->opposite;



    // }

    // template <typename T, size_t dim>
    // auto getConnectedVertices( const std::shared_ptr<HalfEdgeVertex<T, dim>> &vertex) -> std::list< std::shared_ptr< HalfEdgeVertex<T, dim> > >
    // {
    //     std::list< std::shared_ptr< HalfEdgeVertex<T, dim> > > vertices_lst;

    //     auto main_loop = getVertexMainLoop(vertex);
    //     if(main_loop.size())
    //     {
    //         auto it = main_loop.begin();
    //         vertices_lst.push_back(*(it))
    //         std::advance(it,main_loop.size()-1);
    //         vertices_lst.push_back(*(it))
    //     }



    //     // auto ed = vertex->edge;
    //     // auto ed_prev = ed->opposite;
    //     // while (ed_prev)
    //     // {
    //     //     vertices_lst.push_front(ed_prev->vertex);
    //     //     if(ed_prev->previous)
    //     //     {
    //     //         ed_prev = ed_prev->previous->opposite;
    //     //     }
    //     //     else
    //     //     {
    //     //         break;
    //     //     }
    //     // }

    //     // auto ed_next = ed->next;
    //     // while (ed_next)
    //     // {
    //     //     vertices_lst.push_back(ed_next->vertex);
    //     //     if( ed_next->opposite)
    //     //     {
    //     //         ed_next = ed_next->opposite->next;
    //     //     }
    //     //     else
    //     //     {
    //     //         break;
    //     //     }
    //     // }
        
    //     return vertices_lst;
    // }


    template <typename T, size_t dim>
    auto extract_edges_boundary(const auto &faces_lst )
    {
        auto edges_map = getEdgesMap(faces_lst);
        auto bound_start = std::find_if(
            edges_map.begin(), edges_map.end(),
            [](const auto &hed_id){return !hed_id.fisrt->opposite;}
        );

        auto hed = *bound_start;
        do{

        }while (hed != *bound_start);
        
    }

    

    ///////////////////////


    template < typename T, typename _Container>
    auto base_delaunay2d_mesh(const _Container &coords, T tol = 1e-10)
    {
        auto faces_lst = getEncompassingMesh(coords);
        auto vertices = getVerticesVectorFromFaces<T,2>(faces_lst);
        // insert points
        for(const auto &xy : coords)
        {
            boyer_watson<T>(faces_lst, xy, tol);
        }

        // remove external mesh, i.ei faces attached to initial vertices
        for(const auto &vtx : vertices)
        {
            remove_faces(faces_lst, vtx);
        }

        return faces_lst;
    }

    template <typename T>
    bool is_face_edges_inter(const gbs::HalfEdge<T,2> &h_e,const gbs::HalfEdgeFace<T, 2> &h_f )
    {
        assert(h_e.previous);
        for(const auto &h_e_ : gbs::getFaceEdges(h_f) )
        {
            assert(h_e_->previous);
            if(seg_seg_strict_intersection(
                h_e.previous->vertex->coords,
                h_e.vertex->coords,
                h_e_->previous->vertex->coords,
                h_e_->vertex->coords
            ))
            {
                return true;
            }
        }
        return false;
    }

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
                                    is_face_edges_inter(*he1 , *h_face) ||
                                    is_face_edges_inter(*he2 , *h_face) ||
                                    is_face_edges_inter(*he3 , *h_face);
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