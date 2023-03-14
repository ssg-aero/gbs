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
#include "baseGeom.h"

namespace gbs
{

    // TODO: Remove
    template <typename T, size_t dim>
    auto add_vertex(const std::shared_ptr<HalfEdgeFace<T, dim>> &h_f, const std::shared_ptr<HalfEdgeVertex<T, dim>> &h_v)
    {
        return add_vertex(getFaceEdges(h_f), h_v);
    }

    /**
     * @brief returns > 0 if d inside circle formed by the corners of the triangular face return 0. if on
     * 
     * @tparam T 
     * @param pt 
     * @param t 
     * @return T 
     */
    template <typename T>
    T in_circle(const std::array<T,2> &pt, const HalfEdgeFace<T, 2> &t)
    {
        assert(getFaceVertices( t ).size() == 3 );
        auto he= t.edge;
        const auto &a = he->vertex->coords;
        he = he->next;
        const auto &b = he->vertex->coords;
        he = he->next;
        const auto &c = he->vertex->coords;
        return in_circle(a, b, c, pt);
    }
    /**
     * @brief  returns > 0 if d inside circle formed by the corners of the triangular face return 0. if on
     * 
     * @tparam T 
     * @param pt 
     * @param t 
     * @return T 
     */
    template <typename T>
    T in_circle(const std::array<T,2> &pt, const std::shared_ptr<HalfEdgeFace<T, 2>> &t)
    {
        assert(t);
        return in_circle(pt,*t);
    }

    template <typename T>
    T is_ccw(const std::shared_ptr<HalfEdgeFace<T, 2>> &hf)
    {
        auto vertices = getFaceVertices(hf);
        assert(vertices.size()==3);
        return are_ccw(
            (*std::next(vertices.begin(),0))->coords,
            (*std::next(vertices.begin(),1))->coords,
            (*std::next(vertices.begin(),2))->coords
        );
    }

    template <typename T, typename _Container>
    bool is_inside(const std::array<T, 2> &xy, const _Container &h_e_lst)
    {
        size_t count{};
        for (const auto &he : h_e_lst)
        {
            const auto &a = he->previous->vertex->coords;
            const auto &b = he->vertex->coords;   
            if(on_seg(a,b,xy) )
            {
                return true;
            }
            if ( seg_H_strict_end_intersection( a, b, xy))
            {
                count++;
            }
        }
        // When count is odd
        return count % 2 == 1;
    }

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

    template <typename T>
    T is_locally_delaunay(const std::shared_ptr<HalfEdge<T, 2>> &he)
    {
        assert(he);
        if(!he->opposite) return -1;

        auto t1 = he->face;
        auto t2 = he->opposite->face;
        assert(t1);
        assert(t2);
        auto vertices1 = getFaceVertices( t1 );
        auto vertices2 = getFaceVertices( t2 );
        assert(vertices1.size() == 3 );
        assert(vertices2.size() == 3 );

        auto p1 = he->next->vertex->coords;
        auto p2 = he->opposite->next->vertex->coords;
        return std::min( in_circle(p1, t2), in_circle(p2, t1) );
    }

    bool are_face_ccw(const auto &faces_lst)
    {
        for(const auto &hf: faces_lst)
        {
            if(!is_ccw(hf))
            {
                return false;
            }
        }
        return true;
    }

    bool are_edges_2d_ccw(const auto &edges_lst)
    {
        auto sum = std::reduce(
            edges_lst.begin(), edges_lst.end(),
            0.,
            [](auto t , const auto &he){
                auto [x2, y2] = he->vertex->coords;
                auto [x1, y1] = he->previous->vertex->coords;
                return t + (x2-x1)*(y2+y1);
            }
        );
        return sum < 0.;
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
                    // std::cout << in_circle(xy, h_f) <<std::endl;
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
        // // split edge if point is on edge
        // auto it = std::find_if(
        //     h_e_lst.begin(),h_e_lst.end(),
        //     [tol](const auto &he)
        //     {
        //         const auto O = he->previous->vertex->coords;
        //         norm( (he->vertex->coords- O)^(xy- O)  ) < tol
        //     }
        // );
        // if(it != h_e_lst.end())
        // {

        // }
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

    template <typename T>
    inline T sum_area(T Area, const std::shared_ptr<HalfEdgeFace<T, 2>> &h_face)
    {
        auto coords = getFaceCoords(h_face);
        assert(coords.size()==3);
        return Area + tri_area(
                    *std::next(coords.begin(), 0),
                    *std::next(coords.begin(), 1),
                    *std::next(coords.begin(), 2)
                );
    }

    template <typename T>
    inline T sum_area(const std::shared_ptr<HalfEdgeFace<T, 2>> &h_face, T Area)
    {
        return sum_area(Area, h_face);
    }

    template <typename T>
    inline T sum_area( T Area1, T Area2)
    {
        return Area1+Area2;
    }

    template <typename T>
    inline T sum_area(const std::shared_ptr<HalfEdgeFace<T, 2>> &h_face1, const std::shared_ptr<HalfEdgeFace<T, 2>> &h_face2)
    {
        return sum_area(T{},h_face1)+sum_area(T{},h_face2);
    }

    auto getTriangle2dMeshArea(const auto &faces_lst)
    {
        return std::reduce(
            faces_lst.begin(),faces_lst.end(),
            0.,
            [](auto a1, auto a2)
            {
                return sum_area(a1, a2);
            }
        );
    }

    auto getTriangle2dMeshAreaPar(const auto &faces_lst)
    {
        return std::reduce(
            std::execution::par,
            faces_lst.begin(),faces_lst.end(),
            0.,
            [](auto a1, auto a2)
            {
                return sum_area(a1, a2);
            }
        );
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
    struct HalfEdgeMesh
    {
        std::vector<std::shared_ptr<HalfEdge<T, dim>>> edges;
        std::vector<std::shared_ptr<HalfEdgeVertex<T, dim>>> vertices;
        std::vector<std::shared_ptr<HalfEdgeFace<T, dim>>> faces;

        // void addEdge( const std::array<T,dim> &coord )
        // {
        //     edges.push_back( std::make)
        // }

        void addFreeVertices(std::vector<std::array<T, dim>> &free_vertices)
        {
            auto n = vertices.size();
            vertices.insert(vertices.end(),free_vertices.size(), nullptr);

            std::transform(
                free_vertices.begin(), free_vertices.end(),
                std::next(vertices.begin(),n),
                [](const auto &coords)
                {
                    return std::make_shared<HalfEdgeVertex<T, dim>>(
                        HalfEdgeVertex<T, dim>{coords, nullptr}
                    );
                }
            );
        }
        
        /**
         * @brief Must be cw
         * 
         * @param i 
         * @param j 
         * @param k 
         * @return true 
         * @return false 
         */
        bool makeFace( size_t i, size_t j, size_t k)
        {
            // auto vtx1 = vertices[i];
            // auto vtx2 = vertices[j];
            // auto vtx3 = vertices[k];

            // auto p_face = std::make_shared<HalfEdgeFace<T, dim>>()

            // HalfEdge<T, dim> ed1{
            //     .vertex = vtx1,
            //     .face   = p_face,
            // };
            // HalfEdge<T, dim> ed2{
            //     .vertex = vtx2,
            //     .face   = p_face,
            // };
            // HalfEdge<T, dim> ed3{
            //     .vertex = vtx3,
            //     .face   = p_face,
            // };

            // // auto p_ed1 = std::make_shared< HalfEdge<T, dim> > (  HalfEdge<T, dim>{, .} 


            return true;
        }
    };

    template <typename T, size_t dim>
    auto extract_vertices_map_from_faces(const auto &faces_lst ) -> std::map< std::shared_ptr< HalfEdgeVertex<T,dim> >, size_t >
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
    auto extract_vertices_vector_from_faces(const auto &faces_lst ) -> std::vector< std::shared_ptr< HalfEdgeVertex<T,dim> > >
    {
        auto vertices_map = extract_vertices_map_from_faces<T,dim>(faces_lst);

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
    auto extract_edges_map(const auto &faces_lst ) -> std::map< std::shared_ptr< HalfEdgeVertex<T,dim> >, size_t >
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

    template <typename T, size_t dim>
    auto extract_edges_boundary(const auto &faces_lst )
    {
        auto edges_map = extract_edges_map(faces_lst);
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
        auto vertices_map = extract_vertices_map_from_faces<T,2>(faces_lst);
        // insert points
        for(const auto &xy : coords)
        {
            boyer_watson<T>(faces_lst, xy, tol);
        }
        // remove external mesh, i.ei faces attached to initial vertices
        for(const auto &vtx : vertices_map)
        {
            auto it_hf = faces_lst.begin();
            
            while(it_hf != faces_lst.end())
            {
                auto h_e_lst = getFaceEdges(*it_hf);
                auto it = std::find_if( h_e_lst.begin(), h_e_lst.end(),
                        [&vtx](const auto &he){return vtx.first == he->vertex;}
                );
                if( h_e_lst.end() != it )
                {
                    for(auto & h_e : h_e_lst) // no more opposite
                    {
                        if(h_e->opposite)
                        {
                            h_e->opposite->opposite = nullptr;
                        }
                    }
                    it_hf = faces_lst.erase(it_hf);
                }
                else{
                    std::advance(it_hf,1);
                }
            }
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