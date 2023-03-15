#pragma once
#include "halfEdgeMeshData.h"
#include "baseGeom.h"

namespace gbs
{

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
    /**
     * @brief Test if the face is counter clock wise. 
     * Returns > 0 if d inside circle formed by the corners of the triangular face return 0. if on
     * 
     * @tparam T 
     * @param hf 
     * @return T 
     */
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

    template <typename T>
    bool areFacesEdgesIntersect(const gbs::HalfEdge<T,2> &h_e,const gbs::HalfEdgeFace<T, 2> &h_f )
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

} // namespace gbs