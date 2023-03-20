#include <gtest/gtest.h>
#include <memory>
#include <numeric>
#include <exception>
#include <gbs/bscbuild.h>
#include <gbs-render/vtkGbsRender.h>

using Real = float;
const auto tol = std::numeric_limits<Real>::min();
struct Real3
{
    Real x;
    Real y;
    Real z;
};

struct Real2
{
    Real u;
    Real v;
};

struct HalfEdge;
using shr_HalfEdge = std::shared_ptr<HalfEdge>;

struct HalfEdgeVertex
{
    Real3 coords;
    Real2 coordUV;
    shr_HalfEdge p_Hed;
};
using shr_HalfEdgeVertex = std::shared_ptr<HalfEdgeVertex>;

struct HalfEdgeFace
{
    shr_HalfEdge p_Hed;
};
using shr_HalfEdgeFace = std::shared_ptr<HalfEdgeFace>;

struct HalfEdge
{
    shr_HalfEdgeVertex vertex;
    shr_HalfEdgeFace   face;
    shr_HalfEdge       next;
    shr_HalfEdge       previous;
    shr_HalfEdge       opposite;
};


auto connect_edges(shr_HalfEdge &ed1, shr_HalfEdge &ed2)
{
    if( ed1->next || ed2->previous)
    {
        throw std::invalid_argument("HalfEdges already connected.");
    }  
    ed1->next = ed2;
    ed2->previous = ed1;
}

auto make_half_edge(const shr_HalfEdgeVertex &vtx)
{
    auto hEd = std::make_shared<HalfEdge>(
        HalfEdge{
            .vertex = vtx
        }
    );

    return hEd;
}

auto make_half_edge(const Real3 &pt)  -> shr_HalfEdge
{
    auto vtx = std::make_shared<HalfEdgeVertex>(
        HalfEdgeVertex{
            .coords = pt
        }
    );

    return make_half_edge(vtx);
}

auto make_half_edge(const shr_HalfEdgeVertex &vtx, shr_HalfEdge &previous)
{
    auto hEd = std::make_shared<HalfEdge>(
        HalfEdge{
            .vertex = vtx,
        }
    );

    connect_edges(previous, hEd);

    return hEd;
}

auto make_half_edge(const Real3 &pt, shr_HalfEdge &previous)  ->shr_HalfEdge
{
    auto vtx = std::make_shared<HalfEdgeVertex>(
        HalfEdgeVertex{
            .coords = pt
        }
    );
    return make_half_edge(vtx, previous);
}

auto make_face( shr_HalfEdge &first)
{
    auto face = std::make_shared<HalfEdgeFace>(HalfEdgeFace{first});
    auto current = first->next;
    first->face = face;
    while (current != first)
    {
        if( current )
        {
            current->face = face;
        }
        else
        {
            throw std::invalid_argument("Not a closed loop.");
        }
        current = current->next;
    }

    return face;

}

auto get_tail(const shr_HalfEdge &first)
{
    auto tail = first;
    while (tail->next && tail->next != first )
    {
        tail = tail->next;
    }
    return tail;
}

auto get_head(const shr_HalfEdge &first)
{
    auto head = first;
    while (head->previous && head->previous != first )
    {
        head = head->previous;
    }
    return head;
}

auto get_head_tail(const shr_HalfEdge &first)
{
    return std::make_pair(get_head(first), get_tail(first));
}

auto close_to_face(shr_HalfEdge &first)
{
    auto [ head, tail] = get_head_tail(first);
    if( head == first)
    {
        throw std::invalid_argument("First edge belongs to a face.");
    }
    if(!head->opposite)
    {
        throw std::invalid_argument("Vertex is missing.");
    }
    auto hEd = make_half_edge(head->opposite->vertex, tail);
    connect_edges(hEd, head);
    
    return std::make_pair( make_face(first), hEd );
}
/**
 * @brief 
 * 
 * @param pt 
 * @param first 
 * @return auto 
 */
auto close_to_face(const Real3 &pt, shr_HalfEdge &first)
{
    auto [ head, tail] = get_head_tail(first);
    if( head == tail)
    {
        throw std::invalid_argument("First edge belongs to a face.");
    }

    auto hEd = make_half_edge(pt, tail);
    connect_edges(hEd, first);

    return std::make_pair( make_face(first), hEd );
}

auto add_opposite(shr_HalfEdge &hEd)
{
    if( hEd->opposite )
    {
        throw std::invalid_argument("Edge already have an opposite.");
    }
    if(!hEd->previous)
    {
        throw std::invalid_argument("Previous half edge is required.");
    }
    
    auto v_prev = hEd->previous->vertex;
    auto hEdOpp = std::make_shared<HalfEdge>( HalfEdge{
        .vertex   = v_prev,
        .opposite = hEd,
    } );
    hEd->opposite = hEdOpp;

    return hEdOpp;
}

auto is_closed(const shr_HalfEdge & first)
{
    auto [ head, tail] = get_head_tail(first);
    return head->next == tail;
}

TEST(gbs_mesh, HalfEdgeMesh)
{
    auto hEd1 = make_half_edge( {1,0,0} );

    ASSERT_NEAR(hEd1->vertex->coords.x,1.,tol);
    ASSERT_NEAR(hEd1->vertex->coords.y,0.,tol);
    ASSERT_NEAR(hEd1->vertex->coords.z,0.,tol);

    ASSERT_TRUE(hEd1 == get_head(hEd1));
    ASSERT_TRUE(hEd1 == get_tail(hEd1));

    auto hEd2 = make_half_edge( {0.5,1,0}, hEd1 );

    ASSERT_TRUE(hEd1->next == hEd2);
    ASSERT_TRUE(hEd1 == hEd2->previous);

    ASSERT_NEAR(hEd2->vertex->coords.x,0.5,tol);
    ASSERT_NEAR(hEd2->vertex->coords.y,1.,tol);
    ASSERT_NEAR(hEd2->vertex->coords.z,0.,tol);
    ASSERT_EQ(hEd1->next, hEd2);
    ASSERT_EQ(hEd2->previous, hEd1);

    ASSERT_TRUE(hEd1 == get_head(hEd2));
    ASSERT_TRUE(hEd2 == get_tail(hEd1));

    auto [face1, hEd3] = close_to_face({0,0,0}, hEd1);

    ASSERT_TRUE(hEd2->next == hEd3);
    ASSERT_TRUE(hEd2 == hEd3->previous);

    ASSERT_TRUE(is_closed(hEd1));

    ASSERT_FALSE( face1 == nullptr);

    auto hEd4 = add_opposite(hEd2);

    ASSERT_FALSE(  hEd4->opposite == hEd2->opposite );
    ASSERT_FALSE(  hEd2->opposite == hEd4->opposite );

    auto hEd5 =  make_half_edge( {1.5,1,0}, hEd4 );

    auto [face2, hEd6] = close_to_face(hEd5);

    ASSERT_TRUE(is_closed(hEd4));
}

/**
 * @brief Build HalfEdge loop from points
 * @tparam _It 
 * @param points_start 
 * @param points_end 
 * @return std::list<shr_HalfEdge> 
 */
template < typename _It >
auto make_loop( const _It &points_start, const _It &points_end ) -> std::list<shr_HalfEdge>
{
    // init loop
    std::list<shr_HalfEdge> loop;
    auto pt = points_start;
    loop.push_back(
        make_half_edge( {pt->at(0), pt->at(1), pt->at(2)} )
    );
    // add remaining points
    while (pt != points_end)
    {
        loop.push_back(
            make_half_edge( {pt->at(0), pt->at(1), pt->at(2)}, loop.back() )
        );
        pt = std::next(pt);
    }

    return loop;
    
}

TEST(gbs_mesh, AdvancingFront)
{
    auto bound_ext = gbs::build_ellipse<Real,3>(2,1);
    auto bound_int = gbs::build_circle<Real,3>(0.5,{0.5,0,0});

    auto bound_ext_pts = gbs::discretize<Real,3>(bound_ext,10,0.05);
    auto bound_int_pts = gbs::discretize<Real,3>(bound_int,10,0.05);

    std::list<shr_HalfEdge> edges;

    auto loop_ext = make_loop( bound_ext_pts.rbegin(), bound_ext_pts.rend() );
    auto loop_int = make_loop( bound_int_pts.begin(), bound_int_pts.end() );

    std::list< std::list<shr_HalfEdge> > mesh_front{loop_ext, loop_int};


    gbs::plot(
          bound_ext
        , bound_int
        , bound_ext_pts
        , bound_int_pts
    );
}