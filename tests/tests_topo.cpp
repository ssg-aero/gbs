#include <gtest/gtest.h>
#include <topology/vertex.h>
#include <topology/edge.h>
#include <topology/wire.h>
#include <gbs/maths.h>
#include <gbs/bscanalysis.h>

#include <gbs-render/vtkcurvesrender.h>

TEST(tests_topo, vtx)
{
    using namespace gbs;
    using T = double;
    const size_t d{3};

    T precision{1e-5};
    T approximation{1e-4};

    std::array<T, d> pt2{1., 2., 3.};

    Vertex<T, d> vtx1{{}};
    vtx1.setPrecisionTolerance(precision);
    vtx1.setApproximationTolerance(approximation);
    Vertex<T, d> vtx2{pt2, vtx1};

    ASSERT_DOUBLE_EQ(vtx2.precisionTolerance(), precision);
    ASSERT_DOUBLE_EQ(vtx2.approximationTolerance(), approximation);
    ASSERT_NEAR(distance(pt2, vtx2.point()), 0., vtx2.precisionTolerance());
}

TEST(tests_topo, edge_wire)
{

    using namespace gbs;
    using T = double;
    const size_t d{3};

    std::array<T, d> pt1{};
    std::array<T, d> pt2{1., 0., 0.};
    std::array<T, d> pt3{1., 1., 0.};
    std::array<T, d> pt4{0., 1., 0.};

    Edge<T, d> ed1{pt1, pt2};

    auto vtx3 = std::make_shared<gbs::Vertex<T, d>>(pt3);
    Edge<T, d> ed2{ed1.vertex2(), vtx3};

    ASSERT_TRUE(ed1.vertex2()==ed2.vertex1());

    Wire<T,d> w1{ed1};

    ASSERT_TRUE( w1.addEdge(ed2) );

    Edge<T, d> ed3{pt3, pt4};
    Edge<T, d> ed4{pt4, pt1};

    ASSERT_TRUE( w1.addEdge(ed4) );
    ASSERT_TRUE( w1.addEdge(ed3) );

    ASSERT_TRUE( w1.isClosed() );


}




TEST(tests_topo, mesh_wire)
{
    using namespace gbs;
    using T = double;
    const size_t d{3};

    std::array<T, d> pt1{};
    std::array<T, d> pt2{1., 0., 0.};
    std::array<T, d> pt3{1., 1., 0.};
    std::array<T, d> pt4{0., 1., 0.};

    Wire<T,d> w1{{pt1, pt2}};
    w1.addEdge({pt2, pt3});
    w1.addEdge({pt3, pt4});
    w1.addEdge({pt4, pt1});

    T dm{0.1};

    std::vector< std::array<T,d> > bound_vertices;

    auto mesh_edge = [dm](const auto & p_ed)
    {
        auto [u1, u2] = p_ed->bounds();
        const auto &p_crv = p_ed->curve();
        auto l = length(*p_crv, u1, u2);
        size_t n = std::round(l / dm) + 1;
        auto u_lst = uniform_distrib_params(*p_crv, u1, u2, n);
        return make_points( *p_crv, u_lst);
    };

    std::for_each(
        w1.begin(), w1.end(),
        [&bound_vertices, &mesh_edge] (const auto & p_ed)
        {
            auto points = mesh_edge(p_ed);
            bound_vertices.push_back( p_ed->vertex1()->point() );
            bound_vertices.insert( 
                bound_vertices.end(), 
                std::next(points.begin()),
                std::next(points.end(),-1)
            );
        }
    );

    // HalfEdgeMesh<T,d> msh;
    // // Store boundary
    // std::for_each(
    //     bound_vertices.begin(), bound_vertices.end(),
    //     [&msh,&d]( const auto & coords)
    //     {
    //         msh.edges.push_back(
    //             std::make_shared<HalfEdge<T,d>>(coords)
    //         );
    //     }
    // );

    plot(bound_vertices);

}