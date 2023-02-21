#include <gtest/gtest.h>
#include <gbs/maths.h>
#include <gbs-render/vtkcurvesrender.h>
#include <topology/tessellations.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkXMLPolyDataWriter.h>

using namespace gbs;
const bool plot_on = false;

template <typename T, size_t dim>
inline auto mesh_wire_uniform(const Wire<T,dim> &w, T dm)
{
    auto mesh_edge = [dm](const auto & p_ed)
    {
        auto [u1, u2] = p_ed->bounds();
        const auto &p_crv = p_ed->curve();
        auto l = length(*p_crv, u1, u2);
        size_t n = std::round(l / dm) + 1;
        auto u_lst = uniform_distrib_params(*p_crv, u1, u2, n);
        return make_points( *p_crv, u_lst);
    };

    std::vector< std::array<T,dim> > coords;

    std::for_each(
        w.begin(), w.end(),
        [&coords, &mesh_edge] (const auto & p_ed)
        {
            auto points = mesh_edge(p_ed);
            coords.push_back( p_ed->vertex1()->point() );
            coords.insert( 
                coords.end(), 
                std::next(points.begin()),
                std::next(points.end(),-1)
            );
        }
    );

    return coords;
}

template <typename T>
inline auto make_boundary2d_1(T dm = 0.1)
{
    std::array<T, 2> pt1{};
    std::array<T, 2> pt2{1., 0.};
    std::array<T, 2> pt3{1., 1.};
    std::array<T, 2> pt4{0., 1.};

    Wire2d<T> w{{pt1, pt2}};
    w.addEdge({pt2, pt3});
    w.addEdge({pt3, pt4});
    w.addEdge({pt4, pt1});

    return mesh_wire_uniform(w, dm);

}

template <typename T>
inline auto make_boundary2d_2(T dm = 0.1)
{
    auto ell = build_ellipse<T,2>(1.,0.5);
    Wire2d<T> w(Edge<T,2>{std::make_shared<BSCurveRational<T,2>>(ell)});

    return mesh_wire_uniform(w, dm);

}

template <typename T>
inline auto make_boundary2d_3(T dm = 0.1)
{
    point<T,2> P0{0,0};
    point<T,2> P1{1,0};
    point<T,2> P2{1,1};
    point<T,2> P3{2,1};
    point<T,2> P4{2,0};
    point<T,2> P5{3,0};
    point<T,2> P6{3,2};
    point<T,2> P7{2,2};
    point<T,2> P8{2,3};
    point<T,2> P9{1,3};
    point<T,2> P10{0,2};

    Wire2d<T> w{{P0, P1}};
    w.addEdge({P1, P2});
    w.addEdge({P2, P3});
    w.addEdge({P3, P4});
    w.addEdge({P4, P5});
    w.addEdge({P5, P6});
    w.addEdge({P6, P7});
    w.addEdge({P7, P8});
    w.addEdge({P8, P9});
    w.addEdge({P9, P10});
    w.addEdge({P10, P0});

    return mesh_wire_uniform(w, dm);

}

template <typename T>
inline auto add_random_points_grid(std::vector< std::array<T,2> > &coords, size_t n, T x1=0, T x2=1., T ratio=1.)
{
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<double> distrib(x1,1.0);
    for(size_t i{}; i < n; i++)
    {
        coords.push_back(std::array<T,2>{distrib(gen),distrib(gen)});
    }
}

template <typename T>
inline auto add_random_points_ellipse(std::vector< std::array<T,2> > &coords, size_t n, T R=1., T r=0.5)
{
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<double> distribr1(0.0,R);
    std::uniform_real_distribution<double> distribr2(0.0,r);
    std::uniform_real_distribution<double> distribth(0.0,2*std::numbers::pi_v<T>);
 
    for(size_t i{}; i < n; i++)
    {
            auto th = distribth(gen);
            coords.push_back( std::array<T,2>{distribr1(gen)*std::cos(th),distribr2(gen)*std::sin(th)} );
    }
}


TEST(halfEdgeMesh, getVertexMainLoop)
{

    using T = double;
    const size_t d = 2;

    auto he1 = make_shared_h_edge<T,d>({1.0,0.0});
    auto he2 = make_shared_h_edge<T,d>({0.0,1.0});
    auto he3 = make_shared_h_edge<T,d>({0.0,0.0});

    ASSERT_EQ(he1->vertex->edge, he1);
    ASSERT_EQ(he2->vertex->edge, he2);
    ASSERT_EQ(he3->vertex->edge, he3);

    auto lst1 = {he1, he2, he3};
    auto hf1 = make_shared_h_face<T,d>(lst1);

    ASSERT_EQ(he1->face, hf1);
    ASSERT_EQ(he2->face, hf1);
    ASSERT_EQ(he3->face, hf1);


    ASSERT_EQ(he1, getFaceEdge(hf1,he1->vertex));
    ASSERT_EQ(he2, getFaceEdge(hf1,he2->vertex));
    ASSERT_EQ(he3, getFaceEdge(hf1,he3->vertex));

    auto hf2 = add_face(hf1,he2,{1.0,1.0});

    ASSERT_EQ(he2,getCommonEdge(hf1, hf2));

    ASSERT_EQ(he2->opposite,std::get<1>(getCommonEdges(hf1, hf2)));

    {
        auto loop = getFaceEdges(hf2);
        ASSERT_TRUE(loop.size()==3);
        auto it = loop.begin();
        ASSERT_LE(distance((*(it++))->vertex->coords, std::array<T, d>{0.0,1.0}),1e-6);
        ASSERT_LE(distance((*(it++))->vertex->coords, std::array<T, d>{1.0,0.0}),1e-6);
        ASSERT_LE(distance((*(it++))->vertex->coords, std::array<T, d>{1.0,1.0}),1e-6);
    }
    flip(hf1,hf2);
    {
        auto loop = getFaceEdges(hf1);
        ASSERT_TRUE(loop.size()==3);
        auto it = loop.begin();
        ASSERT_LE(distance((*(it++))->vertex->coords, std::array<T, d>{1.0,1.0}),1e-6);
        ASSERT_LE(distance((*(it++))->vertex->coords, std::array<T, d>{0.0,1.0}),1e-6);
        ASSERT_LE(distance((*(it++))->vertex->coords, std::array<T, d>{0.0,0.0}),1e-6);
    }
    {
        auto loop = getFaceEdges(hf2);
        ASSERT_TRUE(loop.size()==3);
        auto it = loop.begin();
        ASSERT_LE(distance((*(it++))->vertex->coords, std::array<T, d>{0.0,0.0}),1e-6);
        ASSERT_LE(distance((*(it++))->vertex->coords, std::array<T, d>{1.0,0.0}),1e-6);
        ASSERT_LE(distance((*(it++))->vertex->coords, std::array<T, d>{1.0,1.0}),1e-6);
    }

    auto faces_lst = {hf1,hf2};

    auto vertices_map = extract_vertices_map<T,d>(faces_lst);

    ASSERT_EQ( vertices_map.size(), 4 );

    flip(hf1,hf2);
    
    auto polyData = make_polydata<T,d>(faces_lst);

    ASSERT_EQ(polyData->GetNumberOfCells(),2);

    // // Create mapper and actor
    // vtkNew<vtkPolyDataMapper> mapper;
    // mapper->SetInputData(polyData);

    // vtkNew<vtkActor> actor;
    // actor->SetMapper(mapper);
    // actor->GetProperty()->SetEdgeVisibility(true);
    // actor->GetProperty()->SetOpacity(0.3);

    // gbs::plot(
    //     actor.Get()
    //     // , ell
    //     );


}


template <typename T>
    const std::array<T,2> &a,
    const std::array<T,2> &b,
{
}

{
}

TEST(halfEdgeMesh, delauney2d_mesh_cloud)
{

}