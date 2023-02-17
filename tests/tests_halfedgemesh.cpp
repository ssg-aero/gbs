#include <gtest/gtest.h>
#include <gbs/maths.h>
#include <gbs-render/vtkcurvesrender.h>
#include <topology/tessellations.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>


TEST(halfEdgeMesh, getVertexMainLoop)
{
    using namespace gbs;
    using T = double;
    const size_t d = 2;

    auto he1 = make_shared_h_edge<T,d>({1.0,0.0});
    auto he2 = make_shared_h_edge<T,d>({0.0,1.0});
    auto he3 = make_shared_h_edge<T,d>({0.0,0.0});

    auto lst1 = {he1, he2, he3};
    auto hf1 = make_shared_h_face<T,d>(lst1);

    ASSERT_EQ(he1->face, hf1);
    ASSERT_EQ(he2->face, hf1);
    ASSERT_EQ(he3->face, hf1);


    ASSERT_EQ(he1, getFaceEdge(hf1,he1->vertex));
    ASSERT_EQ(he2, getFaceEdge(hf1,he2->vertex));
    ASSERT_EQ(he3, getFaceEdge(hf1,he3->vertex));

    auto hf2 = add_face(hf1,he2,{1.0,1.0});

    {
        auto loop = getFaceEdges(hf2);
        ASSERT_TRUE(loop.size()==3);
        auto it = loop.begin();
        ASSERT_LE(distance((*(it++))->vertex->coords, std::array<T, d>{0.0,1.0}),1e-6);
        ASSERT_LE(distance((*(it++))->vertex->coords, std::array<T, d>{1.0,0.0}),1e-6);
        ASSERT_LE(distance((*(it++))->vertex->coords, std::array<T, d>{1.0,1.0}),1e-6);
    }

    auto faces_lst = {hf1,hf2};

    auto vertices_map = extract_vertices_map<T,d>(faces_lst);

    ASSERT_EQ( vertices_map.size(), 4 );
    
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
T orient_2d(
    const std::array<T,2> &a,
    const std::array<T,2> &b,
    const std::array<T,2> &c)
{
    auto [ax,ay] = a;
    auto [bx,by] = b;
    auto [cx,cy] = c;
    return (ax-cx)*(by-cy) - (ay-cy)*(bx-cx);
}

template <typename T>
T are_ccw(
    const std::array<T,2> &a,
    const std::array<T,2> &b,
    const std::array<T,2> &c)
{
    return orient_2d(a,b,c) >= 0;
}

TEST(halfEdgeMesh, delauney2d_mesh_cloud)
{

}