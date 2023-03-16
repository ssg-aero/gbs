#include <gtest/gtest.h>
#include <gbs/maths.h>
#include <gbs/vecop.h>
#include <gbs-render/vtkcurvesrender.h>
#include <topology/tessellations.h>
#include <topology/vertex.h>
#include <topology/edge.h>
#include <topology/wire.h>
#include <topology/baseIntersection.h>
#include <topology/render.h>
#include <topology/halfEdgeMeshQuality.h>
#include <numbers>
#include <random> 
#include <list>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkXMLPolyDataWriter.h>

#include "tests_helpers.h"

using namespace gbs;
const bool plot_on = true;

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

template <typename T, size_t dim>
inline auto mesh_hed_wire_uniform(const Wire<T,dim> &w, T dm)
{
    auto coords = mesh_wire_uniform(w, dm);
    long long n = coords.size();
    std::vector<std::shared_ptr< HalfEdge<T,dim> > > h_edges(n);
    std::transform(
        coords.begin(), coords.end(),
        h_edges.begin(),
        make_shared_h_edge<T,dim>
    );

    auto nm = n-1;
    h_edges.front()->next = h_edges[1];
    h_edges.front()->previous = h_edges.back();
    for( long long i{1}; i < nm; i++)
    {
        h_edges[i]->previous = h_edges[i-1];
        h_edges[i]->next = h_edges[i+1];
    }
    h_edges.back()->next =h_edges.front();
    h_edges.back()->previous =  h_edges[nm];

    return h_edges;
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
inline auto make_boundary2d_2(T dm = 0.1, T R = 1., T r = 0.5)
{
    auto ell = build_ellipse<T, 2>(R, r);
    Wire2d<T> w(Edge<T, 2>{std::make_shared<BSCurveRational<T, 2>>(ell)});

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

    auto vertices_map = getVerticesMapFromFaces<T,d>(faces_lst);

    ASSERT_EQ( vertices_map.size(), 4 );
    
    auto polyData = make_polydata_from_faces<T,d>(faces_lst);

    ASSERT_EQ(polyData->GetNumberOfCells(),2);

    if (plot_on)
    { 
        gbs::plot(
            faces_mesh_actor<T,d>(faces_lst).Get()
        );
    }
}

TEST(halfEdgeMesh, add_vertex)
{

    using T = double;
    const size_t d = 2;

    std::list<std::array<T,d>> coords{{1,1}, {4,-2}, {3,3}};

    std::shared_ptr<HalfEdgeVertex<T,d>> vtx_ref;

    auto A = tri_area(
        *(std::next(coords.begin(), 0)),
        *(std::next(coords.begin(), 1)),
        *(std::next(coords.begin(), 2)));


    auto vertices = make_shared_h_vertices<T,d>(coords);
    ASSERT_EQ(vertices.size(),3);
    for(size_t i{}; i < 3 ; i++)
    {
        ASSERT_TRUE(vertices[i]);
        const std::array<T,d> &XYZ =  *(std::next(coords.begin(), i));
        auto d = distance(vertices[i]->coords, XYZ );
        ASSERT_LE( d, 1e-6);
    }

    auto edges = make_shared_h_edges<T,d>(coords);

    std::list<std::shared_ptr<HalfEdgeFace<T, d>>> faces_lst{make_shared_h_face<T,d>(edges)};

    {
        auto vtx1 = make_shared_h_vertex<T, d>({2., 1});
        auto new_faces = add_vertex(faces_lst.front(), vtx1);
        faces_lst.remove(faces_lst.front());
        faces_lst.insert(faces_lst.end(), new_faces.begin(), new_faces.end());

        ASSERT_EQ(faces_lst.size(), 3);
        ASSERT_NEAR(getTriangle2dMeshArea(faces_lst), A, 1e-9);
    }


    {
        auto vtx = make_shared_h_vertex<T, d>({3., 0});
        auto it_face = std::find_if(
            faces_lst.begin(), faces_lst.end(),
            [pt=vtx->coords](const auto &h_f){
                auto coords = getFaceCoords(h_f);
                auto start = coords.begin();
                return in_triangle(
                    *(std::next(start,0)),
                    *(std::next(start,1)),
                    *(std::next(start,2)),
                    pt
                );
            }
        );
        if(it_face!=faces_lst.end())
        {
            auto new_faces = add_vertex(*it_face, vtx);
            faces_lst.remove(*it_face);
            faces_lst.insert(faces_lst.end(), new_faces.begin(), new_faces.end());
        }

        ASSERT_EQ(faces_lst.size(), 5);
        ASSERT_NEAR(getTriangle2dMeshArea(faces_lst), A, 1e-9);

        vtx_ref = vtx;
    }

    auto local_faces = getFacesAttachedToVertex(vtx_ref);
    ASSERT_EQ(local_faces.size(),3);

    if (plot_on)
    {
        gbs::plot(
            faces_mesh_actor<T,d>(faces_lst).Get(),
            faces_mesh_actor<T,d>(local_faces,{1.,0.,0.}).Get()
        );
    }
}

TEST(halfEdgeMesh, getFacesAttachedToVertex)
{
    using T = double;
    const size_t d = 2;

    auto he1 = make_shared_h_edge<T,d>({1.0,0.0});
    auto he2 = make_shared_h_edge<T,d>({1.0,1.0});
    auto he3 = make_shared_h_edge<T,d>({0.0,1.0});
    auto he4 = make_shared_h_edge<T,d>({0.0,0.0}); 

    auto lst1 = {he1, he2, he3, he4};
    auto hf1 = make_shared_h_face<T,2>(lst1);

    auto hf2 = add_face(he1,{ 0.5,-0.5});
    auto hf3 = add_face(he2,{ 1.5, 0.5});
    auto hf4 = add_face(he3,{ 0.5, 1.5});
    auto hf5 = add_face(he4,{-0.5, 0.5});

    auto neighbors = getNeighboringFaces(hf1);
    ASSERT_EQ(*std::next(neighbors.begin(),0), hf2);
    ASSERT_EQ(*std::next(neighbors.begin(),1), hf3);
    ASSERT_EQ(*std::next(neighbors.begin(),2), hf4);
    ASSERT_EQ(*std::next(neighbors.begin(),3), hf5);

    auto faces_lst = {hf1, hf2, hf3, hf4, hf5};

    if (plot_on)
    {
        gbs::plot(
            faces_mesh_actor<T,d>(faces_lst).Get()
        );
    }

}

TEST(halfEdgeMesh, add_delaunay_point)
{

    using T = double;
    const size_t d = 2;

    auto he1 = make_shared_h_edge<T,d>({1.0,0.0});
    auto he2 = make_shared_h_edge<T,d>({0.0,1.0});
    auto he3 = make_shared_h_edge<T,d>({0.0,0.0});
    auto lst1 = {he1, he2, he3};
    auto hf1 = make_shared_h_face<T,d>(lst1);
    auto hf2 = add_face(hf1,he2,{1.0,1.0});

    {
        auto loop = getFaceEdges(hf2);
        ASSERT_TRUE(loop.size()==3);
        auto it = loop.begin();
        ASSERT_LE(distance((*(it++))->vertex->coords, std::array<T, d>{0.0,1.0}),1e-6);
        ASSERT_LE(distance((*(it++))->vertex->coords, std::array<T, d>{1.0,0.0}),1e-6);
        ASSERT_LE(distance((*(it++))->vertex->coords, std::array<T, d>{1.0,1.0}),1e-6);
    }

    std::list< std::shared_ptr<HalfEdgeFace<T, d>> > faces_lst{hf1,hf2};

    boyerWatson<T>(faces_lst, {0.75,0.25} );
    for(const auto &hf: faces_lst)
    {
        ASSERT_TRUE(is_ccw(hf));
    }
    boyerWatson<T>(faces_lst, {0.25,0.25} );
    for(const auto &hf: faces_lst)
    {
        ASSERT_TRUE(is_ccw(hf));
    }
    boyerWatson<T>(faces_lst, {0.5,1.} );
    for(const auto &hf: faces_lst)
    {
        ASSERT_TRUE(is_ccw(hf));
    }
    boyerWatson<T>(faces_lst, {0.5,0.} );
    for(const auto &hf: faces_lst)
    {
        ASSERT_TRUE(is_ccw(hf));
    }
    boyerWatson<T>(faces_lst, {0.,1./3.} );
    for(const auto &hf: faces_lst)
    {
        ASSERT_TRUE(is_ccw(hf));
    }
    boyerWatson<T>(faces_lst, {0.,2./3.} );
    for(const auto &hf: faces_lst)
    {
        ASSERT_TRUE(is_ccw(hf));
    }
    if (plot_on)
    { 
        gbs::plot(
            faces_mesh_actor<T,d>(faces_lst).Get()
        );
    }
}

TEST(halfEdgeMesh, add_delaunay_points)
{
    using T = double;
    const size_t d = 2;

    auto coords = make_boundary2d_1<T>(0.3);

    auto boundary = make_HalfEdges(coords);

    // add_random_points_grid(coords,11);

    auto faces_lst = delaunay2DBoyerWatson<T>(coords);

    for(const auto &hf: faces_lst)
    {
        ASSERT_TRUE(is_ccw(hf));
    }

    ASSERT_NEAR(getTriangle2dMeshArea(faces_lst), 1., 1e-6);

    if (plot_on)
    { 
        gbs::plot(
            faces_mesh_actor<T,d>(faces_lst).Get(),
            boundary_mesh_actor<T,d>(boundary).Get(),
            coords
        );
    }
}

TEST(halfEdgeMesh, is_inside_boundary)
{
    using T = double;
    const size_t d = 2;

    auto coords = make_boundary2d_2<T>();
    // auto coords = make_boundary2d_1<T>(.1);
    auto boundary = make_HalfEdges(coords);

    ASSERT_TRUE(are_edges_2d_ccw(boundary));
    reverseBoundary(boundary);
    ASSERT_FALSE(are_edges_2d_ccw(boundary));

    add_random_points_ellipse(coords,11);

    // auto faces_lst = delaunay2DBoyerWatson<T>(coords);

    auto faces_lst = getEncompassingMesh(coords);
    auto vertices_map = getVerticesMapFromFaces<T,2>(faces_lst);
    // insert points
    for(const auto &xy : coords)
    {
        boyerWatson<T>(faces_lst, xy);
    }

    ASSERT_TRUE(seg_H_strict_end_intersection<T>({1,-1},{1,1},{0,0}));
    ASSERT_TRUE(seg_H_strict_end_intersection<T>({.1,-1},{.1,1},{0,0}));
    ASSERT_FALSE(seg_H_strict_end_intersection<T>({-1,-1},{-1,1},{0,0}));
    ASSERT_FALSE(seg_H_strict_end_intersection<T>({1,-1},{1,1},{2,1}));
    ASSERT_FALSE(seg_H_strict_end_intersection<T>({1,-1},{1,1},{0,2}));
    ASSERT_TRUE(seg_H_strict_end_intersection<T>({0.5,-1},{1.5,1},{0,0}));
    ASSERT_FALSE(seg_H_strict_end_intersection<T>({0.5,-1},{1.5,1},{1.5,1}));

    ASSERT_TRUE(on_seg<T>({-2,0},{2,1},{0,0.5}));
    ASSERT_FALSE(on_seg<T>({-2,0},{2,1},{0,0.75}));

    ASSERT_TRUE(is_inside<T>({0.,0.}, boundary));

    auto external_faces = takeExternalFaces<T>(faces_lst,boundary);

    for(const auto &hf: faces_lst)
    {
        auto coords = getFaceCoords(hf);
        for(const auto &xy : coords)
        {
            ASSERT_TRUE(is_inside(xy,boundary));
        }
    }
    for(const auto &hf: external_faces)
    {
        auto coords = getFaceCoords(hf);
        bool is_in{true};
        for(const auto &xy : coords)
        {
            is_in = is_in && is_inside(xy,boundary);
        }
        ASSERT_FALSE(is_in);
    }

    if (plot_on)
    { 

        auto polyData = make_polydata_from_faces<T,d>(faces_lst);
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputData(polyData);
        vtkNew<vtkActor> actor;
        actor->SetMapper(mapper);
        actor->GetProperty()->SetEdgeVisibility(true);
        // actor->GetProperty()->SetOpacity(0.3);
        actor->GetProperty()->SetColor(0.,0.,1.);

        auto polyData_external = make_polydata_from_faces<T,d>(external_faces);
        vtkNew<vtkPolyDataMapper> mapper_external;
        mapper_external->SetInputData(polyData_external);
        vtkNew<vtkActor> actor_external;
        actor_external->SetMapper(mapper_external);
        actor_external->GetProperty()->SetEdgeVisibility(true);
        // actor_external->GetProperty()->SetOpacity(0.3);
        actor_external->GetProperty()->SetColor(0.,1.,0.);


        auto polyData_boundary = make_polydata_from_edges_loop<T,d>(
            boundary
        );
        vtkNew<vtkPolyDataMapper> mapper_boundary;
        mapper_boundary->SetInputData(polyData_boundary);

        vtkNew<vtkActor> actor_boundary;
        actor_boundary->SetMapper(mapper_boundary);
        
        actor_boundary->GetProperty()->SetColor(1.,0.,0.);


        gbs::plot(
            faces_mesh_actor<T,d>(faces_lst,{0.,0.,1.}).Get(),
            faces_mesh_actor<T,d>(external_faces,{0.,1.,0.}).Get(),
            boundary_mesh_actor<T,d>(boundary).Get()
        );
    }
}

TEST(halfEdgeMesh, delaunay2d_non_convex_boundary)
{
    using T = double;
    auto coords = make_boundary2d_3<T>(1.);
    auto boundary = make_HalfEdges(coords);
    auto faces_lst = delaunay2DBoyerWatson<T>(coords);
    auto external_faces = takeExternalFaces<T>(faces_lst,boundary);

    ASSERT_NEAR(getTriangle2dMeshArea(faces_lst), 6.5, 1e-6);

    if (plot_on)
    {
        auto polyData = make_polydata_from_faces<T,2>(faces_lst);
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputData(polyData);
        vtkNew<vtkActor> actor;
        actor->SetMapper(mapper);
        actor->GetProperty()->SetEdgeVisibility(true);
        actor->GetProperty()->SetOpacity(0.3);

        auto polyData_boundary = make_polydata_from_edges_loop<T,2>(
            boundary
        );
        vtkNew<vtkPolyDataMapper> mapper_boundary;
        mapper_boundary->SetInputData(polyData_boundary);
        vtkNew<vtkActor> actor_boundary;
        actor_boundary->SetMapper(mapper_boundary);
        actor_boundary->GetProperty()->SetColor(1.,0.,0.);

        gbs::plot(
            faces_mesh_actor<T,2>(faces_lst).Get(),
            boundary_mesh_actor<T,2>(boundary).Get(),
            coords
        );
    }
}

TEST(halfEdgeMesh, delaunay2d_inner_boundary)
{
    using T = double;
    const size_t d{2};
    using namespace gbs;

    T r = 0.3;
    auto coords_outer = make_boundary2d_1<T>(0.1);
    auto coords_inner =  make_boundary2d_2<T>(0.01, r, r);
    for(auto &xy: coords_inner)
    {
        xy = xy + std::array<T,2>{0.35,0.6};
    }
    std::reverse(coords_inner.begin(),coords_inner.end());

    std::vector<std::array<T,d>> coords{coords_outer.begin(), coords_outer.end()};
    coords.insert(coords.end(),coords_inner.begin(), coords_inner.end());

    auto boundary_outer = make_HalfEdges<T>(coords_outer);
    auto boundary_inner = make_HalfEdges<T>(coords_inner);

    auto faces_lst = delaunay2DBoyerWatson<T>(coords);
    auto internal_faces = takeInternalFaces<T>(faces_lst,boundary_inner);

    ASSERT_NEAR(getTriangle2dMeshAreaPar(faces_lst),1-std::numbers::pi_v<T>*r*r,1e-4);
    // ASSERT_NEAR(getTriangle2dMeshArea(faces_lst),1-std::numbers::pi_v<T>*r*r,1e-4);
    if (plot_on)
    {
        gbs::plot(
            faces_mesh_actor<T,d>(faces_lst).Get(),
            boundary_mesh_actor<T,d>(boundary_outer).Get(),
            boundary_mesh_actor<T,d>(boundary_inner).Get(),
            coords);
    }
}    

TEST(halfEdgeMesh, delaunay2d_mesh_cloud)
{
    using T = double;
    const size_t d{2};
    using namespace gbs;

    auto coords_boundary = make_boundary2d_1<T>();
    std::vector< std::array<T,2> > coords_inner;
    add_random_points_grid(coords_inner, 777);

    auto faces_lst = delaunay2DBoyerWatson(coords_boundary, coords_inner, 1e-10);

    std::cout << getTriangle2dMeshArea(faces_lst) << std::endl;

    ASSERT_NEAR(getTriangle2dMeshArea(faces_lst), 1., 1e-6);
    if (plot_on)
    {
        auto polyData = make_polydata_from_faces<T, d>(faces_lst);

        vtkNew<vtkXMLPolyDataWriter> writer;
        writer->SetInputData(polyData);
        writer->SetFileName("toto.vtp");
        writer->Write();

        gbs::plot(
            faces_mesh_actor<T,d>(faces_lst).Get(),
            boundary_mesh_actor<T,d>(getOrientedFacesBoundary(std::list< std::shared_ptr<HalfEdgeFace<T, d>> >(faces_lst))).Get(),
            coords_boundary);
    }
}

template <typename T, size_t dim>
struct DistanceMeshSurface
{
    const Surface<T, dim> &srf;
    DistanceMeshSurface(const auto &srf_) : srf{srf_} {}
    auto operator()(const std::shared_ptr<HalfEdgeFace<T, 2>> &hf)
    {
        auto coordsUV = getFaceCoords(hf);
        auto [u1, u2, v1, v2] = srf.bounds();
        for (auto [u, v] : coordsUV)
        {
            if (u1 > u || u2 < u || v1 > v || v2 < v)
            // if (u1 >= u || u2 <= u || v1 >= v || v2 <= v)
            {
                return std::make_pair(0., std::array<T, 2>{0., 0.}); // ext dummy are perfect
            }
        }
        std::vector<std::array<T, 3>> coords(coordsUV.size());
        std::transform(
            coordsUV.begin(), coordsUV.end(),
            coords.begin(),
            [&srf = this->srf](const auto &UV)
            {
                return srf(UV[0], UV[1]);
            });
        auto G_UV = std::reduce(
                        coordsUV.begin(), coordsUV.end(),
                        std::array<T, 2>{0., 0.},
                        gbs::operator+<T, 2>) /
                    3.;
        auto G = std::reduce(
                     coords.begin(), coords.end(),
                     std::array<T, 3>{0., 0., 0.},
                     gbs::operator+<T, 3>) /
                 3.;

        auto A_UV = 0.5 * (*std::next(coordsUV.begin(), 0) + *std::next(coordsUV.begin(), 1));
        auto A = 0.5 * (*std::next(coords.begin(), 0) + *std::next(coords.begin(), 1));
        auto B_UV = 0.5 * (*std::next(coordsUV.begin(), 1) + *std::next(coordsUV.begin(), 2));
        auto B = 0.5 * (*std::next(coords.begin(), 1) + *std::next(coords.begin(), 2));
        auto C_UV = 0.5 * (*std::next(coordsUV.begin(), 2) + *std::next(coordsUV.begin(), 0));
        auto C = 0.5 * (*std::next(coords.begin(), 2) + *std::next(coords.begin(), 0));

        auto test_lst = {
            std::make_pair(distance(G, srf(G_UV[0], G_UV[1])), G_UV),
            std::make_pair(distance(A, srf(A_UV[0], A_UV[1])), A_UV),
            std::make_pair(distance(B, srf(B_UV[0], B_UV[1])), B_UV),
            std::make_pair(distance(C, srf(C_UV[0], C_UV[1])), C_UV)};

        return *std::max_element(
            test_lst.begin(), test_lst.end(),
            [](const auto &p1, const auto &p2)
            { return p1.first < p2.first; });
    }
};

// template <typename T, size_t dim>
// struct DistanceMeshSurface2
// {
//     const Surface<T, dim> &srf;
//     DistanceMeshSurface2(const auto &srf_) : srf{srf_} {}
//     auto operator()(const std::shared_ptr<HalfEdgeFace<T, 2>> &hf)
//     {
//         auto ed_lst = getFaceEdges(hf);
//         auto n = ed_lst.size();
//         std::vector< std::pair<T, std::array<T,2> > > test_lst(n);

//         std::transform(
//             ed_lst.begin(), ed_lst.end(),
//             test_lst.begin(),
//             [&srf = this->srf](const auto &he){
//                 assert(he && he->vertex && he->previous && he->previous->vertex);
//                 if(!he->opposite)
//                 {
//                     return std::make_pair(0., std::array<T, 2>{0., 0.}); // boundary edges are supposed to respect criteria
//                 }
//                 const auto &p1_uv = he->previous->vertex->coords;
//                 const auto &p2_uv = he->vertex->coords;
//                 auto mid_uv = p1_uv + static_cast<T>(0.5) *(p2_uv-p1_uv);
//                 auto p1 = srf(p1_uv[0], p1_uv[1]);
//                 auto p2 = srf(p2_uv[0], p2_uv[1]);
//                 auto mid = p1 + static_cast<T>(0.5) *(p2-p1);
//                 return std::make_pair(distance(mid, srf(mid_uv[0], mid_uv[1])), mid_uv);
//             }
//         );

//         return *std::max_element(
//             test_lst.begin(), test_lst.end(),
//             [](const auto &p1, const auto &p2)
//             { return p1.first < p2.first; });
//     }
// };

template <typename T, size_t dim>
struct DistanceMeshSurface2
{
    const Surface<T, dim> &srf;
    DistanceMeshSurface2(const Surface<T, dim> &srf_) : srf(srf_) {}

    auto operator()(const std::shared_ptr<HalfEdgeFace<T, 2>> &hf)
    {
        auto ed_lst = getFaceEdges(hf);
        auto n = ed_lst.size();

        std::vector< std::pair<T, std::array<T,2> > > test_lst;
        test_lst.reserve(ed_lst.size());

        std::transform(
            ed_lst.cbegin(), ed_lst.cend(),
            std::back_inserter(test_lst),
            [&srf = this->srf](const auto &he){
                assert(he && he->vertex && he->previous && he->previous->vertex);
                if(!he->opposite)
                {
                    return std::make_pair(0., std::array<T, 2>{0., 0.}); // boundary edges are supposed to respect criteria
                }
                const auto &p1_uv = he->previous->vertex->coords;
                const auto &p2_uv = he->vertex->coords;
                auto mid_uv = p1_uv + static_cast<T>(0.5) *(p2_uv-p1_uv);
                auto mid = srf(mid_uv[0], mid_uv[1]);
                return std::make_pair(distance(mid, static_cast<T>(0.5) *( srf(p1_uv[0], p1_uv[1]) + srf(p2_uv[0], p2_uv[1]) )), mid_uv);
            }
        );

        const auto max_it = std::max_element(
            test_lst.cbegin(), test_lst.cend(),
            [](const auto &p1, const auto &p2) { return p1.first < p2.first; });

        return *max_it;
    }
};

template <typename T, size_t dim>
struct DistanceMeshSurface3
{
    const Surface<T, dim> &srf;
    DistanceMeshSurface3(const auto &srf_) : srf{srf_} {}
    auto operator()(const std::shared_ptr<HalfEdgeFace<T, 2>> &hf)
    {
        auto coords_lst = getFaceCoords(hf);
        auto ed_lst = getFaceEdges(hf);
        auto n = coords_lst.size();

        auto mid_uv = static_cast<T>(1./n) * std::reduce(
            coords_lst.begin(), coords_lst.end(),
            std::array<T,2>{},
            [](const auto &X1, const auto &X2){return X1+X2;}
        );

        // auto mid = static_cast<T>(1./n) * std::reduce(
        //     coords_lst.begin(), coords_lst.end(),
        //     std::array<T,dim>{},
        //     [&srf = this->srf](const auto &X1, const auto &X2)
        //     {
        //         return srf(X1[0],X1[1]) + srf(X2[0],X2[1]) );
        //     }
        // );
        std::array<T,dim> mid{};
        for(const auto &X : coords_lst)
        {
            mid = mid + srf(X[0],X[1]);
        }
        mid = static_cast<T>(1./n) * mid;

        std::vector< std::pair<T, std::array<T,2> > > test_lst(n+1);

        std::transform(
            ed_lst.begin(), ed_lst.end(),
            test_lst.begin(),
            [&srf = this->srf](const auto &he){
                assert(he && he->vertex && he->previous && he->previous->vertex);
                if(!he->opposite)
                {
                    return std::make_pair(0., std::array<T, 2>{0., 0.}); // boundary edges are supposed to respect criteria
                }
                const auto &p1_uv = he->previous->vertex->coords;
                const auto &p2_uv = he->vertex->coords;
                auto mid_uv = p1_uv + static_cast<T>(0.5) *(p2_uv-p1_uv);
                auto p1 = srf(p1_uv[0], p1_uv[1]);
                auto p2 = srf(p2_uv[0], p2_uv[1]);
                auto mid = p1 + static_cast<T>(0.5) *(p2-p1);
                return std::make_pair(distance(mid, srf(mid_uv[0], mid_uv[1])), mid_uv);
            }
        );

        test_lst.back() = std::make_pair(distance(mid, srf(mid_uv[0], mid_uv[1])), mid_uv);
        
        return *std::max_element(
            test_lst.begin(), test_lst.end(),
            [](const auto &p1, const auto &p2)
            { return p1.first < p2.first; });

    }
};

TEST(halfEdgeMesh, delaunay2d_mesh_surface)
{
    using T = double;
    const size_t d{3};
    using namespace gbs;

    BSSurface<T,d> srf{
        {
            {0.,0.,0.},{0.5,0.,0.5},{1.,0.,0.5},
            {0.,1.,0.},{0.5,1.,0.0},{0.,0.,1.}
        },
        {0.,1.},
        {0.,1.},
        {3,3},
        {2,2},
        2,1
    };

    // auto faces_lst = delaunay2DBoyerWatsonSurfaceMesh<T,d,DistanceMeshSurface<T,d>>(std::make_shared<BSSurface<T,d>>(srf), 0.005);
    auto faces_lst = delaunay2DBoyerWatsonSurfaceMesh<T,d,DistanceMeshSurface2<T,d>>(std::make_shared<BSSurface<T,d>>(srf), 0.005);
    // auto faces_lst = delaunay2DBoyerWatsonSurfaceMesh<T,d,DistanceMeshSurface2<T,d>>(std::make_shared<BSSurface<T,d>>(srf), 0.001, 5000, 5, 5, 0.005, 1e-10);
    // auto faces_lst = delaunay2DBoyerWatsonSurfaceMesh<T,d,DistanceMeshSurface3<T,d>>(std::make_shared<BSSurface<T,d>>(srf), 0.001, 500, 5, 5, 0.005 );

    // auto execution_time = measure_execution_time(&delaunay2DBoyerWatsonSurfaceMesh<T,d,DistanceMeshSurface2<T,d>>, std::make_shared<BSSurface<T,d>>(srf), 0.001, 5000, 5, 5, 0.15, 1e-10 );
    // auto execution_time = measure_execution_time<std::chrono::microseconds>(&delaunay2DBoyerWatsonSurfaceMesh<T,d,DistanceMeshSurface2<T,d>>, std::make_shared<BSSurface<T,d>>(srf), 0.001, 5000, 5, 5, 0.005, 1e-10 );
    //     std::cout << "Execution time: " << execution_time / 1000. << " ms" << std::endl;

    if (plot_on)
    {
        gbs::plot(
            surface_mesh_actor<T>(faces_lst, srf, { 51./255.,  161./255.,  201./255.}, true).Get()
        );
    }
}

auto f_curve2d_offset_functor(double r, double h)
{
    auto circle1 = gbs::build_circle<double, 2>(r,{0.5,0.5});
    size_t n_pulse = 10;
    auto f_offset = [n_pulse, h](auto u, size_t d = 0){return d==0 ? (-std::sin(u*2.*std::numbers::pi*n_pulse)*h - h) : -h*2.*std::numbers::pi*n_pulse*std::cos(u*2.*std::numbers::pi*n_pulse);};
    auto p_circle1 = std::make_shared<gbs::BSCurveRational<double, 2>>(circle1);

    gbs::CurveOffset<double, 2,decltype(f_offset)> circle2{
        p_circle1,
        std::make_shared<decltype(f_offset)>( f_offset )
        };
    return std::make_tuple(circle1,circle2,f_offset);
}


TEST(halfEdgeMesh, delaunay2d_mesh_face)
{
    using T = double;
    const size_t d{3};
    using namespace gbs;

    BSSurface<T,d> srf{
        {
            {0.,0.,0.5}  ,{0.33,0.,0.0},  {0.66,0.,0.0}  ,{1.,0.,0.5},
            {0.,0.33,0.},{0.33,0.33,1.5},{0.66,0.33,0.7},{1.,0.33,0.},
            {0.,0.66,0.},{0.33,0.66,0.3},{0.66,0.66,2.5},{1.,0.66,0.},
            {0.,1.,0.5},{0.33,1.,0.0},{0.66,1.,0.0},{1.,1.,0.5},
        },

        {0.,1.},
        {0.,1.},
        {4,4},
        {4,4},
        3,3
    };

    T r = 0.3, h= 0.05;
    auto [circle1,circle2,f_offset] = f_curve2d_offset_functor(r,h); // check if working while build in a factory
    auto coords_inner = discretize(circle2,300);

    auto [circle1o,circle2o,f_offseto] = f_curve2d_offset_functor(r+3*h,h); // check if working while build in a factory
    auto coords_outer = discretize(circle1o,300);

    std::shared_ptr<Surface<T,d>> p_srf = std::make_shared<BSSurface<T,d>>(srf);
    T deviation{0.01};
    T tol{1e-10};
    size_t nu{15}, nv{15};
    auto faces_lst = delaunay2DBoyerWatsonSurfaceBase(p_srf, nu, nv, deviation);

    delaunay2DBoyerWatsonAddInnerBound(faces_lst, coords_inner);
    delaunay2DBoyerWatsonAddOuterBound(faces_lst, coords_outer);

    T crit_max{0.005};
    delaunay2DBoyerWatsonSurfaceMeshRefine<T,d,DistanceMeshSurface2<T,d>>(p_srf, faces_lst, crit_max);

    std::vector<std::array<T,d>> coords_inner_3d(coords_inner.size());
    std::transform(
        coords_inner.begin(), coords_inner.end(),
        coords_inner_3d.begin(),
        [&srf](const auto &uv){
            return srf(uv[0],uv[1]);
        }
    );
    auto boundary_inner_3d = make_HalfEdges<T>(coords_inner_3d);

    std::vector<std::array<T,d>> coords_outer_3d(coords_outer.size());
    std::transform(
        coords_outer.begin(), coords_outer.end(),
        coords_outer_3d.begin(),
        [&srf](const auto &uv){
            return srf(uv[0],uv[1]);
        }
    );
    auto boundary_outer_3d = make_HalfEdges<T>(coords_outer_3d);

    if (plot_on)
    {
        gbs::plot(
            surface_mesh_actor<T>(faces_lst, srf, { 51./255.,  161./255.,  201./255.}, true).Get(),
            boundary_mesh_actor<T,d>(boundary_inner_3d).Get(),
            boundary_mesh_actor<T,d>(boundary_outer_3d,{0.,1.,0.}).Get()
        );
    }
}