#include <gtest/gtest.h>
#include <topology/vertex.h>
#include <topology/edge.h>
#include <topology/wire.h>
#include <gbs-render/vtkcurvesrender.h>
#include <random> 
#include <chrono>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>

template <typename T>
T in_circle(
    const std::array<T,2> &a,
    const std::array<T,2> &b,
    const std::array<T,2> &c,
    const std::array<T,2> &d
    )
{
    auto [ax,ay] = a;
    auto [bx,by] = b;
    auto [cx,cy] = c;
    auto [dx,dy] = d;

    auto A = ax - dx;
    auto B = ay - dy;
    auto C = A*A + B*B;
    
    auto D = bx - dx;
    auto E = by - dy;
    auto F = D*D + E*E;

    auto G = cx - dx;
    auto H = cy - dy;
    auto I = G*G + H*H;

    return A*E*I + D*H*C + B*F*G - ( G*E*C + D*B*I + A*H*F);
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

TEST(tests_mesh_delauney2d, in_circle)
{

    {

        std::array<double,2> a{0.0, 0.0};
        std::array<double,2> b{0.1, 0.0};
        std::array<double,2> c{0.5, 0.5};

        ASSERT_TRUE( are_ccw(a, b, c) );
        ASSERT_FALSE(are_ccw(b, a, c) );

        ASSERT_TRUE( orient_2d(a, b, c) > 0 ); // 
        ASSERT_TRUE( orient_2d(b, a, c) < 0 );

        ASSERT_DOUBLE_EQ( in_circle<double>(a, b, c, a ) , 0);
        ASSERT_DOUBLE_EQ( in_circle<double>(a, b, c, b ) , 0);
        ASSERT_DOUBLE_EQ( in_circle<double>(a, b, c, c ) , 0);

        {
            std::array<double,2> d{0.5,-0.5};

            ASSERT_TRUE( in_circle<double>(a, b, c, d ) < 0);
        }

        {
            std::array<double,2> d{0.0,0.5};

            ASSERT_TRUE( in_circle<double>(a, b, c, d ) < 0);
        }

        {
            std::array<double,2> d{1.0,0.5};

            ASSERT_TRUE( in_circle<double>(a, b, c, d ) < 0);
        }

        {
            std::array<double,2> d{0.5,0.3};

            ASSERT_TRUE( in_circle<double>(a, b, c, d ) > 0);
        }

    }

    // {
    //     size_t n = 500;
    //     std::random_device rd;  // Will be used to obtain a seed for the random number engine
    //     std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    //     std::uniform_real_distribution<double> distrib(-2, 2);
    //     const double sqrt2 = std::sqrt(2);
    //     for (size_t i{}; i < n; i++)
    //     {
    //         std::array<double, 2> d{distrib(gen), distrib(gen)};

    //         auto test = in_circle<double>(
    //             {-1, -1},
    //             {1, -1},
    //             {0., sqrt2},
    //             d);
    //         auto sq_r = d[0] * d[0] + d[1] * d[1];
    //         ASSERT_TRUE((test >= 0) && (sq_r <= 2));
    //     }
    // }
}

TEST(tests_mesh_delauney2d, mesh_cloud)
{

    using T = double;
    const size_t d{2};

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    // std::uniform_real_distribution<double> distrib(0.0,1.0);
 
    // size_t n = 50;

    // std::vector< std::array<T,d> > vertices_cloud(n);

    // std::generate(
    //     vertices_cloud.begin(),
    //     vertices_cloud.end(),
    //     [&](){
    //         return std::array<T,2>{distrib(gen),distrib(gen)};
    //     }
    // );

    std::uniform_real_distribution<double> distribr1(0.0,1.0);
    std::uniform_real_distribution<double> distribr2(0.0,0.5);
    std::uniform_real_distribution<double> distribth(0.0,2*std::acos(-1.));
 
    size_t n = 200;

    std::vector< std::array<T,d> > vertices_cloud(n);

    std::generate(
        vertices_cloud.begin(),
        vertices_cloud.end(),
        [&](){
            auto th = distribth(gen);
            return std::array<T,2>{distribr1(gen)*std::cos(th),distribr2(gen)*std::sin(th)};
        }
    );

    auto points = gbs::make_vtkPoints(vertices_cloud);
    vtkNew<vtkCellArray> triangles;
    const auto t1 = std::chrono::high_resolution_clock::now();
    auto im = n-3;
    auto jm = n-2;
    auto km = n-1;
    // for(size_t i{}; i <= im; i++)
    for(size_t i{}; i <= km; i++)
    {
        auto a = vertices_cloud[i];
        for(size_t j{}; j <= km; j++)
        // for(size_t j{i+1}; j <= jm; j++)
        {
            auto b = vertices_cloud[j];
            // for( size_t k{j+1}; k <= km ; k++)
            for( size_t k{}; k <= km ; k++)
            {
                auto c   = vertices_cloud[k];
                // auto ccw = are_ccw(a, b, c);
                if(are_ccw(a, b, c))
                {
                 
                    auto is_delauney_triangle{true};
                    // for( size_t l{}; l <= km ; l++)
                    // {
                    //     if( l != i && l != j && l != k)
                    //     {
                    //         auto d = vertices_cloud[l];
                    //         if(ccw)
                    //             is_delauney_triangle = in_circle(a, b, c, d) < 0;
                    //         else
                    //             is_delauney_triangle = in_circle(b, a, c, d) < 0;
                    //     }
                    //     if(!is_delauney_triangle)
                    //     {
                    //         break;
                    //     }
                    // }

                    auto it_inside_vtx =  std::find_if(
                        std::execution::seq,
                        // std::execution::par,
                        vertices_cloud.begin(),
                        vertices_cloud.end(),
                        // [&a,&b,&c](const auto &d)
                        [=](const auto &d)
                        {
                            return in_circle(a, b, c, d) > 0;
                        }
                    );
                    is_delauney_triangle = vertices_cloud.end() == it_inside_vtx;

                    if(is_delauney_triangle)
                    {
                        vtkNew<vtkTriangle> triangle;
                        triangle->GetPointIds()->SetId(0, i);
                        triangle->GetPointIds()->SetId(1, j);
                        triangle->GetPointIds()->SetId(2, k);
                        triangles->InsertNextCell(triangle);
                    }
                   
                }
            }
        }
    }

    const std::chrono::duration<double, std::milli> ms_ref = std::chrono::high_resolution_clock::now() - t1;
    printf("Delauney computation time: %.3f ms.\n", ms_ref.count());

    vtkNew<vtkPolyData> trianglePolyData;

    // Add the geometry and topology to the polydata
    trianglePolyData->SetPoints(points);
    trianglePolyData->SetPolys(triangles);

    // Create mapper and actor
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputData(trianglePolyData);

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetEdgeVisibility(true);
    actor->GetProperty()->SetOpacity(0.1);

    gbs::plot(vertices_cloud, actor.Get());
}

TEST(tests_mesh_delauney2d, mesh_wire)
{
    using namespace gbs;
    using T = double;
    const size_t d{2};

    std::array<T, d> pt1{};
    std::array<T, d> pt2{1., 0.};
    std::array<T, d> pt3{1., 1.};
    std::array<T, d> pt4{0., 1.};

    Wire2d<T> w1{{pt1, pt2}};
    w1.addEdge({pt2, pt3});
    w1.addEdge({pt3, pt4});
    w1.addEdge({pt4, pt1});

    // auto ell = build_ellipse<T,d>(1.,0.5);
    // std::shared_ptr<Curve<T,d>> crv = std::make_shared<BSCurveRational<T,d>>(ell);
    // Edge<T,d> ed{crv};

    // Wire2d<T> w1(ed);

    std::vector< std::array<T,d> > vertices_cloud;

    auto mesh_edge = [dm=0.5](const auto & p_ed)
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
        [&vertices_cloud, &mesh_edge] (const auto & p_ed)
        {
            auto points = mesh_edge(p_ed);
            vertices_cloud.push_back( p_ed->vertex1()->point() );
            vertices_cloud.insert( 
                vertices_cloud.end(), 
                std::next(points.begin()),
                std::next(points.end(),-1)
            );
        }
    );

    size_t n = vertices_cloud.size();

    // std::random_device rd;  //Will be used to obtain a seed for the random number engine
    // std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    // std::uniform_real_distribution<double> distrib(0., 1.);

    // for(size_t i{}; i < n; i++)
    // {
    //     vertices_cloud.push_back(std::array<T,2>{distrib(gen),distrib(gen)});
    // }
    // n = vertices_cloud.size();

    auto points = gbs::make_vtkPoints(vertices_cloud);
    vtkNew<vtkCellArray> triangles;
    const auto t1 = std::chrono::high_resolution_clock::now();
    auto im = n-3;
    auto jm = n-2;
    auto km = n-1;
    // size_t i{1};
    for(size_t i{}; i <= im; i++)
    // for(size_t i{}; i <= km; i++)
    {
        auto a = vertices_cloud[i];
        std::vector<size_t> vtx;
        for (size_t j{i + 1}; j <= jm; j++)
        // for(size_t j{}; j <= km; j++)
        {
            auto b = vertices_cloud[j];
            size_t n_tri_per_edge{};
            for (size_t k{j + 1}; k <= km && n_tri_per_edge<2; k++)
            // for (size_t k{j + 1}; k <= km; k++)
            // for( size_t k{}; k <= km ; k++)
            {
                auto c = vertices_cloud[k];
                if (are_ccw(a, b, c))
                {
                    auto is_delauney_triangle{true};

                    bool contains_a_pnt{false};
                    for (size_t l{}; l <= km; l++)
                    {
                        if (l != i && l != j && l != k)
                        {
                            auto d = vertices_cloud[l];
                            contains_a_pnt = in_circle(a, b, c, d) > 0;
                        }
                        if (contains_a_pnt)
                        {
                            is_delauney_triangle = false;
                            break;
                        }
                    }

                    if (is_delauney_triangle && std::find(vtx.begin(),vtx.end(),j) == vtx.end())
                    {
                        vtkNew<vtkTriangle> triangle;
                        triangle->GetPointIds()->SetId(0, i);
                        triangle->GetPointIds()->SetId(1, j);
                        triangle->GetPointIds()->SetId(2, k);
                        triangles->InsertNextCell(triangle);
                        // break;
                        n_tri_per_edge++;
                        vtx.push_back(k);
                    }
                }
            }
        }
    }

    const std::chrono::duration<double, std::milli> ms_ref = std::chrono::high_resolution_clock::now() - t1;
    printf("Delauney computation time: %.3f ms.\n", ms_ref.count());

    vtkNew<vtkPolyData> trianglePolyData;

    // Add the geometry and topology to the polydata
    trianglePolyData->SetPoints(points);
    trianglePolyData->SetPolys(triangles);

    // Create mapper and actor
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputData(trianglePolyData);

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetEdgeVisibility(true);
    actor->GetProperty()->SetOpacity(0.3);

    gbs::plot(
        vertices_cloud
        , actor.Get()
        // , ell
        );

}


TEST(tests_mesh_delauney2d, mesh_wire_half_edge)
{
    using namespace gbs;
    using T = double;
    const size_t d{2};

    std::array<T, d> pt1{};
    std::array<T, d> pt2{1., 0.};
    std::array<T, d> pt3{1., 1.};
    std::array<T, d> pt4{0., 1.};

    Wire2d<T> w1{{pt1, pt2}};
    w1.addEdge({pt2, pt3});
    w1.addEdge({pt3, pt4});
    w1.addEdge({pt4, pt1});

    std::vector< std::array<T,d> > vertices_cloud;

    auto mesh_edge = [dm=0.5](const auto & p_ed)
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
        [&vertices_cloud, &mesh_edge] (const auto & p_ed)
        {
            auto points = mesh_edge(p_ed);
            vertices_cloud.push_back( p_ed->vertex1()->point() );
            vertices_cloud.insert( 
                vertices_cloud.end(), 
                std::next(points.begin()),
                std::next(points.end(),-1)
            );
        }
    );

    // HalfEdgeMesh<T,2> msh;

    // msh.addFreeVertices(vertices_cloud);

    // auto get_vtx_pnt = [&msh](size_t i)
    // {
    //     auto vtx = msh.vertices[i];
    //     auto pnt = vtx->point();
    //     return std::make_pair( vtx, pnt );
    // };

    // auto n = msh.vertices.size();
    // auto im = n-3;
    // auto jm = n-2;
    // auto km = n-1;
    // for(size_t i{}; i <= im; i++)
    // {
    //     auto [vtx_a, a] = get_vtx_pnt(i);
    //     for (size_t j{i + 1}; j <= jm; j++)
    //     {
    //         auto [vtx_b, b] = get_vtx_pnt(j);
    //         for (size_t k{j + 1}; k <= km; k++)
    //         {
    //             auto [vtx_d, c] = get_vtx_pnt(k);

    //             auto is_delauney_triangle{true};
    //             auto contains_a_pnt{false};
    //             auto ccw = are_ccw(a,b,c);
    //             for (size_t l{}; l <= km; l++)
    //             {
    //                 if (l != i && l != j && l != k)
    //                 {
    //                     auto d = vertices_cloud[l];
    //                     if(ccw)
    //                         contains_a_pnt = in_circle(a, b, c, d) > 0;
    //                     else
    //                         contains_a_pnt = in_circle(b, a, c, d) > 0;
    //                 }
    //                 if (contains_a_pnt)
    //                 {
    //                     is_delauney_triangle = false;
    //                     break;
    //                 }
    //             }
    //             if (is_delauney_triangle)
    //             {
    //                 if(ccw)
    //                     msh.makeFace(j,i,k);
    //                 else
    //                     msh.makeFace(i,j,k);
    //             }
    //         }
    //     }
    // }
}