#pragma once
#include "vtkcurvesrender.h"
#include "topology/tessellations.h"
#include "topology/halfEdgeMeshRender.h"
#include "topology/halfEdgeMeshQuality.h"

namespace gbs
{
/**
 * @brief Create a VTK actor to visualize a mesh given by points and triangles.
 *
 * @tparam T Floating-point type for the points.
 * @tparam dim Dimension of the points.
 * @param points_msh Vector of points that form the mesh.
 * @param pts_tri Vector of triples of indices of points that form the triangles.
 * @param col Pointer to an array of 3 doubles representing the RGB color of the mesh.
 * @return vtkSmartPointer<vtkActor> An actor representing the mesh.
 */
    template <std::floating_point T, size_t dim>
    auto make_actor(const points_vector<T, dim> &points_msh, const std::vector<std::array<vtkIdType, 3>> &pts_tri,double *col)
    {
        
        // auto pointActor = gbs::make_actor(pts,5.,true,colors->GetColor4d("Blue").GetData());

        vtkSmartPointer<vtkPolyData> ugrid =
            vtkSmartPointer<vtkPolyData>::New();
        ugrid->Allocate(pts_tri.size());

        std::for_each(pts_tri.begin(), pts_tri.end(),
                      [&ugrid](const auto &tri) { ugrid->InsertNextCell(VTK_TRIANGLE, 3, tri.data()); });

        ugrid->SetPoints( make_vtkPoints( points_msh ) );

        //TODO use surface's normals
        // Generate normals
        vtkSmartPointer<vtkPolyDataNormals> normalGenerator =
            vtkSmartPointer<vtkPolyDataNormals>::New();
        normalGenerator->SetInputData(ugrid);
        normalGenerator->ComputePointNormalsOn();
        normalGenerator->ComputeCellNormalsOff();
        normalGenerator->Update();

        /*
        // Optional settings
        normalGenerator->SetFeatureAngle(0.1);
        normalGenerator->SetSplitting(1);
        normalGenerator->SetConsistency(0);
        normalGenerator->SetAutoOrientNormals(0);
        normalGenerator->SetComputePointNormals(1);
        normalGenerator->SetComputeCellNormals(0);
        normalGenerator->SetFlipNormals(0);
        normalGenerator->SetNonManifoldTraversal(1);
        */

        ugrid = normalGenerator->GetOutput();

        vtkSmartPointer<vtkDataSetMapper> ugridMapper =
            vtkSmartPointer<vtkDataSetMapper>::New();
        ugridMapper->SetInputData(ugrid);

        vtkSmartPointer<vtkActor> ugridActor =
            vtkSmartPointer<vtkActor>::New();
        ugridActor->SetMapper(ugridMapper);
        ugridActor->GetProperty()->SetColor(col);
        // ugridActor->GetProperty()->EdgeVisibilityOn();
        ugridActor->GetProperty()->SetInterpolationToPhong();
        // ugridActor->GetProperty()->SetOpacity(0.7);

        return ugridActor;
    }
/**
 * @brief Create a VTK actor representing a surface mesh with control polygon
 *
 * @tparam T Floating point type
 * @tparam dim Dimension
 * @param points_msh Mesh points
 * @param pts_tri Triangular cells
 * @param poles Control polygon points
 * @param nPolesU Number of poles in the U direction
 * @return vtkSmartPointer<vtkAssembly> VTK assembly representing the surface mesh with control polygon
 */
    template <std::floating_point T, size_t dim>
    auto make_actor(const points_vector<T, dim> &points_msh, const std::vector<std::array<vtkIdType, 3>> &pts_tri, const points_vector<T, dim> &poles,size_t nPolesU) -> vtkSmartPointer<vtkAssembly>
    {
        auto colors = vtkSmartPointer<vtkNamedColors>::New();

        auto srf_actor = vtkSmartPointer<vtkAssembly>::New();

        double Peacock[3] = { 51./255.,  161./255.,  201./255.};
        double Red[3] = { 1.,  0.,  0.};
        double Black[3] = { 0.,  0.,  0.};

        auto srf_msh_actor = make_actor(points_msh, pts_tri, Peacock);
        srf_actor->AddPart(srf_msh_actor);

        auto ctrl_polygon = vtkSmartPointer<vtkAssembly>::New();
        auto polesActor = gbs::make_actor(poles, 20., true, Red);
        polesActor->GetProperty()->SetOpacity(0.3);

        auto ctr_polygon_lines = make_lattice_lines(poles,nPolesU,3.f,0.3,Black);

        ctrl_polygon->AddPart(ctr_polygon_lines);
        ctrl_polygon->AddPart(polesActor);

        srf_actor->AddPart(ctrl_polygon);

        return srf_actor;
    }

/**
 * Create a VTK actor for a surface.
 *
 * @tparam T The type of the surface control points.
 * @tparam dim The dimension of the surface (3 or 4).
 * @param srf The surface to create an actor for.
 * @param col The color of the surface.
 * @param n1 The number of points to use in the first direction when discretizing the surface.
 * @param n2 The number of points to use in the second direction when discretizing the surface.
 * @return A VTK actor representing the surface.
 */
    // template <std::floating_point T, size_t dim>
    // auto make_actor(const Surface<T, dim> &srf, std::array<double,3>  col = { 51./255.,  161./255.,  201./255.},size_t n1 = 5,size_t n2 = 5 )// -> vtkSmartPointer<vtkAssembly>
    // {
    //     auto pts = gbs::discretize(srf, n1, n2); //TODO: improve discretization

    //     std::vector<std::array<vtkIdType, 3>> pts_tri;

    //     vtkIdType nu = n1;
    //     vtkIdType nv = n2;

    //     std::array<vtkIdType, 3> tri;
    //     vtkIdType index;

    //     for (auto j = 0; j < nv - 1; j++)
    //     {
    //         for (auto i = 0; i < nu - 1; i++)
    //         {
    //             index = i + nu * j;
    //             pts_tri.push_back({index, index + 1, index + 1 + nu});
    //             pts_tri.push_back({index + 1 + nu, index + nu, index});
    //         }
    //     }

    //     return make_actor(pts, pts_tri, col.data() );
    // }

    template <std::floating_point T, size_t dim>
    auto make_actor(const Surface<T, dim> &srf, std::array<double,3>  col = { 51./255.,  161./255.,  201./255.},size_t n1 = 5,size_t n2 = 5 )// -> vtkSmartPointer<vtkAssembly>
    {
        auto [u1, u2, v1, v2] = srf.bounds();
        T lRef = 0.5 * distance(srf(u1, v1), srf(u2, v2) );
        T dev = 0.001;
        T dmax=0.005 * lRef;
        auto faces_lst = delaunay2DBoyerWatsonSurfaceMesh<T,dim,DistanceMeshSurface2<T,dim>>(srf, dev, 5000, n1, n2, 0.005, 1e-10);
        return surface_mesh_actor<T>(faces_lst, srf, { 51./255.,  161./255.,  201./255.}, false);
    }

    // template <std::floating_point T, size_t dim>
    // auto make_actor(const Surface<T, dim> &srf, std::array<double,3>  col = { 51./255.,  161./255.,  201./255.},size_t n1 = 5,size_t n2 = 5 )// -> vtkSmartPointer<vtkAssembly>
    // {
    //     auto [u1, u2, v1, v2] = srf.bounds();
    //     T lRef = 0.5 * distance(srf(u1, v1), srf(u2, v2) );
    //     T dev = 0.001;
    //     T dmax=0.005 * lRef;
    //     const auto p_srf = &srf;
    //     auto faces_lst = delaunay2DBoyerWatsonSurfaceMesh<T,dim,DistanceMeshSurface2<T,dim>>(p_srf, dev, 5000, n1, n2, 0.005, 1e-10);
    //     return surface_mesh_actor<T>(faces_lst, srf, { 51./255.,  161./255.,  201./255.}, true);
    // }

// /**
//  * Create a VTK actor to represent a B-Spline surface rational
//  * @tparam T The type of the coordinates of the surface (must be a floating-point type)
//  * @tparam dim The dimension of the surface (must be either 2 or 3)
//  * @param srf The B-Spline surface to be represented
//  * @return A VTK actor representing the surface
//  */
//     template <std::floating_point T, size_t dim>
//     auto make_actor(const BSSurfaceRational<T, dim> &srf) //-> vtkSmartPointer<vtkAssembly>
//     {
//         size_t n1 = fmin(100 * srf.nPolesU(),1000);
//         size_t n2 = fmin(100 * srf.nPolesV(),1000);
//         auto pts = gbs::discretize(srf,n1,n2); //TODO: improve discretization
//         auto poles = srf.polesProjected();
        
//         std::vector<std::array<vtkIdType ,3> > pts_tri;

//         vtkIdType nu = n1;
//         vtkIdType nv = n2;

//         std::array<vtkIdType ,3> tri;
//         vtkIdType index;

//         for (auto j = 0; j < nv - 1; j++)
//         {
//             for (auto i = 0; i < nu - 1; i++)
//             {
//                 index = i + nu * j;
//                 pts_tri.push_back({index, index + 1, index + 1 + nu });
//                 pts_tri.push_back({index + 1 + nu, index + nu, index});
//             }
//         }

//         // auto srf_actor =  make_actor(pts,pts_tri,poles,srf.nPolesU());
//         // return srf_actor;

//         auto colors = vtkSmartPointer<vtkNamedColors>::New();
//         return make_actor(pts, pts_tri, colors->GetColor3d("Peacock").GetData());
//     }
}