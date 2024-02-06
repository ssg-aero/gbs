#pragma once

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkQuad.h>
#include <vtkPolygon.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>
#include <vtkPolyDataMapper.h>

#ifdef GBS_USE_MODULES
    import vecop;
#else
    #include <gbs/vecop.ixx>
#endif

namespace gbs
{
/**
 * @brief Generate VTK points from a vertices map and a surface
 *
 * @tparam T : floating point type
 * @tparam dim : dimension of the surface
 * @param vertices_map : map containing vertices of the faces
 * @param srf : Surface object for generating point coordinates
 * @return vtkSmartPointer<vtkPoints>
 * @requires std::floating_point<T> and (dim < 4)
 */
    template <std::floating_point T, size_t dim>
    requires (dim < 4)
    auto generate_points_from_vertices_map(const auto &vertices_map, const Surface<T, dim> &srf)
    {
        vtkNew<vtkPoints> points;
        points->Allocate(vertices_map.size());
        points->SetNumberOfPoints(vertices_map.size());
        for (const auto &[key, value] : vertices_map)
        {
            auto [u, v] = key->coords;
            points->SetPoint(value, make_vtkPoint(srf(u, v)).data());
        }
        return points;
    }

    auto generate_points_from_vertices_map(const auto &vertices_map)
    {
        vtkNew<vtkPoints> points;
        points->Allocate(vertices_map.size());
        points->SetNumberOfPoints(vertices_map.size());
        for (const auto &[key, value] : vertices_map)
        {
            points->SetPoint(value, make_vtkPoint(key->coords).data());
        }
        return points;
    }

/**
 * @brief Generate points and normals from faces and surface
 *
 * @param vertices_map : map containing vertices of the faces
 * @param srf : surface used to compute points and normals
 * @return std::pair<vtkSmartPointer<vtkPoints>, vtkSmartPointer<vtkFloatArray>>
 */
    template <std::floating_point T>
    auto generate_points_and_normals(const auto &vertices_map, const Surface<T, 3> &srf)
    {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        points->Allocate(vertices_map.size());
        points->SetNumberOfPoints(vertices_map.size());

        vtkSmartPointer<vtkFloatArray> normals = vtkSmartPointer<vtkFloatArray>::New();
        normals->SetNumberOfComponents(3);
        normals->SetNumberOfTuples(vertices_map.size());

        for (const auto &[vtx_ptr, index] : vertices_map)
        {
            auto [u, v] = vtx_ptr->coords;
            points->SetPoint(index, make_vtkPoint(srf(u, v)).data());
            auto N = make_vtkPoint(normalized(cross(srf(u, v, 1, 0), srf(u, v, 0, 1) ) ));
            normals->SetTuple3(index, N[0], N[1], N[2]);
        }

        return std::make_pair(points, normals);
    }

/**
 * @brief Generate cells from a list of faces and a vertices map
 *
 * @param faces_lst : list of faces to be added to the cells
 * @param vertices_map : map containing vertices of the faces
 * @return vtkSmartPointer<vtkCellArray>
 */
    auto generate_cells_from_faces(const auto &faces_lst, const auto &vertices_map)
    {
        vtkNew<vtkCellArray> cells;

        for (const auto &f : faces_lst)
        {
            auto vtx_lst = getFaceVertices(f);
            vtkSmartPointer<vtkCell> cell;
            auto n = vtx_lst.size();
            switch (n)
            {
            case 3:
                cell = vtkSmartPointer<vtkTriangle>::New();
                break;
            case 4:
                cell = vtkSmartPointer<vtkQuad>::New();
                break;
            default:
                cell = vtkSmartPointer<vtkPolygon>::New();
                cell->GetPointIds()->SetNumberOfIds(n);
                break;
            }
            std::ranges::transform(
                vtx_lst, cell->GetPointIds()->begin(),
                [&vertices_map](const auto vtx)
                { return vertices_map.at(vtx); });
            cells->InsertNextCell(cell);
        }

        return cells;
    }

/**
 * @brief Generates a vtkPolyData object using input points and cells.
 *
 * @param points A vtkSmartPointer containing the points of the geometry.
 * @param cells A vtkSmartPointer containing the cells of the topology.
 * @param compute_normals A boolean flag to determine whether to compute normals for the polydata.
 * @return vtkSmartPointer<vtkPolyData> The resulting vtkPolyData object.
 */
    inline auto generate_polydata(const vtkSmartPointer<vtkPoints> &points, const vtkSmartPointer<vtkCellArray> &cells, bool compute_normals=false)
    {
        vtkNew<vtkPolyData> polyData;

        // Add the geometry and topology to the polydata
        polyData->SetPoints(points);
        polyData->SetPolys(cells);
        if(compute_normals)
        { // Build normals from mesh
            vtkSmartPointer<vtkPolyDataNormals> normalGenerator =
                vtkSmartPointer<vtkPolyDataNormals>::New();
            normalGenerator->SetInputData(polyData);
            normalGenerator->ComputePointNormalsOn();
            normalGenerator->ComputeCellNormalsOff();
            normalGenerator->Update();

            polyData->GetPointData()->SetNormals(normalGenerator->GetOutput()->GetPointData()->GetNormals());
        }
        return polyData;
    }

/**
 * @brief Generate polyData from points, normals and cells
 *
 * @param points : vtkPoints to be added to polyData
 * @param normals : vtkFloatArray normals to be added to polyData
 * @param cells : vtkCellArray cells to be added to polyData
 * @return vtkSmartPointer<vtkPolyData>
 */
    inline auto generate_polydata(const vtkSmartPointer<vtkPoints> &points, const vtkSmartPointer<vtkFloatArray> &normals, const vtkSmartPointer<vtkCellArray> &cells)
    {
        vtkNew<vtkPolyData> polyData;

        // Add the geometry and topology to the polydata
        polyData->SetPoints(points);
        polyData->SetPolys(cells);
        polyData->GetPointData()->SetNormals(normals);

        return polyData;
    }

/**
 * @brief Create a vtkPolyData object from a list of faces.
 * 
 * @tparam T Floating point type
 * @tparam dim Dimension (up to 3D)
 * @param faces_lst List of faces to convert into vtkPolyData
 * @return vtkSmartPointer<vtkPolyData> Generated vtkPolyData object
 */
    template <std::floating_point T, size_t dim>
    requires (dim < 4) // Valid only up to 3D
    auto make_polydata_from_faces(const auto &faces_lst)
    {

        auto vertices_map = getVerticesMapFromFaces<T, 2>(faces_lst);

        auto points = generate_points_from_vertices_map(vertices_map);

        auto cells = generate_cells_from_faces(faces_lst, vertices_map);

        vtkNew<vtkPolyData> polyData;
        // Add the geometry and topology to the polydata
        polyData->SetPoints(points);
        polyData->SetPolys(cells);

        return polyData;
    }

/**
 * @brief Create a vtkPolyData object from a list of faces and a surface.
 * 
 * @tparam T Floating point type
 * @tparam dim Dimension (up to 3D)
 * @param faces_lst List of faces to convert into vtkPolyData
 * @param srf Surface object to map the UV coordinates to 3D points
 * @return vtkSmartPointer<vtkPolyData> Generated vtkPolyData object
 */
    template <std::floating_point T, size_t dim>
    requires (dim < 4) // Valid only up to 3D
    auto make_polydata_from_faces(const auto &faces_lst, const Surface<T, dim> &srf)
    {

        auto vertices_map = getVerticesMapFromFaces<T, 2>(faces_lst);

        auto points = generate_points_from_vertices_map(vertices_map, srf);

        auto cells = generate_cells_from_faces(faces_lst, vertices_map);

        return generate_polydata(points, cells, true);
    }

/**
 * @brief Create polyData from face list and surface
 *
 * @param faces_lst : list of faces to be added to polyData
 * @param srf : surface used to compute points and normals
 * @return vtkSmartPointer<vtkPolyData>
 */
    template <std::floating_point T>
    auto make_polydata_from_faces(const auto &faces_lst, const Surface<T, 3> &srf)
    {
        auto vertices_map = getVerticesMapFromFaces<T, 2>(faces_lst);

        auto [points, normals] = generate_points_and_normals<T>(vertices_map, srf);

        auto cells = generate_cells_from_faces(faces_lst, vertices_map);

        return generate_polydata(points, normals, cells);
    }

/**
 * @brief Creates a vtkPolyData object from a list of edges.
 *
 * @tparam T Floating-point type for the input coordinates.
 * @tparam dim Dimension of the input data. Must be less than 4.
 * @param edges_lst A list of half-edge data structures representing the input edges.
 * @return vtkSmartPointer<vtkPolyData> The resulting vtkPolyData object.
 */
    template <std::floating_point T, size_t dim>
    requires (dim < 4)
    auto make_polydata_from_edges_loop(const auto &edges_lst)
    {
        // build vertices map
        std::map<std::shared_ptr<HalfEdgeVertex<T, dim>>, size_t> vertices_map;
        size_t index{};
        for (const auto &e : edges_lst)
        {
            auto vtx = e->vertex;
            if (vertices_map.find(vtx) == vertices_map.end())
            {
                vertices_map[vtx] = index++;
            }
        }

        // store vtk points
        auto points = generate_points_from_vertices_map(vertices_map);

        // make polyline
        vtkNew<vtkPolyLine> polyLine;
        auto n = edges_lst.size();
        polyLine->GetPointIds()->SetNumberOfIds(n + 1);
        size_t i{};
        if (edges_lst.front()->previous)
        {
            polyLine->GetPointIds()->SetId(i++, vertices_map[edges_lst.front()->previous->vertex]);
        }
        for (const auto &e : edges_lst)
        {
            polyLine->GetPointIds()->SetId(i++, vertices_map[e->vertex]);
        }

        // Store cells
        vtkNew<vtkCellArray> cells;
        cells->InsertNextCell(polyLine);

        vtkNew<vtkPolyData> polyData;
        // Add the geometry and topology to the polydata
        polyData->SetPoints(points);
        polyData->SetLines(cells);

        return polyData;
    }

/**
 * @brief Creates a VTK actor for visualizing a mesh represented by a list of faces.
 * @tparam T Floating point type (e.g. float, double)
 * @tparam dim Dimension of the mesh (up to 3D)
 * @param faces_lst List of faces representing the mesh
 * @param color Color of the mesh (default: light gray)
 * @param edges_on Whether to show edges (default: true)
 * @param opacity Opacity of the mesh (default: 1.0)
 * @return VTK actor for rendering the mesh
 */
    template <std::floating_point T, size_t dim>
    auto faces_mesh_actor(const auto &faces_lst, std::array<double, 3> color = {0.9, 0.9, 0.9}, bool edges_on = true, double opacity = 1.)
    {
        auto polyData = make_polydata_from_faces<T, dim>(faces_lst);
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputData(polyData);
        vtkNew<vtkActor> actor;
        actor->SetMapper(mapper);
        actor->GetProperty()->SetEdgeVisibility(edges_on);
        actor->GetProperty()->SetOpacity(opacity);
        actor->GetProperty()->SetColor(color.data());
        return actor;
    }

/**
 * @brief Creates a VTK actor for visualizing a parametric surface mesh represented by a list of faces.
 * @tparam T Floating point type (e.g. float, double)
 * @tparam dim Dimension of the mesh (up to 3D)
 * @param faces_lst List of faces representing the mesh
 * @param srf Parametric surface
 * @param color Color of the mesh (default: light gray)
 * @param edges_on Whether to show edges (default: true)
 * @param opacity Opacity of the mesh (default: 1.0)
 * @return VTK actor for rendering the mesh
 */
    template <std::floating_point T, size_t dim>
    auto surface_mesh_actor(const auto &faces_lst, const Surface<T, dim> &srf, std::array<double, 3> color = {0.9, 0.9, 0.9}, bool edges_on = true, double opacity = 1.)
    {
        auto polyData = make_polydata_from_faces<T, dim>(faces_lst, srf);
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputData(polyData);
        vtkNew<vtkActor> actor;
        actor->SetMapper(mapper);
        actor->GetProperty()->SetEdgeVisibility(edges_on);
        actor->GetProperty()->SetOpacity(opacity);
        actor->GetProperty()->SetColor(color.data());
        return actor;
    }

/**
 * @brief Creates a VTK actor for visualizing a 3D parametric surface mesh represented by a list of faces.
 * @tparam T Floating point type (e.g. float, double)
 * @param faces_lst List of faces representing the mesh
 * @param srf 3D parametric surface
 * @param color Color of the mesh (default: light gray)
 * @param edges_on Whether to show edges (default: true)
 * @param opacity Opacity of the mesh (default: 1.0)
 * @return VTK actor for rendering the mesh
 */
    template <std::floating_point T>
    auto surface_mesh_actor(const auto &faces_lst, const Surface<T, 3> &srf, std::array<double, 3> color = {0.9, 0.9, 0.9}, bool edges_on = true, double opacity = 1.)
    {
        auto polyData = make_polydata_from_faces<T>(faces_lst, srf);
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputData(polyData);
        vtkNew<vtkActor> actor;
        actor->SetMapper(mapper);
        actor->GetProperty()->SetEdgeVisibility(edges_on);
        actor->GetProperty()->SetOpacity(opacity);
        actor->GetProperty()->SetColor(color.data());
        return actor;
    }

/**
 * @brief Creates a VTK actor for visualizing a boundary mesh represented by an edge loop.
 * @tparam T Floating point type (e.g. float, double)
 * @tparam dim Dimension of the mesh (up to 3D)
 * @param boundary Edge loop representing the boundary
 * @param color Color of the boundary (default: red)
 * @param opacity Opacity of the boundary (default: 1.0)
 * @return VTK actor for rendering the boundary
 */
    template <std::floating_point T, size_t dim>
    auto boundary_mesh_actor(const auto &boundary, std::array<double, 3> color = {1., 0., 0.}, double opacity = 1.)
    {
        auto polyData_boundary = make_polydata_from_edges_loop<T, dim>(
            boundary);
        vtkNew<vtkPolyDataMapper> mapper_boundary;
        mapper_boundary->SetInputData(polyData_boundary);

        vtkNew<vtkActor> actor_boundary;
        actor_boundary->SetMapper(mapper_boundary);
        actor_boundary->GetProperty()->SetOpacity(opacity);
        actor_boundary->GetProperty()->SetColor(color.data());
        return actor_boundary;
    }
}