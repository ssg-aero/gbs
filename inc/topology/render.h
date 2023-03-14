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

namespace gbs
{
    template <typename T, size_t dim>
    auto make_polydata_from_faces(const auto &faces_lst)
    {
        // valid only up to 3D
        static_assert(dim < 4);

        // // build vertices map and store vtk points
        vtkNew<vtkPoints> points;
        auto vertices_map = extract_vertices_map_from_faces<T, dim>(faces_lst);
        points->Allocate(vertices_map.size());
        points->SetNumberOfPoints(vertices_map.size());
        for (const auto &vtx : vertices_map)
        {
            points->SetPoint(vtx.second, make_vtkPoint(vtx.first->coords).data());
        }

        // Store cells
        vtkNew<vtkCellArray> cells;

        for (const auto &f : faces_lst)
        {
            auto vtx_lst = getFaceVertices(f);
            vtkSmartPointer<vtkCell> cell;
            auto n = vtx_lst.size();
            switch (n)
            {
            case 3:
            {
                cell = vtkSmartPointer<vtkTriangle>::New();
            }
            break;
            case 4:
            {
                cell = vtkSmartPointer<vtkQuad>::New();
            }
            default:
                cell = vtkSmartPointer<vtkPolygon>::New();
                cell->GetPointIds()->SetNumberOfIds(n);
                break;
            }
            std::transform(
                vtx_lst.begin(), vtx_lst.end(),
                cell->GetPointIds()->begin(),
                [&vertices_map](const auto vtx)
                { return vertices_map[vtx]; });
            cells->InsertNextCell(cell);
        }

        vtkNew<vtkPolyData> polyData;

        // Add the geometry and topology to the polydata
        polyData->SetPoints(points);
        polyData->SetPolys(cells);

        return polyData;
    }

    template <typename T, size_t dim>
    auto make_polydata_from_faces(const auto &faces_lst, const Surface<T, dim> &srf)
    {
        // valid only up to 3D
        static_assert(dim < 4);

        // // build vertices map and store vtk points
        vtkNew<vtkPoints> points;
        auto vertices_map = extract_vertices_map_from_faces<T, 2>(faces_lst);
        points->Allocate(vertices_map.size());
        points->SetNumberOfPoints(vertices_map.size());
        for (const auto &vtx : vertices_map)
        {
            auto [u, v] = vtx.first->coords;
            points->SetPoint(vtx.second, make_vtkPoint(srf(u, v)).data());
        }

        // Store cells
        vtkNew<vtkCellArray> cells;

        for (const auto &f : faces_lst)
        {
            auto vtx_lst = getFaceVertices(f);
            vtkSmartPointer<vtkCell> cell;
            auto n = vtx_lst.size();
            switch (n)
            {
            case 3:
            {
                cell = vtkSmartPointer<vtkTriangle>::New();
            }
            break;
            case 4:
            {
                cell = vtkSmartPointer<vtkQuad>::New();
            }
            default:
                cell = vtkSmartPointer<vtkPolygon>::New();
                cell->GetPointIds()->SetNumberOfIds(n);
                break;
            }
            std::transform(
                vtx_lst.begin(), vtx_lst.end(),
                cell->GetPointIds()->begin(),
                [&vertices_map](const auto vtx)
                { return vertices_map[vtx]; });
            cells->InsertNextCell(cell);
        }

        vtkNew<vtkPolyData> polyData;

        // Add the geometry and topology to the polydata
        polyData->SetPoints(points);
        polyData->SetPolys(cells);

        vtkSmartPointer<vtkPolyDataNormals> normalGenerator =
            vtkSmartPointer<vtkPolyDataNormals>::New();
        normalGenerator->SetInputData(polyData);
        normalGenerator->ComputePointNormalsOn();
        normalGenerator->ComputeCellNormalsOff();
        normalGenerator->Update();

        polyData->GetPointData()->SetNormals(normalGenerator->GetOutput()->GetPointData()->GetNormals());

        return polyData;
    }

    template <typename T>
    auto make_polydata_from_faces(const auto &faces_lst, const Surface<T, 3> &srf)
    {
        // valid only up to 3D
        // static_assert(dim<4);

        // // build vertices map and store vtk points
        auto vertices_map = extract_vertices_map_from_faces<T, 2>(faces_lst);

        vtkNew<vtkPoints> points;
        points->Allocate(vertices_map.size());
        points->SetNumberOfPoints(vertices_map.size());

        vtkNew<vtkFloatArray> normals;
        normals->SetNumberOfComponents(3);
        normals->SetNumberOfTuples(vertices_map.size());

        for (const auto &vtx : vertices_map)
        {
            auto [u, v] = vtx.first->coords;
            points->SetPoint(vtx.second, make_vtkPoint(srf(u, v)).data());
            auto N = make_vtkPoint(normalized(srf(u, v, 1, 0) ^ srf(u, v, 0, 1)));
            normals->SetTuple3(vtx.second, N[0], N[1], N[2]);
        }

        // Store cells
        vtkNew<vtkCellArray> cells;

        for (const auto &f : faces_lst)
        {
            auto vtx_lst = getFaceVertices(f);
            vtkSmartPointer<vtkCell> cell;
            auto n = vtx_lst.size();
            switch (n)
            {
            case 3:
            {
                cell = vtkSmartPointer<vtkTriangle>::New();
            }
            break;
            case 4:
            {
                cell = vtkSmartPointer<vtkQuad>::New();
            }
            default:
                cell = vtkSmartPointer<vtkPolygon>::New();
                cell->GetPointIds()->SetNumberOfIds(n);
                break;
            }
            std::transform(
                vtx_lst.begin(), vtx_lst.end(),
                cell->GetPointIds()->begin(),
                [&vertices_map](const auto vtx)
                { return vertices_map[vtx]; });
            cells->InsertNextCell(cell);
        }

        vtkNew<vtkPolyData> polyData;

        // Add the geometry and topology to the polydata
        polyData->SetPoints(points);
        polyData->SetPolys(cells);
        polyData->GetPointData()->SetNormals(normals);
        std::cout << "using surface normal\n";
        return polyData;
    }

    template <typename T, size_t dim>
    auto make_polydata_from_edges_loop(const auto &edges_lst)
    {
        // valid only up to 3D
        static_assert(dim < 4);

        // build vertices map
        std::map<std::shared_ptr<HalfEdgeVertex<T, dim>>, size_t> vertices_map;
        size_t index{};
        for (const auto &e : edges_lst)
        {
            auto vtx = e->vertex;
            if (!vertices_map.contains(vtx))
            {
                vertices_map[vtx] = index;
                index++;
            }
        }
        // store vtk points
        vtkNew<vtkPoints> points;
        points->Allocate(vertices_map.size());
        points->SetNumberOfPoints(vertices_map.size());
        for (const auto &vtx : vertices_map)
        {
            points->SetPoint(vtx.second, make_vtkPoint(vtx.first->coords).data());
        }
        // make polyline
        vtkNew<vtkPolyLine> polyLine;
        auto n = edges_lst.size();
        polyLine->GetPointIds()->SetNumberOfIds(n + 1);
        size_t i{};
        if (edges_lst.front()->previous)
        {
            polyLine->GetPointIds()->SetId(i, vertices_map[edges_lst.front()->previous->vertex]);
            i++;
        }
        for (const auto &e : edges_lst)
        {
            polyLine->GetPointIds()->SetId(i, vertices_map[e->vertex]);
            i++;
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

    template <typename T, size_t dim>
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

    template <typename T, size_t dim>
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

    template <typename T>
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

    template <typename T, size_t dim>
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