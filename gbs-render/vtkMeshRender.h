#pragma once
/**
 * @file vtkMeshRender.h
 * @brief VTK rendering of a triangle mesh of the native BREP (`gbs-brep/mesh.h`, stage 3).
 *
 * Each corner is its own VTK point carrying the normal of the surface there, so that the
 * shading follows the surfaces and the edges between faces stay sharp.
 */
#include <gbs-brep/mesh.h>

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkSmartPointer.h>

#include <array>

namespace gbs
{
    /// The polydata of a mesh: three points per triangle, with the normals of the corners.
    template <std::floating_point T>
    auto make_polydata(const brep::Mesh<T> &mesh) -> vtkSmartPointer<vtkPolyData>
    {
        auto points = vtkSmartPointer<vtkPoints>::New();
        points->SetDataTypeToDouble();
        points->SetNumberOfPoints(static_cast<vtkIdType>(3 * mesh.triangles.size()));
        auto normals = vtkSmartPointer<vtkDoubleArray>::New();
        normals->SetName("Normals");
        normals->SetNumberOfComponents(3);
        normals->SetNumberOfTuples(static_cast<vtkIdType>(3 * mesh.triangles.size()));
        auto cells = vtkSmartPointer<vtkCellArray>::New();
        cells->AllocateExact(static_cast<vtkIdType>(mesh.triangles.size()), static_cast<vtkIdType>(3 * mesh.triangles.size()));
        for (std::size_t k = 0; k < mesh.triangles.size(); ++k)
        {
            std::array<vtkIdType, 3> ids;
            for (std::size_t c = 0; c < 3; ++c)
            {
                const auto i = static_cast<vtkIdType>(3 * k + c);
                const auto &p = mesh.vertices[mesh.triangles[k][c]];
                const auto &n = mesh.normals[3 * k + c];
                points->SetPoint(i, double(p[0]), double(p[1]), double(p[2]));
                normals->SetTuple3(i, double(n[0]), double(n[1]), double(n[2]));
                ids[c] = i;
            }
            cells->InsertNextCell(3, ids.data());
        }
        auto poly = vtkSmartPointer<vtkPolyData>::New();
        poly->SetPoints(points);
        poly->SetPolys(cells);
        poly->GetPointData()->SetNormals(normals);
        return poly;
    }

    /// An actor of a mesh, of the color col (RGB in [0, 1]).
    template <std::floating_point T>
    auto make_actor(const brep::Mesh<T> &mesh, const std::array<double, 3> &col) -> vtkSmartPointer<vtkActor>
    {
        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(make_polydata(mesh));
        auto actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(col[0], col[1], col[2]);
        actor->GetProperty()->SetInterpolationToPhong();
        return actor;
    }
}
