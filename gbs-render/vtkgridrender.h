#pragma once
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkActor.h>
#include <gbs-render/vtkGbsRender.h>
#include <vtkPointData.h>
namespace gbs
{

    template <typename T>
    auto add_value(vtkStructuredGrid *structuredGrid, const char *name, T default_value)
    {
        vtkSmartPointer<vtkDoubleArray> doubleArray =
            vtkSmartPointer<vtkDoubleArray>::New();

        doubleArray->SetName(name);
        auto dims = structuredGrid->GetDimensions();
        vtkIdType n = dims[0] * dims[1] * dims[2];
        doubleArray->Allocate(n);
        doubleArray->FillValue(default_value);

        structuredGrid->GetPointData()->AddArray(doubleArray);
    }

    template <typename T, size_t dim>
    auto make_structuredgrid(const points_vector<T, dim> &pts, size_t ni, size_t nj) -> vtkSmartPointer<vtkStructuredGrid>
    {
        // Create a grid
        vtkSmartPointer<vtkStructuredGrid> structuredGrid =
            vtkSmartPointer<vtkStructuredGrid>::New();

        vtkSmartPointer<vtkPoints> points =
            vtkSmartPointer<vtkPoints>::New();

        for (size_t j = 0; j < nj; j++)
        {
            for (size_t i = 0; i < ni; i++)
            {
                points->InsertNextPoint(make_vtkPoint<T, dim>(pts[i + ni * j]).data());
            }
        }

        // Specify the dimensions of the grid
        structuredGrid->SetDimensions(ni, nj, 1);
        structuredGrid->SetPoints(points);

        return structuredGrid;
    }

    inline auto make_structuredgrid_actor(const vtkSmartPointer<vtkStructuredGrid> &structuredGrid) -> vtkSmartPointer<vtkActor>
    // inline auto make_structuredgrid_actor(vtkStructuredGrid *structuredGrid) -> vtkSmartPointer<vtkActor>
    {
        vtkSmartPointer<vtkDataSetMapper> mapper =
            vtkSmartPointer<vtkDataSetMapper>::New();
        mapper->SetInputData(structuredGrid);

        vtkSmartPointer<vtkActor> actor =
            vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->EdgeVisibilityOn();

        return actor;
    }

    template <typename T, size_t dim>
    auto make_structuredgrid_actor(const points_vector<T, dim> &pts, size_t ni, size_t nj) -> vtkSmartPointer<vtkActor>
    {

        auto structuredGrid = make_structuredgrid(pts,ni,nj);
        return make_structuredgrid_actor(structuredGrid);

    }
}