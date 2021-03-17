#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
namespace gbs
{
    template <typename T, size_t dim>
    auto make_structuredgrid_actor(const gbs::points_vector<T, dim> &pts, size_t ni, size_t nj) -> vtkSmartPointer<vtkActor>
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

        // Create a mapper and actor
        vtkSmartPointer<vtkDataSetMapper> mapper =
            vtkSmartPointer<vtkDataSetMapper>::New();
        mapper->SetInputData(structuredGrid);

        vtkSmartPointer<vtkActor> actor =
            vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->EdgeVisibilityOn();

        return actor;
    }
}