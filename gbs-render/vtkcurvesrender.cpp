#include <gbs-render/vtkGbsRender.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkNamedColors.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkTexture.h>
#include <vtkNew.h>
#include <vtkSphereSource.h>
#include <vtkCollectionIterator.h>
// #include <vtkOpenGLPolyDataMapper.h>
#include <vtkPolyDataMapper.h>


namespace gbs
{
    auto make_polyline_(vtkPoints *points,std::array<double,3> &color) -> vtkSmartPointer<vtkActor>
    {
        auto n = points->GetNumberOfPoints();
        vtkSmartPointer<vtkPolyLine> polyLine =
            vtkSmartPointer<vtkPolyLine>::New();
        polyLine->GetPointIds()->SetNumberOfIds(n);
        for (unsigned int i = 0; i < n; i++)
        {
            polyLine->GetPointIds()->SetId(i, i);
        }

        // Create a cell array to store the lines in and add the lines to it
        vtkSmartPointer<vtkCellArray> cells =
            vtkSmartPointer<vtkCellArray>::New();
        cells->InsertNextCell(polyLine);

        // Create a polydata to store everything in
        vtkSmartPointer<vtkPolyData> polyData =
            vtkSmartPointer<vtkPolyData>::New();

        // Add the points to the dataset
        polyData->SetPoints(points);

        // Add the lines to the dataset
        polyData->SetLines(cells);

        // Setup actor and mapper
        // vtkSmartPointer<vtkOpenGLPolyDataMapper> mapper =
        //     vtkSmartPointer<vtkOpenGLPolyDataMapper>::New();
        vtkSmartPointer<vtkPolyDataMapper> mapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();


        mapper->SetInputData(polyData);

        vtkSmartPointer<vtkActor> actor =
            vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);

        actor->GetProperty()->SetColor(color.data());


        return actor;
    }

    auto StippledLine(vtkSmartPointer<vtkActor> &actor,
                      int lineStipplePattern,
                      int lineStippleRepeat ) -> void
    {
        vtkSmartPointer<vtkDoubleArray> tcoords =
            vtkSmartPointer<vtkDoubleArray>::New();
        vtkSmartPointer<vtkImageData> image =
            vtkSmartPointer<vtkImageData>::New();
        vtkSmartPointer<vtkTexture> texture =
            vtkSmartPointer<vtkTexture>::New();

        // Create texture
        int dimension = 16 * lineStippleRepeat;

        image->SetDimensions(dimension, 1, 1);
        image->AllocateScalars(VTK_UNSIGNED_CHAR, 4);
        image->SetExtent(0, dimension - 1, 0, 0, 0, 0);
        unsigned char *pixel;
        pixel = static_cast<unsigned char *>(image->GetScalarPointer());
        unsigned char on = 255;
        unsigned char off = 0;
        for (int i = 0; i < 16; ++i)
        {
            unsigned int mask = (1 << i);
            unsigned int bit = (lineStipplePattern & mask) >> i;
            unsigned char value = static_cast<unsigned char>(bit);
            if (value == 0)
            {
                for (int j = 0; j < lineStippleRepeat; ++j)
                {
                    *pixel = on;
                    *(pixel + 1) = on;
                    *(pixel + 2) = on;
                    *(pixel + 3) = off;
                    pixel += 4;
                }
            }
            else
            {
                for (int j = 0; j < lineStippleRepeat; ++j)
                {
                    *pixel = on;
                    *(pixel + 1) = on;
                    *(pixel + 2) = on;
                    *(pixel + 3) = on;
                    pixel += 4;
                }
            }
        }
        vtkPolyData *polyData = dynamic_cast<vtkPolyData *>(actor->GetMapper()->GetInput());

        // Create texture coordinates
        tcoords->SetNumberOfComponents(1);
        tcoords->SetNumberOfTuples(polyData->GetNumberOfPoints());
        for (int i = 0; i < polyData->GetNumberOfPoints(); ++i)
        {
            double value[1] = { static_cast<double>(i) * .5 };
            tcoords->SetTypedTuple(i, value);
        }

        polyData->GetPointData()->SetTCoords(tcoords);
        texture->SetInputData(image);
        texture->InterpolateOff();
        texture->RepeatOn();

        actor->SetTexture(texture);
    }

    auto scale_parts(double s,vtkAssembly *a) ->void
    {
        auto col = vtkPropCollection::New();
        a->GetActors(col);

        auto it = col->NewIterator();
        for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextItem())
        {
            auto actor_ = vtkActor::SafeDownCast(it->GetCurrentObject());
            if (actor_)
                actor_->SetScale(0.05);
        }

        col->Delete();
    }
} // namespace gbs