#include <gtest/gtest.h>

#include <gbslib/bscanalysis.h>
#include <gbslib/bscbuild.h>
#include <gbslib/bscurve.h>

#include <render/vtkcurvesrender.h>

#include <vtkActor.h>
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


TEST(tests_vtk_render, BSC)
{

        auto c = gbs::build_circle<double,3>(1.,{0.,0.,0.});
        c.trim(0.3,0.7);

  vtkSmartPointer<vtkNamedColors> colors =
    vtkSmartPointer<vtkNamedColors>::New();


    auto actor_crv = gbs::make_BSC_actor(c);
    

  // Setup render window, renderer, and interactor
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->SetWindowName("PolyLine");
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderer->AddActor(actor_crv);

    renderer->SetBackground(colors->GetColor3d("White").GetData());

  renderWindow->Render();
  renderWindowInteractor->Start();
}