#include <gtest/gtest.h>



#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkInteractorStyleTrackballActor.h>
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
#include <vtkPointSource.h>
#include <vtkVertexGlyphFilter.h>

#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkInteractorStyleTrackballCamera.h>

TEST(tests_vtkexamples, interctor)
{
  // Sphere 1
  vtkSmartPointer<vtkSphereSource> sphereSource1 = 
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource1->SetCenter(0.0, 0.0, 0.0);
  sphereSource1->SetRadius(4.0);
  sphereSource1->Update();
    
  vtkSmartPointer<vtkPolyDataMapper> mapper1 = 
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper1->SetInputConnection(sphereSource1->GetOutputPort());
  
  vtkSmartPointer<vtkActor> actor1 = 
    vtkSmartPointer<vtkActor>::New();
  actor1->SetMapper(mapper1);
  
  // Sphere 2
  vtkSmartPointer<vtkSphereSource> sphereSource2 = 
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource2->SetCenter(10.0, 0.0, 0.0);
  sphereSource2->SetRadius(3.0);
  sphereSource2->Update();

  // Create a mapper
  vtkSmartPointer<vtkPolyDataMapper> mapper2 = 
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper2->SetInputConnection(sphereSource2->GetOutputPort());

  // Create an actor
  vtkSmartPointer<vtkActor> actor2 = 
    vtkSmartPointer<vtkActor>::New();
  actor2->SetMapper(mapper2);

  // A renderer and render window
  vtkSmartPointer<vtkRenderer> renderer = 
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = 
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  // An interactor
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = 
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  // Add the actors to the scene
  renderer->AddActor(actor1);
  renderer->AddActor(actor2);
  renderer->SetBackground(1,1,1); // Background color white

  // Render an image (lights and cameras are created automatically)
  renderWindow->Render();

  vtkSmartPointer<vtkInteractorStyleTrackballActor> actorStyle =
      vtkSmartPointer<vtkInteractorStyleTrackballActor>::New();
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> cameraStyle =
      vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); //like paraview

  renderWindowInteractor->SetInteractorStyle( actorStyle );
  
  // Begin mouse interaction
  renderWindowInteractor->Start();
  
  
}

#include <vtkActor.h>
#include <vtkAreaPicker.h>
#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkExtractGeometry.h>
#include <vtkGlyph3D.h>
#include <vtkIdFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkInteractorStyleTrackballActor.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPlanes.h>
#include <vtkPointData.h>
#include <vtkPointPicker.h>
#include <vtkPointSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertexGlyphFilter.h>

namespace
{
    // Define interaction style
    class InteractorStyleMoveVertex : public vtkInteractorStyleTrackballActor
    {
    public:
        static InteractorStyleMoveVertex *New();
        vtkTypeMacro(InteractorStyleMoveVertex, vtkInteractorStyleTrackballActor);

        InteractorStyleMoveVertex()
        {

            this->Move = false;
            this->PointPicker = vtkSmartPointer<vtkPointPicker>::New();

            // Setup ghost glyph
            vtkNew<vtkPoints> points;
            points->InsertNextPoint(0, 0, 0);
            this->MovePolyData = vtkSmartPointer<vtkPolyData>::New();
            this->MovePolyData->SetPoints(points);
            this->MoveGlyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
            this->MoveGlyphFilter->SetInputData(this->MovePolyData);
            this->MoveGlyphFilter->Update();

            this->MoveMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            this->MoveMapper->SetInputConnection(
                this->MoveGlyphFilter->GetOutputPort());

            this->MoveActor = vtkSmartPointer<vtkActor>::New();
            this->MoveActor->SetMapper(this->MoveMapper);
            this->MoveActor->VisibilityOff();
            this->MoveActor->GetProperty()->SetPointSize(10);
            this->MoveActor->GetProperty()->SetColor(this->Pink);
        }

        void OnMouseMove() override
        {
            if (!this->Move)
            {
                return;
            }

            vtkInteractorStyleTrackballActor::OnMouseMove();
        }

        void OnMiddleButtonUp() override
        {
            this->EndPan();

            this->Move = false;
            this->MoveActor->VisibilityOff();

            this->Data->GetPoints()->SetPoint(this->SelectedPoint,
                                              this->MoveActor->GetPosition());
            this->Data->Modified();
            /////////////

            vtkNew<vtkVertexGlyphFilter> glyphFilter;
            glyphFilter->SetInputData(this->Data);
            glyphFilter->Update();

            // Create a mapper and actor
            vtkNew<vtkPolyDataMapper> mapper;
            mapper->SetInputConnection(glyphFilter->GetOutputPort());

            this->actor->SetMapper(mapper);
            /////////////
            this->GetCurrentRenderer()->Render();
            this->GetCurrentRenderer()->GetRenderWindow()->Render();
        }
        void OnMiddleButtonDown() override
        {
            // Get the selected point
            int x = this->Interactor->GetEventPosition()[0];
            int y = this->Interactor->GetEventPosition()[1];
            this->FindPokedRenderer(x, y);
            this->PointPicker->Pick(this->Interactor->GetEventPosition()[0],
                                    this->Interactor->GetEventPosition()[1],
                                    0, // always zero.
                                    this->Interactor->GetRenderWindow()
                                        ->GetRenderers()
                                        ->GetFirstRenderer());

            if (this->PointPicker->GetPointId() >= 0)
            {
                this->StartPan();
                this->MoveActor->VisibilityOn();
                this->Move = true;
                this->SelectedPoint = this->PointPicker->GetPointId();

                std::cout << "Dragging point " << this->SelectedPoint << std::endl;

                double p[3];
                this->Data->GetPoint(this->SelectedPoint, p);
                std::cout << "p: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
                this->MoveActor->SetPosition(p);

                this->GetCurrentRenderer()->AddActor(this->MoveActor);
                this->InteractionProp = this->MoveActor;
            }
        }

        // vtkNew<vtkNamedColors> color;
        double Pink[3] = {255./255.,  192./255.,  203./255.};

        vtkUnstructuredGrid *Data;
        vtkPolyData *GlyphData;

        vtkSmartPointer<vtkPolyDataMapper> MoveMapper;
        vtkSmartPointer<vtkActor> MoveActor;
        vtkSmartPointer<vtkActor> actor;
        vtkSmartPointer<vtkPolyData> MovePolyData;
        vtkSmartPointer<vtkVertexGlyphFilter> MoveGlyphFilter;

        vtkSmartPointer<vtkPointPicker> PointPicker;

        bool Move;
        vtkIdType SelectedPoint;
    };
    vtkStandardNewMacro(InteractorStyleMoveVertex);

} // namespace

TEST(tests_vtkexamples, movepoint)
{
    
    vtkNew<vtkPoints> points;
    points->InsertNextPoint(0, 0, 0);
    points->InsertNextPoint(1, 0, 0);
    points->InsertNextPoint(2, 0, 0);

    vtkNew<vtkUnstructuredGrid> input;
    input->SetPoints(points);

    vtkNew<vtkVertexGlyphFilter> glyphFilter;
    glyphFilter->SetInputData(input);
    glyphFilter->Update();

    // Create a mapper and actor
    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(glyphFilter->GetOutputPort());

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetPointSize(10);
    double Tomato[3] = {255./255.,   99./255.,   71./255};
    actor->GetProperty()->SetColor(Tomato);

    // Visualize
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->AddRenderer(renderer);
    renderWindow->SetWindowName("MoveAVertexUnstructuredGrid");

    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    renderWindowInteractor->SetRenderWindow(renderWindow);

    renderer->AddActor(actor);
    renderer->SetBackground(0.95,.95,.95);

    renderWindow->Render();

    vtkNew<InteractorStyleMoveVertex> style;
    renderWindowInteractor->SetInteractorStyle(style);
    style->Data = input;
    style->GlyphData = glyphFilter->GetOutput();
    style->actor = actor;

    renderer->GetActiveCamera()->Zoom(0.9);
    renderWindowInteractor->Start();
}

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkSphereSource.h>
#include <vtkRendererCollection.h>
#include <vtkCellArray.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkPlaneSource.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPropPicker.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>

// Handle mouse events
class MouseInteractorStyle2 : public vtkInteractorStyleTrackballCamera
{
  public:
    static MouseInteractorStyle2* New();
    vtkTypeMacro(MouseInteractorStyle2, vtkInteractorStyleTrackballCamera);

    virtual void OnLeftButtonDown()
    {
      int* clickPos = this->GetInteractor()->GetEventPosition();

      // Pick from this location.
      vtkSmartPointer<vtkPropPicker>  picker =
        vtkSmartPointer<vtkPropPicker>::New();
      picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());

      double* pos = picker->GetPickPosition();
      std::cout << "Pick position (world coordinates) is: "
                << pos[0] << " " << pos[1]
                << " " << pos[2] << std::endl;

      std::cout << "Picked actor: " << picker->GetActor() << std::endl;
      if(picker->GetActor())
      {
          //Create a sphere
          vtkSmartPointer<vtkSphereSource> sphereSource =
              vtkSmartPointer<vtkSphereSource>::New();
          sphereSource->SetCenter(pos[0], pos[1], pos[2]);
          sphereSource->SetRadius(0.1);

          //Create a mapper and actor
          vtkSmartPointer<vtkPolyDataMapper> mapper =
              vtkSmartPointer<vtkPolyDataMapper>::New();
          mapper->SetInputConnection(sphereSource->GetOutputPort());

          vtkSmartPointer<vtkActor> actor =
              vtkSmartPointer<vtkActor>::New();
          actor->SetMapper(mapper);

          //this->GetInteractor()->GetRenderWindow()->GetRenderers()->GetDefaultRenderer()->AddActor(actor);
          this->GetDefaultRenderer()->AddActor(actor);
          this->GetInteractor()->GetRenderWindow()->Render();
      }
      // Forward events
      else
      {
          vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
      }
    }

  private:

};

vtkStandardNewMacro(MouseInteractorStyle2);

// Execute application.
TEST(tests_vtkexamples, Picking)
{
  vtkSmartPointer<vtkPlaneSource> planeSource =
    vtkSmartPointer<vtkPlaneSource>::New();
  planeSource->Update();

  // Create a polydata object
  vtkPolyData* polydata = planeSource->GetOutput();

  // Create a mapper
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  mapper->SetInput ( polydata );
#else
  mapper->SetInputData ( polydata );
#endif

  // Create an actor
  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper ( mapper );

  std::cout << "Actor address: " << actor << std::endl;

  // A renderer and render window
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer ( renderer );

  // An interactor
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow ( renderWindow );

  // Set the custom stype to use for interaction.
  vtkSmartPointer<MouseInteractorStyle2> style =
    vtkSmartPointer<MouseInteractorStyle2>::New();
  style->SetDefaultRenderer(renderer);

  renderWindowInteractor->SetInteractorStyle( style );

  // Add the actors to the scene
  renderer->AddActor ( actor );
  renderer->SetBackground ( 0,0,1 );

  // Render and interact
  renderWindow->Render();
  renderWindowInteractor->Initialize();
  renderWindowInteractor->Start();
}


#include <vtkActor.h>
#include <vtkAreaPicker.h>
#include <vtkCamera.h>
#include <vtkCellPicker.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkExtractGeometry.h>
#include <vtkGlyph3D.h>
#include <vtkIdFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkInteractorStyleTrackballActor.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPlanes.h>
#include <vtkPointData.h>
#include <vtkPointSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertexGlyphFilter.h>

namespace {

// Define interaction style
class InteractorStyleMoveGlyph : public vtkInteractorStyleTrackballActor
{
public:
  static InteractorStyleMoveGlyph* New();
  vtkTypeMacro(InteractorStyleMoveGlyph, vtkInteractorStyleTrackballActor);

  InteractorStyleMoveGlyph()
  {
    this->MoveSphereSource = vtkSmartPointer<vtkSphereSource>::New();
    this->MoveSphereSource->SetRadius(.1);
    this->MoveSphereSource->Update();

    this->MoveMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    this->MoveMapper->SetInputConnection(
        this->MoveSphereSource->GetOutputPort());

    this->MoveActor = vtkSmartPointer<vtkActor>::New();
    this->MoveActor->SetMapper(this->MoveMapper);
    this->MoveActor->GetProperty()->SetColor(
        this->color->GetColor3d("Pink").GetData());
    // this->MoveActor->VisibilityOff();

    this->Move = false;
  }

  void OnMouseMove() override
  {
    if (!this->Move)
    {
      return;
    }

    vtkInteractorStyleTrackballActor::OnMouseMove();
  }

  void OnMiddleButtonUp() override
  {
    // Forward events
    vtkInteractorStyleTrackballActor::OnMiddleButtonUp();
    this->Move = false;
    this->MoveActor->VisibilityOff();

    this->Data->GetPoints()->SetPoint(this->SelectedPoint,
                                      this->MoveActor->GetPosition());
    this->Data->Modified();
    this->GetCurrentRenderer()->Render();
    this->GetCurrentRenderer()->GetRenderWindow()->Render();
  }
  void OnMiddleButtonDown() override
  {
    // Forward events
    vtkInteractorStyleTrackballActor::OnMiddleButtonDown();
    this->MoveActor->VisibilityOn();
    if (static_cast<vtkCellPicker*>(this->InteractionPicker)->GetPointId() >= 0)
    {
      vtkIdType id =
          dynamic_cast<vtkIdTypeArray*>(
              this->GlyphData->GetPointData()->GetArray("InputPointIds"))
              ->GetValue(static_cast<vtkCellPicker*>(this->InteractionPicker)
                             ->GetPointId());
      std::cout << "Id: " << id << std::endl;
      this->Move = true;
      this->SelectedPoint = id;

      double p[3];
      this->Data->GetPoint(id, p);
      std::cout << "p: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
      this->MoveActor->SetPosition(p);
    }

    this->GetCurrentRenderer()->AddActor(this->MoveActor);
    this->InteractionProp = this->MoveActor;
  }
  vtkNew<vtkNamedColors> color;

  vtkPolyData* Data;
  vtkPolyData* GlyphData;

  vtkSmartPointer<vtkPolyDataMapper> MoveMapper;
  vtkSmartPointer<vtkActor> MoveActor;
  vtkSmartPointer<vtkSphereSource> MoveSphereSource;

  bool Move;
  vtkIdType SelectedPoint;
};
vtkStandardNewMacro(InteractorStyleMoveGlyph);

} // namespace

TEST(tests_vtkexamples, moveGlyph)
{
  vtkNew<vtkNamedColors> color;

  vtkNew<vtkPoints> points;
  points->InsertNextPoint(0, 0, 0);
  points->InsertNextPoint(1, 0, 0);
  points->InsertNextPoint(2, 0, 0);

  vtkNew<vtkPolyData> input;
  input->SetPoints(points);

  vtkNew<vtkSphereSource> glyphSource;
  glyphSource->SetRadius(0.1);
  glyphSource->Update();

  vtkNew<vtkGlyph3D> glyph3D;
  glyph3D->GeneratePointIdsOn();
  glyph3D->SetSourceConnection(glyphSource->GetOutputPort());
  glyph3D->SetInputData(input);
  glyph3D->SetScaleModeToDataScalingOff();
  glyph3D->Update();

  // Create a mapper and actor
  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputConnection(glyph3D->GetOutputPort());

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(color->GetColor3d("Tomato").GetData());

  // Visualize
  vtkNew<vtkRenderer> renderer;
  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->AddRenderer(renderer);
  renderWindow->SetWindowName("MoveAGlyph");

  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderer->AddActor(actor);
  renderer->SetBackground(color->GetColor3d("Burlywood").GetData());

  renderWindow->Render();

  vtkNew<InteractorStyleMoveGlyph> style;
  renderWindowInteractor->SetInteractorStyle(style);
  style->Data = input;
  style->GlyphData = glyph3D->GetOutput();

  renderer->GetActiveCamera()->Zoom(0.9);
  renderWindowInteractor->Start();

}
