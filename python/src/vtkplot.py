from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkDataSetMapper,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkRenderer
)
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera
from vtkmodules.vtkCommonColor import vtkNamedColors

def make_sgrid_actor(sgrid, color='Peacock', edges_on=True, grid_only=True):
    sgridMapper = vtkDataSetMapper()
    sgridMapper.SetInputData(sgrid)
    sgridActor = vtkActor()
    sgridActor.SetMapper(sgridMapper)
    sgridActor.GetProperty().SetColor(vtkNamedColors().GetColor3d(color))
    if edges_on and not grid_only: sgridActor.GetProperty().EdgeVisibilityOn()
    if grid_only: 
        sgridActor.GetProperty().SetRepresentationToWireframe()
    return sgridActor

def render_actors(actor_lst, bg_color='White'):
    renderer = vtkRenderer()
    renWin = vtkRenderWindow()
    renWin.AddRenderer(renderer)

    iren = vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    style = vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(style)
    
    for actor in actor_lst:
        renderer.AddActor(actor)
    
    colors = vtkNamedColors()
    renderer.SetBackground(colors.GetColor3d(bg_color))
    renderer.ResetCamera()
    # renWin.SetSize(640, 480)

    # Interact with the data.
    renWin.Render()
    iren.Start()
