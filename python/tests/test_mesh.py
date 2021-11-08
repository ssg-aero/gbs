from pygbs import gbs
import pytest

def test_msh_curves_lattice():
    
    iso_eth_crv1 = gbs.BSCurve2d(
        [[0.,0.],[0.8,0.],[1.2,1.],[2.,1.]],
        [0.,0.5,1.],
        [3,1,3],
        2
    )
    iso_eth_crv2 = gbs.BSCurve2d(
        [[0.,1.0],[0.5,1.2],[1.5,1.5],[2.,1.5]],
        [0.,0.3,1.],
        [3,1,3],
        2
    )
    iso_eth_crv3 = gbs.BSCurve2d(
        [[0.,2.],[0.7,2.],[1.5,2.],[2.,2.]],
        [0.,0.5,1.],
        [3,1,3],
        2
    )

    iso_ksi_crv1 = gbs.BSCurve2d(
        [iso_eth_crv1.begin(),iso_eth_crv3.begin()],
        [0., 1.],
        [2,2],
        1
    )
    iso_ksi_crv2 = gbs.BSCurve2d(
        [iso_eth_crv1.end(),iso_eth_crv3.end()],
        [0., 1.],
        [2,2],
        1
    )

    ksi_i = [0., 0.5, 1.]
    eth_j = [0., 1.]
    
    X_ksi, X_eth, X_ksi_eth, ksi, eth = gbs.msh_curves_lattice(
        [iso_eth_crv1, iso_eth_crv2, iso_eth_crv3],
        [iso_ksi_crv1, iso_ksi_crv2],
        ksi_i,
        eth_j,
        [20],
        [10,10],
    )
    dims = [ len(X_ksi), len(X_eth), 1]

    pts = gbs.tfi_mesh( X_ksi, X_eth, X_ksi_eth, ksi_i, eth_j, ksi, eth )

    assert len(pts) == dims[0] * dims[1] * dims[2]

    assert gbs.dist(iso_eth_crv1.begin(),pts[0]) == pytest.approx(0.)
    assert gbs.dist(iso_eth_crv1.end(),pts[dims[0]-1]) == pytest.approx(0.)
    assert gbs.dist(iso_eth_crv3.begin(),pts[dims[0]*( dims[1] - 1)]) == pytest.approx(0.)
    assert gbs.dist(iso_eth_crv3.end(),pts[dims[0]* dims[1] - 1]) == pytest.approx(0.)

    plot = False

    if plot:
        from vtkmodules.vtkCommonDataModel import vtkStructuredGrid
        from vtkmodules.vtkCommonCore import (
            vtkPoints
        )
        from vtkmodules.vtkRenderingCore import (
            vtkActor,
            vtkDataSetMapper,
            vtkRenderWindow,
            vtkRenderWindowInteractor,
            vtkRenderer
        )
        from vtkmodules.vtkCommonColor import vtkNamedColors
        points = vtkPoints()
        points.Allocate(dims[0] * dims[1] * dims[2])
        offset =0
        for pt in pts:
            if len(pt)==2:
                pt.append(0.)
            points.InsertPoint(offset,pt)
            offset += 1
        sgrid = vtkStructuredGrid()
        sgrid.SetDimensions(dims)
        sgrid.SetPoints(points)

        colors = vtkNamedColors()

        iso_eth_crv1_actor = gbs.make_actor(iso_eth_crv1)
        iso_eth_crv1_actor.GetProperty().SetColor(colors.GetColor3d('Chartreuse'))
        iso_eth_crv2_actor = gbs.make_actor(iso_eth_crv2)
        iso_eth_crv2_actor.GetProperty().SetColor(colors.GetColor3d('Chartreuse'))
        iso_eth_crv3_actor = gbs.make_actor(iso_eth_crv3)
        iso_eth_crv3_actor.GetProperty().SetColor(colors.GetColor3d('Chartreuse'))
        iso_ksi_crv1_actor = gbs.make_actor(iso_ksi_crv1)
        iso_ksi_crv2_actor = gbs.make_actor(iso_ksi_crv2)

        sgridMapper = vtkDataSetMapper()
        sgridMapper.SetInputData(sgrid)
        sgridActor = vtkActor()
        sgridActor.SetMapper(sgridMapper)
        sgridActor.GetProperty().SetColor(colors.GetColor3d('Peacock'))
        sgridActor.GetProperty().EdgeVisibilityOn()

        renderer = vtkRenderer()
        renWin = vtkRenderWindow()
        renWin.AddRenderer(renderer)

        iren = vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)

        renderer.AddActor(sgridActor)
        renderer.AddActor(iso_eth_crv1_actor)
        renderer.AddActor(iso_eth_crv2_actor)
        renderer.AddActor(iso_eth_crv3_actor)
        renderer.AddActor(iso_ksi_crv1_actor)
        renderer.AddActor(iso_ksi_crv2_actor)
        renderer.SetBackground(colors.GetColor3d('White'))
        renderer.ResetCamera()
        # renWin.SetSize(640, 480)

        # Interact with the data.
        renWin.Render()
        iren.Start()
