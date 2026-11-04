from pygbs import gbs
from pygbs import vtkplot as vbs
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

    plot = True

    if plot:
        from vtkmodules.vtkCommonColor import vtkNamedColors

        sgrid = gbs.make_structuredgrid(pts, dims[0], dims[1])

        colors = vtkNamedColors()

        sgridActor = vbs.make_sgrid_actor(sgrid,grid_only=False) # colors.GetColor3d('Peacock')

        iso_eth_crv1_actor = gbs.make_actor(iso_eth_crv1)
        iso_eth_crv1_actor.GetProperty().SetColor(colors.GetColor3d('Chartreuse'))
        iso_eth_crv2_actor = gbs.make_actor(iso_eth_crv2)
        iso_eth_crv2_actor.GetProperty().SetColor(colors.GetColor3d('Chartreuse'))
        iso_eth_crv3_actor = gbs.make_actor(iso_eth_crv3)
        iso_eth_crv3_actor.GetProperty().SetColor(colors.GetColor3d('Chartreuse'))
        iso_ksi_crv1_actor = gbs.make_actor(iso_ksi_crv1)
        iso_ksi_crv2_actor = gbs.make_actor(iso_ksi_crv2)

        vbs.render_actors([
            sgridActor,
            iso_eth_crv1_actor,
            iso_eth_crv2_actor,
            iso_eth_crv3_actor,
            iso_ksi_crv1_actor,
            iso_ksi_crv2_actor,
        ])


def test_mesh_surface():
    srf = gbs.BSSurface3d(
        [
            [0.0,0.0,0.0],[0.5,0.0,0.0],[1.0,0.0,0.0],
            [0.0,0.5,0.0],[0.5,0.5,0.5],[1.0,0.5,0.0],
            [0.0,0.0,1.0],[0.5,0.0,0.0],[1.0,1.0,1.0],
        ],
        [0.,0.,0.,1.,1.,1.],
        [0.,0.,0.,1.,1.,1.],
        2,
        2)
    pts, ni, nj, n_iso_ksi, n_iso_eth = gbs.tfi_mesh(srf,[0., 0.33, 0.66, 1.],[0., 0.33, 0.66, 1.],30,30)
    gbs.project_points(srf, pts)
    sgrid = gbs.make_structuredgrid(pts, ni, nj)
    sgridActor = vbs.make_sgrid_actor(sgrid, color='Black' ,grid_only=True)
    # sgridActor = vbs.make_sgrid_actor(sgrid ,grid_only=False)
    vbs.render_actors([
        gbs.make_actor(srf), 
        gbs.make_actor(pts,render_as_sphere=True), 
        sgridActor,
    ])