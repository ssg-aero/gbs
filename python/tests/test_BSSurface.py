import pygbs.gbs as gbs
from math import sqrt
from pytest import approx
import numpy as np
import pyvista as pv
from pygbs import vistaplot as gbv

tol = 1e-6
def distance(v1,v2):
    return sqrt(sum( (x-y)**2 for x , y in zip(v1,v2)))

def test_ctor():
    srf = gbs.BSSurface3d(
        [[0.,0.,0.],[1.,0.,0.],
         [0.,1.,0.],[1.,1.,0.]],
        [0.,0.,1.,1.],
        [0.,0.,1.,1.],
        1,
        1)

    pt = srf.value(0.5,0.5)

    assert pt[0] == approx(0.5,tol)
    assert pt[1] == approx(0.5,tol)
    assert pt[2] == approx(0.,tol)

    nu, nv = (100, 20)
    u = np.linspace(0, 1, nu)
    v = np.linspace(0, 1, nv)
    u_arr, v_arr = np.meshgrid(u, v)
    u_arr.shape=(-1,)
    v_arr.shape=(-1,)

    pts = srf(u_arr, v_arr)

    assert pts.shape == (len(u_arr), 3)

    for u, v, pt in zip(u_arr, v_arr, pts):
        assert srf(u,v) == approx(pt)


def test_loft():
    crv1 = gbs.BSCurve3d(
        [[0.,0.,0.],[1.,0.,0.]],
        [0.,0.,1.,1.],
        1)
    crv2 = gbs.BSCurve3d(
        [[0.,1.,0.],[1.,1.,0.]],
        [0.,0.,1.,1.],
        1)
    srf = gbs.loft([crv1,crv2])

    crv1 = gbs.BSCurve2d(
        [[0.,0.],[1.,0.]],
        [0.,0.,1.,1.],
        1)
    crv2 = gbs.BSCurve2d(
        [[0.,1.],[1.,1.]],
        [0.,0.,1.,1.],
        1)
    srf = gbs.loft([crv1,crv2])

def test_gordon():

    pts = np.array([
        [[0.,0.,0.], [0.33,0.,0.], [0.66,0.,0.], [1.,0.,0.],],
        [[0.,0.5,0.1], [0.33,0.5,0.1], [0.66,0.5,0.2], [1.,0.5,-0.1],],
        [[0.,1.,0.], [0.33,1.,0.], [0.66,1.,0.], [1.,1.,0.],],
    ])

    pts_t = np.transpose(pts, (1, 0, 2))

    u = [0.,0.33,0.66,1.]
    v = [0.,0.5,1.]
    p = 3
    q = 2

    u_crv_lst =[ gbs.interpolate_cn(pts[i], u, p) for i in range(len(v)) ]
    v_crv_lst =[ gbs.interpolate_cn(pts_t[i], v, q) for i in range(len(u)) ]

    srf = gbs.gordonbs(u_crv_lst, v_crv_lst)

    for u_ , crv in zip(u, v_crv_lst):
        for p1, p2 in zip(srf.isoU(u_).poles(), crv.poles()):
            assert p1 == approx( p2 )

    for v_ , crv in zip(v, u_crv_lst):
        for p1, p2 in zip(srf.isoV(v_).poles(), crv.poles()):
            assert p1 == approx( p2 )

    # plotter = pv.Plotter()
    # gbv.add_surfaces_to_plotter([srf], plotter, use_transparency=True,opacity=0.1)
    # gbv.add_curves_to_plotter(u_crv_lst, plotter)
    # gbv.add_curves_to_plotter(v_crv_lst, plotter, def_col='Blue')
    # plotter.show()

