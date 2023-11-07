import pygbs.gbs as gbs
from math import sqrt
from pytest import approx
import numpy as np

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