import pygbs.gbs as gbs
from math import sqrt
from pytest import approx

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

    x = srf.value(0.5,0.5)

    assert x[0] == approx(0.5,tol)
    assert x[1] == approx(0.5,tol)
    assert x[2] == approx(0.,tol)

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