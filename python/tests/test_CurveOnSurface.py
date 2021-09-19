import pygbs.gbs as gbs
from math import sqrt

tol = 1e-6
def distance(v1,v2):
    return sqrt(sum( (x-y)**2 for x , y in zip(v1,v2)))


def test_ctor():
    crv_2d = gbs.BSCurve2d(
        [[0.,0.5],[1.,0.5]],
        [0.,0.,1.,1.],
        1)
    srf = gbs.BSSurface3d(
        [[0.,0.,0.],[1.,0.,0.],
         [0.,1.,0.],[1.,1.,0.]],
        [0.,0.,1.,1.],
        [0.,0.,1.,1.],
        1,
        1)

    crv_3d = gbs.CurveOnSurface3d(crv_2d,srf)
    u = 0.3
    u_,v_ = crv_2d(u)
    assert distance(crv_3d(u),srf(u_,v_)) < tol

