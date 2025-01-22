import sys
import pytest
import numpy as np
# TODO: check if better using install path
# sys.path.insert(1, 'build/')
import pygbs.gbs as gbs
from math import sqrt

tol = 1e-6
def distance(v1,v2):
    return sqrt(sum( (x-y)**2 for x , y in zip(v1,v2)))
def test_ctor():

    crv = gbs.BSCurve3d(
        [[0.,0.,0.],[1.,0.,0.]],
        [0.,0.,1.,1.],
        1)

    x = crv.value(0.5)

    assert x[0] == pytest.approx(0.5,tol)
    assert x[1] == pytest.approx(0.,tol)
    assert x[2] == pytest.approx(0.,tol)

    crv = gbs.BSCurve2d(
        [[0.,0.],[1.,0.]],
        [0.,0.,1.,1.],
        1)

    x = crv.value(0.5)

    assert x[0] == pytest.approx(0.5,tol)
    assert x[1] == pytest.approx(0.,tol)

    crv = gbs.BSCurve1d(
        [[0.],[1.]],
        [0.,0.,1.,1.],
        1)

    x = crv.value(0.5)

    assert x[0] == pytest.approx(0.5,tol)
    # # simple precision
    # crv = gbs.BSCurve3d_f(
    #     [[0.,0.,0.],[1.,0.,0.]],
    #     [0.,0.,1.,1.],
    #     1)

    # x = crv.value(0.5)

    # assert x[0] == pytest.approx(0.5,tol)
    # assert x[1] == pytest.approx(0.,tol)
    # assert x[2] == pytest.approx(0.,tol)

    # crv = gbs.BSCurve2d_f(
    #     [[0.,0.],[1.,0.]],
    #     [0.,0.,1.,1.],
    #     1)

    # x = crv.value(0.5)

    # assert x[0] == pytest.approx(0.5,tol)
    # assert x[1] == pytest.approx(0.,tol)

    # crv = gbs.BSCurve1d_f(
    #     [[0.],[1.]],
    #     [0.,0.,1.,1.],
    #     1)

    # x = crv.value(0.5)

    # assert x[0] == pytest.approx(0.5,tol)

    # Rational
    crv = gbs.BSCurveRational3d(
        [[0.,0.,0.,1.],[1.,0.,0.,1.]],
        [0.,0.,1.,1.],
        1)

    x = crv.value(0.5)

    assert x[0] == pytest.approx(0.5,tol)
    assert x[1] == pytest.approx(0.,tol)
    assert x[2] == pytest.approx(0.,tol)

    crv = gbs.BSCurve2d(
        [[0.,0.],[1.,0.]],
        [0.,0.,1.,1.],
        1)

    x = crv.value(0.5)

    assert x[0] == pytest.approx(0.5,tol)
    assert x[1] == pytest.approx(0.,tol)

    crv = gbs.BSCurve1d(
        [[0.],[1.]],
        [0.,0.,1.,1.],
        1)

    x = crv.value(0.5)
        
def test_methods():
    crv = gbs.BSCurve3d(
        [
            [0.,0.,0.],
            [0.,1.,0.],
            [1.,1.,0.],
            [1.,1.,1.],
            [1.,1.,2.],
            [3.,1.,1.],
            [0.,4.,1.],
        ],
        [0., 0., 0., 1, 2, 3, 4, 5., 5., 5.],
        2)
    crv_cp = gbs.BSCurve3d( crv )
    crv.increaseDegree()
    p_cp = crv_cp.degree()
    p    = crv.degree()
    assert p_cp == p - 1
    u = np.linspace(0,5,101)
    for u_ in u:
        assert gbs.dist( crv.value( u_ ) , crv_cp.value( u_ ) ) <= tol

    pts = crv(u)

    assert pts.shape == (len(u),3)

    for u_, pt in zip(u, pts):
        assert gbs.dist( crv.value( u_ ) , pt ) <= tol

    writer = gbs.IgesWriter()

    writer.add_geometry(crv)

    writer.write('curves.igs')

def test_to_3d():
    crv_2d = gbs.BSCurve2d(
        [[0.,0.],[1.,0.]],
        [0.,0.,1.,1.],
        1)
    crv_3d = gbs.to_bscurve_3d(crv_2d,1.)

    assert crv_3d(0.5)[0] == pytest.approx(crv_2d(0.5)[0])
    assert crv_3d(0.5)[1] == pytest.approx(crv_2d(0.5)[1])
    assert crv_3d(0.5)[2] == pytest.approx(1.)

def test_offset():
    e = 0.3
    def offset(u,d=0):
        if d == 0:
            return e * u
        elif d == 1:
            return e
        return 0.
    
    crv_2d = gbs.BSCurve2d(
        [[0.,0.],[.5,0.2],[1.,0.]],
        [0.,0., 0., 1.,1.,1.],
    2)

    offset_crv = gbs.CurveOffset2d_func(crv_2d,offset)

    U = gbs.deviation_based_params(offset_crv,np=30)
    assert len(U) >= 30
    for u in U:
        assert gbs.dist( offset_crv(u) , crv_2d(u) ) == pytest.approx( offset(u) )

        
def test_ext():

    crv_2d = gbs.BSCurve2d(
        [[0.,0.],[1.,0.]],
        [0.,0.,1.,1.],
        1)
    crv = gbs.CurveExtended2d(crv_2d)
    
    assert crv( 2.0)[0] == pytest.approx( 2. )
    assert crv(-1.0)[0] == pytest.approx(-1. )

    assert crv( 2.0, 1)[0] == pytest.approx( 1. )
    assert crv(-1.0, 1)[0] == pytest.approx( 1. )

    assert crv( 2.0, 2)[0] == pytest.approx( 0. )
    assert crv(-1.0, 2)[0] == pytest.approx( 0. )