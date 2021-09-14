import sys
import pytest
import numpy as np
# TODO: check if better using install path
sys.path.insert(1, 'build/')
import pygbs as gbs
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
    assert crv_cp.degree() == crv.degree() - 1
    u = np.linspace(0,5,101)
    for u_ in u:
        assert distance( crv.value( u_ ) , crv_cp.value( u_ ) ) <= tol