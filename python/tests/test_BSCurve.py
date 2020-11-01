import sys
import pytest

# TODO: check if better using install path
sys.path.insert(1, 'build/')
import pygbs

tol = 1e-6

def test_ctor():

    crv = pygbs.BSCurve3d_d(
        [[0.,0.,0.],[1.,0.,0.]],
        [0.,0.,1.,1.],
        1)

    x = crv.value(0.5)

    assert x[0] == pytest.approx(0.5,tol)
    assert x[1] == pytest.approx(0.,tol)
    assert x[2] == pytest.approx(0.,tol)

    crv = pygbs.BSCurve2d_d(
        [[0.,0.],[1.,0.]],
        [0.,0.,1.,1.],
        1)

    x = crv.value(0.5)

    assert x[0] == pytest.approx(0.5,tol)
    assert x[1] == pytest.approx(0.,tol)

    crv = pygbs.BSCurve1d_d(
        [[0.],[1.]],
        [0.,0.,1.,1.],
        1)

    x = crv.value(0.5)

    assert x[0] == pytest.approx(0.5,tol)
    # simple precision
    crv = pygbs.BSCurve3d_f(
        [[0.,0.,0.],[1.,0.,0.]],
        [0.,0.,1.,1.],
        1)

    x = crv.value(0.5)

    assert x[0] == pytest.approx(0.5,tol)
    assert x[1] == pytest.approx(0.,tol)
    assert x[2] == pytest.approx(0.,tol)

    crv = pygbs.BSCurve2d_f(
        [[0.,0.],[1.,0.]],
        [0.,0.,1.,1.],
        1)

    x = crv.value(0.5)

    assert x[0] == pytest.approx(0.5,tol)
    assert x[1] == pytest.approx(0.,tol)

    crv = pygbs.BSCurve1d_f(
        [[0.],[1.]],
        [0.,0.,1.,1.],
        1)

    x = crv.value(0.5)

    assert x[0] == pytest.approx(0.5,tol)

    # Rational
    crv = pygbs.BSCurveRational3d_d(
        [[0.,0.,0.,1.],[1.,0.,0.,1.]],
        [0.,0.,1.,1.],
        1)

    x = crv.value(0.5)

    assert x[0] == pytest.approx(0.5,tol)
    assert x[1] == pytest.approx(0.,tol)
    assert x[2] == pytest.approx(0.,tol)

    crv = pygbs.BSCurve2d_d(
        [[0.,0.],[1.,0.]],
        [0.,0.,1.,1.],
        1)

    x = crv.value(0.5)

    assert x[0] == pytest.approx(0.5,tol)
    assert x[1] == pytest.approx(0.,tol)

    crv = pygbs.BSCurve1d_d(
        [[0.],[1.]],
        [0.,0.,1.,1.],
        1)

    x = crv.value(0.5)
        
