import pytest
import sys

# TODO: check if better using install path
sys.path.insert(1, 'build/')
import pygbs

tol = 1e-6


def distance(v1,v2):
    return sum( (x-y)**2 for x , y in zip(v1,v2))
def test_cn():

    pts = [
        [0.,0.,0],
        [0.,0.,1],
        [1.,0.,0.5],
        [1.,1.,1]
    ]
    constrains = []
    for p in pts:
        constrains.append([p])

    crv = pygbs.interpolate_cn_3d_d(
        constrains,
        2,
        pygbs.KnotsCalcMode.CHORD_LENGTH
    )

    assert  distance(crv.begin(),pts[0]) <= tol
    assert  distance(crv.end(),pts[-1]) <= tol

    for p in pts:
        result = pygbs.extrema_PC(crv,p,tol)
        assert result.d <= tol
        assert distance(crv.value(result.u),p) <= tol
    