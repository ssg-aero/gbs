import pytest
import sys

# TODO: check if better using install path
sys.path.insert(1, 'build/')
import pygbs as gbs

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

    crv = gbs.interpolate_cn_3d(
        constrains,
        2,
        gbs.KnotsCalcMode.CHORD_LENGTH
    )

    assert  distance(crv.begin(),pts[0]) <= tol
    assert  distance(crv.end(),pts[-1]) <= tol

    for p in pts:
        result = gbs.extrema_PC_3d(crv,p,tol)
        assert result.d <= tol * 10
        assert distance(crv.value(result.u),p) <= tol
    