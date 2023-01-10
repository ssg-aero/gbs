import pytest
import sys

# TODO: check if better using install path
sys.path.insert(1, 'build/')
import pygbs.gbs as gbs

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

    crv = gbs.interpolate_cn(
        # constraints,
        pts,
        2,
        gbs.KnotsCalcMode.CHORD_LENGTH
    )

    assert  distance(crv.begin(),pts[0]) <= tol
    assert  distance(crv.end(),pts[-1]) <= tol

    for p in pts:
        u, dist = gbs.extrema_curve_point(crv,p,tol)
        assert dist <= tol * 10
        assert distance(crv.value(u),p) <= tol

    gbs.plot_curves([crv])
    
def test_c1_2d():
    pts = [
        [0.,0],
        [0.,1],
        [1.,0.5],
        [1.,1]
    ]
    constraints = []
    for p in pts:
        constraints.append([p,[1.,0.]])

    crv = gbs.interpolate_c1(
        constraints,
        gbs.KnotsCalcMode.CHORD_LENGTH
    )

    assert  distance(crv.begin(),pts[0]) <= tol
    assert  distance(crv.end(),pts[-1]) <= tol

    for p in pts:
        u, dist = gbs.extrema_curve_point(crv,p,tol)
        assert dist <= tol * 10
        assert distance(crv.value(u),p) <= tol
        assert distance(crv.value(u,1),[1.,0.]) <= tol

    u  = [0., 1., 2., 3.]
    crv = gbs.interpolate_c1(
        constraints,
        u
    )
    for p,u_ in zip(pts,u):
        assert distance(crv.value(u_),p) <= tol
        assert distance(crv.value(u_,1),[1.,0.]) <= tol


def test_constrained_2d():
    pt1 = (0.,[0.,0.]) # curve begin
    pt2 = (1.,[1.,0.]) # curve end
    cstr_lst = [
        (0.5,[0.5,0.4],0), # point at param 0.5
        (0.5,[1.2,0.0],1), # tangent at param 0.5
    ]
    p = 3
    crv = gbs.interpolate(pt1,pt2,cstr_lst,p)

    assert  distance(crv(pt1[0]),pt1[1]) <= tol
    assert  distance(crv(pt2[0]),pt2[1]) <= tol
    assert  distance(crv(cstr_lst[0][0],cstr_lst[0][2]),cstr_lst[0][1]) <= tol
    assert  distance(crv(cstr_lst[1][0],cstr_lst[1][2]),cstr_lst[1][1]) <= tol
