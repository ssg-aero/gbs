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


@pytest.mark.parametrize("n, p", [(2, 3), (3, 3), (1, 1), (3, 0)])
def test_interpolate_cn_degree_not_below_points(n, p):
    # issue #228: a ValueError, no longer a segmentation fault
    import numpy as np

    x = list(np.linspace(0.0, 1.0, n))
    with pytest.raises(ValueError):
        gbs.interpolate_cn([0.0] * n, x, p=p)
    with pytest.raises(ValueError):
        gbs.interpolate_cn([[xi, 0.0, 0.0] for xi in x], x, p=p)
    if p >= 1:
        with pytest.raises(ValueError):
            gbs.interpolate_cn([[xi, 0.0] for xi in x], p=p)


def test_interpolate_cn_parameters_checked():
    with pytest.raises(ValueError):
        gbs.interpolate_cn([0.0, 1.0, 0.0], [0.0, 1.0], p=1)  # 2 parameters for 3 values
    with pytest.raises(ValueError):
        gbs.interpolate_cn([0.0, 1.0, 0.0], [0.0, 0.0, 1.0], p=1)  # not strictly increasing
    # n = p + 1 is still an interpolation
    f = gbs.interpolate_cn([0.0, 1.0, 0.0, 2.0], [0.0, 0.3, 0.6, 1.0], p=3)
    assert f(0.3) == pytest.approx(1.0)


def test_interpolate_cn_surface_degree_not_below_points():
    pts = [[i, j, 0.1 * i * j] for j in range(2) for i in range(3)]
    with pytest.raises(ValueError):
        gbs.interpolate_cn(pts, 2, 3, 1)
    with pytest.raises(ValueError):
        gbs.interpolate_cn(pts, 2, 2, 2)
    with pytest.raises(ValueError):
        gbs.interpolate_cn(pts, 0, 2, 1)
    gbs.interpolate_cn(pts, 2, 2, 1)


def test_knots_and_hermite_checked():
    # the uniform knot builder and the Hermite interpolations with parameters (the same out-of-
    # bounds accesses as #228 before: a crash of the interpreter)
    with pytest.raises(ValueError):
        gbs.build_simple_mult_flat_knots(0.0, 1.0, 2, 3)
    assert len(gbs.build_simple_mult_flat_knots(0.0, 1.0, 4, 3)) == 8
    Q = [[[0.0, 0.0], [1.0, 0.0]], [[1.0, 1.0], [1.0, 0.0]], [[2.0, 0.0], [1.0, 0.0]]]
    with pytest.raises(ValueError):
        gbs.interpolate_c1(Q, [0.0, 1.0])  # 2 parameters for 3 points
    with pytest.raises(ValueError):
        gbs.interpolate_c1(Q[:1], [0.0])  # one point
    with pytest.raises(ValueError):
        gbs.interpolate_c1(Q, [0.0, 0.5, 0.5])  # not strictly increasing
    crv = gbs.interpolate_c1(Q, [0.0, 0.5, 1.0])
    assert crv.value(0.5)[1] == pytest.approx(1.0)


def test_interpolate_cn_coincident_points():
    with pytest.raises(ValueError, match="coincide"):
        gbs.interpolate_cn([[0.0, 0.0], [1.0, 1.0], [1.0, 1.0], [2.0, 0.0]], p=2)

