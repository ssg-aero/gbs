"""Defects of the bindings found while writing the documentation of the geometry, each checked against its fix."""
import math

import pytest
import pygbs.gbs as gbs


def _surface():
    return gbs.BSSurface3d(
        [
            [0.0, 0.0, 0.0], [0.5, 0.0, 0.2], [1.0, 0.0, 0.0],
            [0.0, 0.5, 0.1], [0.5, 0.5, 0.5], [1.0, 0.5, 0.1],
            [0.0, 1.0, 0.0], [0.5, 1.0, 0.2], [1.0, 1.0, 0.0],
        ],
        [0., 0., 0., 1., 1., 1.],
        [0., 0., 0., 1., 1., 1.],
        2,
        2,
    )


def test_rotate_and_scale_a_point_return_it():
    # a Python list cannot be changed in place: the in-place forms did nothing
    assert gbs.rotate([1., 0.], math.pi / 2) == pytest.approx([0., 1.], abs=1e-15)
    assert gbs.rotate([2., 1.], math.pi, [1., 1.]) == pytest.approx([0., 1.], abs=1e-15)
    assert gbs.scale([1., 2., 3.], 2.) == pytest.approx([2., 4., 6.])
    assert gbs.scale([1., 2.], 2., [1., 1.]) == pytest.approx([1., 3.])


def test_project_points_returns_them():
    srf = _surface()
    pts = [[0.3, 0.4, 1.0], [0.7, 0.2, -1.0]]
    proj = gbs.project_points(srf, pts)
    for p in proj:
        u, v, d = gbs.extrema_surf_pnt(srf, p, 1e-9)
        assert d == pytest.approx(0., abs=1e-7)


def test_tangential_direction():
    crv = gbs.BSCurve2d([[0., 0.], [1., 1.], [2., 0.]], [0., 0., 0., 1., 1., 1.], 2)
    t = gbs.tangential_direction(crv, 0.)
    assert t == pytest.approx([math.sqrt(0.5), math.sqrt(0.5)])


def test_extended_is_the_whole_surface():
    srf = _surface()
    ext = gbs.extended(srf, 0.2, gbs.SurfaceBound.U_END, True)
    u1, u2, v1, v2 = ext.bounds()
    assert u1 == pytest.approx(0.)
    assert u2 > 1.
    for u, v in [(0.2, 0.3), (0.8, 0.9)]:
        assert gbs.dist(ext(u, v), srf(u, v)) == pytest.approx(0., abs=1e-9)


def test_extreme_coordinates():
    crv = gbs.BSCurve2d([[0., 0.], [1., 3.], [2., -1.], [3., 2.5], [4., 0.5]], [0., 0., 0., 0., 0.5, 1., 1., 1., 1.], 3)
    ys = [crv(i / 10000.)[1] for i in range(10001)]
    u, y = gbs.max_coordinate(crv, 1, 0.)
    assert y == pytest.approx(max(ys), abs=1e-7)
    u, y = gbs.min_coordinate(crv, 1, 0.)
    assert y == pytest.approx(min(ys), abs=1e-7)


def test_join_still_averages():
    # join keeps its documented behaviour (join_curves checks the contact)
    c1 = gbs.BSCurve2d([[0., 0.], [1., 0.]], [0., 0., 1., 1.], 1)
    c2 = gbs.BSCurve2d([[1., 0.], [2., 1.]], [0., 0., 1., 1.], 1)
    j = gbs.join(c1, c2)
    assert j(1.) == pytest.approx([1., 0.])
