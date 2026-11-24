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
    # Xi, unused, may be left out
    assert gbs.max_coordinate(crv, 1)[1] == pytest.approx(max(ys), abs=1e-7)
    # bounds where u_min + (u_max - u_min) passes u_max by an ulp: the last sample was refused
    crv.changeBounds(8.00012375158417, 2829.752226882359)
    u, x = gbs.max_coordinate(crv, 0)
    assert u == 2829.752226882359
    assert x == pytest.approx(4.)


def test_interpolate_c1_to_scale():
    # the end tangents took the chords, not the chords over the parameters: at scale 100 the
    # curve went to x = -3907 for points between 0 and 300
    pts = [[0., 0.], [1., 1.], [2., 1.], [3., 0.]]
    def lerp(a, b, t):  # exact at t = 1, where a + (b - a) may pass b, out of the curve
        return b if t == 1. else a + (b - a) * t

    c1 = gbs.interpolate_c1(pts)
    u1, u2 = c1.bounds()
    for s in (10., 100.):
        cs = gbs.interpolate_c1([[s * x, s * y] for x, y in pts])
        v1, v2 = cs.bounds()
        assert v2 == pytest.approx(s * u2)
        for t in (0., 0.05, 0.5, 0.95, 1.):
            assert gbs.dist(cs(lerp(v1, v2, t)), [s * x for x in c1(lerp(u1, u2, t))]) == pytest.approx(0., abs=1e-9 * s)


def test_bessel_tangents_of_a_parabola():
    u = [0., 1., 3., 4., 7.]
    D = gbs.build_3pt_tg_vec([[x, x * x] for x in u], u)
    for x, d in zip(u, D):
        assert d == pytest.approx([1., 2. * x], abs=1e-12)
    with pytest.raises(ValueError):
        gbs.build_3pt_tg_vec([[0., 0.]], [0.])


def test_poles_set_and_out_of_range():
    crv = gbs.BSCurve2d([[0., 0.], [1., 1.], [2., 0.]], [0., 0., 0., 1., 1., 1.], 2)
    crv.set_pole(1, [1., 2.])
    assert crv.pole(1) == pytest.approx([1., 2.])
    assert crv(0.5) == pytest.approx([1., 1.])
    srf = _surface()
    srf.set_pole(1, 1, [0.5, 0.5, 1.])
    assert srf.pole(1, 1) == pytest.approx([0.5, 0.5, 1.])
    assert srf.pole(2, 0) == pytest.approx([1., 0., 0.])  # i along u
    with pytest.raises(IndexError):
        srf.pole(3, 0)  # past the row: it read the next one


def test_move_to_constraints_delta_returns_the_poles():
    crv = gbs.BSCurve2d([[0., 0.], [1., 1.], [2., 0.], [3., 1.]], [0., 0., 0., 0., 1., 1., 1., 1.], 3)
    k = crv.knotsFlats()
    poles = gbs.move_to_constraints_delta(crv.poles(), [(0.5, [0., 0.25], 0)], k, 3, 0, 3)
    moved = gbs.BSCurve2d(poles, k, 3)
    assert moved(0.5) == pytest.approx([crv(0.5)[0], crv(0.5)[1] + 0.25])


def test_join_tolerance():
    c1 = gbs.BSCurve2d([[0., 0.], [1., 0.]], [0., 0., 1., 1.], 1)
    c2 = gbs.BSCurve2d([[1.1, 0.], [2., 1.]], [0., 0., 1., 1.], 1)
    # without a tolerance: the mean of the ends, as before
    assert gbs.join(c1, c2)(1.) == pytest.approx([1.05, 0.])
    with pytest.raises(ValueError, match="do not meet"):
        gbs.join(c1, c2, tol=1e-3)
    assert gbs.join(c1, c2, 0.2)(1.) == pytest.approx([1.05, 0.])


def test_extended_at_start_keeps_parameters():
    srf = _surface()
    ext = gbs.extended(srf, 0.2, gbs.SurfaceBound.U_START, True)
    u1, u2, v1, v2 = ext.bounds()
    assert u1 < 0.
    assert u2 == pytest.approx(1.)
    for u, v in [(0., 0.3), (0.5, 0.9), (1., 0.)]:
        assert gbs.dist(ext(u, v), srf(u, v)) == pytest.approx(0., abs=1e-9)
