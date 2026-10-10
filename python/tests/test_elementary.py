from math import cos, sin, pi, isclose

import numpy as np
import pytest

import pygbs.gbs as gbs


def test_arcs():
    ax = [[1., 2., 3.], [0., 0., 1.], [1., 0., 0.]]
    c = gbs.build_circle_arc(2., 0.5, 4., ax)
    assert c.bounds() == pytest.approx([0.5, 4.])
    for t in np.linspace(0.5, 4., 15):
        u = gbs.ellipse_arc_parameter(0.5, 4., t)
        assert isclose(gbs.ellipse_arc_angle(0.5, 4., u), t, abs_tol=1e-12)
        assert c.value(u) == pytest.approx([1. + 2. * cos(t), 2. + 2. * sin(t), 3.], abs=1e-12)
    e = gbs.build_ellipse_arc(3., 1., 0., pi)
    assert e.value(gbs.ellipse_arc_parameter(0., pi, 1.)) == pytest.approx([3. * cos(1.), sin(1.), 0.], abs=1e-12)
    c2 = gbs.build_circle_arc(1., 0., pi / 2, [1., 1.])
    assert c2.end() == pytest.approx([1., 2.], abs=1e-12)
    with pytest.raises(ValueError):
        gbs.build_circle_arc(1., 1., 0.)


def test_revolution_extrusion_and_elementary_surfaces():
    seg = gbs.build_segment3d([1., 0., 0.], [2., 0., 1.])
    s = gbs.build_revolution(seg, [[0., 0., 0.], [0., 0., 1.]])
    u = gbs.ellipse_arc_parameter(0., 2. * pi, 1.)
    assert s.value(u, 0.5 * seg.bounds()[1]) == pytest.approx([1.5 * cos(1.), 1.5 * sin(1.), 0.5], abs=1e-12)

    ext = gbs.build_extrusion(gbs.build_circle_arc(1., 0., pi), [0., 0., 2.], 0., 1.)
    assert ext.value(0., 1.) == pytest.approx([1., 0., 2.], abs=1e-12)

    ax = [[0., 0., 0.], [0., 0., 1.], [1., 0., 0.]]
    cyl = gbs.build_cylinder(2., ax, 0., 3.)
    assert cyl.value(u, 1.) == pytest.approx([2. * cos(1.), 2. * sin(1.), 1.], abs=1e-12)
    cone = gbs.build_cone(1., pi / 4, ax, 0., 1.)
    assert cone.value(u, 1.) == pytest.approx([2. * cos(1.), 2. * sin(1.), 1.], abs=1e-12)
    sph = gbs.build_sphere(3.)
    assert np.linalg.norm(sph.value(1.3, 0.4)) == pytest.approx(3., abs=1e-12)
    tor = gbs.build_torus(3., 1.)
    x, y, z = tor.value(2., 1.)
    assert (np.hypot(x, y) - 3.) ** 2 + z ** 2 == pytest.approx(1., abs=1e-12)
