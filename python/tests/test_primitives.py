"""Python bindings of the analytic primitives (gbs/primitives.h)."""
from math import cos, sin, pi, tan

import pytest

import pygbs.gbs as gbs


def test_surfaces():
    cyl = gbs.CylindricalSurface(2., box=[0., 2 * pi, -1., 3.])
    assert cyl.value(1., 2.) == pytest.approx([2 * cos(1.), 2 * sin(1.), 2.])
    assert cyl.parameters([2 * cos(1.), 2 * sin(1.), 2.]) == pytest.approx([1., 2.])
    assert cyl.radius() == 2.
    n = cyl.to_nurbs()
    assert n.value(*cyl.nurbs_parameters(1., 2.)) == pytest.approx(cyl.value(1., 2.), abs=1e-12)
    cone = gbs.ConicalSurface(1., pi / 6)
    assert cone.value(0., 1.) == pytest.approx([1. + tan(pi / 6), 0., 1.])
    sph = gbs.SphericalSurface(3.)
    assert sph.value(0., pi / 2) == pytest.approx([0., 0., 3.], abs=1e-15)
    tor = gbs.ToroidalSurface(3., 1.)
    assert tor.value(0., 0.) == pytest.approx([4., 0., 0.])
    pl = gbs.Plane([[0., 0., 1.], [0., 0., 1.], [1., 0., 0.]], [0., 2., 0., 1.])
    assert pl.value(1., 1.) == pytest.approx([1., 1., 1.])
    assert cyl.iso_u(0.5, -1., 3.).value(1.) == pytest.approx(cyl.value(0.5, 1.))
    with pytest.raises(ValueError):
        gbs.CylindricalSurface(-1.)


def test_curves():
    c = gbs.EllipseCurve(2., 2., [[0., 0., 0.], [0., 0., 1.], [1., 0., 0.]], 5.9, 6.6)
    assert c.is_circle()
    assert c.value(6.2) == pytest.approx([2 * cos(6.2), 2 * sin(6.2), 0.])
    assert c.parameter(c.value(6.5)) == pytest.approx(6.5)
    assert c.reversed().value(-6.2) == pytest.approx(c.value(6.2))
    l = gbs.LineCurve([0., 0., 0.], [0., 2., 0.], -1., 1.)
    assert l.value(0.5) == pytest.approx([0., 1., 0.])
