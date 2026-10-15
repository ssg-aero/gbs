"""Python bindings of the intersections (gbs.intersect)."""
import math

import pygbs.gbs as gbs

inter = gbs.intersect


def segment2d(a, b):
    return gbs.BSCurve2d([a, b], [0., 0., 1., 1.], 1)


def test_curves():
    # two crossing segments; a line through a cylinder
    r = inter.curves(segment2d([0., 0.], [2., 2.]), segment2d([0., 2.], [2., 0.]))
    assert r.complete
    assert len(r.hits) == 1
    assert math.isclose(r.hits[0].p[0], 1., abs_tol=1e-12)
    line = gbs.LineCurve([-2., 0.5, 0.], [1., 0., 0.], 0., 4.)
    cyl = gbs.CylindricalSurface(1., [[0., 0., 0.], [0., 0., 1.], [1., 0., 0.]], [0., 2 * math.pi, -1., 1.])
    rs = inter.curve_surface(line, cyl)
    assert len(rs.hits) == 2
    for h in rs.hits:
        assert math.isclose(abs(h.p[0]), math.sqrt(0.75), abs_tol=1e-9)
    # a range and a 3D curve-curve
    l1 = gbs.LineCurve([0., 0., 0.], [1., 1., 1.], 0., 1.)
    l2 = gbs.LineCurve([1., 0., 0.], [-1., 1., 1.], 0., 1.)
    r3 = inter.curves(l1, l2)
    assert len(r3.hits) == 1
    assert inter.curves(l1, l2, range1=[0., 0.4]).hits == []


def test_surfaces():
    ax = [[0., 0., 0.], [0., 0., 1.], [1., 0., 0.]]
    cyl = gbs.CylindricalSurface(2., ax, [0., 2 * math.pi, -5., 5.])
    plane = gbs.Plane([[0., 0., 1.], [0., 0., 1.], [1., 0., 0.]], [-5., 5., -5., 5.])
    r = inter.surfaces(cyl, plane)
    assert r.complete and not r.coincident
    assert len(r.curves) == 1
    c = r.curves[0]
    assert c.closed and c.deviation <= 1e-6
    # the circle on both surfaces through its pcurves
    for k in range(5):
        t = c.a + (c.b - c.a) * k / 4
        p = c.curve.value(t)
        q1 = c.pc1.value(t)
        q2 = c.pc2.value(t)
        assert math.dist(cyl.value(q1[0], q1[1]), p) < 1e-9
        assert math.dist(plane.value(q2[0], q2[1]), p) < 1e-9
    # the same pair by marching
    m = inter.surfaces(cyl, plane, closed_forms=False)
    assert len(m.curves) >= 1 and not m.singular
