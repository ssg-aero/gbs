"""fill_n_sided (gbs/nsided.h): the octant of the unit sphere filled from its three quarter great
circles, the ribbons' cross directions constant along them."""
import math

import numpy as np
import pygbs.gbs as gbs
import pytest


def octant_ribbons():
    q = math.pi / 2
    frames = [([0, 0, 1], [1, 0, 0]), ([1, 0, 0], [0, 1, 0]), ([0, 1, 0], [0, 0, 1])]
    cross = [[0, 0, 1], [1, 0, 0], [0, 1, 0]]
    ribbons = []
    for (n, x), c in zip(frames, cross):
        arc = gbs.EllipseCurve(1.0, 1.0, [[0.0, 0.0, 0.0], n, x])
        constant = gbs.BSCurve3d([c, c], [0.0, 0.0, q, q], 1)
        ribbons.append(gbs.Ribbon(arc, constant, 0.0, q))
    return ribbons


def inside(polygon, p):
    m = len(polygon)
    for i in range(m):
        a, b = np.array(polygon[i]), np.array(polygon[(i + 1) % m])
        e = b - a
        if e[0] * (p[1] - a[1]) - e[1] * (p[0] - a[0]) <= 0.0:
            return False
    return True


def test_octant():
    r = gbs.fill_n_sided(octant_ribbons())
    assert r.deviation <= 1e-5 and r.angle <= 1e-5
    assert len(r.polygon) == 3
    gap = 0.0
    for u in np.linspace(-1, 1, 31):
        for v in np.linspace(-1, 1, 31):
            if inside(r.polygon, (u, v)):
                gap = max(gap, abs(np.linalg.norm(r.surface.value(u, v)) - 1.0))
    # the surface of least bending energy, not the sphere: within 8.5 % of the radius
    assert gap <= 0.085


def test_guide_sphere():
    # leaning toward the sphere itself: within 1e-3 of it inside
    o = gbs.NSidedOptions()
    o.sphere = gbs.GuideSphere([0.0, 0.0, 0.0], 1.0, 1e4)
    r = gbs.fill_n_sided(octant_ribbons(), o)
    assert r.deviation <= 1e-5 and r.angle <= 1e-5
    gap = 0.0
    for u in np.linspace(-1, 1, 31):
        for v in np.linspace(-1, 1, 31):
            if inside(r.polygon, (u, v)):
                gap = max(gap, abs(np.linalg.norm(r.surface.value(u, v)) - 1.0))
    assert gap <= 1e-3


def test_ribbon_interval_by_default():
    r = octant_ribbons()[0]
    assert (gbs.Ribbon(r.curve, r.cross).t0, gbs.Ribbon(r.curve, r.cross).t1) == (0.0, 1.0)


def test_refusals():
    ribbons = octant_ribbons()
    with pytest.raises(ValueError, match="fewer than three"):
        gbs.fill_n_sided(ribbons[:2])
    with pytest.raises(ValueError, match="do not follow"):
        gbs.fill_n_sided([ribbons[0], ribbons[2], ribbons[1]])
    # the cross direction of a side turned out of the hole: the normals disagree at its corners
    q = math.pi / 2
    out = gbs.BSCurve3d([[-1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]], [0.0, 0.0, q, q], 1)
    with pytest.raises(ValueError, match="normals disagree"):
        gbs.fill_n_sided([ribbons[0], gbs.Ribbon(ribbons[1].curve, out, 0.0, q), ribbons[2]])
    o = gbs.NSidedOptions()
    o.side_weight = 0.0
    with pytest.raises(ValueError, match="not positive"):
        gbs.fill_n_sided(ribbons, o)
    o = gbs.NSidedOptions()
    o.max_spans = 8
    with pytest.raises(RuntimeError, match="the shape moved"):
        gbs.fill_n_sided(ribbons, o)
