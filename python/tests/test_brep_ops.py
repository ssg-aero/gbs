"""Python bindings of stage 2 of the BREP kernel (gbs.brep): classification, section, split,
booleans, sweeps, offset."""
import math
import os

import pytest

import pygbs.gbs as gbs

brep = gbs.brep
STEP_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "in", "step")


def read(m, name):
    root, report = brep.read_step(m, os.path.join(STEP_DIR, name))
    assert report.check.ok()
    return root


def translation(t):
    p = brep.Placement()
    p.t = t
    return p


def volume(m, shape):
    v = 0.
    for s in brep.explore_solids(m, shape):
        sol = m.solid(s)
        v += brep.signed_volume(m, sol.outer, 128)
        for c in sol.voids:
            v += brep.signed_volume(m, c, 128)
    return v


def test_classify_and_section():
    m = brep.Model()
    box = read(m, "box.stp")  # 2 x 3 x 4
    solid = brep.explore_solids(m, box)[0]
    assert brep.classify(m, solid, [1., 1.5, 2.]) == brep.Location.In
    assert brep.classify(m, solid, [0., 1.5, 2.]) == brep.Location.On
    assert brep.classify(m, solid, [5., 5., 5.]) == brep.Location.Out
    face = brep.explore_faces(m, box)[0]
    uv = brep.interior_point(m, face)
    assert brep.classify(m, face, uv) == brep.Location.In
    plane = gbs.Plane([[0., 0., 1.], [0., 0., 1.], [1., 0., 0.]], [-10., 10., -10., 10.])
    sec = brep.section(m, box, plane)
    assert sec.complete and len(sec.edges) == 4 and len(sec.wires) == 1
    assert m.wire(sec.wires[0]).closed


def test_split_and_booleans():
    m = brep.Model()
    a = read(m, "box.stp")
    b = brep.copy_shape(m, a, translation([1., 1., 1.]))
    assert math.isclose(volume(m, brep.fuse(m, a, b).shape), 24. + 24. - 6., rel_tol=2e-3)
    assert math.isclose(volume(m, brep.common(m, a, b).shape), 6., rel_tol=2e-3)
    r = brep.cut(m, a, b)
    assert len(r.solids) == 1
    assert math.isclose(volume(m, r.shape), 18., rel_tol=2e-3)
    assert brep.check(m, r.shape).ok()
    assert {o for _, o in r.faces} == {brep.FaceOrigin.A, brep.FaceOrigin.B}
    # a sphere touching the box: not handled
    sph = read(m, "sphere.stp")
    touching = brep.copy_shape(m, sph, translation([1., 1.5, 7.]))
    with pytest.raises(brep.BRepError):
        brep.fuse(m, a, touching)
    # split by a plane: the four sides in two
    plane = gbs.Plane([[0., 0., 2.], [0.2, 0.1, 1.], [1., 0., 0.]], [-10., 10., -10., 10.])
    s = brep.split(m, a, plane)
    assert len(s.faces) == 4
    assert len(brep.explore_faces(m, a)) == 10


def test_sweeps_and_offset():
    m = brep.Model()
    plate = read(m, "plate_inch.stp")  # 4 x 2 in with a hole, in mm
    face = brep.explore_faces(m, plate)[0]
    prism = brep.make_prism(m, face, [0., 0., 10.])
    area = (8. - math.pi * 0.25) * 25.4 * 25.4
    assert math.isclose(volume(m, prism.shape), area * 10., rel_tol=2e-3)
    ring = brep.make_revolution(m, face, [0., -100., 0.], [1., 0., 0.], math.pi)
    assert brep.check(m, ring.shape).ok()
    cyl = read(m, "cylinder.stp")
    for f in brep.explore_faces(m, cyl):
        g = brep.offset_face(m, f, 1.)
        assert brep.check(m, g).ok()
