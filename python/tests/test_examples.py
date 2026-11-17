"""The examples of docs/sources/examples (cases of the gallery and of the notebooks) stay valid:
every boolean case gives valid, closed results; the known limits are marked xfail, so that their
fix shows as an unexpected pass rather than a broken test."""
import functools
import math
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "docs", "sources", "examples"))
import gbs_examples as ex  # noqa: E402

brep = ex.brep
# known limits, kept visible (xfail) without pinning them: their fix makes the test pass
KNOWN_LIMITS: dict[str, str] = {}


def test_primitives():
    m = brep.Model()
    for s, volume in [
        (ex.box(m, (0, 0, 0), (2, 3, 4)), 24.0),
        (ex.cylinder(m, 1.0, 2.0), 2 * math.pi),
        (ex.cone(m, 1.5, 2.0), math.pi * 1.5**2 * 2 / 3),
        (ex.sphere(m, 1.5), 4 / 3 * math.pi * 1.5**3),
        (ex.torus(m, 3.0, 1.0), 2 * math.pi**2 * 3.0),
    ]:
        assert brep.check(m, s).ok()
        assert brep.signed_volume(m, m.solid(s).outer, 128) == pytest.approx(volume, rel=1e-3)


@functools.cache
def _run(name):
    """The booleans of the case `name`, computed once for the tests that look at them."""
    return ex.run(next(c for c in ex.CASES if c.name == name))


def _cases():
    for c in ex.CASES:
        marks = [pytest.mark.xfail(reason=KNOWN_LIMITS[c.name], strict=False)] if c.name in KNOWN_LIMITS else []
        yield pytest.param(c, id=c.name, marks=marks)


@pytest.mark.parametrize("case", _cases())
def test_boolean_cases(case):
    r = _run(case.name)
    assert not r.errors
    volumes = {}
    for op, shape in r.shapes.items():
        assert brep.check(r.model, shape).ok(), op
        mesh = brep.mesh(r.model, shape)
        assert mesh.open_sides() == 0, op
        assert mesh.folded == [], op
        volumes[op] = mesh.volume()
    # fuse + common = a + b (by the meshes, within their deflections)
    va = brep.mesh(r.model, r.a).volume()
    vb = brep.mesh(r.model, r.b).volume()
    assert volumes["fuse"] + volumes["common"] == pytest.approx(va + vb, rel=5e-3, abs=1e-6)
    assert volumes["cut"] + volumes["common"] == pytest.approx(va, rel=5e-3, abs=1e-6)


def test_nurbs_contact_like_primitives():
    """A sphere touching a box at a point gives the same results whether it is exact or a NURBS
    (the contact certified on the B-spline, brep_palier4.md § 0 quater)."""
    exact = _run("sphere_touching")
    nurbs = _run("nurbs_sphere_touching")
    for op in ex.OPERATIONS:
        e, n = exact.shapes[op], nurbs.shapes[op]
        assert len(brep.explore_solids(nurbs.model, n)) == len(brep.explore_solids(exact.model, e)), op
        assert len(brep.explore_faces(nurbs.model, n)) == len(brep.explore_faces(exact.model, e)), op


# volumes of the sweeps: the area of the profile times the length of the path (Pappus, the
# centroid of the profile on the path)
_HELIX_LENGTH = 2 * math.hypot(2 * math.pi * 2.0, 1.5)
SWEEP_VOLUMES = {
    "pipe_elbow": math.pi * 0.25 * (3.0 + math.pi + 4.0),
    "pipe_helix": 0.16 * _HELIX_LENGTH,
    "pipe_helix_binormal": 0.16 * _HELIX_LENGTH,
    "pipe_stadium": math.pi * 0.25 * (6 * math.pi + 8.0),
}


@pytest.mark.parametrize("case", ex.SWEEPS, ids=lambda c: c.name)
def test_sweeps(case):
    m = brep.Model()
    s, _ = case.build(m)
    assert brep.check(m, s).ok()
    assert brep.signed_volume(m, m.solid(s).outer, 128) == pytest.approx(SWEEP_VOLUMES[case.name], rel=2e-3)
    mesh = brep.mesh(m, s, brep.MeshOptions(deflection=1e-2))
    assert mesh.open_sides() == 0


@pytest.mark.parametrize("case", ex.OFFSETS, ids=lambda c: c.name)
def test_offsets(case):
    m = brep.Model()
    faces, curves, summary = case.build(m)
    for f, _ in faces:
        assert brep.check(m, f).ok()
    assert "accepted" not in summary
    if case.name == "offset_tangent_plane":
        assert summary.endswith("tangent line")
