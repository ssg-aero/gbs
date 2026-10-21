"""The examples of docs/sources/examples (cases of the gallery and of the notebooks) stay valid:
every boolean case gives valid, closed results, except the known limits, which still fail as
documented."""
import math
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "docs", "sources", "examples"))
import gbs_examples as ex  # noqa: E402

brep = ex.brep
KNOWN_LIMITS = {"nurbs_sphere_touching"}


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


@pytest.mark.parametrize("case", ex.CASES, ids=[c.name for c in ex.CASES])
def test_boolean_cases(case):
    r = ex.run(case)
    if case.name in KNOWN_LIMITS:
        assert set(r.errors) == set(ex.OPERATIONS)
        return
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
