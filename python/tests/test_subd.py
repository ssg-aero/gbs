"""Catmull-Clark subdivision from Python (gbs.subd): cages, their refusals, refinement and the exact
limit surface, against the values of OpenSubdiv (tests/in/subd, scripts/subd_reference)."""
import json
import pathlib
import sys

import numpy as np
import pytest

import pygbs.gbs as gbs

subd = gbs.subd
ROOT = pathlib.Path(__file__).resolve().parents[2]
IN = ROOT / "tests" / "in" / "subd"


def cage(name):
    c = json.loads((IN / f"{name}.cage.json").read_text())
    return subd.Cage(c["vertices"], c["faces"], c["creases"], c["corners"])


def test_refusal():
    v = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]
    with pytest.raises(subd.CageError) as e:
        subd.LimitSurface(subd.Cage(v, [[0, 1, 2], [0, 1, 3]]))
    assert isinstance(e.value, ValueError)
    assert e.value.issue == subd.CageIssue.InconsistentOrientation
    assert e.value.index == 1
    with pytest.raises(subd.CageError) as e:
        subd.Topology(subd.Cage(v, []))
    assert e.value.issue == subd.CageIssue.NoFace and e.value.index is None


def test_refine_as_opensubdiv():
    c = cage("ngons")
    ref = json.loads((IN / "ngons.osd.json").read_text())
    r = subd.refine(c, 2)
    ours, theirs = np.array(r.vertices), np.array(ref["level2"]["vertices"])
    assert ours.shape == theirs.shape
    d = np.linalg.norm(ours[:, None, :] - theirs[None, :, :], axis=2)
    assert np.all(d.min(axis=0) < 1e-12)
    assert r.vertices[:len(c.vertices)] != c.vertices  # vertex v stays v, moved by the rules
    assert len(r.faces) == 4 * sum(len(f) for f in c.faces)


def test_limit_as_opensubdiv():
    # on the patches of an n-gon, OpenSubdiv measures the derivatives in the parameters of the
    # n-gon: twice ours per order
    c = cage("hexprism")
    ref = json.loads((IN / "hexprism.osd.json").read_text())
    s = subd.LimitSurface(c)
    assert s.patch_count() == ref["patches"]["count"]
    for sample in ref["patches"]["samples"][:50]:
        p, u, v = sample["patch"], sample["u"], sample["v"]
        k = 2.0 if len(c.faces[s.patch_face(p)[0]]) != 4 else 1.0
        assert np.allclose(s.value(p, u, v), sample["p"], atol=1e-12)
        assert np.allclose(k * np.array(s.value(p, u, v, 1, 0)), sample["du"], atol=1e-11)
        assert np.allclose(k * np.array(s.value(p, u, v, 0, 1)), sample["dv"], atol=1e-11)


def test_extraordinary_vertices():
    # the centres of the hexagons, refined: vertices of valence 6; their limit, a tangent plane,
    # no derivative
    c = subd.refine(cage("hexprism"))
    s = subd.LimitSurface(c)
    t = s.topology()
    ev = [v for v in range(t.n_vertices()) if t.valence(v) == 6]
    assert len(ev) == 2
    for v in ev:
        f = next(i for i, face in enumerate(c.faces) if v in face)
        u, w = [(0, 0), (1, 0), (1, 1), (0, 1)][c.faces[f].index(v)]
        p = s.patch(f)
        assert np.allclose(s.value(p, u, w), s.vertex_limit(v), atol=1e-14)
        n = np.array(s.normal(p, u, w))
        near = np.array(s.normal(p, 1e-6 if u == 0 else 1 - 1e-6, 1e-6 if w == 0 else 1 - 1e-6))
        assert np.linalg.norm(n) == pytest.approx(1.0) and n @ near > 0.999
        with pytest.raises(ValueError):
            s.value(p, u, w, 1, 0)


def test_rules():
    t = subd.Topology(cage("creases"))
    rules = [t.rule(v) for v in range(t.n_vertices())]
    assert rules.count(subd.VertexRule.Dart) == 3 and rules.count(subd.VertexRule.Corner) == 2


def test_bulk_evaluation():
    s = subd.LimitSurface(cage("box"))
    rng = np.random.default_rng(1)
    n = 200
    patches = rng.integers(0, s.patch_count(), n).tolist()
    us, vs = rng.uniform(0, 1, n).tolist(), rng.uniform(0, 1, n).tolist()
    for du, dv in ((0, 0), (1, 0), (1, 1)):
        values = s.values(patches, us, vs, du, dv)
        assert values.shape == (n, 3)
        assert np.allclose(values, [s.value(p, u, v, du, dv) for p, u, v in zip(patches, us, vs)], rtol=0, atol=0)
    normals = s.normals(patches, us, vs)
    assert np.allclose(np.linalg.norm(normals, axis=1), 1.0)
    with pytest.raises(ValueError):
        s.values(patches, us[:-1], vs)
