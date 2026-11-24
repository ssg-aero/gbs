"""Exact files of a BREP model: brep.write_brep / read_brep (gbs-io/brep_json.h)."""
import copy
import json
import pickle

import pytest

import pygbs.gbs as gbs

brep = gbs.brep


def quad(a, b, c, d):
    return gbs.BSSurface3d([a, b, c, d], [0., 0., 1., 1.], [0., 0., 1., 1.], 1, 1)


def box(m, size=1.0):
    corner = lambda i: [size * (i & 1), size * ((i >> 1) & 1), size * ((i >> 2) & 1)]
    faces = []
    for axis in range(3):
        for side in range(2):
            a1, a2 = (axis + 1) % 3, (axis + 2) % 3
            idx = lambda u, v: (side << axis) | (u << a1) | (v << a2)
            faces.append(brep.make_face(m, quad(corner(idx(0, 0)), corner(idx(1, 0)), corner(idx(0, 1)), corner(idx(1, 1)))))
    return brep.make_solid(m, brep.sew(m, faces, tol=1e-6).shells[0])


def same_model(a, b):
    """The same entities, ids and numbers (through the Python copies of the entities)."""
    for t in brep.ShapeType.__members__.values():
        assert a.ids(t) == b.ids(t)
    for v in a.ids(brep.ShapeType.Vertex):
        assert a.vertex(v).pnt == b.vertex(v).pnt and a.vertex(v).tol == b.vertex(v).tol
    for e in a.ids(brep.ShapeType.Edge):
        x, y = a.edge(e), b.edge(e)
        assert (x.u1, x.u2, x.v1, x.v2, x.tol) == (y.u1, y.u2, y.v1, y.v2, y.tol)
    assert a.unit_scale == b.unit_scale


def test_write_read(tmp_path):
    m = brep.Model()
    s = box(m, 2.0)
    m.set_name(s, "box")
    m.set_external_id(brep.FaceId(2), 42)
    dead = brep.make_vertex(m, [5., 5., 5.])
    m.erase(dead)
    m.unit_scale = 25.4
    path = tmp_path / "box.brep"
    brep.write_brep(m, path)
    back = brep.read_brep(path)
    same_model(m, back)
    assert back.name(s) == "box" and back.external_id(brep.FaceId(2)) == 42
    assert not back.alive(dead)  # dead, at the same id
    assert brep.signed_volume(back, back.solid(s).outer) == brep.signed_volume(m, m.solid(s).outer)
    # the text: JSON, versioned, read again into the same text
    text = brep.write_brep_text(m)
    doc = json.loads(text)
    assert doc["format"] == "gbs-brep" and doc["version"] == 1 and doc["vertices"][int(dead)] is None
    assert brep.write_brep_text(brep.read_brep_text(text)) == text
    assert brep.read_brep_text(brep.write_brep_text(m, pretty=True)) is not None
    # pickled exactly, through the same text
    assert brep.write_brep_text(pickle.loads(pickle.dumps(m))) == text
    assert brep.write_brep_text(copy.deepcopy(m)) == text


def test_read_errors(tmp_path):
    m = brep.Model()
    box(m)
    text = brep.write_brep_text(m)
    with pytest.raises(brep.BrepReadError) as e:
        brep.read_brep_text(text[:100])
    assert e.value.code == brep.BrepReadErrc.Syntax and isinstance(e.value, ValueError)
    doc = json.loads(text)
    doc["version"] = 2
    with pytest.raises(brep.BrepReadError) as e:
        brep.read_brep_text(json.dumps(doc))
    assert e.value.code == brep.BrepReadErrc.NewerVersion and e.value.where == "/version" and "version 2" in e.value.reason
    doc["version"] = 1
    doc["edges"][0]["v"][0] = 10_000
    with pytest.raises(brep.BrepReadError) as e:
        brep.read_brep_text(json.dumps(doc))
    assert e.value.code == brep.BrepReadErrc.Invalid and e.value.where == "/edges/0/v/0"
    with pytest.raises(OSError):
        brep.read_brep(tmp_path / "missing.brep")


def test_write_errors():
    # an old Circle2d as a pcurve: its constructor normalizes, it cannot be read back exactly
    m = brep.Model()
    box(m)
    wid = m.face(brep.FaceId(0)).wires[0]
    w = m.wire(wid)
    coedges = list(w.coedges)
    coedges[0].pcurve = gbs.Circle2d(1.0, [0., 0.])
    w.coedges = coedges
    m.set_wire(wid, w)
    with pytest.raises(brep.BrepWriteError) as e:
        brep.write_brep_text(m)
    assert isinstance(e.value, brep.BRepError) and "co-edge 0" in str(e.value)
