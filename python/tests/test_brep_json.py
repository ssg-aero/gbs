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
    assert brep.equal(m, back) and brep.first_difference(m, back) is None
    assert back.name(s) == "box" and back.external_id(brep.FaceId(2)) == 42
    assert not back.alive(dead)  # dead, at the same id
    assert brep.signed_volume(back, back.solid(s).outer) == brep.signed_volume(m, m.solid(s).outer)
    # the text: standard JSON, versioned, read again into the same text
    text = brep.write_brep_text(m)
    doc = json.loads(text)
    assert doc["format"] == "gbs-brep" and doc["version"] == 1 and doc["vertices"][int(dead)] is None
    assert brep.write_brep_text(brep.read_brep_text(text)) == text
    assert brep.equal(brep.read_brep_text(brep.write_brep_text(m, pretty=True)), m)
    # pickled exactly, through the same text
    assert brep.equal(pickle.loads(pickle.dumps(m)), m)
    assert brep.write_brep_text(copy.deepcopy(m)) == text
    # compacted: the dead vertex gone, the remap of compact()
    ctext, remap = brep.write_brep_compact_text(m)
    compacted = copy.deepcopy(m)
    compacted.compact()
    assert brep.equal(brep.read_brep_text(ctext), compacted)
    assert not remap.map(dead).valid() and remap.map(s) == s
    # a difference, described
    other = copy.deepcopy(m)
    other.set_name(s, "other")
    assert brep.first_difference(m, other) == "solid 0: names differ"


def test_read_errors(tmp_path):
    m = brep.Model()
    box(m)
    text = brep.write_brep_text(m)
    with pytest.raises(brep.BRepError) as e:
        brep.read_brep_text(text[:100])
    assert e.value.code == brep.BuildErrc.invalid_file
    doc = json.loads(text)
    doc["version"] = 2
    with pytest.raises(brep.BRepError) as e:
        brep.read_brep_text(json.dumps(doc))
    assert e.value.code == brep.BuildErrc.unsupported_version and "/version: version 2 of the format" in str(e.value)
    doc["version"] = 1
    doc["edges"][0]["v"][0] = 10_000
    with pytest.raises(brep.BRepError) as e:
        brep.read_brep_text(json.dumps(doc))
    assert e.value.code == brep.BuildErrc.invalid_file and ": /edges/0/v/0: " in str(e.value)
    with pytest.raises(brep.BRepError) as e:
        brep.read_brep(tmp_path / "missing.brep")
    assert e.value.code == brep.BuildErrc.invalid_file


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
    with pytest.raises(brep.BRepError) as e:
        brep.write_brep_text(m)
    assert e.value.code == brep.BuildErrc.unsupported_entity and e.value.shapes == (wid,)
    assert "co-edge 0" in str(e.value)
