"""Python bindings of the native BREP core (gbs.brep)."""
import math

import pytest

import pygbs.gbs as gbs

brep = gbs.brep


def quad(a, b, c, d):
    """Bilinear patch S(0,0)=a, S(1,0)=b, S(0,1)=c, S(1,1)=d."""
    return gbs.BSSurface3d([a, b, c, d], [0., 0., 1., 1.], [0., 0., 1., 1.], 1, 1)


def box_faces(m, size=1.0):
    corner = lambda i: [size * (i & 1), size * ((i >> 1) & 1), size * ((i >> 2) & 1)]
    faces = []
    for axis in range(3):
        for side in range(2):
            a1, a2 = (axis + 1) % 3, (axis + 2) % 3
            idx = lambda u, v: (side << axis) | (u << a1) | (v << a2)
            srf = quad(corner(idx(0, 0)), corner(idx(1, 0)), corner(idx(0, 1)), corner(idx(1, 1)))
            faces.append(brep.make_face(m, srf))
    return faces


def test_ids_and_model():
    m = brep.Model()
    v = brep.make_vertex(m, [1., 2., 3.], tol=1e-5)
    assert isinstance(v, brep.VertexId)
    assert v.valid() and int(v) == 0
    assert repr(v) == "VertexId(0)"
    assert not brep.EdgeId().valid()
    assert {v: 1}[brep.VertexId(0)] == 1  # hashable, comparable
    vx = m.vertex(v)
    assert vx.pnt == [1., 2., 3.] and vx.tol == 1e-5
    vx.tol = 1e-4                       # a copy: the model is unchanged…
    assert m.vertex(v).tol == 1e-5
    m.set_vertex(v, vx)                 # …until written back
    assert m.vertex(v).tol == 1e-4
    assert m.count(brep.ShapeType.Vertex) == 1
    assert m.ids(brep.ShapeType.Vertex) == [v]
    assert "vertices=1" in repr(m)


def test_box_to_solid_and_iges(tmp_path):
    m = brep.Model()
    faces = box_faces(m)
    rep = brep.sew(m, faces, tol=1e-6)
    assert len(rep.shells) == 1 and len(rep.merged) == 12
    assert rep.orientable and rep.free_edges == []
    shell = rep.shells[0]
    assert brep.is_closed(m, shell) and brep.is_manifold(m, shell)

    solid = brep.make_solid(m, shell)
    assert brep.shape_type(solid) == brep.ShapeType.Solid
    assert brep.signed_volume(m, shell) == pytest.approx(1.0, abs=1e-12)
    report = brep.check(m, solid)
    assert report.ok() and bool(report)

    assert len(brep.explore_faces(m, solid)) == 6
    assert len(brep.explore_edges(m, solid)) == 12
    assert len(brep.explore_vertices(m, solid)) == 8
    idx = brep.TopologyIndex(m, solid)
    for e in brep.explore_edges(m, solid):
        assert len(idx.faces_of(e)) == 2
        assert len(idx.coedges_of(e)) == 2
    lo, hi = brep.bounding_box(m, solid)
    assert lo == pytest.approx([0., 0., 0.], abs=1e-5) and hi == pytest.approx([1., 1., 1.], abs=1e-5)

    w = gbs.IgesWriter()
    out = w.add_brep(m, solid, "box")
    assert out.faces == 6 and out.approximated_surfaces == 0
    f = tmp_path / "box.igs"
    w.write(str(f))
    de_types = [int(line[:8]) for i, line in enumerate(l for l in f.read_text().splitlines() if len(l) >= 73 and l[72] == "D") if i % 2 == 0]
    assert de_types.count(144) == 6 and de_types.count(142) == 6


def test_face_from_wire_and_pcurves():
    m = brep.Model()
    srf = gbs.BSSurface3d([[-2., -2., 0.], [2., -2., 0.], [-2., 2., 0.], [2., 2., 0.]], [-2., -2., 2., 2.], [-2., -2., 2., 2.], 1, 1)
    pts = [[-1., -1., 0.], [1., -1., 0.], [1., 1., 0.], [-1., 1., 0.]]
    edges = [brep.make_edge(m, pts[i], pts[(i + 1) % 4]) for i in (0, 2, 1, 3)]  # any order
    wire = brep.make_wire(m, edges)
    face = brep.make_face(m, srf, wire)
    assert brep.uv_signed_area(m, wire) == pytest.approx(4.0, abs=1e-9)
    assert all(ce.pcurve is not None for ce in m.wire(wire).coedges)
    assert brep.check(m, face).ok()

    seg = lambda a, b: gbs.BSCurve2d([a, b], [0., 0., 1., 1.], 1)
    tri = [seg([0., .5], [.5, -.5]), seg([.5, -.5], [-.5, -.5]), seg([-.5, -.5], [0., .5])]
    sq = [seg([-1.5, -1.], [1.5, -1.]), seg([1.5, -1.], [1.5, 1.]), seg([1.5, 1.], [-1.5, 1.]), seg([-1.5, 1.], [-1.5, -1.])]
    holed = brep.make_face_from_pcurves(m, srf, sq, [tri])
    wires = m.face(holed).wires
    assert len(wires) == 2
    assert brep.uv_signed_area(m, wires[1]) < 0  # hole turned clockwise
    assert brep.check(m, holed).ok()


def test_natural_faces_and_closure():
    m = brep.Model()
    plane = quad([0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.])
    c = brep.surface_closure(plane, 1e-6)
    assert not (c.closed_u or c.closed_v)
    f = brep.make_face(m, plane)
    assert m.face(f).natural_bounds
    assert len(brep.explore_edges(m, f)) == 4


def test_errors_raise():
    m = brep.Model()
    with pytest.raises(brep.BRepError, match="invalid tolerance"):
        brep.make_vertex(m, [0., 0., 0.], tol=0.)
    e = brep.make_edge(m, [0., 0., 0.], [1., 0., 0.])
    with pytest.raises(brep.BRepError, match="wire not closed"):
        srf = quad([0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.])
        brep.make_face(m, srf, brep.make_wire(m, [e]))
    with pytest.raises(brep.BRepError, match="empty input"):
        brep.sew(m, [])
    # the C++ accessor of a dead id raises too
    with pytest.raises(brep.BRepError):
        m.edge(brep.EdgeId(99))


def test_check_reports_issues():
    m = brep.Model()
    shell = brep.sew(m, box_faces(m)).shells[0]
    solid = brep.make_solid(m, shell)
    # turn the shell inwards: getters return copies, so rebuild the face uses and write back
    sh = m.shell(shell)
    sh.faces = [brep.FaceUse(fu.face, brep.reverse(fu.orient)) for fu in sh.faces]
    m.set_shell(shell, sh)
    report = brep.check(m, solid)
    assert not report.ok()
    assert report.has(brep.Issue.solid_not_outward, solid)
    assert report.count(brep.Issue.solid_not_outward) == 1
    assert "solid not outward" in repr(report.entries[0])


def test_attributes():
    m = brep.Model()
    faces = box_faces(m)
    m.set_name(faces[0], "bottom")
    m.set_external_id(faces[0], 1234)
    assert m.name(faces[0]) == "bottom" and m.name(faces[1]) == ""
    assert m.external_id(faces[0]) == 1234 and m.external_id(faces[1]) is None
    m.set_external_id(faces[0], None)
    assert m.external_id(faces[0]) is None
    assert m.unit_scale == 1.0
    m.unit_scale = 25.4
    assert m.unit_scale == 25.4
    with pytest.raises(brep.BRepError):
        m.unit_scale = 0.0
    with pytest.raises(brep.BRepError):
        m.set_name(brep.FaceId(99), "x")
