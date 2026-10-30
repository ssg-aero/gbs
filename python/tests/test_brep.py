"""Python bindings of the native BREP core (gbs.brep)."""
import copy
import math
import pickle

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


def test_ordered_wire_and_face_use():
    from math import pi
    m = brep.Model()
    R, h = 1.5, 2.
    ax = [[0., 0., 0.], [0., 0., 1.], [1., 0., 0.]]
    cyl = gbs.build_cylinder(R, ax, 0., h)
    vb = brep.make_vertex(m, [R, 0., 0.])
    vt = brep.make_vertex(m, [R, 0., h])

    def circle(z, v):
        c = gbs.build_circle_arc(R, 0., 2 * pi, [[0., 0., z], [0., 0., 1.], [1., 0., 0.]])
        return brep.make_edge(m, c, 0., 2 * pi, v, v)

    bottom, top = circle(0., vb), circle(h, vt)
    seam = brep.make_edge(m, vb, vt)
    F, Rv = brep.Orientation.Forward, brep.Orientation.Reversed
    w = brep.make_wire_ordered(m, [(bottom, F), (seam, F), (top, Rv), (seam, Rv)])
    use = brep.make_face_use(m, cyl, w)
    assert use.orient == F
    assert brep.check(m, use.face).ok()
    with pytest.raises(brep.BRepError):
        brep.make_wire_ordered(m, [(seam, F), (seam, F)])


def test_split_edge():
    m = brep.Model()
    a = brep.make_vertex(m, [0., 0., 0.])
    b = brep.make_vertex(m, [2., 0., 0.])
    e = brep.make_edge(m, a, b)
    w = brep.make_wire_ordered(m, [(e, brep.Orientation.Forward)])
    first, second, v = brep.split_edge(m, e, 1.)
    assert first == e
    assert m.vertex(v).pnt == pytest.approx([1., 0., 0.])
    assert len(m.wire(w).coedges) == 2
    with pytest.raises(brep.BRepError):
        brep.split_edge(m, e, 5.)


def box_solid(m):
    return brep.make_solid(m, brep.sew(m, box_faces(m)).shells[0])


def test_errors_carry_code_and_shapes():
    m = brep.Model()
    a = brep.make_edge(m, [0., 0., 0.], [1., 0., 0.])
    b = brep.make_edge(m, [5., 0., 0.], [6., 0., 0.])
    with pytest.raises(brep.BRepError) as err:
        brep.make_wire(m, [a, b])
    assert err.value.code == brep.BuildErrc.disconnected
    with pytest.raises(brep.BRepError) as err:
        m.edge(brep.EdgeId(99))
    assert err.value.code == brep.BuildErrc.invalid_id and err.value.shapes == [brep.EdgeId(99)]
    # an operation keeps the reason of the error it catches: an operand with an erased vertex
    s1, s2 = box_solid(m), box_solid(m)
    v = brep.explore_vertices(m, s2)[0]
    m.erase(v)
    with pytest.raises(brep.BRepError) as err:
        brep.fuse(m, s1, s2)
    assert err.value.code == brep.BuildErrc.invalid_id and err.value.shapes == [v]
    # a BRepError raised from Python has the class attributes
    assert brep.BRepError("x").code is None and brep.BRepError("x").shapes == ()


def test_ids_compare_pickle_and_copy():
    m = brep.Model()
    s = box_solid(m)
    v = brep.VertexId(0)
    assert v != None and not (v == brep.EdgeId(0)) and v != brep.EdgeId(0)  # noqa: E711
    assert brep.EdgeId(0) in [brep.VertexId(0), brep.EdgeId(0)]
    with pytest.raises(TypeError):
        brep.VertexId(1) < brep.EdgeId(0)
    # entries on edges and wires compared with a solid
    e = brep.explore_edges(m, s)[0]
    x = m.edge(e)
    x.tol = -1.
    m.set_edge(e, x)
    entries = brep.check(m, s).entries
    assert entries and [en for en in entries if en.shape == s] == []
    assert pickle.loads(pickle.dumps(brep.FaceId(3))) == brep.FaceId(3)
    assert copy.copy(brep.ShellId(2)) == brep.ShellId(2)


def test_setters_refuse_missing_geometry():
    # an edge without curve or a face without surface crashed the operations (mesh, section, booleans...)
    m = brep.Model()
    s = box_solid(m)
    e = brep.explore_edges(m, s)[0]
    x = m.edge(e)
    x.curve = None
    with pytest.raises(brep.BRepError) as err:
        m.set_edge(e, x)
    assert err.value.code == brep.BuildErrc.null_curve and err.value.shapes == [e]
    assert m.edge(e).curve is not None
    f = brep.explore_faces(m, s)[0]
    y = m.face(f)
    y.surface = None
    with pytest.raises(brep.BRepError) as err:
        m.set_face(f, y)
    assert err.value.code == brep.BuildErrc.null_surface
    assert brep.mesh(m, s).volume() == pytest.approx(1.0)
    # a degenerate edge has no curve
    d = brep.make_degenerate_edge(m, brep.make_vertex(m, [0., 0., 0.]))
    m.set_edge(d, m.edge(d))


def test_topology_index_skips_dead_vertices():
    m = brep.Model()
    s = box_solid(m)
    e = brep.explore_edges(m, s)[0]
    x = m.edge(e)
    x.v1 = brep.VertexId(50_000_000)
    m.set_edge(e, x)
    idx = brep.TopologyIndex(m, s)
    assert len(idx.faces_of(e)) == 2 and idx.edges_of(brep.VertexId(50_000_000)) == []


def test_face_options_checked():
    m = brep.Model()
    srf = quad([0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.])
    pts = [[.2, .2, 0.], [.8, .2, 0.], [.8, .8, 0.]]
    w = brep.make_wire(m, [brep.make_edge(m, pts[i], pts[(i + 1) % 3]) for i in range(3)])
    with pytest.raises(brep.BRepError) as err:
        brep.make_face(m, srf, w, options=brep.MakeFaceOptions(pcurve_degree=0))
    assert err.value.code == brep.BuildErrc.invalid_options
    assert m.count(brep.ShapeType.Face) == 0
    # the default options read as Python in the signatures
    assert "options: pygbs.gbs.brep.MakeFaceOptions = MakeFaceOptions()" in brep.make_face.__doc__


def test_iges_coedge_without_pcurve():
    m = brep.Model()
    s = box_solid(m)
    f = brep.explore_faces(m, s)[3]
    w = m.face(f).wires[0]
    x = m.wire(w)
    c = x.coedges[1]
    c.pcurve = None
    x.coedges = x.coedges[:1] + [c] + x.coedges[2:]
    m.set_wire(w, x)
    with pytest.raises(brep.BRepError) as err:
        gbs.IgesWriter().add_brep(m, s)
    assert err.value.code == brep.BuildErrc.missing_pcurve and err.value.shapes == [f]


def test_copies_and_repr():
    m = brep.Model()
    s = box_solid(m)
    c = copy.deepcopy(m)
    c.erase(s)
    assert m.alive(s) and not c.alive(s)
    assert copy.copy(m).count(brep.ShapeType.Face) == 6
    e = brep.explore_edges(m, s)[0]
    x = copy.copy(m.edge(e))
    assert x.u2 == m.edge(e).u2
    assert repr(m.edge(e)).startswith("Edge(curve=")
    assert repr(m.vertex(brep.VertexId(0))) == "Vertex(pnt=[0.0, 0.0, 0.0], tol=1e-06)"
    assert repr(m.shell(m.solid(s).outer)) == "Shell(faces=6, closed=True)"
    assert repr(brep.FaceUse(brep.FaceId(2))) == "FaceUse(face=FaceId(2), orient=Forward)"


def test_long_operations_release_the_gil():
    import threading
    import time

    m = brep.Model()
    s = box_solid(m)
    count = [0]
    stop = threading.Event()

    def spin():
        while not stop.is_set():
            count[0] += 1
            time.sleep(0)

    t = threading.Thread(target=spin)
    t.start()
    try:
        time.sleep(0.01)
        before = count[0]
        brep.mesh(m, s, brep.MeshOptions(max_edge=0.02))  # about 0.1 s, 75 000 triangles
        during = count[0] - before
    finally:
        stop.set()
        t.join()
    assert during > 1000  # about 21 when the call keeps the GIL, 600 000 when it releases it
