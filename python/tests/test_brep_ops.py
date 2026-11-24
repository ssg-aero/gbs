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
    # a sphere touching the box at a point: two solids
    sph = read(m, "sphere.stp")
    touching = brep.copy_shape(m, sph, brep.Placement.frame([1., 1.5, 7.], [1., 0., 0.], [0., 0., -1.], [0., 1., 0.]))
    f = brep.fuse(m, a, touching)
    assert len(f.solids) == 2
    assert math.isclose(volume(m, f.shape), 24. + 36. * math.pi, rel_tol=2e-3)
    # a box resting on an edge on the top of the other: two solids touching along it
    c = math.sqrt(0.5)
    tilted = brep.copy_shape(m, a, brep.Placement.frame([0., 1.5, 4.], [1., 0., 0.], [0., c, c], [0., -c, c]))
    assert brep.section(m, a, tilted).edges_on_faces
    f = brep.fuse(m, a, tilted)
    assert len(f.solids) == 2
    assert math.isclose(volume(m, f.shape), 48., rel_tol=2e-3)
    assert brep.check(m, f.shape).ok()
    # stacked boxes: the sides of both unified into one face each (a box 2 x 3 x 8), or kept apart
    up = brep.copy_shape(m, a, translation([0., 0., 4.]))
    assert len(brep.explore_faces(m, brep.fuse(m, a, up).shape)) == 6
    raw = brep.fuse(m, a, up, unify=False)
    assert len(brep.explore_faces(m, raw.shape)) == 10
    # unified afterwards: the same box, its shell updated
    faces, rep = brep.unify_faces(m, m.shell(brep.explore_shells(m, raw.shape)[0]).faces)
    assert len(faces) == 6 and rep.faces_merged == 4 and len(rep.history) == 4
    assert len(brep.explore_faces(m, raw.shape)) == 6
    assert brep.check(m, raw.shape).ok()
    # pieces of one B-spline surface, each on a surface of its own: one face on their join
    mb = brep.Model()
    ku, kv = [0., 0., 0., 0., 1., 2., 3., 3., 3., 3.], [0., 0., 0., 0., 1., 2., 2., 2., 2.]
    wavy = gbs.BSSurface3d([[float(i), float(j), 0.3 * math.sin(1.1 * i) * math.cos(0.7 * j)] for j in range(5) for i in range(6)], ku, kv, 3, 3)
    left, right = gbs.split_u(wavy, 1.5)
    sewn = brep.sew(mb, [brep.make_face(mb, left), brep.make_face(mb, right)], tol=1e-6)
    faces, rep = brep.unify_faces(mb, mb.shell(sewn.shells[0]).faces)
    assert len(faces) == 1 and rep.surfaces_joined == 1 and rep.faces_merged == 1
    assert brep.check(mb, sewn.shells[0]).ok()
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


def test_pipe():
    m = brep.Model()
    # a disk r 1 around the circle R 4: a torus 2 pi² R r²
    ax = [[4., 0., 0.], [0., 1., 0.], [1., 0., 0.]]
    e = brep.make_edge(m, gbs.EllipseCurve(1., 1., ax))
    disk = brep.make_face(m, gbs.Plane(ax, [-2., 2., -2., 2.]), brep.make_wire(m, [e]))
    path = brep.make_edge(m, gbs.EllipseCurve(4., 4., [[0., 0., 0.], [0., 0., 1.], [1., 0., 0.]]))
    t = brep.make_pipe(m, disk, path)
    assert t.approximated == 0
    assert brep.check(m, t.shape).ok()
    assert math.isclose(volume(m, t.shape), 8. * math.pi**2, rel_tol=2e-3)
    o = brep.PipeOptions(frame=brep.PipeFrame.Frenet, approx_tol=1e-6)
    assert o.frame == brep.PipeFrame.Frenet and o.approx_tol == 1e-6
    assert brep.check(m, brep.make_pipe(m, disk, path, o).shape).ok()
    # Frenet along a segment: no frame
    seg = brep.make_edge(m, [4., 0., 0.], [4., 3., 0.])
    try:
        brep.make_pipe(m, disk, seg, o)
        assert False
    except brep.BRepError:
        pass


def test_offset_solid_and_thick():
    m = brep.Model()
    pts = [[0., 0., 0.], [2., 0., 0.], [2., 1.5, 0.], [0., 1.5, 0.]]
    edges = [brep.make_edge(m, pts[i], pts[(i + 1) % 4]) for i in range(4)]
    base = brep.make_face(m, gbs.Plane([[0., 0., 0.], [0., 0., 1.], [1., 0., 0.]], [-5., 5., -5., 5.]), brep.make_wire(m, edges))
    box = brep.explore_solids(m, brep.make_prism(m, base, [0., 0., 1.]).shape)[0]
    big = brep.offset_solid(m, box, 0.25)
    assert brep.check(m, big).ok()
    assert math.isclose(brep.signed_volume(m, m.solid(big).outer, 64), 2.5 * 2. * 1.5, rel_tol=1e-9)
    top = [f for f in brep.explore_faces(m, box) if brep.mesh(m, f).vertices[:, 2].min() > 1. - 1e-12]
    assert len(top) == 1
    cup = brep.make_thick(m, box, top, 0.1, brep.OffsetOptions(tol=1e-7))
    assert brep.check(m, cup).ok()
    assert math.isclose(brep.mesh(m, cup).volume(), 3. - 1.8 * 1.3 * 0.9, rel_tol=1e-9)
    try:
        brep.offset_solid(m, box, -0.8)
        assert False
    except brep.BRepError:
        pass


def test_fillet_and_chamfer_geometry():
    m = brep.Model()
    pts = [[0., 0., 0.], [2., 0., 0.], [2., 1.5, 0.], [0., 1.5, 0.]]
    edges = [brep.make_edge(m, pts[i], pts[(i + 1) % 4]) for i in range(4)]
    base = brep.make_face(m, gbs.Plane([[0., 0., 0.], [0., 0., 1.], [1., 0., 0.]], [-5., 5., -5., 5.]), brep.make_wire(m, edges))
    box = brep.explore_solids(m, brep.make_prism(m, base, [0., 0., 1.]).shape)[0]
    # the top edge along y at x = 2
    edge = next(e for e in brep.explore_edges(m, box)
                if all(abs(m.edge(e).curve.value(t)[0] - 2.) < 1e-12 and abs(m.edge(e).curve.value(t)[2] - 1.) < 1e-12
                       for t in (m.edge(e).u1, m.edge(e).u2)))
    r = 0.3
    f = brep.fillet_geometry(m, box, edge, r)
    assert f.convex and f.exact and isinstance(f.surface, gbs.CylindricalSurface)
    c = f.spine.value(f.t1)
    assert math.isclose(c[0], 2. - r, abs_tol=1e-9) and math.isclose(c[2], 1. - r, abs_tol=1e-9)
    for ct in f.contacts:
        assert math.dist(ct.curve.value(f.t1), c) == pytest.approx(r, abs=1e-9)
    # 0.2 on the top face, 0.3 on the side
    top = next(f for f in brep.explore_faces(m, box)
               if all(abs(m.vertex(v).pnt[2] - 1.) < 1e-12 for v in brep.explore_vertices(m, f)))
    ch = brep.chamfer_geometry(m, box, edge, top, 0.2, 0.3, brep.FilletOptions(extension=0.))
    assert isinstance(ch.surface, gbs.Plane) and ch.exact
    assert ch.contacts[0].face == top
    assert (ch.a, ch.b) == (ch.t1, ch.t2)
    p = ch.spine.value(ch.t1)
    assert [math.dist(ct.curve.value(ch.t1), p) for ct in ch.contacts] == pytest.approx([0.2, 0.3], abs=1e-9)
    assert ch.contacts[0].curve.value(ch.t1)[2] == pytest.approx(1., abs=1e-12)
    with pytest.raises(brep.BRepError):
        brep.fillet_geometry(m, box, edge, 1.2)  # wider than the side face
    with pytest.raises(brep.BRepError):
        brep.fillet_geometry(m, box, edge, 0.1, brep.FilletOptions(extension=-1.))


def test_fillet_and_chamfer():
    m = brep.Model()
    pts = [[0., 0., 0.], [2., 0., 0.], [2., 1.5, 0.], [0., 1.5, 0.]]
    edges = [brep.make_edge(m, pts[i], pts[(i + 1) % 4]) for i in range(4)]
    base = brep.make_face(m, gbs.Plane([[0., 0., 0.], [0., 0., 1.], [1., 0., 0.]], [-5., 5., -5., 5.]), brep.make_wire(m, edges))
    box = brep.explore_solids(m, brep.make_prism(m, base, [0., 0., 1.]).shape)[0]

    def top_edge(x):
        return next(e for e in brep.explore_edges(m, box)
                    if all(abs(m.edge(e).curve.value(t)[0] - x) < 1e-12 and abs(m.edge(e).curve.value(t)[2] - 1.) < 1e-12
                           for t in (m.edge(e).u1, m.edge(e).u2)))

    def volume(s):
        return brep.signed_volume(m, m.solid(s).outer, 64)

    r = 0.3
    # the two top edges along y: each takes off the corner (1 - pi / 4) r^2 along 1.5
    f = brep.fillet(m, box, [top_edge(0.), top_edge(2.)], r)
    assert f != box and brep.check(m, f).ok()
    assert len(brep.explore_faces(m, f)) == 8
    assert volume(f) == pytest.approx(3. - 2. * (1. - math.pi / 4.) * r * r * 1.5, abs=1e-9)
    assert volume(box) == pytest.approx(3., abs=1e-12)  # unchanged
    c = brep.chamfer(m, box, [top_edge(2.)], 0.2)
    assert volume(c) == pytest.approx(3. - 0.02 * 1.5, abs=1e-9)
    top = next(f for f in brep.explore_faces(m, box)
               if all(abs(m.vertex(v).pnt[2] - 1.) < 1e-12 for v in brep.explore_vertices(m, f)))
    c2 = brep.chamfer(m, box, [top_edge(2.)], top, 0.2, 0.3)
    assert volume(c2) == pytest.approx(3. - 0.03 * 1.5, abs=1e-9)
    with pytest.raises(brep.BRepError):
        brep.fillet(m, box, [top_edge(2.)], 1.2)  # wider than the side face
    # variable chamfers (stage 7): d1 linear from 0.1 to 0.4 on the top, d2 = 0.2 on the side: the
    # triangle d1 d2 / 2 taken off at each abscissa; the same distance on both faces; distance
    # and angle (d2 = d1 tan(angle) between planes at right angles)
    e = top_edge(2.)
    c3 = brep.chamfer(m, box, [e], top, brep.BlendLaw.linear(e, 0.1, 0.4), brep.BlendLaw.constant(e, 0.2))
    assert volume(c3) == pytest.approx(3. - 0.25 * 0.2 / 2. * 1.5, abs=1e-9)
    c4 = brep.chamfer(m, box, [e], brep.BlendLaw([brep.BlendPoint(e, 0.0, 0.1), brep.BlendPoint(e, 1.0, 0.3)], brep.BlendVariation.Linear))
    assert volume(c4) == pytest.approx(3. - 1.5 * (0.01 + 0.1 * 0.2 + 0.04 / 3.) / 2., abs=1e-9)
    c5 = brep.chamfer_angle(m, box, [e], top, 0.3, math.pi / 6.)
    assert volume(c5) == pytest.approx(3. - 1.5 * 0.3 * 0.3 * math.tan(math.pi / 6.) / 2., abs=1e-9)
    g = brep.chamfer_geometry(m, box, e, top, brep.BlendLaw([brep.BlendPoint(e, 0., 0.2), brep.BlendPoint(e, 0.5, 0.4), brep.BlendPoint(e, 1., 0.25)],
                                                            brep.RadiusVariation.Linear), brep.BlendLaw.constant(e, 0.2))
    assert len(g.pieces) == 2 and len(g.creases) == 1
    assert brep.RadiusLaw is brep.BlendLaw and brep.RadiusPoint(e, 0.5, 0.2).radius == 0.2
    # the keywords of stage 5 still taken
    assert brep.RadiusPoint(e, at=0.5, radius=0.2).value == 0.2
    assert [p.value for p in brep.RadiusLaw.linear(e, r1=0.1, r2=0.3).points] == [0.1, 0.3]
    assert [p.value for p in brep.BlendLaw.linear(e, v1=0.1, v2=0.3).points] == [0.1, 0.3]
    # conic sections (stage 7, PR 2): rho of the circle gives the circle; rho = 0.5 the parabola,
    # its area r^2 / 6 between it and the edge; a chord c at right angles: the radius c / sqrt(2)
    o = brep.FilletOptions()
    o.section = brep.BlendSection.conic(math.cos(math.pi / 4) / (1 + math.cos(math.pi / 4)))
    assert volume(brep.fillet(m, box, [e], 0.3, o)) == pytest.approx(volume(brep.fillet(m, box, [e], 0.3)), abs=1e-9)
    o.section = brep.BlendSection.conic(0.5)
    assert volume(brep.fillet(m, box, [e], 0.3, o)) == pytest.approx(3. - 1.5 * 0.09 / 6., abs=1e-9)
    assert not brep.fillet_geometry(m, box, e, 0.3, o).circular
    c6 = brep.fillet_chordal(m, box, [e], 0.4)
    # between two planes: the exact cylinder of radius c / (2 sin 45) (it was a rational B-spline on it)
    g6 = brep.fillet_chordal_geometry(m, box, e, 0.4)
    assert g6.exact and not g6.variable and type(g6.surface).__name__ == "CylindricalSurface"
    assert volume(c6) == pytest.approx(3. - 1.5 * (1. - math.pi / 4.) * 0.08, abs=1e-7)
    # the G2 section (stage 7, PR 3): the contacts at 1.305 r from the edge, the area of the
    # quartic 0.232084 r^2 (the arc: 1 - pi / 4 = 0.2146); not with a constant chord
    o.section = brep.BlendSection.g2()
    assert o.section.kind == brep.SectionKind.G2
    g2 = brep.fillet_geometry(m, box, e, 0.3, o)
    assert g2.exact and not g2.circular
    assert volume(brep.fillet(m, box, [e], 0.3, o)) == pytest.approx(3. - 1.5 * 0.09 * 0.23208374949396718, abs=1e-8)
    with pytest.raises(brep.BRepError, match="constant chord with a G2 section"):
        brep.fillet_chordal(m, box, [e], 0.4, o)


def test_fillet_faces():
    # stage 8, PR 1: the box 2 x 1.5 x 1, its top edge at x = 2 chamfered at 0.1; the fillet
    # between the top and the side over the chamfer is that of the edge before the chamfer
    m = brep.Model()
    pts = [[0., 0., 0.], [2., 0., 0.], [2., 1.5, 0.], [0., 1.5, 0.]]
    edges = [brep.make_edge(m, pts[i], pts[(i + 1) % 4]) for i in range(4)]
    base = brep.make_face(m, gbs.Plane([[0., 0., 0.], [0., 0., 1.], [1., 0., 0.]], [-5., 5., -5., 5.]), brep.make_wire(m, edges))
    box = brep.explore_solids(m, brep.make_prism(m, base, [0., 0., 1.]).shape)[0]
    e = next(e for e in brep.explore_edges(m, box)
             if all(abs(m.edge(e).curve.value(t)[0] - 2.) < 1e-12 and abs(m.edge(e).curve.value(t)[2] - 1.) < 1e-12 for t in (m.edge(e).u1, m.edge(e).u2)))
    ch = brep.chamfer(m, box, [e], 0.1)

    def face(keep):
        return next(f for f in brep.explore_faces(m, ch) if all(keep(m.vertex(v).pnt) for v in brep.explore_vertices(m, f)))

    top = face(lambda p: abs(p[2] - 1.) < 1e-12)
    side = face(lambda p: abs(p[0] - 2.) < 1e-12)
    r = 0.3
    f = brep.fillet_faces(m, ch, [top], [side], r)
    assert brep.check(m, f).ok()
    assert brep.signed_volume(m, m.solid(f).outer, 64) == pytest.approx(3. - 1.5 * (1. - math.pi / 4.) * r * r, abs=1e-9)
    with pytest.raises(brep.BRepError, match="spills over"):
        brep.fillet_faces(m, ch, [top], [side], 0.13)
    with pytest.raises(brep.BRepError, match="on both sides"):
        brep.fillet_faces(m, ch, [top], [top], 0.03)


def _planar(m, pts):
    """The planar face bounded by the polygon pts (counter-clockwise about its normal)."""
    a, b, c = (pts[0], pts[1], pts[2])
    u = [b[i] - a[i] for i in range(3)]
    v = [c[i] - a[i] for i in range(3)]
    n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
    ln, lu = math.sqrt(sum(x * x for x in n)), math.sqrt(sum(x * x for x in u))
    plane = gbs.Plane([list(a), [x / ln for x in n], [x / lu for x in u]], [-10., 10., -10., 10.])
    edges = [brep.make_edge(m, list(pts[i]), list(pts[(i + 1) % len(pts)])) for i in range(len(pts))]
    return brep.make_face(m, plane, brep.make_wire(m, edges))


def _u_channel(m):
    """The U of stage 8, PR 2: wall y = -0.6 and floor z = 0, along x over 2."""
    yz = [(-1., -0.6), (1., -0.6), (1., 1.), (0.6, 1.), (0.6, 0.), (-0.6, 0.), (-0.6, 1.), (-1., 1.)]
    return brep.explore_solids(m, brep.make_prism(m, _planar(m, [(0., y, z) for y, z in yz]), [2., 0., 0.]).shape)[0]


def _tapered_u(m):
    """The U with a face tapering from 0.12 wide at x = 0 to a point on the corner at x = 1."""
    a, b, c, d = (0., -0.65, -0.05), (0., -0.65, 0.17), (0., -0.43, -0.05), (0.22 / 0.12, -0.65, -0.05)
    faces = [_planar(m, q) for q in ((a, c, b), (a, b, d), (a, d, c), (b, c, d))]
    tet = brep.make_solid(m, brep.sew(m, faces, 1e-9, 16).shells[0], [])
    return brep.explore_solids(m, brep.fuse(m, _u_channel(m), tet).shape)[0]


def _u_faces(m, s):
    def face(keep):
        return next(f for f in brep.explore_faces(m, s) if all(keep(m.vertex(v).pnt) for v in brep.explore_vertices(m, f)))
    wall = face(lambda p: abs(p[1] + 0.6) < 1e-12 and p[2] > -1e-12)
    floor = face(lambda p: abs(p[2]) < 1e-12 and abs(p[1]) < 0.6 + 1e-12)
    return wall, floor


def _corner(m, s):
    return next(e for e in brep.explore_edges(m, s)
                if all(abs(m.edge(e).curve.value(t)[1] + 0.6) < 1e-12 and abs(m.edge(e).curve.value(t)[2]) < 1e-12 for t in (m.edge(e).u1, m.edge(e).u2)))


def test_fillet_faces_partial_strip():
    # stage 8, PR 2: a U with a face tapering to a point on the corner of its wall and floor, which
    # meet along the rest of it; the face is under the continuation of the corner (it touches both,
    # its edges with them end on its tip, a vertex of the corner): swallowed, the fillet that of the
    # corner of the sharp U; at r = 0.2 it spills over the face (0.085 wide from the corner, the arc
    # 0.083 from it); fillet of the corner alone does not cross the face
    m = brep.Model()
    s = _tapered_u(m)
    assert brep.check(m, s).ok() and len(brep.explore_faces(m, s)) == 11
    wall, floor = _u_faces(m, s)
    sharp = _u_channel(m)
    for r in (0.3, 0.5):
        f = brep.fillet_faces(m, s, [wall], [floor], r)
        assert brep.check(m, f).ok()
        ref = brep.fillet(m, sharp, [_corner(m, sharp)], r)
        assert len(brep.explore_faces(m, f)) == len(brep.explore_faces(m, ref))
        assert brep.signed_volume(m, m.solid(f).outer, 512) == pytest.approx(brep.signed_volume(m, m.solid(ref).outer, 512), abs=1e-9)
    assert brep.signed_volume(m, m.solid(f).outer, 512) == pytest.approx(4.107300918992928, abs=1e-8)
    with pytest.raises(brep.BRepError, match="spills over"):
        brep.fillet_faces(m, s, [wall], [floor], 0.2)
    with pytest.raises(brep.BRepError, match="fillet_faces swallows them"):
        brep.fillet(m, s, [_corner(m, s)], 0.3)
    assert "under the common edge" in brep.fillet_faces.__doc__


def test_fillet_faces_law_by_position():
    # stage 8, PR 2: a radius from 0.2 to 0.4 along the virtual edge of the tapered U, its points
    # given by position: the solid of the same law on the corner of the sharp U
    m = brep.Model()
    s = _tapered_u(m)
    wall, floor = _u_faces(m, s)
    p = brep.BlendPoint.at_point([0., -0.55, 0.05], 0.4)
    assert p.where is not None and p.value == 0.4
    law = brep.RadiusLaw([p, brep.BlendPoint.at_point([2.1, -0.6, 0.], 0.5)], brep.BlendVariation.Linear)
    f = brep.fillet_faces(m, s, [wall], [floor], law)
    assert brep.check(m, f).ok()
    sharp = _u_channel(m)
    e = _corner(m, sharp)
    at0 = 0. if m.edge(e).curve.value(m.edge(e).u1)[0] < 1. else 1.
    ref = brep.fillet(m, sharp, [e], brep.RadiusLaw([brep.BlendPoint(e, at0, 0.4), brep.BlendPoint(e, 1. - at0, 0.5)], brep.BlendVariation.Linear))
    assert brep.signed_volume(m, m.solid(f).outer, 512) == pytest.approx(brep.signed_volume(m, m.solid(ref).outer, 512), abs=1e-9)
    with pytest.raises(brep.BRepError, match="away from the edges"):
        brep.fillet_faces(m, s, [wall], [floor], brep.RadiusLaw([brep.BlendPoint.at_point([1., 20., 20.], 0.4)]))
    # the geometry of the corner alone with a law of one value given by position: that of the value
    g = brep.fillet_geometry(m, sharp, e, brep.RadiusLaw([brep.BlendPoint.at_point([1., -0.55, 0.05], 0.3)]))
    g0 = brep.fillet_geometry(m, sharp, e, 0.3)
    for k in range(2):
        for w in (0.1, 0.5, 0.9):
            t = g.a + (g.b - g.a) * w
            p, q = g.contacts[k].curve.value(t), g0.contacts[k].curve.value(t)
            assert math.dist(p, q) < 1e-12


def test_fillet_setbacks_option():
    # a tuple, set by assigning a list: appending to it fails instead of changing a copy
    o = brep.FilletOptions()
    assert o.setbacks == ()
    o.setbacks = [brep.Setback(brep.VertexId(), brep.EdgeId(), 0.4)]
    assert len(o.setbacks) == 1 and o.setbacks[0].distance == 0.4
    with pytest.raises(AttributeError):
        o.setbacks.append(brep.Setback(brep.VertexId(), brep.EdgeId(), 0.5))
    assert len(o.setbacks) == 1


def test_variable_fillet():
    import sys

    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "docs", "sources", "examples"))
    import gbs_examples as ex

    # a box edge from 0.1 to 0.3 (linear), then with a point in between: a crease (peak) or a
    # lune (hollow) of the ball where the law kinks
    m = brep.Model()
    s = ex.box(m, size=(2.0, 1.0, 1.0))
    e = ex.edge_near(m, s, (1.0, 0.0, 1.0))
    g = brep.fillet_geometry(m, s, e, brep.RadiusLaw.linear(e, 0.1, 0.3))
    assert not g.exact and not g.pieces and g.surface is not None
    peak = brep.RadiusLaw([brep.RadiusPoint(e, 0.0, 0.1), brep.RadiusPoint(e, 0.5, 0.3), brep.RadiusPoint(e, 1.0, 0.1)], brep.RadiusVariation.Linear)
    g = brep.fillet_geometry(m, s, e, peak)
    assert len(g.pieces) == 2 and len(g.creases) == 1 and not g.lunes and g.surface is None
    hollow = brep.RadiusLaw([brep.RadiusPoint(e, 0.0, 0.3), brep.RadiusPoint(e, 0.5, 0.1), brep.RadiusPoint(e, 1.0, 0.3)], brep.RadiusVariation.Linear)
    g = brep.fillet_geometry(m, s, e, hollow)
    assert len(g.pieces) == 2 and len(g.lunes) == 1 and not g.creases
    r = brep.fillet(m, s, [e], hollow)
    assert brep.check(m, r).ok()
    assert len(brep.explore_faces(m, r)) == 9
    # a point on an edge not filleted
    with pytest.raises(brep.BRepError):
        brep.fillet(m, s, [e], brep.RadiusLaw([brep.RadiusPoint(ex.edge_near(m, s, (1.0, 1.0, 1.0)), 0.0, 0.1)]))
