"""Shared helpers of the GBS examples: primitive solids of the native BREP, the cases of the
booleans, their meshes and their PyVista rendering.

The primitives are built with the BREP API itself (planar faces, prisms, revolutions), placed by
rigid placements; each face of a result gets its own color, so that the splits of the faces are
seen, and the edges of the BREP are drawn as lines: the sharp ones, not the seams where a periodic
surface closes on itself (they part no two surfaces) nor the degenerate edges of the poles.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

import numpy as np
import pygbs.gbs as gbs

brep = gbs.brep
F, R = brep.Orientation.Forward, brep.Orientation.Reversed


# ---------------------------------------------------------------------------- primitives
def _axes(origin, z):
    """A direct frame at `origin` whose z axis is along `z`: (origin, x, y, z) as lists, x taken
    normal to z and to the first of the x and y axes that is not almost along z."""
    z = np.asarray(z, float)
    z = z / np.linalg.norm(z)
    a = np.array([1.0, 0.0, 0.0]) if abs(z[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    x = np.cross(a, z)
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    return list(map(float, origin)), list(x), list(y), list(z)


def _frame(origin, z):
    """A placement putting the z axis along `z`, at `origin`."""
    return brep.Placement.frame(*_axes(origin, z))


def _ax2(origin, z):
    """The same frame as [origin, z, x], the form of the elementary curves and surfaces."""
    o, x, _, z = _axes(origin, z)
    return [o, z, x]


def _placed(m, shape, origin, axis):
    return brep.copy_shape(m, shape, _frame(origin, axis))


def planar_face(m, points):
    """The planar face bounded by the polygon `points` (counter-clockwise around its normal)."""
    p = [np.asarray(q, float) for q in points]
    n = np.cross(p[1] - p[0], p[2] - p[0])
    n /= np.linalg.norm(n)
    x = (p[1] - p[0]) / np.linalg.norm(p[1] - p[0])
    size = 2.0 * max(np.linalg.norm(q - p[0]) for q in p)
    plane = gbs.Plane([list(p[0]), list(n), list(x)], [-size, size, -size, size])
    v = [brep.make_vertex(m, list(q)) for q in p]
    edges = [brep.make_edge(m, v[i], v[(i + 1) % len(v)]) for i in range(len(v))]
    return brep.make_face(m, plane, brep.make_wire_ordered(m, [(e, F) for e in edges]))


def circle(radius, ax, t1=0.0, t2=2 * math.pi, exact=True):
    """The arc of circle of `radius` in the frame `ax` = [center, normal, x direction], from the
    angle t1 to t2: exact (`EllipseCurve`, swept into exact cylinders, cones, spheres, tori) or
    its NURBS (swept into NURBS surfaces)."""
    if exact:
        return gbs.EllipseCurve(radius, radius, ax, t1, t2)
    return gbs.build_circle_arc(radius, t1, t2, ax)


def disk_face(m, radius, center=(0.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0), exact=True):
    """The planar disk of `radius` at `center`, normal to `normal`."""
    ax = _ax2(center, normal)
    c = circle(radius, ax, exact=exact)
    v = brep.make_vertex(m, list(np.asarray(c.value(0.0))))
    e = brep.make_edge(m, c, 0.0, 2 * math.pi, v, v)
    plane = gbs.Plane(ax, [-2 * radius, 2 * radius, -2 * radius, 2 * radius])
    return brep.make_face(m, plane, brep.make_wire_ordered(m, [(e, F)]))


def solid_of(m, sweep):
    return brep.explore_solids(m, sweep.shape)[0]


def box(m, origin=(0.0, 0.0, 0.0), size=(1.0, 1.0, 1.0)):
    """The box [origin, origin + size]."""
    ox, oy, oz = origin
    dx, dy, dz = size
    f = planar_face(m, [[ox, oy, oz], [ox + dx, oy, oz], [ox + dx, oy + dy, oz], [ox, oy + dy, oz]])
    return solid_of(m, brep.make_prism(m, f, [0.0, 0.0, dz]))


def cylinder(m, radius, height, base=(0.0, 0.0, 0.0), axis=(0.0, 0.0, 1.0), exact=True):
    """The cylinder of `radius` and `height` from `base` along `axis`."""
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    f = disk_face(m, radius, base, axis, exact)
    return solid_of(m, brep.make_prism(m, f, list(axis * height)))


def _revolved(m, profile_xz, origin, axis):
    """The solid of revolution of the polygon `profile_xz` (in the (x, z) half plane) around z,
    then placed at `origin` with z along `axis`."""
    pts = [[x, 0.0, z] for x, z in profile_xz]
    f = planar_face(m, pts)
    s = solid_of(m, brep.make_revolution(m, f, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]))
    return brep.explore_solids(m, _placed(m, s, origin, axis))[0]


def cone(m, radius, height, base=(0.0, 0.0, 0.0), axis=(0.0, 0.0, 1.0)):
    """The cone of base `radius`, apex at `height` along `axis`."""
    return _revolved(m, [(0.0, 0.0), (radius, 0.0), (0.0, height)], base, axis)


def sphere(m, radius, center=(0.0, 0.0, 0.0), axis=(0.0, 0.0, 1.0), exact=True):
    """The sphere of `radius` at `center`, its poles along `axis` (`exact`: a spherical surface,
    else its NURBS)."""
    # a half disk in the (x, z) plane: the diameter on the axis, the half circle
    ax = [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]  # angle 0 at +z, pi / 2 at +x, pi at -z
    arc = circle(radius, ax, 0.0, math.pi, exact)
    top, bottom = brep.make_vertex(m, [0.0, 0.0, radius]), brep.make_vertex(m, [0.0, 0.0, -radius])
    e_arc = brep.make_edge(m, arc, 0.0, math.pi, top, bottom)
    e_axis = brep.make_edge(m, bottom, top)
    plane = gbs.Plane([[0.0, 0.0, 0.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], [-2 * radius, 2 * radius, -2 * radius, 2 * radius])
    f = brep.make_face(m, plane, brep.make_wire_ordered(m, [(e_arc, F), (e_axis, F)]))
    s = solid_of(m, brep.make_revolution(m, f, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]))
    return brep.explore_solids(m, _placed(m, s, center, axis))[0]


def torus(m, major, minor, center=(0.0, 0.0, 0.0), axis=(0.0, 0.0, 1.0), exact=True):
    """The torus of radii `major` and `minor` at `center` around `axis`."""
    f = disk_face(m, minor, (major, 0.0, 0.0), (0.0, 1.0, 0.0), exact)
    s = solid_of(m, brep.make_revolution(m, f, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]))
    return brep.explore_solids(m, _placed(m, s, center, axis))[0]


# ---------------------------------------------------------------------------- boolean cases
@dataclass
class Case:
    """Two operands built in a model, and what is known of the results."""

    name: str
    title: str
    build: Callable  # (model) -> (a, b)
    note: str = ""


def _box_cylinder(m):
    return box(m, (0, 0, 0), (2, 3, 4)), cylinder(m, 0.8, 6, (1, 1.5, -1))


def _box_sphere(m):
    return box(m, (0, 0, 0), (2, 3, 4)), sphere(m, 1.5, (2, 3, 4), (1, 1, 1))


def _crossed_cylinders(m):
    return cylinder(m, 1.0, 6, (0, 0, -3)), cylinder(m, 0.6, 6, (-3, 0, 0), (1, 0, 0))


def _steinmetz(m):
    return cylinder(m, 1.0, 6, (0, 0, -3)), cylinder(m, 1.0, 6, (-3, 0, 0), (1, 0, 0))


def _boxes(m):
    return box(m, (0, 0, 0), (2, 3, 4)), box(m, (1, 1, 1), (2, 3, 4))


def _coplanar_boxes(m):
    return box(m, (0, 0, 0), (2, 3, 4)), box(m, (1, 1, 0), (2, 3, 4))


def _stacked_boxes(m):
    return box(m, (0, 0, 0), (2, 3, 4)), box(m, (1, 0, 4), (2, 3, 4))


def _coaxial_cylinders(m):
    return cylinder(m, 1.0, 2.0), cylinder(m, 1.0, 2.0, (0, 0, 1))


def _sphere_touching(m):
    return box(m, (0, 0, 0), (2, 3, 4)), sphere(m, 1.5, (1, 1.5, 5.5), (1, 1, 1))


def _cylinder_lying(m):
    return box(m, (0, 0, 0), (2, 3, 4)), cylinder(m, 1.0, 5, (1, -1, 5), (0, 1, 0))


def _cone_through_apex(m):
    return box(m, (-2, -2, 1), (4, 4, 4)), cone(m, 1.5, 2.0)


def _torus_box(m):
    return box(m, (-1, -5, -2), (2, 10, 4)), torus(m, 3.0, 1.0)


def _nurbs_sphere_touching(m):
    return box(m, (0, 0, 0), (2, 3, 4)), sphere(m, 1.5, (1, 1.5, 5.5), (1, 1, 1), exact=False)


def _sphere_cylinder(m):
    return cylinder(m, 1.0, 4, (0, 0, -2)), sphere(m, 1.2, (1.0, 0, 0), (0, 1, 0))


CASES = [
    Case("boxes", "Two boxes", _boxes),
    Case("box_cylinder", "Box and cylinder", _box_cylinder),
    Case("box_sphere", "Box and sphere at its corner", _box_sphere),
    Case("crossed_cylinders", "Crossed cylinders", _crossed_cylinders),
    Case("steinmetz", "Equal crossed cylinders (branch points)", _steinmetz),
    Case("coplanar_boxes", "Coplanar boxes", _coplanar_boxes, "coincident faces"),
    Case("stacked_boxes", "Stacked boxes", _stacked_boxes, "coincident faces"),
    Case("coaxial_cylinders", "Coaxial cylinders", _coaxial_cylinders, "coincident cylinders"),
    Case("sphere_touching", "Sphere touching a box", _sphere_touching, "contact at a point"),
    Case("cylinder_lying", "Cylinder lying on a box", _cylinder_lying, "contact along a line"),
    Case("cone_apex", "Cone cut through its apex", _cone_through_apex, "apex on a face"),
    Case("torus_box", "Torus and box", _torus_box),
    Case("sphere_cylinder", "Sphere centred on a cylinder", _sphere_cylinder, "poles on the section"),
    Case("nurbs_sphere_touching", "NURBS sphere touching a box", _nurbs_sphere_touching,
         "contact at a point on a B-spline surface"),
]

OPERATIONS = {"fuse": brep.fuse, "cut": brep.cut, "common": brep.common}


# ---------------------------------------------------------------------------- sweeps
@dataclass
class SweepCase:
    """A profile swept along a path (`make_pipe`): `build(m)` returns (the swept shape, the path)."""

    name: str
    title: str
    build: Callable
    note: str = ""


def _path(m, edges):
    return brep.make_wire_ordered(m, [(e, F) for e in edges])


def _elbow(m):
    """A disk along a segment, a quarter of circle and a segment: two cylinders and a torus."""
    e1 = brep.make_edge(m, [0.0, 0.0, 0.0], [0.0, 0.0, 3.0])
    arc = gbs.EllipseCurve(2.0, 2.0, [[2.0, 0.0, 3.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], 0.0, math.pi / 2)
    e2 = brep.make_edge(m, arc)
    e3 = brep.make_edge(m, list(arc.value(math.pi / 2)), [6.0, 0.0, 5.0])
    path = _path(m, [e1, e2, e3])
    return solid_of(m, brep.make_pipe(m, disk_face(m, 0.5), path)), path


def _helix(m, turns=2.0, radius=2.0, pitch=1.5, n=96):
    """A helix along z from (radius, 0, 0), as an interpolating B-spline."""
    pts = []
    for i in range(n + 1):
        t = 2 * math.pi * turns * i / n
        pts.append([radius * math.cos(t), radius * math.sin(t), pitch * t / (2 * math.pi)])
    return brep.make_edge(m, gbs.interpolate_cn(pts, 3, gbs.KnotsCalcMode.CHORD_LENGTH))


def _square_on_helix(m, frame):
    """A square 0.4 x 0.4 normal to the helix at its start, carried by `frame`."""
    path = _helix(m)
    t = np.array([0.0, 2.0, 1.5 / (2 * math.pi)])  # tangent at the start
    t /= np.linalg.norm(t)
    x = np.array([1.0, 0.0, 0.0])
    y = np.cross(t, x)
    c = np.array([2.0, 0.0, 0.0])
    sq = planar_face(m, [c - 0.2 * x - 0.2 * y, c + 0.2 * x - 0.2 * y, c + 0.2 * x + 0.2 * y, c - 0.2 * x + 0.2 * y])
    o = brep.PipeOptions(frame=frame, binormal=[0.0, 0.0, 1.0])
    return solid_of(m, brep.make_pipe(m, sq, path, o)), path


def _stadium(m):
    """A disk around a stadium (two half circles joined by segments): closed, no caps."""
    right = gbs.EllipseCurve(3.0, 3.0, [[4.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], 0.0, math.pi)
    left = gbs.EllipseCurve(3.0, 3.0, [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], 0.0, math.pi)
    v = [brep.make_vertex(m, p) for p in ([4.0, -3.0, 0.0], [4.0, 3.0, 0.0], [0.0, 3.0, 0.0], [0.0, -3.0, 0.0])]
    edges = [brep.make_edge(m, right, 0.0, math.pi, v[0], v[1]), brep.make_edge(m, v[1], v[2]),
             brep.make_edge(m, left, 0.0, math.pi, v[2], v[3]), brep.make_edge(m, v[3], v[0])]
    path = _path(m, edges)
    return solid_of(m, brep.make_pipe(m, disk_face(m, 0.5, (4.0, -3.0, 0.0), (1.0, 0.0, 0.0)), path)), path


SWEEPS = [
    SweepCase("pipe_elbow", "Disk along an elbow", _elbow, "exact cylinders and torus"),
    SweepCase("pipe_helix", "Square along a helix, rotation-minimizing frame",
              lambda m: _square_on_helix(m, brep.PipeFrame.RotationMinimizing), "B-spline faces within approx_tol"),
    SweepCase("pipe_helix_binormal", "Square along a helix, binormal along z",
              lambda m: _square_on_helix(m, brep.PipeFrame.FixedBinormal), "the square stays level"),
    SweepCase("pipe_stadium", "Disk around a closed path", _stadium, "no caps"),
]


# ---------------------------------------------------------------------------- offset solids, thick solids
@dataclass
class ThickCase:
    """A solid offset (`offset_solid`) or hollowed out (`make_thick`): `build(m)` returns (the
    result, the solid it comes from, the exact volume of the result or None)."""

    name: str
    title: str
    build: Callable
    note: str = ""


def _top_face(m, solid):
    """The face of `solid` highest along z (all its points at the top)."""
    faces = brep.explore_faces(m, solid)
    return max(faces, key=lambda f: brep.mesh(m, f).vertices[:, 2].min())


def _stadium_slab(m, L=3.0, R=1.0, H=0.8):
    """A slab of height H on a stadium: planes and half cylinders tangent along vertical edges."""
    v = [brep.make_vertex(m, p) for p in ([0.0, -R, 0.0], [L, -R, 0.0], [L, R, 0.0], [0.0, R, 0.0])]
    right = gbs.EllipseCurve(R, R, [[L, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], -math.pi / 2, math.pi / 2)
    left = gbs.EllipseCurve(R, R, [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], math.pi / 2, 3 * math.pi / 2)
    edges = [brep.make_edge(m, v[0], v[1]), brep.make_edge(m, right, -math.pi / 2, math.pi / 2, v[1], v[2]),
             brep.make_edge(m, v[2], v[3]), brep.make_edge(m, left, math.pi / 2, 3 * math.pi / 2, v[3], v[0])]
    f = brep.make_face(m, gbs.Plane([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], [-10.0, 10.0, -10.0, 10.0]), brep.make_wire(m, edges))
    return solid_of(m, brep.make_prism(m, f, [0.0, 0.0, H]))


def _offset_box(m):
    s = box(m, size=(2.0, 1.5, 1.0))
    return brep.offset_solid(m, s, 0.25), s, 2.5 * 2.0 * 1.5


def _offset_slab(m):
    s = _stadium_slab(m)
    d = 0.2
    return brep.offset_solid(m, s, d), s, (2 * (1.0 + d) * 3.0 + math.pi * (1.0 + d) ** 2) * (0.8 + 2 * d)


def _offset_helix(m):
    s, _ = _square_on_helix(m, brep.PipeFrame.RotationMinimizing)
    return brep.offset_solid(m, s, 0.05), s, None


def _thick_box(m):
    s = box(m, size=(2.0, 1.5, 1.0))
    return brep.make_thick(m, s, [_top_face(m, s)], 0.1), s, 3.0 - 1.8 * 1.3 * 0.9


def _thick_cup(m):
    s = cylinder(m, 1.0, 2.0)
    return brep.make_thick(m, s, [_top_face(m, s)], 0.1), s, math.pi * (2.0 - 0.81 * 1.9)


THICK = [
    ThickCase("offset_box", "Box offset by 0.25", _offset_box, "planes moved, edges and corners where they meet"),
    ThickCase("offset_slab", "Rounded slab offset by 0.2", _offset_slab, "tangent edges between planes and half cylinders"),
    ThickCase("offset_helix", "Square swept along a helix, offset by 0.05", _offset_helix,
              "sharp edges between offset B-spline surfaces, by marching"),
    ThickCase("thick_box", "Box hollowed out to 0.1, open at its top", _thick_box, "make_thick"),
    ThickCase("thick_cup", "Cup: a cylinder hollowed out to 0.1, open at its top", _thick_cup, "make_thick"),
]


# ---------------------------------------------------------------------------- fillets and chamfers: geometry
@dataclass
class BlendCase:
    """Blends along edges of a solid, their geometry only (`fillet_geometry`, `chamfer_geometry`):
    `build(m)` returns the solid and its blends (`EdgeBlend`), continued past the vertices of
    their edges for the cut by the faces at their ends."""

    name: str
    title: str
    build: Callable
    note: str = ""


def edge_near(m, shape, point):
    """The edge of `shape` passing nearest `point` (17 points along each)."""
    p = np.asarray(point, float)
    best, gap = None, math.inf
    for e in brep.explore_edges(m, shape):
        ed = m.edge(e)
        if ed.degenerate or ed.curve is None:
            continue
        for k in range(17):
            d = np.linalg.norm(np.asarray(ed.curve.value(ed.u1 + (ed.u2 - ed.u1) * k / 16)) - p)
            if d < gap:
                best, gap = e, d
    return best


def _l_block(m):
    """An L-shaped prism: convex edges outside, a concave one in the inner corner."""
    f = planar_face(m, [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 2.0, 0.0], [0.0, 2.0, 0.0]])
    s = solid_of(m, brep.make_prism(m, f, [0.0, 0.0, 1.0]))
    return s, [brep.fillet_geometry(m, s, edge_near(m, s, (2.0, 0.5, 1.0)), 0.3),
               brep.fillet_geometry(m, s, edge_near(m, s, (1.0, 1.0, 0.5)), 0.25),
               brep.chamfer_geometry(m, s, edge_near(m, s, (0.5, 2.0, 1.0)), 0.2, 0.3)]


def _cylinder_blends(m):
    """The top edge of a cylinder filleted (a torus), its bottom edge chamfered (a cone)."""
    s = cylinder(m, 1.0, 1.5)
    a = (math.cos(1.0), math.sin(1.0))
    return s, [brep.fillet_geometry(m, s, edge_near(m, s, (a[0], a[1], 1.5)), 0.25),
               brep.chamfer_geometry(m, s, edge_near(m, s, (a[0], a[1], 0.0)), 0.2, 0.15)]


def _crossed_blend(m):
    """A cylinder standing on a lying one: the saddle where they meet, concave."""
    lying = cylinder(m, 1.0, 4.0, (-2.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    s = brep.explore_solids(m, brep.fuse(m, lying, cylinder(m, 0.4, 1.5)).shape)[0]
    a = 0.4 * math.cos(1.0), 0.4 * math.sin(1.0)
    return s, [brep.fillet_geometry(m, s, edge_near(m, s, (a[0], a[1], math.sqrt(1.0 - a[1] ** 2))), 0.15)]


def _bump_blends(m):
    """The prism of a B-spline bump: its top edge filleted, its bottom edge chamfered."""
    bump = gbs.interpolate_cn([[0.0, 0.0, 0.0], [0.5, -0.3, 0.0], [1.0, 0.2, 0.0], [1.5, -0.1, 0.0], [2.0, 0.0, 0.0]], 3,
                              gbs.KnotsCalcMode.CHORD_LENGTH)
    a, b = bump.bounds()
    v = [brep.make_vertex(m, p) for p in ([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [0.0, 1.0, 0.0])]
    edges = [brep.make_edge(m, bump, a, b, v[0], v[1]), brep.make_edge(m, v[1], v[2]), brep.make_edge(m, v[2], v[3]),
             brep.make_edge(m, v[3], v[0])]
    f = brep.make_face(m, gbs.Plane([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], [-10.0, 10.0, -10.0, 10.0]), brep.make_wire(m, edges))
    s = solid_of(m, brep.make_prism(m, f, [0.0, 0.0, 0.8]))
    mid = np.asarray(bump.value((a + b) / 2))
    return s, [brep.fillet_geometry(m, s, edge_near(m, s, mid + [0.0, 0.0, 0.8]), 0.1),
               brep.chamfer_geometry(m, s, edge_near(m, s, mid), 0.08, 0.05)]


BLENDS = [
    BlendCase("blend_l_block", "Fillets and a chamfer on an L-shaped block", _l_block,
              "two planes: a cylinder around the centre line (convex, then concave), a plane"),
    BlendCase("blend_cylinder", "Fillet and chamfer of the edges of a cylinder", _cylinder_blends,
              "a plane and a cylinder: a torus around a circular centre line, a cone"),
    BlendCase("blend_crossed", "Fillet where two cylinders cross", _crossed_blend,
              "a marched centre line, arcs of circle skinned into a rational B-spline"),
    BlendCase("blend_bump", "Fillet and chamfer along a B-spline face", _bump_blends,
              "a plane and a B-spline face: approximated within approx_tol"),
]


# ---------------------------------------------------------------------------- offset surfaces
@dataclass
class OffsetCase:
    """Faces and their offsets (`offset_face`, `gbs.OffsetSurface`): `build(m)` returns the faces
    with their colors, the curves to draw in red (gbs curves on [a, b]) and a summary."""

    name: str
    title: str
    build: Callable
    note: str = ""


def paraboloid_surface(a=1.0):
    """z = (x² + y²) / 2 on [-a, a]², exact as a biquadratic Bézier patch: curvature 1 at the
    origin, toward the normal (0, 0, 1) there."""
    x, c = [-a, 0.0, a], [a * a / 2, -a * a / 2, a * a / 2]
    poles = [[x[i], x[j], c[i] + c[j]] for j in range(3) for i in range(3)]
    k = [-a, -a, -a, a, a, a]
    return gbs.BSSurface3d(poles, k, k, 2, 2)


def wavy_surface(scale=1.0, shift=(0.0, 0.0, 0.0)):
    """A bicubic B-spline with interior knots, 5 x 5."""
    poles = []
    for j in range(6):
        for i in range(6):
            z = 0.3 * math.sin(1.3 * i) * math.cos(0.9 * j) + 0.05 * i * j
            poles.append([scale * i + shift[0], scale * j + shift[1], z + shift[2]])
    k = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0]
    return gbs.BSSurface3d(poles, k, k, 3, 3)


def _wavy_offsets(m):
    f = brep.make_face(m, wavy_surface())
    faces = [(f, "#bbbbbb")] + [(brep.offset_face(m, f, d), c) for d, c in ((0.4, "#1f77b4"), (-0.4, "#ff7f0e"))]
    return faces, [], "the face (grey) and its offsets at +0.4 (blue) and -0.4 (orange)"


def _paraboloid_fold(m):
    f = brep.make_face(m, paraboloid_surface())
    faces = [(f, "#bbbbbb")] + [(brep.offset_face(m, f, d), c) for d, c in ((0.5, "#1f77b4"), (0.9, "#2ca02c"))]
    margin, _ = gbs.offset_margin(m.face(f).surface, 0.9)
    try:
        brep.offset_face(m, f, 1.0)
        refused = "**1.0 accepted**"
    except brep.BRepError:
        refused = "1.0 refused"
    return faces, [], f"offsets 0.5 and 0.9 (margin {margin:.2g} at the origin), {refused}"


def _tangent_plane(m):
    cyl = gbs.CylindricalSurface(1.0, box=[0.0, 2 * math.pi, 0.0, 2.0]).to_nurbs()
    off = gbs.OffsetSurface(cyl, 0.5)
    plane = gbs.Plane([[1.5, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], [-1.0, 1.0, -0.5, 2.5])
    r = gbs.intersect.surfaces(off, plane)
    faces = [(brep.make_face(m, off), "#1f77b4"), (brep.make_face(m, plane), "#bbbbbb")]
    curves = [(c.curve, c.a, c.b) for c in r.curves]
    kinds = ", ".join("tangent line" if c.tangent else "crossing curve" for c in r.curves) or "none"
    return faces, curves, f"a NURBS cylinder r 1 offset by 0.5 and a plane tangent to it: {kinds}"


OFFSETS = [
    OffsetCase("offset_wavy", "Offsets of a B-spline face", _wavy_offsets, "exact derivatives of any order"),
    OffsetCase("offset_paraboloid", "Near the fold of a paraboloid", _paraboloid_fold,
               "the offset toward the centre of curvature folds at d = 1, refused"),
    OffsetCase("offset_tangent_plane", "Plane tangent to an offset cylinder", _tangent_plane,
               "the intersection, in red, is a tangent line: second derivatives of the offset"),
]


@dataclass
class Result:
    """The operands and the results of a case: shapes, or the error raised."""

    model: object
    a: object
    b: object
    shapes: dict = field(default_factory=dict)  # operation -> ShapeId
    errors: dict = field(default_factory=dict)  # operation -> message


def run(case: Case) -> Result:
    """Builds the operands of `case` and computes the three booleans."""
    m = brep.Model()
    a, b = case.build(m)
    r = Result(m, a, b)
    for op, f in OPERATIONS.items():
        try:
            r.shapes[op] = f(m, a, b).shape
        except brep.BRepError as e:
            r.errors[op] = str(e)
    return r


# ---------------------------------------------------------------------------- rendering
PALETTE = [
    "#4e79a7", "#f28e2b", "#59a14f", "#e15759", "#76b7b2", "#edc948",
    "#b07aa1", "#ff9da7", "#9c755f", "#bab0ac", "#86bcb6", "#d37295",
]


def seam_edges(m, shape):
    """The edges of `shape` that are no sharp edge: where a periodic surface (cylinder, cone,
    sphere, torus) closes on itself. Either the boundary of one face runs the edge twice (once in
    each sense: the seam of the face), or the edge parts two faces on the same surface (a region
    across the seam split into two faces there, left so by the booleans with `unify=False`)."""
    seams = set()
    surfaces = {}  # edge -> the surfaces of the faces using it
    for f in brep.explore_faces(m, shape):
        srf = m.face(f).surface_id  # (the identity of the surface object, shared by its faces)
        seen = set()
        for w in m.face(f).wires:
            for ce in m.wire(w).coedges:
                (seams if ce.edge.index in seen else seen).add(ce.edge.index)
                surfaces.setdefault(ce.edge.index, []).append(srf)
    for e, srfs in surfaces.items():
        if len(srfs) == 2 and srfs[0] == srfs[1]:
            seams.add(e)
    return seams


def edge_lines(m, shape, samples=48, seams=False):
    """The edges of the BREP of `shape` as a PyVista polyline set: the sharp edges, without the
    seams of the periodic surfaces (unless `seams`) and the degenerate edges of their poles."""
    import pyvista as pv

    hidden = set() if seams else seam_edges(m, shape)
    points, lines = [], []
    for e in brep.explore_edges(m, shape):
        ed = m.edge(e)
        if ed.degenerate or ed.curve is None or e.index in hidden:
            continue
        start = len(points)
        for k in range(samples + 1):
            points.append(np.asarray(ed.curve.value(ed.u1 + (ed.u2 - ed.u1) * k / samples)))
        lines.append([samples + 1, *range(start, start + samples + 1)])
    if not points:
        return None
    return pv.PolyData(np.array(points), lines=np.hstack(lines))


def face_polydata(m, shape, options=None):
    """The mesh of `shape` as PyVista polydata, with a `face` cell array (index of the face of the
    model of each triangle) and the normals of the surfaces at the corners."""
    import pyvista as pv

    mesh = brep.mesh(m, shape, options or brep.MeshOptions())
    t = mesh.triangles
    if len(t) == 0:  # (an empty result: nothing in common)
        return None, mesh
    poly = pv.PolyData(mesh.vertices[t.reshape(-1)], np.hstack([np.full((len(t), 1), 3), np.arange(3 * len(t)).reshape(-1, 3)]).ravel())
    poly.point_data["Normals"] = mesh.normals
    face = np.empty(len(t), dtype=int)
    for k, fr in enumerate(mesh.faces):
        face[fr.first:fr.first + fr.count] = k
    poly.cell_data["face"] = face
    return poly, mesh


def add_shape(plotter, m, shape, opacity=1.0, color=None, edges=True, options=None, seams=False):
    """Adds `shape` to `plotter`: one color per face (or `color`), its edges as lines (the seams of
    the periodic surfaces too when `seams`)."""
    poly, mesh = face_polydata(m, shape, options)
    if poly is not None:
        if color is None:
            colors = np.array([PALETTE[i % len(PALETTE)] for i in poly.cell_data["face"]])
            rgb = np.array([[int(c[1:3], 16), int(c[3:5], 16), int(c[5:7], 16)] for c in colors], dtype=np.uint8)
            poly.cell_data["rgb"] = rgb
            plotter.add_mesh(poly, scalars="rgb", rgb=True, opacity=opacity, smooth_shading=True)
        else:
            plotter.add_mesh(poly, color=color, opacity=opacity, smooth_shading=True)
    if edges:
        lines = edge_lines(m, shape, seams=seams)
        if lines is not None:
            plotter.add_mesh(lines, color="#222222", line_width=2)
    return mesh


def plot_case(result: Result, case: Case, off_screen=False, window_size=(1600, 450)):
    """The operands and the three results of a case side by side; returns the plotter and the
    meshes of the results (operation -> mesh)."""
    import pyvista as pv

    p = pv.Plotter(shape=(1, 4), off_screen=off_screen, window_size=window_size, border=False)
    p.subplot(0, 0)
    add_shape(p, result.model, result.a, opacity=0.6, color="#4e79a7")
    add_shape(p, result.model, result.b, opacity=0.6, color="#f28e2b")
    p.add_text(case.title, font_size=9)
    meshes = {}
    for k, op in enumerate(OPERATIONS, start=1):
        p.subplot(0, k)
        if op in result.shapes:
            mesh = meshes[op] = add_shape(p, result.model, result.shapes[op])
            n_faces = len(brep.explore_faces(result.model, result.shapes[op]))
            closed = "closed" if mesh.open_sides() == 0 else f"{mesh.open_sides()} open sides"
            p.add_text(f"{op}: {n_faces} faces, volume {mesh.volume():.4g}, {closed}", font_size=8)
        else:
            p.add_text(f"{op}: {result.errors.get(op, '')}"[:80], font_size=8, color="red")
    p.link_views()
    p.view_isometric()
    return p, meshes
