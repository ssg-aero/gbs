"""Shared helpers of the GBS examples: primitive solids of the native BREP, the cases of the
booleans, their meshes and their PyVista rendering.

The primitives are built with the BREP API itself (planar faces, prisms, revolutions), placed by
rigid placements; each face of a result gets its own color, so that the splits of the faces are
seen, and the edges of the BREP are drawn as lines: the sharp ones, not the seams where a periodic
surface closes on itself (they bound no region) nor the degenerate edges of the poles.
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
def _frame(origin, z):
    """A placement putting the z axis along `z`, at `origin`."""
    z = np.asarray(z, float)
    z = z / np.linalg.norm(z)
    a = np.array([1.0, 0.0, 0.0]) if abs(z[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    x = np.cross(a, z)
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    return brep.Placement.frame(list(map(float, origin)), list(x), list(y), list(z))


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


def _ax2(origin, z):
    z = np.asarray(z, float)
    z = z / np.linalg.norm(z)
    a = np.array([1.0, 0.0, 0.0]) if abs(z[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    x = np.cross(a, z)
    x /= np.linalg.norm(x)
    return [list(map(float, origin)), list(z), list(x)]


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
         "known limit: a point contact with a NURBS surface stays singular (NonGeneric)"),
]

OPERATIONS = {"fuse": brep.fuse, "cut": brep.cut, "common": brep.common}


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
    """The seams of the faces of `shape`: the edges that the boundary of one face runs twice (once in
    each sense), where a periodic surface (cylinder, cone, sphere, torus) closes on itself."""
    seams = set()
    for f in brep.explore_faces(m, shape):
        seen = set()
        for w in m.face(f).wires:
            for ce in m.wire(w).coedges:
                (seams if ce.edge.index in seen else seen).add(ce.edge.index)
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
    """The operands and the three results of a case side by side."""
    import pyvista as pv

    p = pv.Plotter(shape=(1, 4), off_screen=off_screen, window_size=window_size, border=False)
    p.subplot(0, 0)
    add_shape(p, result.model, result.a, opacity=0.6, color="#4e79a7")
    add_shape(p, result.model, result.b, opacity=0.6, color="#f28e2b")
    p.add_text(case.title, font_size=9)
    for k, op in enumerate(OPERATIONS, start=1):
        p.subplot(0, k)
        if op in result.shapes:
            mesh = add_shape(p, result.model, result.shapes[op])
            n_faces = len(brep.explore_faces(result.model, result.shapes[op]))
            closed = "closed" if mesh.open_sides() == 0 else f"{mesh.open_sides()} open sides"
            p.add_text(f"{op}: {n_faces} faces, volume {mesh.volume():.4g}, {closed}", font_size=8)
        else:
            p.add_text(f"{op}: {result.errors.get(op, '')}"[:80], font_size=8, color="red")
    p.link_views()
    p.view_isometric()
    return p
