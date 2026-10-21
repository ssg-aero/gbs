"""Python bindings of stage 3 of the BREP kernel (gbs.brep): triangle meshes, STL and OBJ."""
import math
import os

import numpy as np
import pytest

import pygbs.gbs as gbs

brep = gbs.brep
STEP_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "in", "step")


def read(m, name):
    root, report = brep.read_step(m, os.path.join(STEP_DIR, name))
    assert report.check.ok()
    return root


def test_mesh_of_a_cylinder():
    # the cylinder R 5 h 10: a closed mesh facing outwards, within the deflection
    m = brep.Model()
    root = read(m, "cylinder.stp")
    mesh = brep.mesh(m, root, brep.MeshOptions(deflection=1e-2))
    v, t = mesh.vertices, mesh.triangles
    assert v.shape[1] == 3 and t.shape[1] == 3
    assert t.dtype == np.uint32 and t.max() < len(v)
    assert mesh.uv.shape == (3 * len(t), 2)
    assert mesh.normals.shape == (3 * len(t), 3)
    assert np.allclose(np.linalg.norm(mesh.normals, axis=1), 1.0)
    assert mesh.open_sides() == 0
    assert mesh.incomplete == []
    assert mesh.folded == []
    assert sum(f.count for f in mesh.faces) == len(t)
    exact = math.pi * 25.0 * 10.0
    assert mesh.volume() == pytest.approx(exact, abs=mesh.area() * 1e-2)
    assert mesh.volume() < exact  # (inscribed)
    # the volume again from the arrays
    a, b, c = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
    assert np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0 == pytest.approx(mesh.volume())


def test_options_and_errors():
    o = brep.MeshOptions()
    assert o.deflection == 0.0 and o.angle == 0.5 and o.max_edge == 0.0 and o.max_points == 200000
    assert o.min_angle == pytest.approx(math.radians(20.0))
    m = brep.Model()
    root = read(m, "box.stp")
    mesh = brep.mesh(m, root, brep.MeshOptions(max_edge=0.5))
    assert len(mesh.triangles) > 12
    sides = mesh.vertices[mesh.triangles[:, [1, 2, 0]]] - mesh.vertices[mesh.triangles]
    assert np.linalg.norm(sides, axis=2).max() <= 0.5 + 1e-12
    with pytest.raises(brep.BRepError):
        brep.write_stl(mesh, os.path.join("no_such_directory", "x.stl"))
    # the limit of points: the face reported
    root = read(m, "cone.stp")
    assert len(brep.mesh(m, root, brep.MeshOptions(deflection=1e-4, max_points=100)).incomplete) > 0


def test_stl_and_obj(tmp_path):
    m = brep.Model()
    root = read(m, "sphere.stp")
    mesh = brep.mesh(m, root)
    n = len(mesh.triangles)
    stl = tmp_path / "sphere.stl"
    brep.write_stl(mesh, str(stl))
    assert stl.stat().st_size == 84 + 50 * n
    data = np.frombuffer(stl.read_bytes()[84:], dtype=np.dtype([("n", "<f4", 3), ("p", "<f4", (3, 3)), ("a", "<u2")]))
    assert np.allclose(data["p"][:, 0], mesh.vertices[mesh.triangles[:, 0]], atol=1e-5)
    ascii_stl = tmp_path / "sphere_ascii.stl"
    brep.write_stl(mesh, str(ascii_stl), binary=False)
    assert ascii_stl.read_text().count("facet normal") == n
    obj = tmp_path / "sphere.obj"
    brep.write_obj(mesh, str(obj))
    lines = obj.read_text().splitlines()
    assert sum(l.startswith("v ") for l in lines) == len(mesh.vertices)
    assert sum(l.startswith("f ") for l in lines) == n


def test_vtk_and_pyvista():
    m = brep.Model()
    mesh = brep.mesh(m, read(m, "sphere.stp"), brep.MeshOptions(deflection=1e-2))
    actor = gbs.make_actor(mesh)
    poly = actor.GetMapper().GetInput()
    assert poly.GetNumberOfCells() == len(mesh.triangles)
    assert poly.GetNumberOfPoints() == 3 * len(mesh.triangles)
    pytest.importorskip("pyvista")
    from pygbs.vistaplot import mesh_to_polydata

    poly = mesh_to_polydata(mesh)
    assert poly.n_cells == len(mesh.triangles)
    assert poly.volume == pytest.approx(mesh.volume(), rel=1e-9)
