"""Python bindings of the STEP reader (gbs.brep.read_step)."""
import math
import os
import pathlib

import pytest

import pygbs.gbs as gbs

brep = gbs.brep
STEP_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "in", "step")


def step_file(name):
    return os.path.join(STEP_DIR, name)


def volume(m, root):
    return sum(brep.signed_volume(m, m.solid(s).outer, 256) for s in brep.explore_solids(m, root))


def test_read_box():
    m = brep.Model()
    root, report = brep.read_step(m, step_file("box.stp"))
    assert isinstance(root, brep.SolidId)
    assert report.solids == 1 and report.faces == 6 and report.edges == 12 and report.vertices == 8
    assert report.check.ok()
    assert not report.failed and not report.unsupported
    assert report.length_unit == "MILLI METRE"
    assert report.file_tolerance == pytest.approx(1e-7)
    assert m.name(root) == "BOX"
    assert m.external_id(root) is not None
    assert volume(m, root) == pytest.approx(24., abs=1e-9)
    assert "solids=1" in repr(report)


def test_units_and_text():
    m = brep.Model()
    root, report = brep.read_step(m, step_file("plate_inch.stp"), target_unit_mm=25.4)
    assert m.unit_scale == 25.4
    assert report.length_unit == "INCH"
    assert report.length_factor == pytest.approx(1.)
    with open(step_file("cylinder.stp")) as f:
        text = f.read()
    m2 = brep.Model()
    root2, report2 = brep.read_step_text(m2, text)
    assert volume(m2, root2) == pytest.approx(math.pi * 25 * 10, abs=0.2)


def test_assembly_and_copies():
    m = brep.Model()
    root, report = brep.read_step(m, step_file("assembly.stp"))
    assert isinstance(root, brep.CompoundId)
    assert m.name(root) == "ASM"
    assert report.solids == 3
    shapes = m.compound(root).shapes
    assert [m.name(s) for s in shapes] == ["BOX", "BOX", "SUB"]
    # a further copy, moved by a placement
    p = brep.Placement.frame([0., 0., 100.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.])
    c = brep.copy_shape(m, shapes[0], p)
    (lo0, _), (lo1, _) = brep.bounding_box(m, shapes[0]), brep.bounding_box(m, c)
    assert lo1[2] - lo0[2] == pytest.approx(100., abs=1e-6)
    assert (p * p.inverse()).is_identity()


def test_partial_strict_and_errors():
    m = brep.Model()
    root, report = brep.read_step(m, step_file("box_unsupported_face.stp"))
    assert isinstance(root, brep.ShellId)
    assert report.unsupported[0].type == "SURFACE_REPLICA"
    assert len(report.failed) == 2
    ms = brep.Model()
    with pytest.raises(brep.BRepError):
        brep.read_step(ms, step_file("box_unsupported_face.stp"), strict=True)
    assert ms.count(brep.ShapeType.Face) == 0
    with pytest.raises(brep.BRepError):
        brep.read_step(ms, step_file("no_such_file.stp"))
    with pytest.raises(brep.BRepError):
        brep.read_step_text(ms, "not a step file")


def test_assembly_provenance():
    m = brep.Model()
    root, report = brep.read_step(m, step_file("assembly.stp"))
    s = m.compound(root).shapes
    i1, i2 = brep.instance(m, s[0]), brep.instance(m, s[1])
    assert i1.part == i2.part
    assert i1.place.apply([0., 0., 0.]) == pytest.approx([10., 0., 0.])
    k = brep.new_part(m)
    brep.set_instance(m, s[2], brep.InstanceOf(k))
    assert brep.instance(m, s[2]).part == k
    brep.set_instance(m, s[2], None)
    assert brep.instance(m, s[2]) is None


def test_write_and_read_back(tmp_path):
    m = brep.Model()
    root, _ = brep.read_step(m, step_file("cylinder.stp"))
    path = str(tmp_path / "cyl.stp")
    report = brep.write_step(m, root, path)
    assert report.parts == 1 and report.solids == 1
    assert report.approximated_surfaces == 0
    m2 = brep.Model()
    root2, r2 = brep.read_step(m2, path)
    assert r2.check.ok()
    assert volume(m2, root2) == pytest.approx(volume(m, root), rel=1e-6)
    text, report = brep.write_step_text(m, [root], schema="AP242", product_name="tube")
    assert "AP242" in text
    m3 = brep.Model()
    root3, r3 = brep.read_step_text(m3, text)
    assert r3.solids == 1


def test_write_assembly_and_errors(tmp_path):
    m = brep.Model()
    root, _ = brep.read_step(m, step_file("assembly.stp"))
    text, report = brep.write_step_text(m, root)
    assert report.assemblies == 2 and report.solids == 2 and report.occurrences == 4
    text, report = brep.write_step_text(m, root, flatten=True)
    assert report.parts == 1 and report.assemblies == 0
    with pytest.raises(brep.BRepError):
        brep.write_step_text(m, root, schema="AP203")
    with pytest.raises(brep.BRepError):
        brep.write_step(m, root, str(tmp_path / "no" / "such" / "dir.stp"))


def test_paths_and_iterables(tmp_path):
    # os.PathLike paths were refused (TypeError); shapes may be any iterable
    m = brep.Model()
    root, _ = brep.read_step(m, pathlib.Path(step_file("box.stp")))
    rep = brep.write_step(m, {root}, tmp_path / "box.stp")
    assert rep.solids == 1
    m2 = brep.Model()
    root2, report = brep.read_step(m2, tmp_path / "box.stp")
    assert report.solids == 1
    text, rep = brep.write_step_text(m, (s for s in [root]))
    assert rep.solids == 1
    with pytest.raises(TypeError):
        brep.write_step_text(m, 3)
    with pytest.raises(TypeError):
        brep.write_step_text(m, [root, "x"])
    with pytest.raises(brep.BRepError) as err:
        brep.read_step(brep.Model(), step_file("box.stp"), face_options=brep.MakeFaceOptions(pcurve_degree=0))
    assert err.value.code == brep.BuildErrc.invalid_options
