"""Blades as shells from Python: contours cut in pieces skinned one by one (loft_pieces), extended
together (extrapolated_v_by_lengths), sewn into one shell (brep.sewn_shell), meshed with their
edges (Mesh.edges)."""
import numpy as np
import pytest

import pygbs.gbs as gbs

brep = gbs.brep
ZS, V = [0.0, 0.6, 1.5, 2.0], [0.0, 0.3, 0.75, 1.0]


def contour_pieces(zs):
    # a cubic side from (0, 0) to (1 + 0.3 z, 0) and a straight edge back, closing it
    k = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0]
    side = [gbs.BSCurve3d([[0.0, 0.0, z], [0.3, 0.2 + 0.1 * z, z], [0.6, 0.25, z], [0.9 + 0.2 * z, 0.1, z], [1.0 + 0.3 * z, 0.0, z]], k, 3)
            for z in zs]
    edge = [gbs.BSCurve3d([[1.0 + 0.3 * z, 0.0, z], [0.0, 0.0, z]], [0.0, 0.0, 1.0, 1.0], 1) for z in zs]
    return [side, edge]


def test_pieces_extended_sewn_meshed():
    pieces = contour_pieces(ZS)
    lofted = gbs.loft_pieces(pieces, V)
    assert len(lofted) == 2
    for i, v in enumerate(V):
        assert np.allclose(lofted[0].value(0.5, v), pieces[0][i].value(0.5), atol=1e-12)
    surfaces = gbs.extrapolated_v_by_lengths(lofted, [0.1, 0.1])
    v1, v2 = surfaces[0].bounds()[2:]
    assert v1 < 0.0 and v2 > 1.0
    for v in np.linspace(v1, v2, 11):
        assert np.allclose(surfaces[0].value(1.0, v), surfaces[1].value(0.0, v), atol=1e-14)

    m = brep.Model()
    sewn = brep.sewn_shell(m, surfaces)
    assert len(sewn.faces) == 2
    assert len(sewn.report.merged) == 2
    mesh = brep.mesh(m, sewn.shell)
    assert len(mesh.edges) == 6  # the two joins, the four ends
    for edge, vertices in mesh.edges:
        e = m.edge(edge)
        assert vertices.dtype == np.uint32 and len(vertices) >= 2
        assert np.allclose(mesh.vertices[vertices[0]], m.vertex(e.v1).pnt, atol=1e-9)
        assert np.allclose(mesh.vertices[vertices[-1]], m.vertex(e.v2).pnt, atol=1e-9)


def test_refusals():
    pieces = contour_pieces(ZS)
    with pytest.raises(ValueError):
        gbs.loft_pieces(pieces, [0.0, 0.75, 0.3, 1.0])
    with pytest.raises(ValueError):
        gbs.loft_pieces([pieces[0], pieces[1][:3]], V)
    lofted = gbs.loft_pieces(pieces, V)
    with pytest.raises(ValueError):
        gbs.extrapolated_v_by_lengths(lofted, [float("nan"), 0.1])
    # surfaces apart: no single shell, the model left empty
    far = gbs.loft_pieces([[gbs.BSCurve3d([[x + 5.0 for x in p] for p in c.poles()], c.knotsFlats(), c.degree()) for c in piece]
                           for piece in pieces], V)
    m = brep.Model()
    with pytest.raises(brep.BRepError, match="disconnected"):
        brep.sewn_shell(m, [lofted[0], far[1]])
    assert m.empty()
