"""The notebooks of docs/sources/examples run without error against the current API (images only,
off screen): a renamed option or function shows here, not in front of a reader."""
import os
import pathlib
import sys

import pytest

nbformat = pytest.importorskip("nbformat")
nbclient = pytest.importorskip("nbclient")

EXAMPLES = pathlib.Path(__file__).parents[2] / "docs" / "sources" / "examples"


# the notebooks render with VTK: the Windows runners of the CI have no OpenGL context (the kernel
# dies at the first render: "failed to get wglChoosePixelFormatARB"); Linux (xvfb) and macOS run them
@pytest.mark.skipif(sys.platform == "win32" and bool(os.environ.get("CI")), reason="no OpenGL context on the Windows runners of the CI")
@pytest.mark.parametrize("name", ["booleans.ipynb", "meshing.ipynb"])
def test_notebook_runs(name, monkeypatch):
    monkeypatch.setenv("GBS_NOTEBOOK_BACKEND", "static")
    monkeypatch.setenv("PYVISTA_OFF_SCREEN", "true")
    nb = nbformat.read(EXAMPLES / name, as_version=4)
    client = nbclient.NotebookClient(nb, timeout=300, kernel_name="python3", resources={"metadata": {"path": str(EXAMPLES)}})
    client.execute()
    errors = [o for c in nb.cells if c.cell_type == "code" for o in c.get("outputs", []) if o.get("output_type") == "error"]
    assert not errors
