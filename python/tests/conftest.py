"""Pytest configuration for the gbs Python tests.

Several tests end with an interactive VTK render call:
  * ``gbs.plot_curves`` / ``plot_curves_2d`` / ``plot_surfaces`` (gbs bindings)
  * ``pygbs.vtkplot.render_actors`` (gbs VTK helper)
  * pyvista's ``Plotter.show()``

These spin a blocking ``vtkRenderWindowInteractor`` event loop that never
returns without a user: a local run waits for each window to be closed by hand,
and a headless CI runner hangs (the interactor holds the GIL, so even
pytest-timeout cannot interrupt it).

The tests therefore display nothing unless the environment variable
GBS_TEST_PLOT_ON is set (to anything but "", "0" or "false"), as TEST_PLOT_ON
for the C++ tests: the blocking entry points are neutralized here, for every
test. The numerical assertions still run, and the actors and plotters are still
built. In a headless / CI environment the display stays off whatever
GBS_TEST_PLOT_ON says, and pyvista renders off-screen.
"""
import os
import sys

_HEADLESS = bool(os.environ.get("CI") or os.environ.get("PYVISTA_OFF_SCREEN"))
_PLOT_ON = os.environ.get("GBS_TEST_PLOT_ON", "").strip().lower() not in ("", "0", "false") and not _HEADLESS


def _nothing(*args, **kwargs):
    return None


try:
    import pyvista

    if _HEADLESS:
        # render off-screen (test_json.py calls Plotter.show())
        pyvista.OFF_SCREEN = True
    # On the headless Windows runner, even the off-screen render window of
    # VTK 9.6 crashes (access violation in Plotter.render, no GL context),
    # killing the whole pytest process. Skip the render there; the plotter
    # is still built, so the geometry/mesh code feeding it is exercised.
    # Locally, without GBS_TEST_PLOT_ON, no window either.
    if (_HEADLESS and sys.platform == "win32") or (not _HEADLESS and not _PLOT_ON):
        pyvista.Plotter.show = _nothing
except Exception:
    pass

if not _PLOT_ON:
    # gbs bindings: the plot_* helpers open a blocking interactive window.
    try:
        import pygbs.gbs as _gbs

        for _name in ("plot_curves", "plot_curves_2d", "plot_surfaces"):
            if hasattr(_gbs, _name):
                setattr(_gbs, _name, _nothing)
    except Exception:
        pass

    # gbs VTK helper used by the mesh tests.
    try:
        from pygbs import vtkplot as _vtkplot

        if hasattr(_vtkplot, "render_actors"):
            _vtkplot.render_actors = _nothing
    except Exception:
        pass
