"""Pytest configuration for the gbs Python tests.

Several tests end with an interactive VTK render call:
  * ``gbs.plot_curves`` / ``plot_curves_2d`` / ``plot_surfaces`` (gbs bindings)
  * ``pygbs.vtkplot.render_actors`` (gbs VTK helper)
  * pyvista's ``Plotter.show()``

These spin a blocking ``vtkRenderWindowInteractor`` event loop that never
returns without a user, which hangs the whole suite on a headless CI runner
(the interactor holds the GIL, so even pytest-timeout cannot interrupt it).

In a headless / CI environment we therefore force off-screen rendering and
neutralize the blocking entry points. The numerical assertions in those tests
still run; only the interactive display is short-circuited. Local interactive
runs (no CI / PYVISTA_OFF_SCREEN set) are left untouched.
"""
import os

_HEADLESS = bool(os.environ.get("CI") or os.environ.get("PYVISTA_OFF_SCREEN"))

if _HEADLESS:
    # pyvista: render off-screen (test_json.py calls Plotter.show()).
    try:
        import pyvista

        pyvista.OFF_SCREEN = True
    except Exception:
        pass

    # gbs bindings: the plot_* helpers open a blocking interactive window.
    try:
        import pygbs.gbs as _gbs

        for _name in ("plot_curves", "plot_curves_2d", "plot_surfaces"):
            if hasattr(_gbs, _name):
                setattr(_gbs, _name, lambda *args, **kwargs: None)
    except Exception:
        pass

    # gbs VTK helper used by the mesh tests.
    try:
        from pygbs import vtkplot as _vtkplot

        if hasattr(_vtkplot, "render_actors"):
            _vtkplot.render_actors = lambda *args, **kwargs: None
    except Exception:
        pass
