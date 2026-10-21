Examples
========

Visual examples of the native BREP (``gbs.brep``), in ``docs/sources/examples``:

``booleans.ipynb``
   The booleans (``fuse``, ``cut``, ``common``) of primitive solids in 3D: the operands in
   transparency, each result with one color per face (so that the splits of the faces show) and
   the edges of the BREP in black. The cases are those of the comparison with OpenCascade:
   generic intersections, coincident faces, tangent contacts, branch points, poles.

``meshing.ipynb``
   Triangle meshes of shapes: the deflection, the angle and the smallest angle, a boolean result,
   STL and OBJ export.

``gallery.py``
   Renders every case off screen to PNG images and writes ``gallery.md``, which gathers them with
   the number of faces, the volume, the validity and the closure of each result::

      python docs/sources/examples/gallery.py --out gallery

   On Linux without a display, run it under ``xvfb-run -a``. The CI renders the gallery on every
   build and publishes it as the ``gallery`` artifact.

``gbs_examples.py``
   The shared module: primitives built with the BREP API (``box``, ``cylinder``, ``cone``,
   ``sphere``, ``torus``, exact or NURBS), the boolean cases, the meshes and their PyVista
   rendering. ``python/tests/test_examples.py`` checks that every case stays valid.

The notebooks need ``pygbs`` and ``pyvista`` (with ``trame`` for interactive views); they are
not run by the CI, the gallery is.
