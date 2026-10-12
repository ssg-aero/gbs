.. GBS documentation master file, created by
   sphinx-quickstart on Thu Dec 17 12:02:18 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to GBS's documentation!
===============================

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Docs
====

.. doxygenfile:: bscbuild.h
   :project: GBS

.. doxygenfile:: bssurf.h
   :project: GBS

.. doxygenfile:: bscurve.h
   :project: GBS

Native BREP core (gbs-brep)
===========================

Design: ``docs/sources/design/brep_core.md`` and the architecture notes
``docs/sources/design/brep_pr01_model.md`` to ``brep_pr10_python.md``.

.. doxygenfile:: model.h
   :project: GBS

.. doxygenfile:: explore.h
   :project: GBS

.. doxygenfile:: builders.h
   :project: GBS

.. doxygenfile:: pcurve.h
   :project: GBS

.. doxygenfile:: closure.h
   :project: GBS

.. doxygenfile:: sew.h
   :project: GBS

.. doxygenfile:: check.h
   :project: GBS

.. doxygenfile:: edit.h
   :project: GBS

.. doxygenfile:: placement.h
   :project: GBS

Exact elementary geometry
=========================

.. doxygenfile:: bselementary.h
   :project: GBS

.. doxygenfile:: primitives.h
   :project: GBS

STEP reader (gbs-io/step)
=========================

Design: ``docs/sources/design/step_reader.md`` and the architecture notes
``docs/sources/design/step_pr01_p21.md`` to ``step_pr10_python.md``.
Python: ``gbs.brep.read_step(model, path)`` returns ``(root, report)``.

.. doxygenfile:: reader.h
   :project: GBS

.. doxygenfile:: geometry.h
   :project: GBS

.. doxygenfile:: units.h
   :project: GBS

.. doxygenfile:: p21.h
   :project: GBS

.. doxygenfile:: p21_writer.h
   :project: GBS
