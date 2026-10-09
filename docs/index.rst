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
