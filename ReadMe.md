GBS is header library for computing bspline curves and surface.

GBS is compatible with OpenCASCADE, curves and surfaces can be converted to OpenCASCADE's objects.

![Screencast](foilApproximation.png)

For now, this lib is to be used inside a conda environment with the following package installed:
* nlopt 
* eigen3

The optional module occt-utils requires the additional package:
* occt >=7.4.0

The optional module render requires the additional package:
* vtk >=9.0

The test library needs:
* gtest
* occt>=7.4.0
* sundials

**Warning tests relative to performances evaluation should be run in release mode**

As GBS base is a header library it doesn’t need compilation.

If one needs to compile the optional module occt-utils, -DUSE_OCCT_UTILS:BOOL=TRUE shall be added to cmake command.
If one needs to compile the optional module occt-utils, -DUSE_RENDER:BOOL=TRUE shall be added to cmake command.

The full test suite, which require the optional module occt-utils, please add -DGBS_BUILD_TESTS:BOOL=TRUE to the cmake command.
