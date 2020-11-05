mkdir build
cd build

cmake .. ^
-D CMAKE_BUILD_TYPE:STRING="Release" ^
-D USE_OCCT_UTILS:BOOL=TRUE ^
-D GBS_BUILD_TESTS:BOOL=TRUE ^
-D USE_RENDER:BOOL=FALSE ^
-D USE_PYTHON_BINDINGS=TRUE ^
-D CMAKE_INSTALL_PREFIX=%CONDA_PREFIX%/Library ^
-G "Ninja" ^
-Wno-dev

ninja install

cd ..
