mkdir build-conda
cd build-conda

cmake .. ^
-D CMAKE_BUILD_TYPE:STRING="Release" ^
-D USE_OCCT_UTILS:BOOL=TRUE ^
-D GBS_BUILD_TESTS:BOOL=TRUE ^
-D USE_RENDER:BOOL=FALSE ^
-D USE_PYTHON_BINDINGS=TRUE ^
-D BUILD_DOC=TRUE ^
-D CMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX% ^
-G "Ninja" ^
-Wno-dev

ninja install

cd ..
