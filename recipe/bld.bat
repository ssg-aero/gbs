mkdir build-conda
cd build-conda
cmake .. ^
-DCMAKE_BUILD_TYPE:STRING="Release" ^
-DUSE_OCCT_UTILS:BOOL=TRUE ^
-DGBS_BUILD_TESTS:BOOL=TRUE ^
-DUSE_RENDER:BOOL=TRUE ^
-DUSE_PYTHON_BIDING=TRUE ^
-G"Ninja" ^
-Wno-dev

@REM ninja install
ninja