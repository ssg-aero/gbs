
mkdir build-conda
cd build-conda

cmake .. ^
-D CMAKE_BUILD_TYPE:STRING="Release" ^
-D GBS_USE_OCCT_UTILS=OFF ^
-D GBS_USE_CUDA=OFF ^
-D GBS_BUILD_TESTS=OFF ^
-D GBS_USE_RENDER=ON ^
-D GBS_USE_PYTHON_BINDINGS=ON ^
-D GBS_BUILD_DOC=OFF ^
-D CMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX% ^
-G "Ninja" ^
--log-level=DEBUG ^
-Wno-dev

ninja install

cd ..
