mkdir build-conda
cd build-conda

cmake .. ^
-D CMAKE_BUILD_TYPE:STRING="Release" ^
-D BUILD_GBS=ON ^
-D BUILD_GBS_MESH=ON ^
-D BUILD_GBS_RENDER=ON ^
-D BUILD_PYTHON_BINDINGS=ON ^
-D CMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX% ^
-G "Ninja" ^
-Wno-dev

ninja install

cd ../python
%PYTHON% -m pip install ../python --no-deps -vv

cd ..
