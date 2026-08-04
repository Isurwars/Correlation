@echo on

mkdir build-conda
cd build-conda

cmake -G "Ninja" ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%" ^
  -DCMAKE_PREFIX_PATH="%LIBRARY_PREFIX%" ^
  -DBUILD_GUI=OFF ^
  -DBUILD_WITH_HDF5=ON ^
  -DBUILD_WITH_ARROW=OFF ^
  -DPYTHON_EXECUTABLE="%PYTHON%" ^
  "%SRC_DIR%"
if errorlevel 1 exit 1

cmake --build . --config Release --parallel
if errorlevel 1 exit 1

cmake --install .
if errorlevel 1 exit 1
