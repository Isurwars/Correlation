@echo on

:: 1. Build and install C++ core engine and CLI binary
cmake -G "Ninja" ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%" ^
  -DCMAKE_PREFIX_PATH="%LIBRARY_PREFIX%" ^
  -DBUILD_GUI=OFF ^
  -DBUILD_TESTING=OFF ^
  -DBUILD_WITH_HDF5=ON ^
  -DBUILD_WITH_ARROW=OFF ^
  -DBUILD_PYTHON_BINDINGS=OFF ^
  -B build-conda ^
  -S "%SRC_DIR%"
if errorlevel 1 exit 1

cmake --build build-conda --config Release --parallel
if errorlevel 1 exit 1

cmake --install build-conda
if errorlevel 1 exit 1

:: 2. Build and install Python extension and packages
"%PYTHON%" -m pip install "%SRC_DIR%" --no-deps --no-build-isolation -vv
if errorlevel 1 exit 1
