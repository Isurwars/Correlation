#!/usr/bin/env bash
set -ex

mkdir -p build-conda
cd build-conda

cmake -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_PREFIX_PATH="${PREFIX}" \
  -DBUILD_GUI=OFF \
  -DBUILD_WITH_HDF5=ON \
  -DBUILD_WITH_ARROW=OFF \
  -DPYTHON_EXECUTABLE="${PYTHON}" \
  ${SRC_DIR}

cmake --build . --config Release --parallel
cmake --install .
