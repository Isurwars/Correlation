# SPDX-License-Identifier: AGPL-3.0-only
# Patch miniz/CMakeLists.txt to replace /Zi with /Z7 on MSVC.
# Hardcoded /Zi causes CMake to emit /Fd<TARGET_COMPILE_PDB>, which conflicts
# with sccache and CMAKE_MSVC_DEBUG_INFORMATION_FORMAT=Embedded (-Z7).
if(EXISTS "${MINIZ_DIR}/CMakeLists.txt")
  file(READ "${MINIZ_DIR}/CMakeLists.txt" _miniz_cmake)
  string(REPLACE "/Zi" "/Z7" _miniz_cmake "${_miniz_cmake}")
  file(WRITE "${MINIZ_DIR}/CMakeLists.txt" "${_miniz_cmake}")
  message(STATUS "Patched miniz/CMakeLists.txt (/Zi -> /Z7 for sccache compatibility)")
endif()

