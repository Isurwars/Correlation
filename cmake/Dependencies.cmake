# -----------------------------------------------------------
# Dependencies
# -----------------------------------------------------------
include(FetchContent)

# 1. TBB
message(STATUS "Downloading TBB from GitHub...")
FetchContent_Declare(
  TBB
  GIT_REPOSITORY https://github.com/uxlfoundation/oneTBB.git
  GIT_TAG v2022.3.0  
)
FetchContent_MakeAvailable(TBB)

# 2. Slint
message(STATUS "Downloading Slint from GitHub...")
FetchContent_Declare(
  Slint
  GIT_REPOSITORY https://github.com/slint-ui/slint.git
  GIT_TAG v1.14.1
  SOURCE_SUBDIR api/cpp
)
FetchContent_MakeAvailable(Slint)

# 3. HDF5
message(STATUS "Downloading HDF5 from GitHub...")
# HDF5 options for FetchContent
set(HDF5_BUILD_CPP_LIB ON CACHE BOOL "Build HDF5 C++ Library" FORCE)
set(HDF5_BUILD_EXAMPLES OFF CACHE BOOL "Build HDF5 Examples" FORCE)
set(HDF5_BUILD_HL_LIB ON CACHE BOOL "Build HDF5 High Level Library" FORCE)
set(HDF5_BUILD_TOOLS OFF CACHE BOOL "Build HDF5 Tools" FORCE)
set(HDF5_BUILD_FORTRAN OFF CACHE BOOL "Build HDF5 Fortran" FORCE)
set(HDF5_BUILD_JAVA OFF CACHE BOOL "Build HDF5 Java" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "Build HDF5 Tests" FORCE)
set(HDF5_PACK_EXAMPLES OFF CACHE BOOL "Pack HDF5 Examples" FORCE)

FetchContent_Declare(
  HDF5
  GIT_REPOSITORY https://github.com/HDFGroup/hdf5.git
  GIT_TAG hdf5_2.0.0
)
FetchContent_MakeAvailable(HDF5)

# Alias logic for HighFive if we built HDF5 ourselves
if(NOT TARGET HDF5::HDF5)
  if(TARGET hdf5_cpp-shared)
    add_library(HDF5::HDF5 ALIAS hdf5_cpp-shared)
  elseif(TARGET hdf5_cpp-static)
    add_library(HDF5::HDF5 ALIAS hdf5_cpp-static)
  endif()
  if(TARGET hdf5_hl-shared)
    add_library(HDF5::HDF5_HL ALIAS hdf5_hl-shared)
  elseif(TARGET hdf5_hl-static)
    add_library(HDF5::HDF5_HL ALIAS hdf5_hl-static)
  endif()
endif()

# Ensure HDF5::HDF5_HL is available for linkage
if(NOT TARGET HDF5::HDF5_HL)
  message(STATUS "Attempting to find HDF5 HL target to alias...")
  if(TARGET hdf5::hdf5_hl-static)
     add_library(HDF5::HDF5_HL ALIAS hdf5::hdf5_hl-static)
  elseif(TARGET hdf5::hdf5_hl-shared)
     add_library(HDF5::HDF5_HL ALIAS hdf5::hdf5_hl-shared)
  elseif(TARGET hdf5::hdf5_hl)
     add_library(HDF5::HDF5_HL ALIAS hdf5::hdf5_hl)
  elseif(TARGET hdf5_hl-static)
     add_library(HDF5::HDF5_HL ALIAS hdf5_hl-static)
  elseif(TARGET hdf5_hl-shared)
     add_library(HDF5::HDF5_HL ALIAS hdf5_hl-shared)
  endif()
endif()

# 4. HighFive
set(HIGHFIVE_USE_BOOST OFF CACHE BOOL "Do not use Boost")
set(HIGHFIVE_EXAMPLES OFF CACHE BOOL "Do not build examples")
set(HIGHFIVE_UNIT_TESTS OFF CACHE BOOL "Do not build tests")
set(HIGHFIVE_USE_INSTALL OFF CACHE BOOL "Do not use install")
set(HIGHFIVE_FIND_HDF5 OFF CACHE BOOL "Do not let HighFive find HDF5" FORCE)

FetchContent_Declare(
  HighFive
  GIT_REPOSITORY https://github.com/highfive-devs/highfive.git
  GIT_TAG v3.3.0
)
FetchContent_MakeAvailable(HighFive)

# 5. GoogleTest (only fetch, but don't enable testing here)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.16.0
)
FetchContent_MakeAvailable(googletest)
