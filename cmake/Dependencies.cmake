# -----------------------------------------------------------
# Dependencies
# -----------------------------------------------------------
include(FetchContent)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "Force static libraries" FORCE)
set(TBB_USE_STATIC_LIBS ON CACHE BOOL "Use static TBB libraries" FORCE)

# 1. TBB
find_package(TBB QUIET)
if (TBB_FOUND)
  message(STATUS "Found TBB: ${TBB_DIR} (Version: ${TBB_VERSION})")
else()
  message(STATUS "TBB not found. Downloading TBB from GitHub...")
  FetchContent_Declare(
    TBB
    GIT_REPOSITORY https://github.com/uxlfoundation/oneTBB.git
    GIT_TAG v2022.3.0  
  )
  FetchContent_MakeAvailable(TBB)
endif()

# 2. Slint
find_package(Slint QUIET)
if(Slint_FOUND)
  message(STATUS "Found Slint: ${Slint_DIR} (Version: ${Slint_VERSION})")
else()
  set(SLINT_STYLE "material-dark" CACHE STRING "Slint style to use")
  message(STATUS "Slint not found. Downloading Slint from GitHub...")
  FetchContent_Declare(
    Slint
    GIT_REPOSITORY https://github.com/slint-ui/slint.git
    GIT_TAG v1.15.0
    SOURCE_SUBDIR api/cpp
  )
  FetchContent_MakeAvailable(Slint)
endif()

# 3. HDF5
find_package(HDF5 COMPONENTS C CXX HL QUIET)
if(HDF5_FOUND)
  message(STATUS "Found HDF5: ${HDF5_DIR} (Version: ${HDF5_VERSION})")
else()
  message(STATUS "HDF5 not found. Downloading HDF5 from GitHub...")
  # HDF5 options for FetchContent
  set(HDF5_BUILD_CPP_LIB ON CACHE BOOL "Build HDF5 C++ Library" FORCE)
  set(HDF5_BUILD_EXAMPLES OFF CACHE BOOL "Build HDF5 Examples" FORCE)
  set(HDF5_BUILD_HL_LIB ON CACHE BOOL "Build HDF5 High Level Library" FORCE)
  set(HDF5_BUILD_TOOLS OFF CACHE BOOL "Build HDF5 Tools" FORCE)
  set(HDF5_BUILD_FORTRAN OFF CACHE BOOL "Build HDF5 Fortran" FORCE)
  set(HDF5_BUILD_JAVA OFF CACHE BOOL "Build HDF5 Java" FORCE)
  set(BUILD_TESTING OFF CACHE BOOL "Build HDF5 Tests" FORCE)
  set(HDF5_PACK_EXAMPLES OFF CACHE BOOL "Pack HDF5 Examples" FORCE)
  set(H5_HAVE_C99_COMPLEX_NUMBERS OFF CACHE BOOL "Disable C99 Complex numbers for MSVC" FORCE)
  set(H5_HAVE_COMPLEX_NUMBERS OFF CACHE BOOL "Disable Complex numbers for MSVC" FORCE)
  set(H5_HAVE_COMPLEX_H OFF CACHE BOOL "Disable complex.h for MSVC" FORCE)
  set(HDF5_ENABLE_NONSTANDARD_FEATURE_FLOAT16 OFF CACHE BOOL "Disable _Float16 support" FORCE)
  set(HDF5_ENABLE_NONSTANDARD_FEATURE_COMPLEX OFF CACHE BOOL "Disable _Complex support" FORCE)
  set(HDF5_ENABLE_NONSTANDARD_FEATURE_COMPLEX_H OFF CACHE BOOL "Disable _Complex_H support" FORCE)
  # Static specific options
  set(HDF5_BUILD_SHARED_LIBS OFF CACHE BOOL "Build HDF5 Shared Library" FORCE)
  set(HDF5_USE_STATIC_LIBRARIES ON CACHE BOOL "Use HDF5 Static Libraries" FORCE)
  set(HDF5_BUILD_STATIC_TOOLS OFF CACHE BOOL "Build HDF5 Static Tools" FORCE)


  FetchContent_Declare(
    HDF5
    GIT_REPOSITORY https://github.com/HDFGroup/hdf5.git
    GIT_TAG hdf5_2.0.0
  )
  FetchContent_MakeAvailable(HDF5)
endif()

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
find_package(HighFive QUIET)
if(HighFive_FOUND)
  message(STATUS "Found HighFive: ${HighFive_DIR} (Version: ${HighFive_VERSION})")
else()
  message(STATUS "HighFive not found. Downloading HighFive from GitHub...")
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
endif()

# 5. GoogleTest (only fetch, but don't enable testing here)
find_package(GTest QUIET)
if(GTest_FOUND)
  message(STATUS "Found GTest: ${GTest_DIR} (Version: ${GTest_VERSION})")
else()
  message(STATUS "GTest not found. Downloading GoogleTest from GitHub...")
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.16.0
  )
  FetchContent_MakeAvailable(googletest)

  # Check if targets exist and create aliases if needed
  if(NOT TARGET GTest::gtest)
    if(TARGET gtest)
      add_library(GTest::gtest ALIAS gtest)
    else()
      # Fallback if gtest target has a different name (e.g. GTest) in some versions, though 'gtest' is standard for FetchContent
      message(WARNING "Target 'gtest' not found after FetchContent. Adjust alias logic.")
    endif()
  endif()

  if(NOT TARGET GTest::gtest_main)
     if(TARGET gtest_main)
       add_library(GTest::gtest_main ALIAS gtest_main)
     endif()
  endif()
endif()
