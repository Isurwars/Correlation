# -----------------------------------------------------------
# Dependencies
# -----------------------------------------------------------
include(FetchContent)

set(BUILD_SHARED_LIBS ON CACHE BOOL "Force shared libraries")

# Save original BUILD_TESTING cache state if it exists
get_property(BUILD_TESTING_EXISTS CACHE BUILD_TESTING PROPERTY VALUE SET)
if(BUILD_TESTING_EXISTS)
  get_property(ORIG_BUILD_TESTING CACHE BUILD_TESTING PROPERTY VALUE)
  get_property(ORIG_BUILD_TESTING_TYPE CACHE BUILD_TESTING PROPERTY TYPE)
  get_property(ORIG_BUILD_TESTING_HELP CACHE BUILD_TESTING PROPERTY HELPSTRING)
endif()

# Force BUILD_TESTING to OFF for all dependencies to avoid building their tests
set(BUILD_TESTING OFF CACHE BOOL "Disable testing for dependencies" FORCE)

# 1. TBB
find_package(TBB QUIET)
if (TBB_FOUND)
  message(STATUS "Found TBB: ${TBB_DIR} (Version: ${TBB_VERSION})")
else()
  message(STATUS "TBB not found. Downloading TBB from GitHub...")
  FetchContent_Declare(
    TBB
    GIT_REPOSITORY https://github.com/uxlfoundation/oneTBB.git
    GIT_TAG v2023.0.0  
  )
  set(TBB_TEST OFF CACHE BOOL "Disable TBB tests" FORCE)
  FetchContent_MakeAvailable(TBB)
endif()

# 2. Slint
if(BUILD_GUI)
  find_package(Slint QUIET)
  if(Slint_FOUND)
    message(STATUS "Found Slint: ${Slint_DIR} (Version: ${Slint_VERSION})")
  else()
    set(SLINT_STYLE "material-dark" CACHE STRING "Slint style to use")
    message(STATUS "Slint not found. Downloading Slint from GitHub...")

    # Temporarily disable BUILD_SHARED_LIBS so Slint is built statically
    set(TEMP_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Force shared libraries" FORCE)

    FetchContent_Declare(
      Slint
      GIT_REPOSITORY https://github.com/slint-ui/slint.git
      GIT_TAG v1.16.1
      SOURCE_SUBDIR api/cpp
    )
    set(SLINT_FEATURE_JEMALLOC OFF CACHE BOOL "Disable jemalloc on macOS" FORCE)
    FetchContent_MakeAvailable(Slint)

    # Restore BUILD_SHARED_LIBS
    set(BUILD_SHARED_LIBS ${TEMP_BUILD_SHARED_LIBS} CACHE BOOL "Force shared libraries" FORCE)
  endif()

  # Transitively propagate platform dependencies for statically built Slint
  foreach(_target IN ITEMS Slint::Slint Slint slint_cpp slint_cpp-static)
    if(TARGET ${_target})
      set(_real_target "${_target}")
      get_target_property(_aliased_target ${_target} ALIASED_TARGET)
      if(_aliased_target)
        set(_real_target "${_aliased_target}")
      endif()

      get_target_property(_target_type ${_real_target} TYPE)
      if(NOT "${_target_type}" STREQUAL "UTILITY")
        if(APPLE)
          foreach(_fw IN ITEMS
              OpenGL CoreVideo Cocoa Carbon IOKit QuartzCore
              AppKit CoreGraphics Metal CoreFoundation Foundation Security)
            find_library(_${_fw}_FW ${_fw} REQUIRED)
            set_property(TARGET ${_real_target} APPEND PROPERTY INTERFACE_LINK_LIBRARIES "${_${_fw}_FW}")
            unset(_${_fw}_FW CACHE)
          endforeach()
        elseif(WIN32)
          set_property(TARGET ${_real_target} APPEND PROPERTY INTERFACE_LINK_LIBRARIES
            Imm32 Comctl32 dwrite propsys opengl32 dwmapi uxtheme
          )
        else()
          find_package(Fontconfig REQUIRED)
          set_property(TARGET ${_real_target} APPEND PROPERTY INTERFACE_LINK_LIBRARIES "${Fontconfig_LIBRARIES}")
        endif()
      endif()
    endif()
  endforeach()
endif()

# 3. HDF5
if(BUILD_WITH_HDF5)
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
    set(HDF5_PACK_EXAMPLES OFF CACHE BOOL "Pack HDF5 Examples" FORCE)
    set(H5_HAVE_C99_COMPLEX_NUMBERS OFF CACHE BOOL "Disable C99 Complex numbers for MSVC" FORCE)
    set(H5_HAVE_COMPLEX_NUMBERS OFF CACHE BOOL "Disable Complex numbers for MSVC" FORCE)
    set(H5_HAVE_COMPLEX_H OFF CACHE BOOL "Disable complex.h for MSVC" FORCE)
    set(HDF5_ENABLE_NONSTANDARD_FEATURE_FLOAT16 OFF CACHE BOOL "Disable _Float16 support" FORCE)
    set(HDF5_ENABLE_NONSTANDARD_FEATURE_COMPLEX OFF CACHE BOOL "Disable _Complex support" FORCE)
    set(HDF5_ENABLE_NONSTANDARD_FEATURE_COMPLEX_H OFF CACHE BOOL "Disable _Complex_H support" FORCE)
    # Shared specific options
    set(HDF5_BUILD_SHARED_LIBS ON CACHE BOOL "Build HDF5 Shared Library" FORCE)
    set(HDF5_USE_STATIC_LIBRARIES OFF CACHE BOOL "Use HDF5 Static Libraries" FORCE)
    set(HDF5_BUILD_STATIC_TOOLS OFF CACHE BOOL "Build HDF5 Static Tools" FORCE)


    FetchContent_Declare(
      HDF5
      GIT_REPOSITORY https://github.com/HDFGroup/hdf5.git
      GIT_TAG 2.1.1
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
    if(TARGET hdf5::hdf5_hl-shared)
       add_library(HDF5::HDF5_HL ALIAS hdf5::hdf5_hl-shared)
    elseif(TARGET hdf5_hl-shared)
       add_library(HDF5::HDF5_HL ALIAS hdf5_hl-shared)
    elseif(TARGET hdf5::hdf5_hl)
       add_library(HDF5::HDF5_HL ALIAS hdf5::hdf5_hl)
    elseif(TARGET hdf5::hdf5_hl-static)
       add_library(HDF5::HDF5_HL ALIAS hdf5::hdf5_hl-static)
    elseif(TARGET hdf5_hl-static)
       add_library(HDF5::HDF5_HL ALIAS hdf5_hl-static)
    endif()
  endif()
endif()

# 4. HighFive
if(BUILD_WITH_HDF5)
  set(HIGHFIVE_USE_BOOST OFF CACHE BOOL "Do not use Boost")
  set(HIGHFIVE_EXAMPLES OFF CACHE BOOL "Do not build examples")
  set(HIGHFIVE_UNIT_TESTS OFF CACHE BOOL "Do not build tests")
  set(HIGHFIVE_USE_INSTALL OFF CACHE BOOL "Do not use install")
  set(HIGHFIVE_FIND_HDF5 OFF CACHE BOOL "Do not let HighFive find HDF5" FORCE)

  find_package(HighFive QUIET)
  if(HighFive_FOUND)
    message(STATUS "Found HighFive: ${HighFive_DIR} (Version: ${HighFive_VERSION})")
  else()
    message(STATUS "HighFive not found. Downloading HighFive from GitHub...")
    set(HIGHFIVE_FIND_HDF5 OFF CACHE BOOL "Do not let HighFive find HDF5" FORCE)

    FetchContent_Declare(
      HighFive
      GIT_REPOSITORY https://github.com/highfive-devs/highfive.git
      GIT_TAG v3.3.0
    )
    FetchContent_MakeAvailable(HighFive)
  endif()
endif()

# 5. GoogleTest (only fetch, but don't enable testing here)
if(ORIG_BUILD_TESTING)
  find_package(GTest QUIET)
  if(GTest_FOUND)
    message(STATUS "Found GTest: ${GTest_DIR} (Version: ${GTest_VERSION})")
  else()
    message(STATUS "GTest not found. Downloading GoogleTest from GitHub...")
    FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest.git
      GIT_TAG v1.17.0
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
endif()

# 6. Arrow/Parquet
if(BUILD_WITH_ARROW)
  find_package(Arrow QUIET)
  find_package(Parquet QUIET)
  if (Arrow_FOUND AND Parquet_FOUND)
    message(STATUS "Found Arrow: ${Arrow_DIR} (Version: ${Arrow_VERSION})")
    message(STATUS "Found Parquet: ${Parquet_DIR} (Version: ${Parquet_VERSION})")

    # Create ALIAS targets so the rest of the project can just link 'arrow_shared' and 'parquet_shared'
    if (NOT TARGET arrow_shared)
      if (TARGET Arrow::arrow_shared)
        add_library(arrow_shared ALIAS Arrow::arrow_shared)
      elseif (TARGET arrow_static)
        add_library(arrow_shared ALIAS arrow_static)
      elseif (TARGET Arrow::arrow_static)
        add_library(arrow_shared ALIAS Arrow::arrow_static)
      endif()
    endif()

    if (NOT TARGET parquet_shared)
      if (TARGET Parquet::parquet_shared)
        add_library(parquet_shared ALIAS Parquet::parquet_shared)
      elseif (TARGET parquet_static)
        add_library(parquet_shared ALIAS parquet_static)
      elseif (TARGET Parquet::parquet_static)
        add_library(parquet_shared ALIAS Parquet::parquet_static)
      endif()
    endif()
  else()
    message(STATUS "Arrow/Parquet not found. Downloading Arrow from GitHub...")
    FetchContent_Declare(
      arrow
      GIT_REPOSITORY https://github.com/apache/arrow.git
      GIT_TAG        apache-arrow-23.0.1 # Use the latest stable version
      SOURCE_SUBDIR  cpp                 # We only need the C++ source
    )

    # Configure Arrow Build Options (Crucial for speed)
    # We disable things you don't need for a simple exporter to save compile time.
    # AUTO mode lets Arrow find system ZLIB (always available on Linux/macOS) without
    # conflicting with HDF5 FetchContent which also defines ZLIB::ZLIB.
    # ARROW_WITH_BOOST=OFF is the explicit kill switch that prevents Boost use.
    set(ARROW_DEPENDENCY_SOURCE "AUTO" CACHE INTERNAL "")
    set(Thrift_SOURCE "BUNDLED" CACHE STRING "Build bundled Thrift" FORCE)
    set(xsimd_SOURCE "BUNDLED" CACHE STRING "Build bundled xsimd" FORCE)
    set(RapidJSON_SOURCE "BUNDLED" CACHE STRING "Build bundled RapidJSON" FORCE)
    set(Boost_SOURCE "BUNDLED" CACHE STRING "Build bundled Boost" FORCE)
    set(ARROW_WITH_BOOST OFF CACHE INTERNAL "") # Never use Boost
    set(ARROW_PARQUET ON CACHE INTERNAL "")
    set(ARROW_WITH_SNAPPY OFF CACHE INTERNAL "")
    set(ARROW_WITH_ZSTD OFF CACHE INTERNAL "")
    set(ARROW_WITH_LZ4 OFF CACHE INTERNAL "")
    set(ARROW_WITH_BROTLI OFF CACHE INTERNAL "")
    set(ARROW_WITH_BZ2 OFF CACHE INTERNAL "")
    set(ARROW_BUILD_STATIC OFF CACHE INTERNAL "")
    set(ARROW_BUILD_SHARED ON CACHE INTERNAL "")
    set(ARROW_COMPUTE OFF CACHE INTERNAL "")
    set(ARROW_CSV OFF CACHE INTERNAL "")
    set(ARROW_JSON OFF CACHE INTERNAL "")
    set(ARROW_DATASET OFF CACHE INTERNAL "")
    set(ARROW_FILESYSTEM ON CACHE INTERNAL "")
    set(ARROW_IPC ON CACHE INTERNAL "")
    set(ARROW_BUILD_TESTS OFF CACHE INTERNAL "")
    set(ARROW_BUILD_BENCHMARKS OFF CACHE INTERNAL "")
    set(ARROW_SIMD_LEVEL "NONE" CACHE STRING "Arrow SIMD Level" FORCE)

    # Arrow pulls in rapidjson which has a broken CMake < 3.5 minimum version check for CMake 4.0+.
    # We enforce a policy version minimum before making it available to avoid errors.
    set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "" FORCE)

    FetchContent_MakeAvailable(arrow)

    # Create alias targets for the shared/static libraries built from source
    if (NOT TARGET arrow_shared)
      if (TARGET Arrow::arrow_shared)
        add_library(arrow_shared ALIAS Arrow::arrow_shared)
      elseif (TARGET arrow_static)
        add_library(arrow_shared ALIAS arrow_static)
      elseif (TARGET Arrow::arrow_static)
        add_library(arrow_shared ALIAS Arrow::arrow_static)
      endif()
    endif()

    if (NOT TARGET parquet_shared)
      if (TARGET Parquet::parquet_shared)
        add_library(parquet_shared ALIAS Parquet::parquet_shared)
      elseif (TARGET parquet_static)
        add_library(parquet_shared ALIAS parquet_static)
      elseif (TARGET Parquet::parquet_static)
        add_library(parquet_shared ALIAS Parquet::parquet_static)
      endif()
    endif()

    # Arrow targets built from source don't set the correct INCLUDE directories by default
    # We manually expose source and generated header folders.
    if (TARGET arrow_static)
      target_include_directories(arrow_static INTERFACE 
        $<BUILD_INTERFACE:${arrow_SOURCE_DIR}/cpp/src>
        $<BUILD_INTERFACE:${arrow_BINARY_DIR}/src>
      )
    elseif (TARGET arrow_shared)
      target_include_directories(arrow_shared INTERFACE 
        $<BUILD_INTERFACE:${arrow_SOURCE_DIR}/cpp/src>
        $<BUILD_INTERFACE:${arrow_BINARY_DIR}/src>
      )
    endif()

    if (TARGET parquet_static)
      target_include_directories(parquet_static INTERFACE 
        $<BUILD_INTERFACE:${arrow_SOURCE_DIR}/cpp/src>
        $<BUILD_INTERFACE:${arrow_BINARY_DIR}/src>
      )
    elseif (TARGET parquet_shared)
      target_include_directories(parquet_shared INTERFACE 
        $<BUILD_INTERFACE:${arrow_SOURCE_DIR}/cpp/src>
        $<BUILD_INTERFACE:${arrow_BINARY_DIR}/src>
      )
    endif()
  endif()
endif()

# 7. pybind11
if(BUILD_PYTHON_BINDINGS)
  find_package(pybind11 QUIET)
  if (pybind11_FOUND)
    message(STATUS "Found pybind11: ${pybind11_DIR} (Version: ${pybind11_VERSION})")
  else()
    message(STATUS "pybind11 not found. Downloading pybind11 from GitHub...")
    set(PYBIND11_FINDPYTHON NEW CACHE BOOL "Use new FindPython discovery")
    FetchContent_Declare(
      pybind11
      GIT_REPOSITORY https://github.com/pybind/pybind11.git
      GIT_TAG        v3.0.4
    )
    FetchContent_MakeAvailable(pybind11)
  endif()
endif()

# 8. CLI11
find_package(CLI11 QUIET)
if (CLI11_FOUND)
  message(STATUS "Found CLI11: ${CLI11_DIR}")
else()
  message(STATUS "CLI11 not found. Downloading CLI11 from GitHub...")
  FetchContent_Declare(
    cli11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG        v2.4.2
  )
  FetchContent_MakeAvailable(cli11)
endif()

# 9. nativefiledialog-extended
if(BUILD_GUI)
  message(STATUS "Downloading nativefiledialog-extended from GitHub...")
  set(NFD_BUILD_TESTS OFF CACHE BOOL "Disable NFD tests" FORCE)

  FetchContent_Declare(
    nfd
    GIT_REPOSITORY https://github.com/btzy/nativefiledialog-extended.git
    GIT_TAG        v1.2.1
  )
  FetchContent_MakeAvailable(nfd)
endif()

# 10. FFT Library Configuration
set(CORRELATION_FFT_BACKEND "Auto" CACHE STRING "FFT Backend to use (Auto, FFTW3, MKL, Fallback)")
set_property(CACHE CORRELATION_FFT_BACKEND PROPERTY STRINGS "Auto" "FFTW3" "MKL" "Fallback")

set(CORRELATION_USE_FFTW3 OFF)
set(CORRELATION_USE_MKL OFF)

if(CORRELATION_FFT_BACKEND STREQUAL "Auto")
  # By default, try to find FFTW3. If not found, fall back silently.
  find_package(FFTW3 QUIET)
  if(NOT FFTW3_FOUND)
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
      pkg_check_modules(FFTW3 QUIET IMPORTED_TARGET fftw3)
    endif()
  endif()
  if(FFTW3_FOUND OR FFTW3_LIBRARIES)
    set(CORRELATION_USE_FFTW3 ON)
    message(STATUS "FFT Backend: FFTW3 (Auto)")
  else()
    message(STATUS "FFT Backend: Fallback (Cooley-Tukey)")
  endif()

elseif(CORRELATION_FFT_BACKEND STREQUAL "FFTW3")
  find_package(FFTW3 QUIET)
  if(NOT FFTW3_FOUND)
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
      pkg_check_modules(FFTW3 QUIET IMPORTED_TARGET fftw3)
    endif()
  endif()
  if(FFTW3_FOUND OR FFTW3_LIBRARIES)
    set(CORRELATION_USE_FFTW3 ON)
    message(STATUS "FFT Backend: FFTW3")
  else()
    message(WARNING "FFTW3 requested but not found. Falling back to Cooley-Tukey.")
  endif()

elseif(CORRELATION_FFT_BACKEND STREQUAL "MKL")
  find_package(MKL QUIET)
  if(NOT MKL_FOUND)
    find_package(oneMKL QUIET)
  endif()
  if(NOT MKL_FOUND AND NOT oneMKL_FOUND)
    find_path(MKL_INCLUDE_DIR mkl_dfti.h)
    find_library(MKL_CORE_LIBRARY mkl_core)
    if(MKL_INCLUDE_DIR AND MKL_CORE_LIBRARY)
      set(MKL_FOUND TRUE)
      set(MKL_LIBRARIES ${MKL_CORE_LIBRARY})
      set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
    endif()
  endif()

  if(MKL_FOUND OR oneMKL_FOUND)
    set(CORRELATION_USE_MKL ON)
    message(STATUS "FFT Backend: Intel MKL")
  else()
    message(WARNING "Intel MKL requested but not found. Falling back to Cooley-Tukey.")
  endif()

else()
  message(STATUS "FFT Backend: Fallback (Cooley-Tukey)")
endif()

# Create correlation_fft interface target
add_library(correlation_fft INTERFACE)
if(CORRELATION_USE_FFTW3)
  target_compile_definitions(correlation_fft INTERFACE CORRELATION_USE_FFTW3)
  if(TARGET PkgConfig::FFTW3)
    target_link_libraries(correlation_fft INTERFACE PkgConfig::FFTW3)
  else()
    if(FFTW3_INCLUDE_DIRS)
      target_include_directories(correlation_fft INTERFACE ${FFTW3_INCLUDE_DIRS})
    elseif(FFTW3_INCLUDE_DIR)
      target_include_directories(correlation_fft INTERFACE ${FFTW3_INCLUDE_DIR})
    endif()
    target_link_libraries(correlation_fft INTERFACE ${FFTW3_LIBRARIES})
    if(FFTW3_LIBRARY_DIRS)
      target_link_directories(correlation_fft INTERFACE ${FFTW3_LIBRARY_DIRS})
    endif()
  endif()
elseif(CORRELATION_USE_MKL)
  target_compile_definitions(correlation_fft INTERFACE CORRELATION_USE_MKL)
  if(TARGET MKL::MKL)
    target_link_libraries(correlation_fft INTERFACE MKL::MKL)
  elseif(TARGET oneMKL::oneMKL)
    target_link_libraries(correlation_fft INTERFACE oneMKL::oneMKL)
  else()
    if(MKL_INCLUDE_DIRS)
      target_include_directories(correlation_fft INTERFACE ${MKL_INCLUDE_DIRS})
    endif()
    target_link_libraries(correlation_fft INTERFACE ${MKL_LIBRARIES})
  endif()
endif()


# 11. voro++
message(STATUS "Downloading voro++ from GitHub...")
set(VORO_BUILD_EXAMPLES OFF CACHE BOOL "Disable voro++ examples" FORCE)
set(VORO_BUILD_CMD_LINE OFF CACHE BOOL "Disable voro++ command line" FORCE)
set(VORO_ENABLE_DOXYGEN OFF CACHE BOOL "Disable voro++ doxygen" FORCE)

# voro++ does not export symbols and cannot be built as a DLL under MSVC.
# We build it statically on Windows, but let it build as a shared library on macOS/Linux.
if(WIN32)
  set(TEMP_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "Force shared libraries" FORCE)
endif()

FetchContent_Declare(
  voro
  GIT_REPOSITORY https://github.com/chr1shr/voro.git
  GIT_TAG        b0dac575a47af0f90b5b100e6dc199a493c7cb83
)
FetchContent_MakeAvailable(voro)

if(WIN32)
  set(BUILD_SHARED_LIBS ${TEMP_BUILD_SHARED_LIBS} CACHE BOOL "Force shared libraries" FORCE)
endif()

# Ensure voro++ is built with position-independent code (PIC) since it might be linked into shared libraries/modules
if(TARGET voro++)
  set_target_properties(voro++ PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

# Restore original BUILD_TESTING cache state
if(BUILD_TESTING_EXISTS)
  set(BUILD_TESTING "${ORIG_BUILD_TESTING}" CACHE ${ORIG_BUILD_TESTING_TYPE} "${ORIG_BUILD_TESTING_HELP}" FORCE)
else()
  unset(BUILD_TESTING CACHE)
endif()