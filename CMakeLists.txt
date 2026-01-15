cmake_minimum_required(VERSION 3.21)
project(Correlation VERSION 1.7 LANGUAGES CXX)

# Set the C++ standard as a target property, not a global variable.
set(CXX_STANDARD 23)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set default install path
set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation directory")

# -----------------------------------------------------------
# Dependencies
# -----------------------------------------------------------
# Find dependencies
find_package(TBB QUIET COMPONENTS tbb tbbmalloc)
if (NOT TBB_FOUND)
  message("TBB could not be located. Downloading it from GitHub and building it locally.")
  include(FetchContent)
  FetchContent_Declare(
    TBB
    GIT_REPOSITORY https://github.com/uxlfoundation/oneTBB.git
    GIT_TAG v2022.2.0
  )
  FetchContent_MakeAvailable(TBB)
endif()

find_package(Slint QUIET)
if (NOT Slint_FOUND)
  message("Slint could not be located in the CMake module search path. Downloading it from Git and building it locally")
  include(FetchContent)
  FetchContent_Declare(
    Slint
    GIT_REPOSITORY https://github.com/slint-ui/slint.git
    GIT_TAG v1.13.1
    SOURCE_SUBDIR api/cpp
  )
  FetchContent_MakeAvailable(Slint)
endif (NOT Slint_FOUND)


# -----------------------------------------------------------
# Helper Function for Bundling DLLs on Windows
# -----------------------------------------------------------
function(bundle_windows_dlls target_name)
    if(WIN32)
        add_custom_command(
            TARGET ${target_name} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_RUNTIME_DLLS:${target_name}>
            $<TARGET_FILE_DIR:${target_name}>
            COMMAND_EXPAND_LISTS
            COMMENT "Copying runtime DLLs for ${target_name}"
        )
    endif()
endfunction()


#-----------------------------------------------------------
# Library setup
#-----------------------------------------------------------

# Add sources to library
add_library(correlation_lib STATIC
    src/Cell.cpp
    src/DistributionFunctions.cpp
    src/FileIO.cpp
    src/FileWriter.cpp
    src/StructureAnalyzer.cpp
)

# Include directories
target_include_directories(correlation_lib PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(correlation_lib PUBLIC cxx_std_${CXX_STANDARD})
target_compile_definitions(correlation_lib PUBLIC _PSTL_PAR_BACKEND_TBB)
target_link_libraries(correlation_lib PRIVATE TBB::tbb TBB::tbbmalloc)

#------------------------------------------------------------------
# Executable Target (correlation)
#------------------------------------------------------------------
set(CORRELATION_SOURCES
    src/main.cpp
    src/AppBackend.cpp
    src/AppController.cpp
)

if (WIN32)
    list(APPEND CORRELATION_SOURCES resource.rc)
    add_executable(correlation WIN32 ${CORRELATION_SOURCES})
else()
    add_executable(correlation ${CORRELATION_SOURCES})
endif()


# Link correlation_lib to correlation
target_link_libraries(correlation
  correlation_lib
  Slint::Slint
)

# Add the include directory for the generated slint header
target_include_directories(correlation PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}
)

slint_target_sources(correlation
  ui/AppWindow.slint
  LIBRARY_PATHS material=${CMAKE_CURRENT_SOURCE_DIR}/material-1.0/material.slint
)

# Call the helper function to bundle DLLs on Windows
bundle_windows_dlls(correlation)

#-----------------------------------------------------------
# Testing setup (using GoogleTest)
#-----------------------------------------------------------
enable_testing()

# Modern GoogleTest fetch
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.16.0
)
FetchContent_MakeAvailable(googletest)

# Explicitly list test sources.
add_executable(correlation_tests
    tests/Atom_test.cpp
    tests/Cell_test.cpp
    tests/Distribution_Functions_tests.cpp
    tests/FileIO_test.cpp
    tests/FileWriter_test.cpp
    tests/StructureAnalyzer_test.cpp
    tests/AngleTests.cpp
)

# The test executable also links to the library and inherits its properties.
target_link_libraries(correlation_tests PRIVATE
    correlation_lib
    gtest_main
)

# Call the helper function to bundle DLLs on Windows
bundle_windows_dlls(correlation_tests)

include(GoogleTest)
gtest_discover_tests(correlation_tests DISCOVERY_TIMEOUT 60)

#-----------------------------------------------------------
# Install setup
#-----------------------------------------------------------

# Install the executable and the library.
install(TARGETS correlation
    RUNTIME DESTINATION bin
)

# -----------------------------------------------------------
# Install setup for Linux Desktop Integration
# -----------------------------------------------------------
if (NOT WIN32)
    # 1. Install the Icon (adjust path if you use a different icon directory)
    install(FILES "${PROJECT_SOURCE_DIR}/Images/Logo.png"
            DESTINATION share/icons/hicolor/scalable/apps
            RENAME correlation.png)

    # 2. Configure and Install the .desktop file
    configure_file("${PROJECT_SOURCE_DIR}/correlation.desktop"
                   "${PROJECT_BINARY_DIR}/correlation.desktop"
                   COPYONLY)
    install(FILES "${PROJECT_BINARY_DIR}/correlation.desktop"
            DESTINATION share/applications)
endif()
