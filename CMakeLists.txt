cmake_minimum_required(VERSION 3.21)
project(Correlation VERSION 1.3 LANGUAGES CXX)

# Set the C++ standard as a target property, not a global variable.
set(CXX_STANDARD 23)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set default install path
set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation directory")


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
    GIT_TAG release/1
    SOURCE_SUBDIR api/cpp
  )
  FetchContent_MakeAvailable(Slint)
endif (NOT Slint_FOUND)

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
if (WIN32)
add_executable(correlation WIN32
  src/main.cpp
  src/AppBackend.cpp
  src/AppController.cpp
  resource.rc
)
else()
add_executable(correlation
  src/main.cpp
  src/AppBackend.cpp
  src/AppController.cpp
)
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
  ui/app_window.slint
)

if (WIN32)
    add_custom_command(
      TARGET correlation POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_RUNTIME_DLLS:correlation>
      $<TARGET_FILE_DIR:correlation> COMMAND_EXPAND_LISTS)
endif()
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
)

# The test executable also links to the library and inherits its properties.
target_link_libraries(correlation_tests PRIVATE
    correlation_lib
    gtest_main
)

if (WIN32)
    add_custom_command(
      TARGET correlation_tests POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_RUNTIME_DLLS:correlation_tests>
      $<TARGET_FILE_DIR:correlation_tests> COMMAND_EXPAND_LISTS)
endif()

include(GoogleTest)
gtest_discover_tests(correlation_tests TEST_DISCOVERY_TIMEOUT 300)

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
                   "${PROJECT_BINARY_DIR}/correlation.desktop" COPYONLY)
    install(FILES "${PROJECT_BINARY_DIR}/correlation.desktop"
            DESTINATION share/applications)
endif()
