cmake_minimum_required(VERSION 3.20)
project(Correlation VERSION 1.2 LANGUAGES CXX)

# Set the C++ standard as a target property, not a global variable.
set(CXX_STANDARD 23)

# Set default install path
set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation directory")

# Find dependencies
find_package(TBB REQUIRED COMPONENTS tbb tbbmalloc)
find_package(Slint REQUIRED)

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
    src/ProgramOptions.cpp
    src/Smoothing.cpp
    src/StructureFactor.cpp
)

# Include directories
target_include_directories(correlation_lib PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include> # Path when installed
)

target_compile_features(correlation_lib PUBLIC cxx_std_${CXX_STANDARD})
target_compile_definitions(correlation_lib PUBLIC _PSTL_PAR_BACKEND_TBB)

# Link dependencies to libraries
target_link_libraries(correlation_lib PRIVATE
    TBB::tbb
    TBB::tbbmalloc
)

#------------------------------------------------------------------
# Executable Target (correlation)
#------------------------------------------------------------------
add_executable(correlation src/main.cpp)

# Link correlation_lib to correlation
target_link_libraries(correlation PRIVATE correlation_lib)

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

include(GoogleTest)
gtest_discover_tests(correlation_tests)

#-----------------------------------------------------------
# Install setup
#-----------------------------------------------------------

# Install the executable and the library.
install(TARGETS correlation correlation_lib
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

# Install the public headers so that other projects can use this library.
install(DIRECTORY include/ DESTINATION include)
