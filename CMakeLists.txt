cmake_minimum_required(VERSION 3.20)
project(Correlation VERSION 1.1 LANGUAGES CXX)

# Set default install path (customize this path)
set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation directory")

# Enable C++23
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find TBB (critical for PSTL)
find_package(TBB REQUIRED COMPONENTS tbb tbbmalloc)
# Find Slint
find_package(Slint REQUIRED)  # Add Slint package

#-----------------------------------------------------------
# Library and executable setup
#-----------------------------------------------------------

# Collect source files
file(GLOB_RECURSE ALL_SOURCES CONFIGURE_DEPENDS "src/*.cpp")
# Filter out main.cpp for the library
list(FILTER ALL_SOURCES EXCLUDE REGEX ".*main\\.cpp$")

# Create STATIC library
add_library(correlation_lib STATIC ${ALL_SOURCES})

# Link libraries to executable
target_link_libraries(correlation_lib PRIVATE
    TBB::tbb
    TBB::tbbmalloc
)

# Define executable with just main.cpp
add_executable(correlation src/main.cpp)

# Link libraries to executable
target_link_libraries(correlation PRIVATE
    correlation_lib
)

# Set PSTL backend
target_compile_definitions(correlation_lib PRIVATE _PSTL_PAR_BACKEND_TBB)
target_compile_definitions(correlation PRIVATE _PSTL_PAR_BACKEND_TBB)

# Include directories
target_include_directories(correlation_lib PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)
target_include_directories(correlation PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

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

# Collect test sources
file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS "tests/*.cpp")


# Create single test executable
add_executable(correlation_tests ${TEST_SOURCES})

# Link with dependencies
target_link_libraries(correlation_tests PRIVATE
    correlation_lib
    gtest_main
)

# Add include directories
target_include_directories(correlation_tests PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

# Add PSTL definition
target_compile_definitions(correlation_tests PRIVATE _PSTL_PAR_BACKEND_TBB)

# Discover and add all tests
include(GoogleTest)
gtest_discover_tests(correlation_tests
    EXTRA_ARGS --gtest_output=xml:${CMAKE_BINARY_DIR}/test_results/
)

#-----------------------------------------------------------
# Install setup
#-----------------------------------------------------------

# Install rules (executable)
install(TARGETS correlation
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)
