cmake_minimum_required(VERSION 3.14)
project(Correlation VERSION 1.1 LANGUAGES CXX)

# Set default install path (customize this path)
set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation directory")

# Enable C++23
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find TBB (critical for PSTL)
find_package(TBB REQUIRED COMPONENTS tbb tbbmalloc)

# Include directory for headers
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add src subdirectory
add_subdirectory(src)


#-----------------------------------------------------------
# Testing setup (using doctest)
#-----------------------------------------------------------
# Enable testing framework
enable_testing()

include(FetchContent)
FetchContent_Declare(
  GoogleTest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.16.0
)

# Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(GoogleTest)

add_subdirectory(tests)


#-----------------------------------------------------------
# Install setup
#-----------------------------------------------------------

# Install rules (executable)
install(TARGETS correlation
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)
