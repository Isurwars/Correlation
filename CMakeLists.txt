cmake_minimum_required(VERSION 3.20)
project(Correlation)

# Enable C++23
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find and link TBB (critical for PSTL parallel execution)
find_package(TBB REQUIRED COMPONENTS tbb tbbmalloc)  # Include both TBB components

# Define PSTL's TBB backend explicitly (required for some compilers)
add_compile_definitions(_PSTL_PAR_BACKEND_TBB)  # Force PSTL to use TBB

include_directories(src)
add_subdirectory(src)

# Link TBB to your target
target_link_libraries(correlation PRIVATE
    TBB::tbb
    TBB::tbbmalloc  # Include if using scalable memory allocation
)

install(TARGETS correlation DESTINATION bin)
