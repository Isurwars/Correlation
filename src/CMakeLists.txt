# Create a STATIC library from source files
add_library(correlation_lib STATIC
    Atom.cpp
    Cell.cpp
    Constants.cpp
    ReadFiles.cpp
    Smoothing.cpp
    StructureFactor.cpp
    WriteFiles.cpp
)

# Define executable
add_executable(correlation
    Atom.cpp
    Cell.cpp
    Constants.cpp
    ReadFiles.cpp
    Smoothing.cpp
    StructureFactor.cpp
    WriteFiles.cpp
    main.cpp
)

# PSTL TBB backend definition (target-specific)
target_compile_definitions(correlation PRIVATE _PSTL_PAR_BACKEND_TBB)

# Link TBB libraries
target_link_libraries(correlation PRIVATE correlation_lib
    TBB::tbb
    TBB::tbbmalloc
)

# Modern include directory handling (if needed)
target_include_directories(correlation PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)
