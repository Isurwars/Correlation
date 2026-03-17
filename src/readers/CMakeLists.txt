# Correlation - Liquid and Amorphous Solid Analysis Tool

# Readers as an OBJECT library to ensure static registration triggers
add_library(readers_obj OBJECT
    ArcReader.cpp
    CarReader.cpp
    CellReader.cpp
    CastepMdReader.cpp
    CifReader.cpp
    LammpsDumpReader.cpp
    OnetepDatReader.cpp
    OutmolReader.cpp
    ReaderFactory.cpp
)

target_include_directories(readers_obj PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

# Apply SIMD flags to translation units that use SIMD kernels if needed
# (Currently readers don't seem to use them, but we could add them here if they did)
