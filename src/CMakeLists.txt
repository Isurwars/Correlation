# Correlation - Liquid and Amorphous Solid Analysis Tool

# Readers as an OBJECT library to ensure static registration triggers
add_library(readers_obj OBJECT
    readers/ArcReader.cpp
    readers/CarReader.cpp
    readers/CellReader.cpp
    readers/CastepMdReader.cpp
    readers/CifReader.cpp
    readers/LammpsDumpReader.cpp
    readers/OnetepDatReader.cpp
    readers/OutmolReader.cpp
    readers/FileReader.cpp
    readers/ReaderFactory.cpp
)
target_include_directories(readers_obj PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

if (CORRELATION_SIMD_FLAGS)
    # Convert string of flags to a list for target_compile_options
    separate_arguments(SIMD_FLAGS_LIST NATIVE_COMMAND "${CORRELATION_SIMD_FLAGS}")
    target_compile_options(readers_obj PRIVATE ${SIMD_FLAGS_LIST})
endif()

# Calculators as an OBJECT library to ensure static registration triggers
add_library(calculators_obj OBJECT
    calculators/AngleCalculator.cpp
    calculators/DihedralCalculator.cpp
    calculators/MotifFinder.cpp
    calculators/CNCalculator.cpp
    calculators/DistanceCalculator.cpp
    calculators/RDFCalculator.cpp
    calculators/PADCalculator.cpp
    calculators/DADCalculator.cpp
    calculators/RDCalculator.cpp
    calculators/VACFCalculator.cpp
    calculators/VDOSCalculator.cpp
    calculators/StructureFactorCalculator.cpp
    calculators/XRDCalculator.cpp
    calculators/SteinhardtCalculator.cpp
    calculators/MSDCalculator.cpp
    calculators/CalculatorFactory.cpp
)
target_include_directories(calculators_obj PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(calculators_obj PRIVATE TBB::tbb)
target_compile_definitions(calculators_obj PRIVATE _PSTL_PAR_BACKEND_TBB)
if (CORRELATION_SIMD_FLAGS)
    separate_arguments(SIMD_FLAGS_LIST NATIVE_COMMAND "${CORRELATION_SIMD_FLAGS}")
    target_compile_options(calculators_obj PRIVATE ${SIMD_FLAGS_LIST})
endif()

# Writers as an OBJECT library to ensure static registration triggers
add_library(writers_obj OBJECT
    writers/ArrowWriter.cpp
    writers/CSVWriter.cpp
    writers/HDF5Writer.cpp
    writers/FileWriter.cpp
    writers/WriterFactory.cpp
)
target_include_directories(writers_obj PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(writers_obj PRIVATE HighFive HDF5::HDF5 HDF5::HDF5_HL arrow_shared parquet_shared)
if (CORRELATION_SIMD_FLAGS)
    separate_arguments(SIMD_FLAGS_LIST NATIVE_COMMAND "${CORRELATION_SIMD_FLAGS}")
    target_compile_options(writers_obj PRIVATE ${SIMD_FLAGS_LIST})
endif()

# Add sources to library
add_library(correlation_lib STATIC
    core/Cell.cpp
    analysis/DistributionFunctions.cpp
    analysis/DynamicsAnalyzer.cpp
    core/NeighborGraph.cpp
    analysis/StructureAnalyzer.cpp
    core/Trajectory.cpp
    analysis/TrajectoryAnalyzer.cpp
    $<TARGET_OBJECTS:readers_obj>
    $<TARGET_OBJECTS:calculators_obj>
    $<TARGET_OBJECTS:writers_obj>
)

# Include directories
target_include_directories(correlation_lib PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_definitions(correlation_lib PUBLIC _PSTL_PAR_BACKEND_TBB)

if (CORRELATION_SIMD_FLAGS)
  # Ensure the entire library uses consistent SIMD flags for ABI safety
  separate_arguments(SIMD_FLAGS_LIST NATIVE_COMMAND "${CORRELATION_SIMD_FLAGS}")
  target_compile_options(correlation_lib PUBLIC ${SIMD_FLAGS_LIST})
endif()

target_link_libraries(correlation_lib PRIVATE TBB::tbb TBB::tbbmalloc)
target_link_libraries(correlation_lib PUBLIC HighFive HDF5::HDF5 HDF5::HDF5_HL arrow_shared parquet_shared)

#------------------------------------------------------------------
# Executable Target (correlation)
#------------------------------------------------------------------
set(CORRELATION_SOURCES
    main.cpp
    app/AppBackend.cpp
    app/AppController.cpp
)

if (WIN32)
    list(APPEND CORRELATION_SOURCES ../packaging/resource.rc)
    add_executable(correlation WIN32 ${CORRELATION_SOURCES})
else()
    add_executable(correlation ${CORRELATION_SOURCES})
endif()

# Link correlation_lib to correlation
set(CORRELATION_LIBS correlation_lib Slint::Slint)

if (CORRELATION_SIMD_FLAGS)
    # Main executable must also use SIMD flags to avoid ABI issues with Vector3
    separate_arguments(SIMD_FLAGS_LIST NATIVE_COMMAND "${CORRELATION_SIMD_FLAGS}")
    target_compile_options(correlation PRIVATE ${SIMD_FLAGS_LIST})
endif()

if (WIN32)
    list(APPEND CORRELATION_LIBS Imm32 Comctl32 dwrite propsys opengl32 dwmapi uxtheme)
endif()

# We link readers_obj directly to ensure static registration
target_link_libraries(correlation PRIVATE ${CORRELATION_LIBS} readers_obj writers_obj calculators_obj)

if (WIN32)
    target_include_directories(correlation PRIVATE ${PROJECT_SOURCE_DIR})
endif()

# Add the include directory for the generated slint header
target_include_directories(correlation PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}
)

slint_target_sources(correlation
  ../ui/AppWindow.slint
  LIBRARY_PATHS material=${CMAKE_CURRENT_SOURCE_DIR}/../material-1.0/material.slint
)

set_property(TARGET correlation PROPERTY SLINT_EMBED_RESOURCES "embed-files")

# Copy DLLs to the executable directory on Windows
bundle_windows_dlls(correlation)
install_windows_dlls(correlation ${CMAKE_INSTALL_BINDIR})

if (NOT WIN32)
    set_target_properties(correlation PROPERTIES
        INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}"
    )
endif()
