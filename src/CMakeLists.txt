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
    readers/ReaderFactory.cpp
)
target_include_directories(readers_obj PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

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
    calculators/SQCalculator.cpp
    calculators/XRDCalculator.cpp
    calculators/CalculatorFactory.cpp
)
target_include_directories(calculators_obj PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(calculators_obj PRIVATE TBB::tbb)
target_compile_definitions(calculators_obj PRIVATE _PSTL_PAR_BACKEND_TBB)

# Add sources to library
add_library(correlation_lib STATIC
    Cell.cpp
    DistributionFunctions.cpp
    DynamicsAnalyzer.cpp
    FileReader.cpp
    FileWriter.cpp
    writers/ArrowWriter.cpp
    writers/CSVWriter.cpp
    writers/HDF5Writer.cpp
    NeighborGraph.cpp
    StructureAnalyzer.cpp
    Trajectory.cpp
    TrajectoryAnalyzer.cpp
    $<TARGET_OBJECTS:readers_obj>
    $<TARGET_OBJECTS:calculators_obj>
)

# Include directories
target_include_directories(correlation_lib PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# target_compile_features(correlation_lib PUBLIC cxx_std_${CXX_STANDARD})
target_compile_definitions(correlation_lib PUBLIC _PSTL_PAR_BACKEND_TBB)

if (CORRELATION_SIMD_FLAGS)
  # Apply SIMD flags to translation units that use SIMD kernels
  set_source_files_properties(
      calculators/DistanceCalculator.cpp
      calculators/SQCalculator.cpp
      calculators/XRDCalculator.cpp
      calculators/RDFCalculator.cpp
      calculators/AngleCalculator.cpp
      PROPERTIES COMPILE_FLAGS "${CORRELATION_SIMD_FLAGS}"
  )
endif()

target_link_libraries(correlation_lib PRIVATE TBB::tbb TBB::tbbmalloc)
target_link_libraries(correlation_lib PUBLIC HighFive HDF5::HDF5 HDF5::HDF5_HL arrow_shared parquet_shared)

#------------------------------------------------------------------
# Executable Target (correlation)
#------------------------------------------------------------------
set(CORRELATION_SOURCES
    main.cpp
    AppBackend.cpp
    AppController.cpp
)

if (WIN32)
    list(APPEND CORRELATION_SOURCES ../resource.rc)
    add_executable(correlation WIN32 ${CORRELATION_SOURCES})
else()
    add_executable(correlation ${CORRELATION_SOURCES})
endif()

# Link correlation_lib to correlation
set(CORRELATION_LIBS correlation_lib Slint::Slint)

if (WIN32)
    list(APPEND CORRELATION_LIBS Imm32 Comctl32 dwrite propsys opengl32 dwmapi uxtheme)
endif()

# We link readers_obj directly to ensure static registration
target_link_libraries(correlation PRIVATE ${CORRELATION_LIBS} readers_obj)

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
if (WIN32)
    add_custom_command(TARGET correlation POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:correlation>" "$<TARGET_RUNTIME_DLLS:correlation>" "$<TARGET_FILE_DIR:correlation>"
        COMMAND_EXPAND_LISTS
    )
endif()

# Installation logic for the executable
install(TARGETS correlation
    RUNTIME DESTINATION bin
)

if (NOT WIN32)
    set_target_properties(correlation PROPERTIES
        INSTALL_RPATH "$ORIGIN/../lib"
    )
endif()
