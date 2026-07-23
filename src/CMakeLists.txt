# Correlation - Liquid and Amorphous Solid Analysis Tool

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Readers as an OBJECT library to ensure static registration triggers
add_library(readers_obj OBJECT
    readers/ArcReader.cpp
    readers/CarReader.cpp
    readers/CellReader.cpp
    readers/CastepMdReader.cpp
    readers/CifReader.cpp
    readers/GromacsReader.cpp
    readers/LammpsDumpReader.cpp
    readers/OnetepDatReader.cpp
    readers/OutmolReader.cpp
    readers/PdbReader.cpp
    readers/VaspReader.cpp
    readers/XdatcarReader.cpp
    readers/FileReader.cpp
    readers/ReaderFactory.cpp
    readers/QEReader.cpp
    readers/CP2KReader.cpp
    readers/XYZReader.cpp
)
target_include_directories(readers_obj PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(readers_obj PRIVATE correlation_simd_options)

# Calculators as an OBJECT library to ensure static registration triggers
add_library(calculators_obj OBJECT
    calculators/AngleCalculator.cpp
    calculators/CNACalculator.cpp
    calculators/ChiralityCalculator.cpp
    calculators/DihedralCalculator.cpp
    calculators/HBondCalculator.cpp
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
    calculators/SDFCalculator.cpp
    calculators/ClusterCalculator.cpp
    calculators/LocalEntropyCalculator.cpp
    calculators/VoronoiCalculator.cpp
    calculators/HyperuniformityCalculator.cpp
    calculators/MLIPCalculator.cpp
)

if(NOT (BUILD_WITH_CUDA OR BUILD_WITH_HIP))
    set_source_files_properties(
        calculators/GPUSQCalculator.cu
        calculators/GPUDistanceCalculator.cu
        calculators/GPUSteinhardtCalculator.cu
        PROPERTIES LANGUAGE CXX
    )
    target_sources(calculators_obj PRIVATE
        calculators/GPUSQCalculator.cu
        calculators/GPUDistanceCalculator.cu
        calculators/GPUSteinhardtCalculator.cu
    )
endif()
target_sources(calculators_obj PRIVATE
    calculators/gpu/SYCLUtils.cpp
    calculators/gpu/SYCLDistanceCalculator.cpp
    calculators/gpu/SYCLSQCalculator.cpp
)
target_include_directories(calculators_obj PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(calculators_obj PRIVATE TBB::tbb voro++)
target_compile_definitions(calculators_obj PRIVATE _PSTL_PAR_BACKEND_TBB)
target_link_libraries(calculators_obj PRIVATE correlation_simd_options)

# Optional GPU-accelerated calculators (CUDA / HIP)
if(BUILD_WITH_CUDA OR BUILD_WITH_HIP)
    add_library(gpu_calculators_obj OBJECT
        calculators/GPUSQCalculator.cu
        calculators/GPUDistanceCalculator.cu
        calculators/GPUSteinhardtCalculator.cu
    )
    target_include_directories(gpu_calculators_obj PRIVATE
        ${PROJECT_SOURCE_DIR}/include
    )
    target_link_libraries(gpu_calculators_obj PRIVATE
        TBB::tbb
    )
    if(BUILD_WITH_CUDA)
        target_link_libraries(gpu_calculators_obj PRIVATE CUDA::cudart)
        set_target_properties(gpu_calculators_obj PROPERTIES
            CUDA_STANDARD 20
            POSITION_INDEPENDENT_CODE ON
        )
    elseif(BUILD_WITH_HIP)
        set_target_properties(gpu_calculators_obj PROPERTIES
            HIP_STANDARD 20
            POSITION_INDEPENDENT_CODE ON
        )
    endif()
endif()

# Writers as an OBJECT library to ensure static registration triggers
set(WRITERS_SOURCES
    writers/CSVWriter.cpp
    writers/FileWriter.cpp
    writers/WriterFactory.cpp
)
if(BUILD_WITH_ARROW)
    list(APPEND WRITERS_SOURCES writers/ArrowWriter.cpp)
endif()
if(BUILD_WITH_HDF5)
    list(APPEND WRITERS_SOURCES writers/HDF5Writer.cpp)
endif()

add_library(writers_obj OBJECT ${WRITERS_SOURCES})
target_include_directories(writers_obj PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(writers_obj PRIVATE correlation_simd_options)
if(BUILD_WITH_HDF5)
    target_link_libraries(writers_obj PRIVATE HighFive HDF5::HDF5 HDF5::HDF5_HL)
endif()
if(BUILD_WITH_ARROW)
    target_link_libraries(writers_obj PRIVATE arrow_shared parquet_shared)
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
    plotters/pdfgen.c
    plotters/PathFont.cpp
    plotters/RobotoData_ASCII_1.cpp
    plotters/RobotoData_ASCII_2.cpp
    plotters/RobotoData_Extended.cpp
    $<TARGET_OBJECTS:readers_obj>
    $<TARGET_OBJECTS:calculators_obj>
    $<TARGET_OBJECTS:writers_obj>
    $<$<OR:$<BOOL:${BUILD_WITH_CUDA}>,$<BOOL:${BUILD_WITH_HIP}>>:$<TARGET_OBJECTS:gpu_calculators_obj>>
)

target_link_libraries(correlation_lib PRIVATE readers_obj calculators_obj writers_obj)
if(BUILD_WITH_CUDA OR BUILD_WITH_HIP)
    target_link_libraries(correlation_lib PRIVATE gpu_calculators_obj)
endif()

# Include directories
target_include_directories(correlation_lib PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_definitions(correlation_lib PUBLIC _PSTL_PAR_BACKEND_TBB)
if(BUILD_WITH_HDF5)
    target_compile_definitions(correlation_lib PUBLIC CORRELATION_USE_HDF5)
endif()
if(BUILD_WITH_ARROW)
    target_compile_definitions(correlation_lib PUBLIC CORRELATION_USE_ARROW)
endif()

target_link_libraries(correlation_lib PUBLIC correlation_simd_options correlation_fft voro++)

# Pybind11 creates a shared module and requires PIC (handled globally via CMAKE_POSITION_INDEPENDENT_CODE ON)

target_link_libraries(correlation_lib PRIVATE TBB::tbb TBB::tbbmalloc)
if(BUILD_WITH_HDF5)
    target_link_libraries(correlation_lib PUBLIC HighFive HDF5::HDF5 HDF5::HDF5_HL)
endif()
if(BUILD_WITH_ARROW)
    target_link_libraries(correlation_lib PUBLIC arrow_shared parquet_shared)
endif()
if(BUILD_WITH_CUDA)
    target_link_libraries(correlation_lib PRIVATE CUDA::cudart)
    set_target_properties(correlation_lib PROPERTIES CUDA_STANDARD 20)
endif()
if(BUILD_WITH_HIP)
    set_target_properties(correlation_lib PROPERTIES HIP_STANDARD 20)
endif()

#------------------------------------------------------------------
# Executable Target (correlation)
#------------------------------------------------------------------
if(BUILD_GUI)
    set(CORRELATION_SOURCES
        main.cpp
        app/AnalysisRunner.cpp
        app/AppBackend.cpp
        app/AppController.cpp
        app/FileIOHandler.cpp
        app/InputValidator.cpp
        app/PlotController.cpp
        app/PresetController.cpp
        app/PresetManager.cpp
    )

    if (WIN32)
        list(APPEND CORRELATION_SOURCES ../packaging/resource.rc)
        add_executable(correlation WIN32 ${CORRELATION_SOURCES})
    else()
        add_executable(correlation ${CORRELATION_SOURCES})
    endif()

    # Link correlation_lib to correlation
    set(CORRELATION_LIBS correlation_lib Slint::Slint nfd)

    target_link_libraries(correlation PRIVATE correlation_simd_options)

    if (WIN32)
        list(APPEND CORRELATION_LIBS Imm32 Comctl32 dwrite propsys opengl32 dwmapi uxtheme)
    endif()

    target_link_libraries(correlation PRIVATE ${CORRELATION_LIBS})

    # Link platform-specific deps required by a statically-built Slint
    # Must come AFTER Slint::Slint so fontconfig/frameworks follow the
    # static lib that references them (GNU ld resolves left-to-right).
    link_slint_platform_deps(correlation)

    if (WIN32)
        target_include_directories(correlation PRIVATE ${PROJECT_SOURCE_DIR})
    endif()

    # Add the include directory for the generated slint header
    target_include_directories(correlation PRIVATE
      ${CMAKE_CURRENT_BINARY_DIR}
    )

    slint_target_sources(correlation
      ../ui/AppWindow.slint
      LIBRARY_PATHS material=${CMAKE_CURRENT_SOURCE_DIR}/../ui/material/material.slint
    )

    set_property(TARGET correlation PROPERTY SLINT_EMBED_RESOURCES "embed-files")

    # Copy DLLs to the executable directory on Windows
    bundle_windows_dlls(correlation)
    install_windows_dlls(correlation ${CMAKE_INSTALL_BINDIR})

    if (NOT WIN32)
        if(APPLE)
            set_target_properties(correlation PROPERTIES
                INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}"
            )
        else()
            set_target_properties(correlation PROPERTIES
                INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}"
            )
        endif()
    endif()
    if(BUILD_WITH_CUDA)
        set_target_properties(correlation PROPERTIES CUDA_STANDARD 20)
    endif()
    if(BUILD_WITH_HIP)
        set_target_properties(correlation PROPERTIES HIP_STANDARD 20)
    endif()
endif()

#------------------------------------------------------------------
# CLI Parser (OBJECT library — shared with tests)
#------------------------------------------------------------------
add_library(cli_parser_obj OBJECT
    cli/CliParser.cpp
)
target_include_directories(cli_parser_obj PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(cli_parser_obj PUBLIC correlation_simd_options CLI11::CLI11)

#------------------------------------------------------------------
# CLI Executable Target (correlation-cli) — no GUI dependency
#------------------------------------------------------------------
add_executable(correlation-cli
    cli/main_cli.cpp
    app/AppBackend.cpp
)
target_link_libraries(correlation-cli PRIVATE
    correlation_lib readers_obj writers_obj calculators_obj cli_parser_obj
)
target_include_directories(correlation-cli PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)
target_compile_definitions(correlation-cli PRIVATE _USE_MATH_DEFINES)
target_link_libraries(correlation-cli PRIVATE correlation_simd_options)
bundle_windows_dlls(correlation-cli)
install(TARGETS correlation-cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT correlation)
if(BUILD_WITH_CUDA)
    set_target_properties(correlation-cli PROPERTIES CUDA_STANDARD 20)
endif()
if(BUILD_WITH_HIP)
    set_target_properties(correlation-cli PROPERTIES HIP_STANDARD 20)
endif()

#------------------------------------------------------------------
# Python Bindings
#------------------------------------------------------------------
option(BUILD_PYTHON_BINDINGS "Build pybind11 bindings for Correlation" OFF)
if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(bindings)
endif()
