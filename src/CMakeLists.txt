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
    readers/OrcaReader.cpp
    readers/GpawReader.cpp
    readers/AbinitReader.cpp
    readers/DftbReader.cpp
    readers/MaceReader.cpp
    readers/ChgnetReader.cpp
    readers/GapReader.cpp
    readers/NequipReader.cpp
)
target_include_directories(readers_obj PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(readers_obj PRIVATE correlation_simd_options)

# Calculators as an OBJECT library to ensure static registration triggers
add_library(calculators_obj OBJECT
    # Base and factory
    calculators/CalculatorFactory.cpp

    # Machine learning calculators
    calculators/ml/MLIPCalculator.cpp

    # Spatial calculators
    calculators/spatial/AngleCalculator.cpp
    calculators/spatial/ClusterCalculator.cpp
    calculators/spatial/CNACalculator.cpp
    calculators/spatial/CNCalculator.cpp
    calculators/spatial/DADCalculator.cpp
    calculators/spatial/DihedralCalculator.cpp
    calculators/spatial/DistanceCalculator.cpp
    calculators/spatial/HBondCalculator.cpp
    calculators/spatial/MotifFinder.cpp
    calculators/spatial/PADCalculator.cpp
    calculators/spatial/RDCalculator.cpp
    calculators/spatial/RDFCalculator.cpp
    calculators/spatial/SDFCalculator.cpp
    calculators/spatial/VoronoiCalculator.cpp

    # Order and symmetry calculators
    calculators/order/ChiralityCalculator.cpp
    calculators/order/HyperuniformityCalculator.cpp
    calculators/order/LocalEntropyCalculator.cpp
    calculators/order/SteinhardtCalculator.cpp

    # Temporal / dynamic calculators
    calculators/temporal/MSDCalculator.cpp
    calculators/temporal/VACFCalculator.cpp
    calculators/temporal/VDOSCalculator.cpp

    # Scattering and reciprocal space calculators
    calculators/scattering/StructureFactorCalculator.cpp
    calculators/scattering/TDOSCalculator.cpp
    calculators/scattering/XRDCalculator.cpp
)

if(NOT (BUILD_WITH_CUDA OR BUILD_WITH_HIP))
    set_source_files_properties(
        calculators/gpu/GPUSQCalculator.cu
        calculators/gpu/GPUDistanceCalculator.cu
        calculators/gpu/GPUSteinhardtCalculator.cu
        calculators/gpu/GPUXRDCalculator.cu
        PROPERTIES LANGUAGE CXX
    )
    target_sources(calculators_obj PRIVATE
        calculators/gpu/GPUSQCalculator.cu
        calculators/gpu/GPUDistanceCalculator.cu
        calculators/gpu/GPUSteinhardtCalculator.cu
        calculators/gpu/GPUXRDCalculator.cu
    )
endif()
target_sources(calculators_obj PRIVATE
    calculators/gpu/SYCLUtils.cpp
    calculators/gpu/SYCLDistanceCalculator.cpp
    calculators/gpu/SYCLAngleCalculator.cpp
    calculators/gpu/SYCLSQCalculator.cpp
    calculators/gpu/SYCLXRDCalculator.cpp
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
        calculators/gpu/GPUSQCalculator.cu
        calculators/gpu/GPUDistanceCalculator.cu
        calculators/gpu/GPUSteinhardtCalculator.cu
        calculators/gpu/GPUXRDCalculator.cu
    )
    target_include_directories(gpu_calculators_obj PRIVATE
        ${PROJECT_SOURCE_DIR}/include
    )
    target_link_libraries(gpu_calculators_obj PRIVATE
        TBB::tbb
        correlation_simd_options
    )
    if(BUILD_WITH_CUDA)
        target_link_libraries(gpu_calculators_obj PRIVATE CUDA::cudart)
        target_compile_options(gpu_calculators_obj PRIVATE $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CUDA_COMPILER_ID:NVIDIA>>:--expt-relaxed-constexpr>)
        set_target_properties(gpu_calculators_obj PROPERTIES
            CUDA_STANDARD 23
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
    writers/ExportBundleWriter.cpp
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
target_link_libraries(writers_obj PRIVATE correlation_simd_options miniz)
if(BUILD_WITH_HDF5)
    target_link_libraries(writers_obj PRIVATE HighFive HDF5::HDF5 HDF5::HDF5_HL)
endif()
if(BUILD_WITH_ARROW)
    target_link_libraries(writers_obj PRIVATE arrow_shared parquet_shared)
endif()

# MLIP as an OBJECT library
set(MLIP_SOURCES
    mlip/PeriodicGraphBuilder.cpp
    mlip/TorchGNNModel.cpp
    mlip/GraphDescriptors.cpp
)

add_library(mlip_obj OBJECT ${MLIP_SOURCES})
target_include_directories(mlip_obj PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(mlip_obj PRIVATE TBB::tbb correlation_simd_options)
if(CORRELATION_ENABLE_LIBTORCH AND Torch_FOUND)
    target_link_libraries(mlip_obj PRIVATE "${TORCH_LIBRARIES}")
    target_include_directories(mlip_obj PRIVATE ${TORCH_INCLUDE_DIRS})
endif()

# Vendor: pdfgen C library
add_library(correlation_pdfgen STATIC
    plotters/vendor/pdfgen.c
)
set_target_properties(correlation_pdfgen PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    C_STANDARD 11
)
target_include_directories(correlation_pdfgen PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Plotters as an OBJECT library
set(PLOTTERS_SOURCES
    plotters/PathFont.cpp
    plotters/PngPlotter.cpp
    plotters/fonts/RobotoData_ASCII_1.cpp
    plotters/fonts/RobotoData_ASCII_2.cpp
    plotters/fonts/RobotoData_Extended.cpp
)

add_library(plotters_obj OBJECT ${PLOTTERS_SOURCES})
target_include_directories(plotters_obj PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(plotters_obj PRIVATE correlation_simd_options correlation_pdfgen)

# Add sources to library
add_library(correlation_lib STATIC
    core/Cell.cpp
    analysis/DistributionFunctions.cpp
    analysis/DynamicsAnalyzer.cpp
    analysis/CorrelationEngine.cpp
    analysis/StructuralElectronicCorrelation.cpp
    core/NeighborGraph.cpp
    analysis/StructureAnalyzer.cpp
    core/Trajectory.cpp
    analysis/TrajectoryAnalyzer.cpp
    utils/ZipArchiveWriter.cpp
    $<TARGET_OBJECTS:readers_obj>
    $<TARGET_OBJECTS:calculators_obj>
    $<TARGET_OBJECTS:writers_obj>
    $<TARGET_OBJECTS:mlip_obj>
    $<TARGET_OBJECTS:plotters_obj>
    $<$<OR:$<BOOL:${BUILD_WITH_CUDA}>,$<BOOL:${BUILD_WITH_HIP}>>:$<TARGET_OBJECTS:gpu_calculators_obj>>
)

target_link_libraries(correlation_lib PRIVATE readers_obj calculators_obj writers_obj mlip_obj plotters_obj miniz)
target_link_libraries(correlation_lib PUBLIC correlation_pdfgen)
if(CORRELATION_ENABLE_LIBTORCH AND Torch_FOUND)
    target_link_libraries(correlation_lib PRIVATE "${TORCH_LIBRARIES}")
endif()
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
    target_link_libraries(correlation_lib PUBLIC CUDA::cudart)
    set_target_properties(correlation_lib PROPERTIES CUDA_STANDARD 20)
endif()
if(BUILD_WITH_HIP)
    set_target_properties(correlation_lib PROPERTIES HIP_STANDARD 20)
endif()

#------------------------------------------------------------------
# Application Services (Headless & GUI)
#------------------------------------------------------------------
# Headless Application Services as an OBJECT library (shared across GUI, CLI, and tests)
add_library(app_services_obj OBJECT
    app/services/AnalysisDispatcher.cpp
    app/services/BondCutoffService.cpp
    app/services/PhysicsService.cpp
    app/services/TrajectoryLoader.cpp
    app/services/PresetManager.cpp
    app/services/UpdateChecker.cpp
    app/services/ValidationRuleService.cpp
    app/core/SettingsManager.cpp
    app/formatters/BondCutoffMapper.cpp
    app/formatters/PlotExportService.cpp
)
target_include_directories(app_services_obj PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(app_services_obj PUBLIC
    correlation_lib
    correlation_simd_options
)

#------------------------------------------------------------------
# Executable Target (correlation)
#------------------------------------------------------------------
if(BUILD_GUI)
    # GUI ViewModels & Controllers compiled once
    add_library(app_gui_obj OBJECT
        app/core/AppController.cpp
        app/viewmodel/BondCutoffController.cpp
        app/viewmodel/PlotController.cpp
        app/viewmodel/PresetController.cpp
        app/services/AnalysisRunner.cpp
        app/services/FileIOHandler.cpp
        app/services/InputValidator.cpp
        app/services/OptionsResetService.cpp
        app/services/OptionsSyncService.cpp
        app/formatters/PlotSeriesManager.cpp
        app/formatters/PlotTableFormatter.cpp
    )
    target_include_directories(app_gui_obj PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_BINARY_DIR}
    )
    target_link_libraries(app_gui_obj PUBLIC
        correlation_lib
        app_services_obj
        Slint::Slint
        nfd
        correlation_simd_options
    )
    set_property(TARGET app_gui_obj PROPERTY SLINT_EMBED_RESOURCES "embed-files")
    slint_target_sources(app_gui_obj
        ../ui/AppWindow.slint
        LIBRARY_PATHS material=${CMAKE_CURRENT_SOURCE_DIR}/../ui/material/material.slint
    )

    set(CORRELATION_SOURCES
        main.cpp
    )

    if (WIN32)
        list(APPEND CORRELATION_SOURCES "${PROJECT_BINARY_DIR}/packaging/resource.rc")
        add_executable(correlation WIN32 ${CORRELATION_SOURCES})
    else()
        add_executable(correlation ${CORRELATION_SOURCES})
    endif()

    set(CORRELATION_LIBS
        "$<LINK_LIBRARY:WHOLE_ARCHIVE,correlation_lib>"
        app_services_obj
        app_gui_obj
        Slint::Slint
        nfd
    )

    target_link_libraries(correlation PRIVATE correlation_simd_options)

    if (WIN32)
        list(APPEND CORRELATION_LIBS Imm32 Comctl32 dwrite propsys opengl32 dwmapi uxtheme wininet Shell32)
    endif()

    target_link_libraries(correlation PRIVATE ${CORRELATION_LIBS})

    # Link platform-specific deps required by a statically-built Slint
    link_slint_platform_deps(correlation)

    if (WIN32)
        target_include_directories(correlation PRIVATE ${PROJECT_SOURCE_DIR})
    endif()

    # Copy DLLs to the executable directory on Windows
    bundle_windows_dlls(correlation)
    install_windows_dlls(correlation ${CMAKE_INSTALL_BINDIR})

    if(BUILD_WITH_CUDA)
        set_target_properties(correlation PROPERTIES CUDA_STANDARD 20)
    endif()
    if(BUILD_WITH_HIP)
        set_target_properties(correlation PROPERTIES HIP_STANDARD 20)
    endif()
endif()

#------------------------------------------------------------------
# CLI Parser & Executable Target (disabled for WebAssembly builds)
#------------------------------------------------------------------
if(NOT (EMSCRIPTEN OR BUILD_WASM))
    # CLI Parser (OBJECT library — shared with tests)
    add_library(cli_parser_obj OBJECT
        cli/CliParser.cpp
    )
    target_include_directories(cli_parser_obj PUBLIC
        ${PROJECT_SOURCE_DIR}/include
    )
    target_link_libraries(cli_parser_obj PUBLIC correlation_simd_options CLI11::CLI11)

    # CLI Executable Target (correlation-cli) — no GUI dependency
    add_executable(correlation-cli
        cli/main_cli.cpp
    )
    target_link_libraries(correlation-cli PRIVATE
        "$<LINK_LIBRARY:WHOLE_ARCHIVE,correlation_lib>"
        cli_parser_obj
        app_services_obj
        correlation_simd_options
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
endif()

#------------------------------------------------------------------
# Language & Platform Bindings (Python / WebAssembly)
#------------------------------------------------------------------
if(BUILD_PYTHON_BINDINGS OR EMSCRIPTEN OR BUILD_WASM)
    add_subdirectory(bindings)
endif()
