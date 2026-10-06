# Correlation - Liquid and Amorphous Solid Analysis Tool
# Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
# SPDX-License-Identifier: AGPL-3.0-only
# Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#-----------------------------------------------------------
# Testing setup (using GoogleTest)
#-----------------------------------------------------------

include(GoogleTest)

#-----------------------------------------------------------
# 1. Unit Tests (Unitary logic and limit cases)
#-----------------------------------------------------------
add_executable(correlation_unit_tests
    unit/core/AtomTests.cpp
    unit/core/CellTests.cpp
    unit/core/TrajectoryTests.cpp
    unit/core/NeighborGraphTests.cpp
    unit/core/MappedFileTests.cpp
    unit/core/GPUPortabilityTests.cpp
    unit/readers/FileReaderTests.cpp
    unit/writers/WriterTests.cpp
    unit/readers/VaspReaderTests.cpp
    unit/readers/XdatcarReaderTests.cpp
    unit/readers/XYZReaderTests.cpp
    unit/readers/MaceReaderTests.cpp
    unit/readers/ChgnetReaderTests.cpp
    unit/readers/GapReaderTests.cpp
    unit/readers/NequipReaderTests.cpp
    unit/readers/QEReaderTests.cpp
    unit/readers/CP2KReaderTests.cpp
    unit/readers/OrcaReaderTests.cpp
    unit/readers/GpawReaderTests.cpp
    unit/readers/AbinitReaderTests.cpp
    unit/readers/DftbReaderTests.cpp
    unit/readers/GromacsReaderTests.cpp
    unit/readers/PdbReaderTests.cpp
    unit/readers/ArcReaderTests.cpp
    unit/readers/CarReaderTests.cpp
    unit/readers/CellReaderTests.cpp
    unit/readers/CifReaderTests.cpp
    unit/readers/CastepMdReaderTests.cpp
    unit/readers/LammpsDumpReaderTests.cpp
    unit/readers/OnetepDatReaderTests.cpp
    unit/readers/OutmolReaderTests.cpp
    unit/calculators/DihedralCalculatorTests.cpp
    unit/calculators/MotifFinderTests.cpp
    unit/calculators/RDFCalculatorTests.cpp
    unit/calculators/XRDCalculatorTests.cpp
    unit/calculators/PADCalculatorTests.cpp
    unit/calculators/VACFCalculatorTests.cpp
    unit/calculators/VDOSCalculatorTests.cpp
    unit/calculators/DADCalculatorTests.cpp
    unit/calculators/RDCalculatorTests.cpp
    unit/calculators/RingCellTests.cpp
    unit/calculators/StructureFactorCalculatorTests.cpp
    unit/calculators/SteinhardtCalculatorTests.cpp
    unit/calculators/CNACalculatorTests.cpp
    unit/calculators/HBondCalculatorTests.cpp
    unit/calculators/SDFCalculatorTests.cpp
    unit/calculators/ClusterCalculatorTests.cpp
    unit/calculators/LocalEntropyCalculatorTests.cpp
    unit/calculators/VoronoiCalculatorTests.cpp
    unit/calculators/ChiralityCalculatorTests.cpp
    unit/calculators/HyperuniformityCalculatorTests.cpp
    unit/calculators/MLIPCalculatorTests.cpp
    unit/calculators/TDOSCalculatorTests.cpp
    unit/mlip/PeriodicGraphBuilderTests.cpp
    unit/mlip/TorchGNNModelTests.cpp
    unit/mlip/GraphDescriptorsTests.cpp
    unit/calculators/GPUSteinhardtCalculatorTests.cpp
    unit/calculators/GPUDistanceCalculatorTests.cpp
    unit/calculators/GPUSQCalculatorTests.cpp
    unit/calculators/GPUXRDCalculatorTests.cpp
    unit/calculators/SYCLCalculatorTests.cpp
    unit/physics/PhysicalDataTests.cpp
    unit/calculators/AngleCalculatorTests.cpp
    unit/calculators/CNCalculatorTests.cpp
    unit/calculators/DistanceCalculatorTests.cpp
    unit/calculators/MSDCalculatorTests.cpp
    unit/calculators/CalculatorFactoryTests.cpp
    unit/readers/ReaderFactoryTests.cpp
    unit/writers/WriterFactoryTests.cpp
    unit/writers/ArrowWriterTests.cpp
    unit/plotters/SvgPlotterTests.cpp
    unit/plotters/PdfPlotterTests.cpp
    unit/plotters/PathFontTests.cpp
    unit/utils/SIMDUtilsTests.cpp
    unit/utils/ZipArchiveWriterTests.cpp
    unit/math/Vector3Tests.cpp
    unit/math/Matrix3Tests.cpp
    unit/math/PrecisionTests.cpp
    unit/math/SpecialFunctionsTests.cpp
    unit/math/SmoothingTests.cpp
    unit/math/FFTUtilsTests.cpp
    unit/math/ConstantsTests.cpp
    unit/analysis/CorrelationEngineTests.cpp
    unit/analysis/StructuralElectronicCorrelationTests.cpp
    unit/app/AnalysisDispatcherTests.cpp
    unit/app/BondCutoffMapperTests.cpp
    unit/app/BondCutoffServiceTests.cpp
    unit/app/PhysicsServiceTests.cpp
    unit/app/PlotExportServiceTests.cpp
    unit/app/PresetManagerTests.cpp
    unit/app/SettingsManagerTests.cpp
    unit/app/TrajectoryLoaderTests.cpp
    unit/app/UpdateCheckerTests.cpp
    unit/app/ValidationRuleServiceTests.cpp
    unit/cli/CliParserTests.cpp
)

target_link_libraries(correlation_unit_tests PRIVATE
    "$<LINK_LIBRARY:WHOLE_ARCHIVE,correlation_lib>"
    cli_parser_obj
    app_services_obj
    miniz
    GTest::gtest
    GTest::gmock
    GTest::gtest_main
)
if(WIN32)
    target_link_libraries(correlation_unit_tests PRIVATE wininet Shell32)
endif()
add_dependencies(correlation_unit_tests correlation_lib cli_parser_obj app_services_obj)
# (Definitions and discovery are configured at the end of the file in a loop)

#-----------------------------------------------------------
# 2. Functional Tests (Real-world scenarios)
#-----------------------------------------------------------
add_executable(correlation_functional_tests
    functional/core/CellTests.cpp
    functional/core/AtomTests.cpp
    functional/core/NeighborGraphTests.cpp
    functional/core/TrajectoryTests.cpp
    functional/core/MappedFileTests.cpp
    functional/analysis/StructureAnalyzerTests.cpp
    functional/analysis/TrajectoryAnalyzerTests.cpp
    functional/analysis/DynamicsAnalyzerTests.cpp
    functional/analysis/NumericalParityTests.cpp
    functional/analysis/ThreadDeterminismTests.cpp
    functional/app/ServiceIntegrationTests.cpp
)

target_link_libraries(correlation_functional_tests PRIVATE
    "$<LINK_LIBRARY:WHOLE_ARCHIVE,correlation_lib>"
    app_services_obj
    TBB::tbb
    GTest::gtest
    GTest::gmock
    GTest::gtest_main
)
add_dependencies(correlation_functional_tests correlation_lib app_services_obj)


#-----------------------------------------------------------
# 3. GUI tests (requires Slint — runs on platforms with a display)
#-----------------------------------------------------------
if(BUILD_GUI)
    add_executable(correlation_gui_tests
        functional/app/AppControllerTests.cpp
        functional/app/SlintSvgIntegrationTests.cpp
        unit/app/BondCutoffControllerTests.cpp
        unit/app/OptionsResetServiceTests.cpp
        unit/app/OptionsSyncServiceTests.cpp
        unit/app/PlotSeriesManagerTests.cpp
        unit/app/PlotTableFormatterTests.cpp
    )

    target_link_libraries(correlation_gui_tests PRIVATE
        "$<LINK_LIBRARY:WHOLE_ARCHIVE,correlation_lib>"
        app_services_obj
        app_gui_obj
        GTest::gtest
        GTest::gmock
        GTest::gtest_main
        Slint::Slint
        nfd
    )
    add_dependencies(correlation_gui_tests correlation_lib app_services_obj app_gui_obj)

    # Link platform-specific deps required by a statically-built Slint
    link_slint_platform_deps(correlation_gui_tests)
endif()

# Configure common test target settings and discovery
set(TEST_TARGETS correlation_unit_tests correlation_functional_tests)
if(BUILD_GUI)
    list(APPEND TEST_TARGETS correlation_gui_tests)
endif()

foreach(tgt IN LISTS TEST_TARGETS)
    target_include_directories(${tgt} PRIVATE ${CMAKE_SOURCE_DIR}/include)
    target_compile_definitions(${tgt} PRIVATE _USE_MATH_DEFINES)
    if(MSVC)
        target_compile_options(${tgt} PRIVATE /FI "${CMAKE_CURRENT_SOURCE_DIR}/PrecisionTestHelper.hpp")
    else()
        target_compile_options(${tgt} PRIVATE -include "${CMAKE_CURRENT_SOURCE_DIR}/PrecisionTestHelper.hpp")
    endif()
    bundle_windows_dlls(${tgt})
    if(BUILD_WITH_CUDA)
        set_target_properties(${tgt} PROPERTIES CUDA_STANDARD 20)
    endif()
    if(BUILD_WITH_HIP)
        set_target_properties(${tgt} PROPERTIES HIP_STANDARD 20)
    endif()
    gtest_discover_tests(${tgt} DISCOVERY_MODE PRE_TEST DISCOVERY_TIMEOUT 120 PROPERTIES TIMEOUT 180)
endforeach()

#-----------------------------------------------------------
# Python binding smoke test (only when bindings are built)
#-----------------------------------------------------------
if(BUILD_PYTHON_BINDINGS)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)

    # The shared module is placed in the bindings build directory.
    set(CORRELATION_PY_DIR "${CMAKE_BINARY_DIR}/src/bindings")

    add_test(
        NAME PythonBindings
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_bindings.py"
    )

    # Make the .so importable without installing.
    set_tests_properties(PythonBindings PROPERTIES
        ENVIRONMENT  "PYTHONPATH=${CORRELATION_PY_DIR}"
        DEPENDS      _correlation
        LABELS       "python;bindings"
    )
endif()
