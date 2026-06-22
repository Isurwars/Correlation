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
    unit/readers/FileReaderTests.cpp
    unit/writers/WriterTests.cpp
    unit/readers/VaspReaderTests.cpp
    unit/readers/XdatcarReaderTests.cpp
    unit/readers/XYZReaderTests.cpp
    unit/readers/QEReaderTests.cpp
    unit/readers/CP2KReaderTests.cpp
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
    unit/calculators/RDFTests.cpp
    unit/calculators/XRDTests.cpp
    unit/calculators/PADTests.cpp
    unit/calculators/VACFTests.cpp
    unit/calculators/vDoSTests.cpp
    unit/calculators/DADTests.cpp
    unit/calculators/RDTests.cpp
    unit/calculators/StructureFactorCalculatorTests.cpp
    unit/calculators/SteinhardtCalculatorTests.cpp
    unit/calculators/CNACalculatorTests.cpp
    unit/calculators/HBondCalculatorTests.cpp
    unit/calculators/SDFCalculatorTests.cpp
    unit/calculators/ClusterCalculatorTests.cpp
    unit/physics/PhysicalDataTests.cpp
    unit/calculators/AngleCalculatorTests.cpp
    unit/calculators/CNCalculatorTests.cpp
    unit/calculators/DistanceCalculatorTests.cpp
    unit/calculators/MSDCalculatorTests.cpp
    unit/calculators/CalculatorFactoryTests.cpp
    unit/readers/ReaderFactoryTests.cpp
    unit/writers/WriterFactoryTests.cpp
    unit/plotters/SvgPlotterTests.cpp
    unit/plotters/PdfPlotterTests.cpp
    unit/utils/SIMDUtilsTests.cpp
    unit/math/LinearAlgebraTests.cpp
    unit/math/SpecialFunctionsTests.cpp
    unit/math/SmoothingTests.cpp
    unit/math/FFTUtilsTests.cpp
    unit/app/AppBackendTests.cpp
    unit/app/PresetManagerTests.cpp
    unit/cli/CliParserTests.cpp
    ../src/app/AppBackend.cpp
    ../src/app/PresetManager.cpp
)

target_link_libraries(correlation_unit_tests PRIVATE
    correlation_lib
    readers_obj
    writers_obj
    calculators_obj
    cli_parser_obj
    GTest::gtest
    GTest::gtest_main
)
# (Definitions and discovery are configured at the end of the file in a loop)

#-----------------------------------------------------------
# 2. Functional Tests (Real-world scenarios)
#-----------------------------------------------------------
add_executable(correlation_functional_tests
    functional/core/CellTests.cpp
    functional/analysis/StructureAnalyzerTests.cpp
    functional/analysis/TrajectoryAnalyzerTests.cpp
    functional/analysis/DynamicsAnalyzerTests.cpp
)

target_link_libraries(correlation_functional_tests PRIVATE
    correlation_lib
    readers_obj
    writers_obj
    calculators_obj
    GTest::gtest
    GTest::gtest_main
)

#-----------------------------------------------------------
# 2b. CLI End-to-End Tests (subprocess-based)
#-----------------------------------------------------------
add_executable(correlation_cli_e2e_tests
    functional/cli/CliEndToEndTests.cpp
)

target_link_libraries(correlation_cli_e2e_tests PRIVATE
    GTest::gtest
    GTest::gtest_main
)

# Tell the test binary where the CLI executable lives
target_compile_definitions(correlation_cli_e2e_tests PRIVATE
    CORRELATION_CLI_PATH="$<TARGET_FILE:correlation-cli>"
)

# The E2E tests depend on the CLI binary being built
add_dependencies(correlation_cli_e2e_tests correlation-cli)
# (Definitions and discovery are configured at the end of the file in a loop)

#-----------------------------------------------------------
# 3. GUI tests (requires Slint — runs on platforms with a display)
#-----------------------------------------------------------
if(BUILD_GUI)
    add_executable(correlation_gui_tests
        functional/app/AppControllerTests.cpp
        functional/app/SlintSvgIntegrationTests.cpp
        ../src/app/AnalysisRunner.cpp
        ../src/app/AppBackend.cpp
        ../src/app/AppController.cpp
        ../src/app/FileIOHandler.cpp
        ../src/app/InputValidator.cpp
        ../src/app/PlotController.cpp
        ../src/app/PresetController.cpp
        ../src/app/PresetManager.cpp
    )

    target_link_libraries(correlation_gui_tests PRIVATE
        correlation_lib
        readers_obj
        writers_obj
        calculators_obj
        GTest::gtest
        GTest::gtest_main
        Slint::Slint
        nfd
    )

    # Link platform-specific deps required by a statically-built Slint
    link_slint_platform_deps(correlation_gui_tests)

    # Add slint target sources for AppWindow
    slint_target_sources(correlation_gui_tests
        ../ui/AppWindow.slint
        LIBRARY_PATHS material=${CMAKE_CURRENT_SOURCE_DIR}/../ui/material/material.slint
    )
    target_include_directories(correlation_gui_tests PRIVATE
      ${CMAKE_CURRENT_BINARY_DIR}
    )
endif()

# Configure common test target settings and discovery
set(TEST_TARGETS correlation_unit_tests correlation_functional_tests correlation_cli_e2e_tests)
if(BUILD_GUI)
    list(APPEND TEST_TARGETS correlation_gui_tests)
endif()

foreach(tgt IN LISTS TEST_TARGETS)
    target_compile_definitions(${tgt} PRIVATE _USE_MATH_DEFINES)
    bundle_windows_dlls(${tgt})
    gtest_discover_tests(${tgt} DISCOVERY_TIMEOUT 60)
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
