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
    unit/io/FileReaderTests.cpp
    unit/io/WriterTests.cpp
    unit/io/VaspReaderTests.cpp
    unit/io/XdatcarReaderTests.cpp
    unit/io/XYZReaderTests.cpp
    unit/io/QEReaderTests.cpp
    unit/io/CP2KReaderTests.cpp
    unit/io/GromacsReaderTests.cpp
    unit/io/PdbReaderTests.cpp
    unit/io/ArcReaderTests.cpp
    unit/io/CarReaderTests.cpp
    unit/io/CellReaderTests.cpp
    unit/io/CifReaderTests.cpp
    unit/io/CastepMdReaderTests.cpp
    unit/io/LammpsDumpReaderTests.cpp
    unit/io/OnetepDatReaderTests.cpp
    unit/io/OutmolReaderTests.cpp
    unit/analysis/DihedralCalculatorTests.cpp
    unit/analysis/MotifFinderTests.cpp
    unit/analysis/RDFTests.cpp
    unit/analysis/XRDTests.cpp
    unit/analysis/PADTests.cpp
    unit/analysis/VACFTests.cpp
    unit/analysis/vDoSTests.cpp
    unit/analysis/DADTests.cpp
    unit/analysis/RDTests.cpp
    unit/analysis/StructureFactorCalculatorTests.cpp
    unit/analysis/SteinhardtCalculatorTests.cpp
    unit/analysis/CNACalculatorTests.cpp
    unit/analysis/HBondCalculatorTests.cpp
    unit/analysis/SDFCalculatorTests.cpp
    unit/analysis/ClusterCalculatorTests.cpp
    unit/physics/PhysicalDataTests.cpp
    unit/analysis/AngleCalculatorTests.cpp
    unit/analysis/CNCalculatorTests.cpp
    unit/analysis/DistanceCalculatorTests.cpp
    unit/analysis/MSDCalculatorTests.cpp
    unit/analysis/CalculatorFactoryTests.cpp
    unit/io/ReaderFactoryTests.cpp
    unit/io/WriterFactoryTests.cpp
    unit/plotters/SvgPlotterTests.cpp
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
        ../src/app/AppBackend.cpp
        ../src/app/AppController.cpp
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
