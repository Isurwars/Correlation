# Correlation - Liquid and Amorphous Solid Analysis Tool
# Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
# SPDX-License-Identifier: MIT
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
    unit/io/FileReaderTests.cpp
    unit/io/WriterTests.cpp
    unit/io/VaspReaderTests.cpp
    unit/io/XdatcarReaderTests.cpp
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
    unit/utils/SIMDUtilsTests.cpp
    unit/app/AppBackendTests.cpp
    ../src/app/AppBackend.cpp
)

target_link_libraries(correlation_unit_tests PRIVATE
    correlation_lib
    readers_obj
    writers_obj
    calculators_obj
    GTest::gtest
    GTest::gtest_main
)
target_compile_definitions(correlation_unit_tests PRIVATE _USE_MATH_DEFINES)
bundle_windows_dlls(correlation_unit_tests)
gtest_discover_tests(correlation_unit_tests DISCOVERY_TIMEOUT 60)

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
target_compile_definitions(correlation_functional_tests PRIVATE _USE_MATH_DEFINES)
bundle_windows_dlls(correlation_functional_tests)
gtest_discover_tests(correlation_functional_tests DISCOVERY_TIMEOUT 60)

#-----------------------------------------------------------
# 3. GUI tests (requires Slint — runs on platforms with a display)
#-----------------------------------------------------------
add_executable(correlation_gui_tests
    functional/app/AppControllerTests.cpp
    ../src/app/AppBackend.cpp
    ../src/app/AppController.cpp
)

target_link_libraries(correlation_gui_tests PRIVATE
    correlation_lib
    readers_obj
    writers_obj
    calculators_obj
    GTest::gtest
    GTest::gtest_main
    Slint::Slint
)
target_compile_definitions(correlation_gui_tests PRIVATE _USE_MATH_DEFINES)

# Add slint target sources for AppWindow
slint_target_sources(correlation_gui_tests
    ../ui/AppWindow.slint
    LIBRARY_PATHS material=${CMAKE_CURRENT_SOURCE_DIR}/../ui/material/material.slint
)
target_include_directories(correlation_gui_tests PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}
)

# Call the helper function to bundle DLLs on Windows
bundle_windows_dlls(correlation_gui_tests)

gtest_discover_tests(correlation_gui_tests DISCOVERY_TIMEOUT 60)

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
