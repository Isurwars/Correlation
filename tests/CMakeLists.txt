# Correlation - Liquid and Amorphous Solid Analysis Tool
# Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
# SPDX-License-Identifier: MIT
# Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#-----------------------------------------------------------
# Testing setup (using GoogleTest)
#-----------------------------------------------------------

include(GoogleTest)

#-----------------------------------------------------------
# 1. Headless tests (no Slint dependency — runs on ALL platforms)
#-----------------------------------------------------------
add_executable(correlation_tests
    _01_Atom_Tests.cpp
    _02_Cell_Tests.cpp
    _03_Trajectory_Tests.cpp
    _04_FileReader_Tests.cpp
    _05_Writer_Tests.cpp
    _06_DihedralCalculator_Tests.cpp
    _07_MotifFinder_Tests.cpp
    _08_StructureAnalyzer_Tests.cpp
    _09_TrajectoryAnalyzer_Tests.cpp
    _10_DynamicsAnalyzer_Tests.cpp
    _11_RDF_Tests.cpp
    _12_XRD_Tests.cpp
    _13_PAD_Tests.cpp
    _14_VACF_Tests.cpp
    _15_vDoS_Tests.cpp
    _16_AppBackend_Tests.cpp
    _18_DAD_Tests.cpp
    _19_RD_Tests.cpp
    _20_SIMDUtils_Tests.cpp
    _21_StructureFactorCalculator_Tests.cpp
    _22_SteinhardtCalculator_Tests.cpp
    _23_VaspReader_Tests.cpp
    _24_XdatcarReader_Tests.cpp
    ../src/app/AppBackend.cpp
)

target_link_libraries(correlation_tests PRIVATE
    correlation_lib
    readers_obj
    writers_obj
    calculators_obj
    GTest::gtest
    GTest::gtest_main
)
target_compile_definitions(correlation_tests PRIVATE _USE_MATH_DEFINES)

# Call the helper function to bundle DLLs on Windows
bundle_windows_dlls(correlation_tests)

gtest_discover_tests(correlation_tests DISCOVERY_TIMEOUT 60)

#-----------------------------------------------------------
# 2. GUI tests (requires Slint — runs on platforms with a display)
#-----------------------------------------------------------
add_executable(correlation_gui_tests
    _17_AppController_Tests.cpp
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
