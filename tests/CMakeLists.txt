# Correlation - Liquid and Amorphous Solid Analysis Tool
# Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
# SPDX-License-Identifier: MIT
# Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#-----------------------------------------------------------
# Testing setup (using GoogleTest)
#-----------------------------------------------------------

include(GoogleTest)

# Explicitly list test sources.
add_executable(correlation_tests
    _01_Atom_Tests.cpp
    _02_Cell_Tests.cpp
    _03_Trajectory_Tests.cpp
    _04_FileReader_Tests.cpp
    _05_FileWriter_Tests.cpp
    _06_DihedralCalculator_Tests.cpp
    _07_MotifFinder_Tests.cpp
    _08_StructureAnalyzer_Tests.cpp
    _09_TrajectoryAnalyzer_Tests.cpp
    _10_DynamicsAnalyzer_Tests.cpp
    _11_RDF_Tests.cpp
    _12_SQ_Tests.cpp
    _13_XRD_Tests.cpp
    _14_PAD_Tests.cpp
    _15_VACF_Tests.cpp
    _16_vDoS_Tests.cpp
    _17_AppBackend_Tests.cpp
    _18_AppController_Tests.cpp
    _19_DAD_Tests.cpp
    _20_MD_Tests.cpp
    ../src/AppBackend.cpp
    ../src/AppController.cpp
)

# The test executable also links to the library and inherits its properties.
target_link_libraries(correlation_tests PRIVATE
    correlation_lib
    GTest::gtest
    GTest::gtest_main
    Slint::Slint
)
target_compile_definitions(correlation_tests PRIVATE _USE_MATH_DEFINES)

# Add slint target sources for AppWindow
slint_target_sources(correlation_tests
    ../ui/AppWindow.slint
    LIBRARY_PATHS material=${CMAKE_CURRENT_SOURCE_DIR}/../material-1.0/material.slint
)
target_include_directories(correlation_tests PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}
)

# Call the helper function to bundle DLLs on Windows
bundle_windows_dlls(correlation_tests)

gtest_discover_tests(correlation_tests DISCOVERY_TIMEOUT 60)
