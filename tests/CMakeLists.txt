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
    Test01_Atom.cpp
    Test02_Cell.cpp
    Test03_Trajectory.cpp
    Test04_StructureAnalyzer.cpp
    Test05_TrajectoryAnalyzer.cpp
    Test06_RDF.cpp
    Test07_SQ.cpp
    Test08_XRD.cpp
    Test09_PAD.cpp
    Test10_VACF.cpp
    FileWriter_test.cpp
    FileIO_test.cpp
    DynamicsAnalyzer_test.cpp
)

# The test executable also links to the library and inherits its properties.
target_link_libraries(correlation_tests PRIVATE
    correlation_lib
    GTest::gtest
    GTest::gtest_main
)
target_compile_definitions(correlation_tests PRIVATE _USE_MATH_DEFINES)

# Call the helper function to bundle DLLs on Windows
bundle_windows_dlls(correlation_tests)

gtest_discover_tests(correlation_tests DISCOVERY_TIMEOUT 60)
