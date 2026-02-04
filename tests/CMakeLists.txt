#-----------------------------------------------------------
# Testing setup (using GoogleTest)
#-----------------------------------------------------------

include(GoogleTest)

# Explicitly list test sources.
add_executable(correlation_tests
    Atom_test.cpp
    Cell_test.cpp
    Distribution_Functions_tests.cpp
    FileIO_test.cpp
    FileWriter_test.cpp
    DynamicsAnalyzer_test.cpp
    StructureAnalyzer_test.cpp
    AngleTests.cpp
    TrajectoryAnalyzer_test.cpp
    TrajectoryRevision_test.cpp
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
