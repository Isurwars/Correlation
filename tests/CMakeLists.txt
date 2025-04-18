add_executable(atom_test
  Atom_test.cpp
  ../src/Atom.cpp
)

target_link_libraries(atom_test PRIVATE correlation_lib GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(atom_test)
