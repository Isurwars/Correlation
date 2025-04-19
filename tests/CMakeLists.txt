# Find testing framework
find_package(GTest REQUIRED)

# Find all test source files
file(GLOB TEST_SOURCES "*.cpp")

# Create and register tests
foreach(test_source ${TEST_SOURCES})
    get_filename_component(test_name ${test_source} NAME_WE)

    add_executable(${test_name} ${test_source})

    target_link_libraries(${test_name} PRIVATE
	correlation_lib
	GTest::GTest
	GTest::Main
	TBB::tbb
	TBB::tbbmalloc
    )

    target_include_directories(${test_name} PRIVATE
	${PROJECT_SOURCE_DIR}/include
    )

    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
