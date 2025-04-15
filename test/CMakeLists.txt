cmake
cmake_minimum_required(VERSION 3.10)

project(Tests)

find_package(Math REQUIRED)

add_executable(test_1)
target_link_libraries(test_1 PRIVATE m)
add_test(NAME test_1 COMMAND test_1)

add_executable(test_2)
target_link_libraries(test_2 PRIVATE m)
add_test(NAME test_2 COMMAND test_2)

add_executable(test_3)
target_link_libraries(test_3 PRIVATE m)
add_test(NAME test_3 COMMAND test_3)

add_executable(test_4)
target_link_libraries(test_4 PRIVATE m)
add_test(NAME test_4 COMMAND test_4)

add_executable(test_5)
target_link_libraries(test_5 PRIVATE m)
add_test(NAME test_5 COMMAND test_5)

add_executable(test_6)
target_link_libraries(test_6 PRIVATE m)
add_test(NAME test_6 COMMAND test_6)

add_executable(test_7)
target_link_libraries(test_7 PRIVATE m)
add_test(NAME test_7 COMMAND test_7)

add_executable(test_8)
target_link_libraries(test_8 PRIVATE m)
add_test(NAME test_8 COMMAND test_8)