cmake
cmake_minimum_required(VERSION 3.10)
project(Correlation)

include_directories(src test)

add_subdirectory(test)
add_subdirectory(src)

install(TARGETS correlation DESTINATION bin)

