cmake_minimum_required(VERSION 3.10)

project(correlation)

include_directories(../include)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(correlation
    Atom.cpp
    Cell.cpp
    Constants.cpp
    ReadFiles.cpp
    Smoothing.cpp
    StructureFactor.cpp
    WriteFiles.cpp
    main.cpp
)

target_link_libraries(correlation m)