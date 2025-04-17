cmake_minimum_required(VERSION 3.20)

project(correlation)

include_directories(../include)

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
