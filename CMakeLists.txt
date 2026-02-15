# Correlation - Liquid and Amorphous Solid Analysis Tool
# Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
# SPDX-License-Identifier: MIT
# Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

cmake_minimum_required(VERSION 3.21)
project(Correlation VERSION 2.0.0 LANGUAGES CXX)

# Set the C++ standard as a target property, not a global variable.
set(CXX_STANDARD 23)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set default install path
set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation directory")

# -----------------------------------------------------------
# Utils & Dependencies
# -----------------------------------------------------------
include(cmake/Utils.cmake)
include(cmake/Dependencies.cmake)

#-----------------------------------------------------------
# Library setup
#-----------------------------------------------------------

# Add sources to library
add_library(correlation_lib STATIC
    src/Cell.cpp
    src/DistributionFunctions.cpp
    src/DynamicsAnalyzer.cpp
    src/FileIO.cpp
    src/FileWriter.cpp
    src/StructureAnalyzer.cpp
    src/Trajectory.cpp
    src/TrajectoryAnalyzer.cpp
)

# Include directories
target_include_directories(correlation_lib PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(correlation_lib PUBLIC cxx_std_${CXX_STANDARD})
target_compile_definitions(correlation_lib PUBLIC _PSTL_PAR_BACKEND_TBB)
target_link_libraries(correlation_lib PRIVATE TBB::tbb TBB::tbbmalloc)
target_link_libraries(correlation_lib PUBLIC HighFive HDF5::HDF5 HDF5::HDF5_HL)

#------------------------------------------------------------------
# Executable Target (correlation)
#------------------------------------------------------------------
set(CORRELATION_SOURCES
    src/main.cpp
    src/AppBackend.cpp
    src/AppController.cpp
)

if (WIN32)
    list(APPEND CORRELATION_SOURCES resource.rc)
    add_executable(correlation WIN32 ${CORRELATION_SOURCES})
else()
    add_executable(correlation ${CORRELATION_SOURCES})
endif()


# Link correlation_lib to correlation
set(CORRELATION_LIBS correlation_lib Slint::Slint)

if (WIN32)
    list(APPEND CORRELATION_LIBS Imm32 Comctl32 dwrite propsys opengl32 dwmapi uxtheme)
endif()

target_link_libraries(correlation PRIVATE ${CORRELATION_LIBS})


# Add the include directory for the generated slint header
target_include_directories(correlation PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}
)

slint_target_sources(correlation
  ui/AppWindow.slint
  LIBRARY_PATHS material=${CMAKE_CURRENT_SOURCE_DIR}/material-1.0/material.slint
)

#-----------------------------------------------------------
# Testing
#-----------------------------------------------------------
enable_testing()
add_subdirectory(tests)

#-----------------------------------------------------------
# Install setup
#-----------------------------------------------------------

# Install the executable and the library.
install(TARGETS correlation
    RUNTIME DESTINATION bin
)

# -----------------------------------------------------------
# Install setup for Linux Desktop Integration
# -----------------------------------------------------------
if (NOT WIN32)
    # 1. Install the Icon (adjust path if you use a different icon directory)
    install(FILES "${PROJECT_SOURCE_DIR}/Images/Logo.png"
            DESTINATION share/icons/hicolor/scalable/apps
            RENAME correlation.png)

    # 2. Configure and Install the .desktop file
    configure_file("${PROJECT_SOURCE_DIR}/correlation.desktop"
                   "${PROJECT_BINARY_DIR}/correlation.desktop"
                   COPYONLY)
    install(FILES "${PROJECT_BINARY_DIR}/correlation.desktop"
            DESTINATION share/applications)
endif()
