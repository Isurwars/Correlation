# Correlation - Liquid and Amorphous Solid Analysis Tool
# Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
# SPDX-License-Identifier: MIT
# Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

cmake_minimum_required(VERSION 3.21)
project(Correlation VERSION 2.0.0 LANGUAGES CXX)
set(CORRELATION_VERSION "2.0.0")

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
    src/calculators/AngleCalculator.cpp
    src/calculators/DihedralCalculator.cpp
    src/calculators/MotifFinder.cpp
    src/calculators/CNCalculator.cpp
    src/calculators/DistanceCalculator.cpp
    src/calculators/RDFCalculator.cpp
    src/calculators/PADCalculator.cpp
    src/calculators/DADCalculator.cpp
    src/calculators/MotifFinder.cpp
    src/calculators/RDCalculator.cpp
    src/calculators/VACFCalculator.cpp
    src/calculators/VDOSCalculator.cpp
    src/calculators/SQCalculator.cpp
    src/calculators/XRDCalculator.cpp
    src/DynamicsAnalyzer.cpp
    src/FileReader.cpp
    src/readers/ArcReader.cpp
    src/readers/CarReader.cpp
    src/readers/CellReader.cpp
    src/readers/CastepMdReader.cpp
    src/readers/CifReader.cpp
    src/readers/LammpsDumpReader.cpp
    src/readers/OnetepDatReader.cpp
    src/readers/OutmolReader.cpp
    src/FileWriter.cpp
    src/writers/ArrowWriter.cpp
    src/writers/CSVWriter.cpp
    src/writers/HDF5Writer.cpp
    src/NeighborGraph.cpp
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
target_link_libraries(correlation_lib PUBLIC HighFive HDF5::HDF5 HDF5::HDF5_HL arrow_shared parquet_shared)

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

# Copy DLLs to the executable directory so it can be run without adding to PATH
if (WIN32)
    add_custom_command(TARGET correlation POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:correlation> $<TARGET_FILE_DIR:correlation>
        COMMAND_EXPAND_LISTS
    )
endif()


# Add the include directory for the generated slint header
target_include_directories(correlation PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}
)

slint_target_sources(correlation
  ui/AppWindow.slint
  LIBRARY_PATHS material=${CMAKE_CURRENT_SOURCE_DIR}/material-1.0/material.slint
)

set_property(TARGET correlation PROPERTY SLINT_EMBED_RESOURCES "embed-files")

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
# Install setup for Linux shared libraries
# -----------------------------------------------------------
if (NOT WIN32)
    # Ensure correlation finds shared libraries in <prefix>/lib
    set_target_properties(correlation PROPERTIES
        INSTALL_RPATH "$ORIGIN/../lib"
    )
    # If Arrow/Parquet were built from source, install them
    if (TARGET arrow_shared AND TARGET parquet_shared AND NOT Arrow_FOUND)
        install(TARGETS arrow_shared parquet_shared
            LIBRARY DESTINATION lib
            RUNTIME DESTINATION bin
        )
    endif()
endif()

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

# -----------------------------------------------------------
# CPack Packaging
# -----------------------------------------------------------
set(CPACK_PACKAGE_NAME "Correlation")
set(CPACK_PACKAGE_VENDOR "Isaías Rodríguez")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Analysis tool for liquid and amorphous solid structures")
set(CPACK_PACKAGE_VERSION "${CORRELATION_VERSION}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "Correlation")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.md")

if (WIN32)
    set(CPACK_GENERATOR "NSIS")
    set(CPACK_NSIS_PACKAGE_NAME "Correlation Analysis Tool")
    set(CPACK_NSIS_MUI_ICON "${PROJECT_SOURCE_DIR}/Images/AppIcon.ico") # Corrected icon name
    set(CPACK_NSIS_HELP_LINK "https://github.com/Isurwars/Correlation")
    set(CPACK_NSIS_URL_INFO_ABOUT "https://github.com/Isurwars/Correlation")
    set(CPACK_NSIS_CONTACT "isurwars@gmail.com")
    set(CPACK_NSIS_MODIFY_PATH ON)
elseif (APPLE)
    set(CPACK_GENERATOR "DragNDrop")
    set(CPACK_DMG_VOLUME_NAME "Correlation")
    set(CPACK_SYSTEM_NAME "macOS")
else ()
    # Linux
    set(CPACK_GENERATOR "DEB;TGZ")
    set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Isaías Rodríguez <isurwars@gmail.com>")
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "libtbb12 (>= 2021.5), libhdf5-103-1 | libhdf5-103") # Approximated dependencies
    set(CPACK_DEBIAN_PACKAGE_SECTION "science")
    set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
    set(CPACK_STRIP_FILES TRUE)
endif ()

include(CPack)
