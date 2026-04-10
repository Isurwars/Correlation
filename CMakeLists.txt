# Correlation - Liquid and Amorphous Solid Analysis Tool
# Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
# SPDX-License-Identifier: MIT
# Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

cmake_minimum_required(VERSION 3.21)
project(Correlation VERSION 2.6.3 LANGUAGES CXX)

# -----------------------------------------------------------
# Auto-generate packaging files from templates.
# correlation.iss is the single versioned output — edit
# packaging/correlation.iss.in, never the .iss directly.
# -----------------------------------------------------------
configure_file(
    "${PROJECT_SOURCE_DIR}/packaging/correlation.iss.in"
    "${PROJECT_SOURCE_DIR}/packaging/correlation.iss"
    @ONLY
)

include(GNUInstallDirs)

# Set the C++ standard globally.
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Suppress MSVC warning C4624: 'destructor was implicitly defined as deleted'
# This is frequently triggered by Slint's internal headers.
if (MSVC)
    add_compile_options(/wd4624)
endif()

# Set default install path
set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation directory")

# -----------------------------------------------------------
# Utils & Dependencies
# -----------------------------------------------------------
include(cmake/Utils.cmake)
include(cmake/Dependencies.cmake)

# -----------------------------------------------------------
# SIMD Feature Detection
# -----------------------------------------------------------
include(cmake/SIMD.cmake)

#-----------------------------------------------------------
# Targets setup
#-----------------------------------------------------------
add_subdirectory(src)

#-----------------------------------------------------------
# Testing
#-----------------------------------------------------------
include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

#-------------------------------------------
# Documentation
#-------------------------------------------
find_package(Doxygen)
if (DOXYGEN_FOUND)
    set(DOXYGEN_IN ${PROJECT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${PROJECT_BINARY_DIR}/docs/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    add_custom_target(doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
endif()

#-----------------------------------------------------------
# Install setup
#-----------------------------------------------------------

# Install the executable and the library.
install(TARGETS correlation
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT correlation
)

# -----------------------------------------------------------
# Install setup for Linux shared libraries
# -----------------------------------------------------------
if (NOT WIN32)
    # Ensure correlation finds shared libraries in <prefix>/lib
    set_target_properties(correlation PROPERTIES
        INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}"
    )
    # If Arrow/Parquet were built from source, install them
    if (TARGET arrow_shared AND TARGET parquet_shared AND NOT Arrow_FOUND)
        install(TARGETS arrow_shared parquet_shared
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            COMPONENT correlation
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
            RENAME correlation.png
            COMPONENT correlation)

    # 2. Configure and Install the .desktop file
    configure_file("${PROJECT_SOURCE_DIR}/packaging/correlation.desktop"
                   "${PROJECT_BINARY_DIR}/correlation.desktop"
                   COPYONLY)
    install(FILES "${PROJECT_BINARY_DIR}/correlation.desktop"
            DESTINATION share/applications
            COMPONENT correlation)
endif()

# -----------------------------------------------------------
# CPack Packaging
# -----------------------------------------------------------
set(CPACK_PACKAGE_NAME "correlation")
set(CPACK_PACKAGE_VENDOR "Isaías Rodríguez")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Analysis tool for liquid and amorphous solid structures")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "correlation")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.md")
set(CPACK_COMPONENTS_ALL correlation)

if (WIN32)
    set(CPACK_GENERATOR "NSIS")
    set(CPACK_PACKAGE_FILE_NAME "correlation-${PROJECT_VERSION}-win64")
    set(CPACK_NSIS_PACKAGE_NAME "Correlation Analysis Tool")
    set(CPACK_NSIS_DISPLAY_NAME "Correlation Analysis Tool")
    set(CPACK_PACKAGE_EXECUTABLES "correlation" "Correlation Analysis Tool")
    set(CPACK_CREATE_DESKTOP_LINKS "correlation")
    # Icon shown in Add/Remove Programs — relative to $INSTDIR (no leading backslash)
    set(CPACK_NSIS_INSTALLED_ICON_NAME "correlation.exe")
    string(REPLACE "/" "\\\\" CPACK_NSIS_MUI_ICON "${PROJECT_SOURCE_DIR}/Images/AppIcon.ico")
    string(REPLACE "/" "\\\\" CPACK_NSIS_MUI_UNIICON "${PROJECT_SOURCE_DIR}/Images/AppIcon.ico")
    # Branded installer graphics (24-bit BMP required by NSIS)
    string(REPLACE "/" "\\\\" CPACK_NSIS_MUI_WELCOMEFINISHPAGE_BITMAP "${PROJECT_SOURCE_DIR}/packaging/nsis_welcome.bmp")
    string(REPLACE "/" "\\\\" CPACK_NSIS_MUI_HEADERIMAGE "${PROJECT_SOURCE_DIR}/packaging/nsis_header.bmp")
    set(CPACK_NSIS_MUI_HEADERIMAGE_RIGHT ON)
    # Use default zlib compression instead of lzma to avoid false-positive virus detections
    # set(CPACK_NSIS_COMPRESSOR "lzma")
    set(CPACK_NSIS_HELP_LINK "https://github.com/Isurwars/Correlation")
    set(CPACK_NSIS_URL_INFO_ABOUT "https://github.com/Isurwars/Correlation")
    set(CPACK_NSIS_CONTACT "isurwars@gmail.com")
    set(CPACK_NSIS_MODIFY_PATH ON)
elseif (APPLE)
    set(CPACK_GENERATOR "DragNDrop")
    set(CPACK_DMG_VOLUME_NAME "correlation")
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

set(CPACK_INSTALL_CMAKE_PROJECTS "${PROJECT_BINARY_DIR};Correlation;correlation;/")
# Unset the internal CPack included guard to force CPack to overwrite
# the CPackConfig with our top-level settings, overriding fetched dependencies.
unset(CPack_CMake_INCLUDED)
include(CPack)
