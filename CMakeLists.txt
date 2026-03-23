# Correlation - Liquid and Amorphous Solid Analysis Tool
# Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
# SPDX-License-Identifier: MIT
# Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

cmake_minimum_required(VERSION 3.21)
project(Correlation VERSION 2.2.0 LANGUAGES CXX)
set(CORRELATION_VERSION "2.2.0")

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
    configure_file("${PROJECT_SOURCE_DIR}/packaging/correlation.desktop"
                   "${PROJECT_BINARY_DIR}/correlation.desktop"
                   COPYONLY)
    install(FILES "${PROJECT_BINARY_DIR}/correlation.desktop"
            DESTINATION share/applications)
endif()

# -----------------------------------------------------------
# CPack Packaging
# -----------------------------------------------------------
set(CPACK_PACKAGE_NAME "correlation")
set(CPACK_PACKAGE_VENDOR "Isaías Rodríguez")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Analysis tool for liquid and amorphous solid structures")
set(CPACK_PACKAGE_VERSION "${CORRELATION_VERSION}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "correlation")
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

include(CPack)
