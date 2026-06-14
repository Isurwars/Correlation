# Correlation - Liquid and Amorphous Solid Analysis Tool
# Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
# SPDX-License-Identifier: AGPL-3.0-only
# Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#
# Fuzz testing targets using libFuzzer (Clang-only).
# Build with: CC=clang CXX=clang++ cmake -B build-fuzz -G Ninja \
#             -DBUILD_FUZZ_TESTS=ON -DBUILD_GUI=OFF
#             cmake --build build-fuzz

if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(STATUS "Fuzz tests skipped: Clang compiler (libFuzzer) is required.")
    return()
endif()

# ---------------------------------------------------------------
# Sanitizer flags for all fuzz targets
# ---------------------------------------------------------------
set(FUZZ_SANITIZER_FLAGS -fsanitize=fuzzer,address,undefined)
set(FUZZ_COMPILE_FLAGS   ${FUZZ_SANITIZER_FLAGS} -g -O1 -fno-omit-frame-pointer)
set(FUZZ_LINK_FLAGS      ${FUZZ_SANITIZER_FLAGS})

# Common libraries every fuzz target needs
set(FUZZ_LINK_LIBS
    correlation_lib
    readers_obj
    writers_obj
    calculators_obj
)

# ---------------------------------------------------------------
# Helper macro: define a fuzz target with its corpus directory
# ---------------------------------------------------------------
macro(add_fuzz_target TARGET_NAME CORPUS_DIR)
    add_executable(${TARGET_NAME} ${TARGET_NAME}.cpp)

    target_compile_options(${TARGET_NAME} PRIVATE ${FUZZ_COMPILE_FLAGS})
    target_link_options(${TARGET_NAME}    PRIVATE ${FUZZ_LINK_FLAGS})

    target_include_directories(${TARGET_NAME} PRIVATE
        ${PROJECT_SOURCE_DIR}/include
    )

    target_link_libraries(${TARGET_NAME} PRIVATE ${FUZZ_LINK_LIBS})
    target_link_libraries(${TARGET_NAME} PRIVATE correlation_simd_options)

    # Register a CTest smoke test: run for 10 seconds via wrapper script to auto-cleanup.
    add_test(
        NAME ${TARGET_NAME}
        COMMAND ${CMAKE_COMMAND}
            "-DFUZZ_TARGET=$<TARGET_FILE:${TARGET_NAME}>"
            "-DSEED_DIR=${CMAKE_CURRENT_SOURCE_DIR}/corpus/${CORPUS_DIR}"
            "-DOUT_DIR=${CMAKE_CURRENT_BINARY_DIR}/corpus_out_${CORPUS_DIR}"
            -P "${PROJECT_SOURCE_DIR}/cmake/RunFuzzTest.cmake"
    )
    set_tests_properties(${TARGET_NAME} PROPERTIES
        LABELS "fuzz"
        TIMEOUT 30
    )
endmacro()

# ---------------------------------------------------------------
# Fuzz targets — one per reader
# ---------------------------------------------------------------

# In-memory readers (no temp file needed — fastest)
add_fuzz_target(fuzz_lammps_dump  lammps)
add_fuzz_target(fuzz_gromacs      gromacs)

# File-based readers
add_fuzz_target(fuzz_vasp         vasp)
add_fuzz_target(fuzz_xdatcar      xdatcar)
add_fuzz_target(fuzz_pdb          pdb)
add_fuzz_target(fuzz_cif          cif)
add_fuzz_target(fuzz_xyz          xyz)
add_fuzz_target(fuzz_arc          arc)
add_fuzz_target(fuzz_car          car)
add_fuzz_target(fuzz_cell         cell)
add_fuzz_target(fuzz_castep_md    castep_md)
add_fuzz_target(fuzz_cp2k         cp2k)
add_fuzz_target(fuzz_onetep_dat   onetep_dat)
add_fuzz_target(fuzz_outmol       outmol)
add_fuzz_target(fuzz_qe           qe)
