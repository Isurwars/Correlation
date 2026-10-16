# Correlation - Liquid and Amorphous Solid Analysis Tool
# Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
# SPDX-License-Identifier: AGPL-3.0-only
# Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

add_executable(correlation_benchmarks
    main.cpp
    BM_DistanceCalculator.cpp
    BM_MotifFinder.cpp
    BM_GraphDescriptors.cpp
)

target_include_directories(correlation_benchmarks PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_compile_definitions(correlation_benchmarks PRIVATE _USE_MATH_DEFINES)

target_link_libraries(correlation_benchmarks PRIVATE
    correlation_lib
    benchmark::benchmark
)

if(BUILD_WITH_CUDA)
    set_target_properties(correlation_benchmarks PROPERTIES CUDA_STANDARD 20)
endif()
if(BUILD_WITH_HIP)
    set_target_properties(correlation_benchmarks PROPERTIES HIP_STANDARD 20)
endif()
