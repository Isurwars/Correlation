if(BUILD_PYTHON_BINDINGS)
    pybind11_add_module(_correlation
        module.cpp
        math_bindings.cpp
        core_bindings.cpp
        io_bindings.cpp
        analysis_bindings.cpp
        calculator_bindings.cpp
        writer_bindings.cpp
        mlip_bindings.cpp
    )

    # Keep backward compatibility: the library file is named _correlation.so/.pyd
    # The Python package python/correlation/__init__.py re-exports via:
    #   from correlation._correlation import *

    target_link_libraries(_correlation PRIVATE correlation_lib readers_obj writers_obj calculators_obj)

    # Match the C++ standard used by correlation_lib
    set_target_properties(_correlation PROPERTIES CXX_STANDARD 23)

    # If SIMD is enabled, compile the module with the same SIMD flags to match ABI
    if (CORRELATION_SIMD_FLAGS)
        separate_arguments(SIMD_FLAGS_LIST NATIVE_COMMAND "${CORRELATION_SIMD_FLAGS}")
        target_compile_options(_correlation PRIVATE ${SIMD_FLAGS_LIST})
    endif()

    # Install rule for scikit-build-core packaging
    install(TARGETS _correlation DESTINATION correlation)
endif()

if(EMSCRIPTEN OR BUILD_WASM)
    add_executable(correlation_wasm
        wasm_bindings.cpp
    )

    target_include_directories(correlation_wasm PRIVATE
        ${PROJECT_SOURCE_DIR}/include
    )

    target_link_libraries(correlation_wasm PRIVATE
        correlation_lib
        readers_obj
        calculators_obj
        writers_obj
    )

    set_target_properties(correlation_wasm PROPERTIES
        CXX_STANDARD 23
        OUTPUT_NAME "correlation_wasm"
        SUFFIX ".js"
    )

    target_link_options(correlation_wasm PRIVATE
        -O3
        -msimd128
        --bind
        -s WASM=1
        -s MODULARIZE=1
        -s EXPORT_NAME="createCorrelationModule"
        -s ALLOW_MEMORY_GROWTH=1
        -s MAXIMUM_MEMORY=2GB
        -s FORCE_FILESYSTEM=1
        -s NO_DISABLE_EXCEPTION_CATCHING
        -s ENVIRONMENT=web,worker
    )
endif()

