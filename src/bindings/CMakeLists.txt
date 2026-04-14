pybind11_add_module(correlation_py
    module.cpp
    math_bindings.cpp
    core_bindings.cpp
    io_bindings.cpp
    analysis_bindings.cpp
    calculator_bindings.cpp
    writer_bindings.cpp
)

set_target_properties(correlation_py PROPERTIES OUTPUT_NAME "correlation")

target_link_libraries(correlation_py PRIVATE correlation_lib)

# Match the C++ standard used by correlation_lib
set_target_properties(correlation_py PROPERTIES CXX_STANDARD 23)

# If SIMD is enabled, compile the module with the same SIMD flags to match ABI
if (CORRELATION_SIMD_FLAGS)
    separate_arguments(SIMD_FLAGS_LIST NATIVE_COMMAND "${CORRELATION_SIMD_FLAGS}")
    target_compile_options(correlation_py PRIVATE ${SIMD_FLAGS_LIST})
endif()
