pybind11_add_module(correlation_py
    module.cpp
    core_bindings.cpp
    io_bindings.cpp
)

set_target_properties(correlation_py PROPERTIES OUTPUT_NAME "correlation")

target_link_libraries(correlation_py PRIVATE correlation_lib)

# It seems correlation_lib is built with CXX23, we should ensure pybind target is the same
set_target_properties(correlation_py PROPERTIES CXX_STANDARD 23)

# If SIMD is enabled, we need to compile the module with SIMD flags as well to match ABI
if (CORRELATION_SIMD_FLAGS)
    separate_arguments(SIMD_FLAGS_LIST NATIVE_COMMAND "${CORRELATION_SIMD_FLAGS}")
    target_compile_options(correlation_py PRIVATE ${SIMD_FLAGS_LIST})
endif()
