pybind11_add_module(_correlation
    module.cpp
    math_bindings.cpp
    core_bindings.cpp
    io_bindings.cpp
    analysis_bindings.cpp
    calculator_bindings.cpp
    writer_bindings.cpp
)

# Keep backward compatibility: the library file is named _correlation.so/.pyd
# The Python package python/correlation/__init__.py re-exports via:
#   from correlation._correlation import *

target_link_libraries(_correlation PRIVATE correlation_lib)

# Match the C++ standard used by correlation_lib
set_target_properties(_correlation PROPERTIES CXX_STANDARD 23)

# If SIMD is enabled, compile the module with the same SIMD flags to match ABI
if (CORRELATION_SIMD_FLAGS)
    separate_arguments(SIMD_FLAGS_LIST NATIVE_COMMAND "${CORRELATION_SIMD_FLAGS}")
    target_compile_options(_correlation PRIVATE ${SIMD_FLAGS_LIST})
endif()

# Install rule for scikit-build-core packaging
install(TARGETS _correlation DESTINATION correlation)

