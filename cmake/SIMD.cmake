# SIMD Feature Detection
# Checks both the host CPU runtime capabilities AND compiler
# flag support to avoid building ISA code the CPU cannot run.

include(CheckCXXCompilerFlag)

option(CORRELATION_OPTIMIZE_NATIVE "Enable architecture-specific optimizations (AVX2/AVX512)" OFF)

set(CORRELATION_SIMD_FLAGS "")
set(HOST_HAS_AVX512 FALSE)
set(HOST_HAS_AVX2 FALSE)

if (CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
  execute_process(
    COMMAND grep -m1 flags /proc/cpuinfo
    OUTPUT_VARIABLE CPU_FLAGS_LINE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if (CPU_FLAGS_LINE MATCHES "avx512f")
    set(HOST_HAS_AVX512 TRUE)
  endif()
  if (CPU_FLAGS_LINE MATCHES "avx2")
    set(HOST_HAS_AVX2 TRUE)
  endif()
elseif (CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin")
  execute_process(
    COMMAND sysctl -n machdep.cpu.features
    OUTPUT_VARIABLE CPU_FLAGS_LINE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  string(TOLOWER "${CPU_FLAGS_LINE}" CPU_FLAGS_LINE)
  if (CPU_FLAGS_LINE MATCHES "avx512f")
    set(HOST_HAS_AVX512 TRUE)
  endif()
  if (CPU_FLAGS_LINE MATCHES "avx2")
    set(HOST_HAS_AVX2 TRUE)
  endif()
elseif (CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows")
  # On Windows, we only enable SIMD if CORRELATION_OPTIMIZE_NATIVE is ON.
  # This prevents 'Illegal Instruction' crashes on older CPUs when the compiler
  # capability is misaligned with the host hardware.
  if (CORRELATION_OPTIMIZE_NATIVE)
    check_cxx_compiler_flag("/arch:AVX512" COMPILER_SUPPORTS_AVX512)
    if (COMPILER_SUPPORTS_AVX512)
      set(HOST_HAS_AVX512 TRUE)
    endif()
    check_cxx_compiler_flag("/arch:AVX2" COMPILER_SUPPORTS_AVX2)
    if (COMPILER_SUPPORTS_AVX2)
      set(HOST_HAS_AVX2 TRUE)
    endif()
  endif()
endif()

if (NOT MSVC)
  if (HOST_HAS_AVX512)
    check_cxx_compiler_flag("-mavx512f" COMPILER_SUPPORTS_AVX512)
    if (COMPILER_SUPPORTS_AVX512)
      set(CORRELATION_SIMD_FLAGS "-mavx512f -mfma")
      message(STATUS "SIMD: AVX-512 enabled (-mavx512f -mfma)")
    endif()
  endif()
  if (NOT CORRELATION_SIMD_FLAGS AND HOST_HAS_AVX2)
    check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
    if (COMPILER_SUPPORTS_AVX2)
      set(CORRELATION_SIMD_FLAGS "-mavx2 -mfma")
      message(STATUS "SIMD: AVX2 enabled (-mavx2 -mfma)")
    endif()
  endif()
  if (NOT CORRELATION_SIMD_FLAGS)
    message(STATUS "SIMD: No AVX support detected – using scalar fallback")
  endif()
else()
  if (HOST_HAS_AVX512 AND COMPILER_SUPPORTS_AVX512)
    set(CORRELATION_SIMD_FLAGS "/arch:AVX512")
    message(STATUS "SIMD: AVX-512 enabled (/arch:AVX512)")
  elseif (HOST_HAS_AVX2 AND COMPILER_SUPPORTS_AVX2)
    set(CORRELATION_SIMD_FLAGS "/arch:AVX2")
    message(STATUS "SIMD: AVX2 enabled (/arch:AVX2)")
  else()
    if (CORRELATION_OPTIMIZE_NATIVE)
      message(STATUS "SIMD: Requested optimization but no compatible AVX level found")
    else()
      message(STATUS "SIMD: Optimization disabled – using scalar fallback")
    endif()
  endif()
endif()
