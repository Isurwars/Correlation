/**
 * @file GPUPortability.hpp
 * @brief Unified GPU portability layer between CUDA and HIP with host fallback stubs.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/CompilerPortability.hpp" // IWYU pragma: export
#include <cstddef>                      // IWYU pragma: keep
#include <utility>                      // IWYU pragma: keep

#ifdef CORRELATION_USE_HIP
#include <hip/hip_runtime.h>
#elif defined(CORRELATION_USE_CUDA)
#include <cuda_runtime.h>

// Map HIP API names to CUDA equivalents
using hipError_t = cudaError_t;
#define hipSuccess cudaSuccess
#define hipGetDeviceCount cudaGetDeviceCount
#define hipMalloc cudaMalloc
#define hipFree cudaFree
#define hipMemcpy cudaMemcpy
#define hipMemcpyHostToDevice cudaMemcpyHostToDevice
#define hipMemcpyDeviceToHost cudaMemcpyDeviceToHost
#define hipDeviceSynchronize cudaDeviceSynchronize
#define hipGetErrorString cudaGetErrorString

#if defined(__CUDACC__)
/**
 * @brief Portable kernel launch function mapping to CUDA chevron syntax.
 * @tparam K Kernel function type.
 * @tparam Args Kernel argument types.
 * @param[in] kernel Kernel function.
 * @param[in] grid Grid dimensions.
 * @param[in] block Block dimensions.
 * @param[in] shared Dynamic shared memory in bytes.
 * @param[in] stream CUDA execution stream.
 * @param[in] args Arguments forwarded to the kernel.
 */
template <typename K, typename... Args>
inline void hipLaunchKernelGGL(K kernel, dim3 grid, dim3 block, std::size_t shared,
                               cudaStream_t stream, Args &&...args) {
  kernel<<<grid, block, shared, stream>>>(std::forward<Args>(args)...);
}
#endif

#else

/// Return error type for host fallback emulation.
using hipError_t = int;
/// Successful execution code.
constexpr int hipSuccess = 0;
/// Device unavailable error code.
constexpr int hipErrorNoDevice = 100;
/// Operation not supported in host fallback mode.
constexpr int hipErrorNotSupported = 801;

/**
 * @brief Returns the number of compute-capable GPU devices.
 * @param[out] count Pointer to integer receiving device count (always 0 in host fallback).
 * @return hipSuccess.
 */
inline hipError_t hipGetDeviceCount(int *count) noexcept {
  if (count != nullptr) {
    *count = 0;
  }
  return hipSuccess;
}

/**
 * @brief Host fallback allocation stub returning failure when GPU support is disabled.
 * @tparam T Pointer element type.
 * @param[out] ptr Output pointer set to nullptr.
 * @param[in] size Requested byte size.
 * @return hipErrorNoDevice indicating no GPU acceleration is compiled in.
 */
template <typename T>
inline hipError_t hipMalloc(T **ptr, [[maybe_unused]] std::size_t size) noexcept {
  if (ptr != nullptr) {
    *ptr = nullptr;
  }
  return hipErrorNoDevice;
}

/**
 * @brief Host fallback deallocation stub.
 * @param[in] ptr Pointer to free (ignored in fallback mode).
 * @return hipSuccess.
 */
inline hipError_t hipFree([[maybe_unused]] void *ptr) noexcept { return hipSuccess; }

/**
 * @brief Memory copy parameter struct for host fallback.
 */
struct MemcpyParams {
  void *dst{nullptr};
  const void *src{nullptr};
  std::size_t count{0};
  int kind{0};
};

/**
 * @brief Host fallback memory copy stub.
 * @param[in] params Memory copy parameters.
 * @return hipErrorNoDevice.
 */
inline hipError_t hipMemcpy([[maybe_unused]] MemcpyParams params) noexcept {
  return hipErrorNoDevice;
}

/**
 * @brief Host fallback memory copy overload.
 * @param[out] dst Destination buffer pointer.
 * @param[in] src Source buffer pointer.
 * @param[in] count Byte count.
 * @param[in] kind Transfer direction.
 * @return hipErrorNoDevice.
 */
inline hipError_t hipMemcpy(void *dst, const void *src, std::size_t count, int kind) noexcept {
  return hipMemcpy(MemcpyParams{
      .dst = dst,
      .src = src,
      .count = count,
      .kind = kind,
  });
}

/// Host-to-device memory copy transfer direction flag.
constexpr int hipMemcpyHostToDevice = 0;
/// Device-to-host memory copy transfer direction flag.
constexpr int hipMemcpyDeviceToHost = 0;

/**
 * @brief Host fallback device synchronization stub.
 * @return hipErrorNoDevice.
 */
inline hipError_t hipDeviceSynchronize() noexcept { return hipErrorNoDevice; }

/**
 * @brief Translates an error code to a human-readable description.
 * @param[in] err Error code.
 * @return Static string description of the error code.
 */
inline const char *hipGetErrorString(hipError_t err) noexcept {
  switch (err) {
  case hipSuccess:
    return "Success";
  case hipErrorNoDevice:
    return "No GPU device available (compiled without CUDA/HIP support)";
  case hipErrorNotSupported:
    return "GPU operation not supported in host fallback mode";
  default:
    return "Unknown GPU error";
  }
}

/**
 * @brief Host fallback atomic addition.
 * @tparam T Address data type.
 * @tparam U Increment value type.
 * @param[in,out] addr Target memory location.
 * @param[in] val Value to add.
 * @return Old value before addition.
 */
template <typename T, typename U> inline T atomicAdd(T *addr, U val) noexcept {
  if (addr == nullptr) {
    return T{};
  }
  T const old = *addr;
  *addr += static_cast<T>(val);
  return old;
}

/**
 * @brief 3D vector dimension descriptor for CUDA/HIP emulation.
 */
struct dim3 {
  unsigned int x{1};
  unsigned int y{1};
  unsigned int z{1};
};

/// Emulated thread index inside thread block.
inline constexpr dim3 threadIdx{
    .x = 0,
    .y = 0,
    .z = 0,
};
/// Emulated block index inside grid.
inline constexpr dim3 blockIdx{
    .x = 0,
    .y = 0,
    .z = 0,
};
/// Emulated block dimension.
inline constexpr dim3 blockDim{
    .x = 1,
    .y = 1,
    .z = 1,
};

/**
 * @brief Host fallback stub for kernel launch when GPU acceleration is disabled.
 * @tparam K Kernel function type.
 * @tparam Grid Grid dimension type.
 * @tparam Block Block dimension type.
 * @tparam Shared Shared memory size type.
 * @tparam Stream Stream type.
 * @tparam Args Kernel argument types.
 * @param[in] kernel Kernel function pointer or callable.
 * @param[in] grid Grid dimensions.
 * @param[in] block Block dimensions.
 * @param[in] shared Shared memory allocation size in bytes.
 * @param[in] stream Execution stream.
 * @param[in] args Kernel launch arguments.
 */
template <typename K, typename Grid, typename Block, typename Shared, typename Stream,
          typename... Args>
constexpr void hipLaunchKernelGGL([[maybe_unused]] K &&kernel, [[maybe_unused]] Grid &&grid,
                                  [[maybe_unused]] Block &&block, [[maybe_unused]] Shared &&shared,
                                  [[maybe_unused]] Stream &&stream,
                                  [[maybe_unused]] Args &&...args) noexcept {}

#endif
