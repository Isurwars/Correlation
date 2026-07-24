/**
 * @file GPUPortability.hpp
 * @brief Unified GPU portability layer between CUDA and HIP with host fallback stubs.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/CompilerPortability.hpp" // IWYU pragma: export
#include <cstddef>
#include <utility> // IWYU pragma: keep

#if defined(CORRELATION_USE_HIP)
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

#if defined(__CUDACC__)
// Portable kernel launch function
template <typename K, typename... Args>
inline void hipLaunchKernelGGL(K kernel, dim3 grid, dim3 block, std::size_t shared, cudaStream_t stream,
                               Args &&...args) {
  kernel<<<grid, block, shared, stream>>>(std::forward<Args>(args)...);
}
#endif

#else

using hipError_t = int;
constexpr int hipSuccess = 0;

inline hipError_t hipGetDeviceCount(int *count) {
  if (count != nullptr) {
    *count = 0;
  }
  return hipSuccess;
}

template <typename T> inline hipError_t hipMalloc(T **ptr, [[maybe_unused]] std::size_t size) {
  if (ptr) {
    *ptr = nullptr;
  }
  return 0;
}

inline hipError_t hipFree([[maybe_unused]] void *ptr) { return 0; }

struct MemcpyParams {
  void *dst{nullptr};
  const void *src{nullptr};
  std::size_t count{0};
  int kind{0};
};

inline hipError_t hipMemcpy([[maybe_unused]] MemcpyParams params) noexcept { return 0; }

inline hipError_t hipMemcpy(void *dst, const void *src, std::size_t count, int kind) noexcept {
  return hipMemcpy(MemcpyParams{
      .dst = dst,
      .src = src,
      .count = count,
      .kind = kind,
  });
}

constexpr int hipMemcpyHostToDevice = 0;
constexpr int hipMemcpyDeviceToHost = 0;
inline hipError_t hipDeviceSynchronize() { return 0; }

template <typename T, typename U> inline T atomicAdd(T *addr, U val) {
  if (addr != nullptr) {
    return T{};
  }
  T old = *addr;
  *addr += static_cast<T>(val);
  return old;
}

struct dim3 {
  unsigned int x{1};
  unsigned int y{1};
  unsigned int z{1};
};

inline constexpr dim3 threadIdx{
    .x = 0,
    .y = 0,
    .z = 0,
};
inline constexpr dim3 blockIdx{
    .x = 0,
    .y = 0,
    .z = 0,
};
inline constexpr dim3 blockDim{
    .x = 1,
    .y = 1,
    .z = 1,
};

#define hipLaunchKernelGGL(kernel, grid, block, shared, stream, ...)

#endif
