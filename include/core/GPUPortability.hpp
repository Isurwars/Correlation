/**
 * @file GPUPortability.hpp
 * @brief Unified GPU portability layer between CUDA and HIP with host fallback stubs.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

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
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
inline hipError_t hipMemcpy([[maybe_unused]] void *dst, [[maybe_unused]] const void *src,
                            [[maybe_unused]] std::size_t count,
                            [[maybe_unused]] int kind) noexcept {
  return 0;
}
// NOLINTEND(bugprone-easily-swappable-parameters)
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
  unsigned int x{1}, y{1}, z{1};
  constexpr dim3(unsigned int x_dim = 1, unsigned int y_dim = 1, unsigned int z_dim = 1) noexcept  // NOLINT(bugprone-easily-swappable-parameters)
      : x(x_dim), y(y_dim), z(z_dim) {}
};

inline constexpr dim3 threadIdx{0, 0, 0};
inline constexpr dim3 blockIdx{0, 0, 0};
inline constexpr dim3 blockDim{1, 1, 1};

// NOLINTBEGIN(bugprone-reserved-identifier, cert-dcl37-c, cert-dcl51-cpp)
#ifndef __device__
#define __device__
#endif
#ifndef __global__
#define __global__
#endif
#ifndef __restrict__
#define __restrict__
#endif
#ifndef __forceinline__
#define __forceinline__ inline
#endif
// NOLINTEND(bugprone-reserved-identifier, cert-dcl37-c, cert-dcl51-cpp)

#define hipLaunchKernelGGL(kernel, grid, block, shared, stream, ...)

#endif
