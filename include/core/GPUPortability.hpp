/**
 * @file GPUPortability.hpp
 * @brief Unified GPU portability layer between CUDA and HIP.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include <utility>

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

// Portable kernel launch function
template <typename K, typename... Args>
inline void hipLaunchKernelGGL(K kernel, dim3 grid, dim3 block, size_t shared, cudaStream_t stream, Args&&... args) {
    kernel<<<grid, block, shared, stream>>>(std::forward<Args>(args)...);
}

#endif

