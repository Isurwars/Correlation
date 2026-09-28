// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

/**
 * @file GPUPortabilityTests.cpp
 * @brief Unit tests for GPU portability layer, error checking, and host fallback safety.
 */

#include "core/DeviceBuffer.hpp"
#include "core/GPUErrorCheck.hpp"
#include "core/GPUPortability.hpp"

#include <gtest/gtest.h>
#include <string>

namespace correlation::core::gpu::testing {

namespace {

TEST(GPUPortabilityTests, HostFallbackDeviceCountReturnsSuccess) {
  int count = -1;
  hipError_t const err = hipGetDeviceCount(&count);
  EXPECT_EQ(err, hipSuccess);
#if !defined(CORRELATION_USE_CUDA) && !defined(CORRELATION_USE_HIP)
  EXPECT_EQ(count, 0);
#endif
}

TEST(GPUPortabilityTests, HostFallbackMallocFailsSafe) {
#if !defined(CORRELATION_USE_CUDA) && !defined(CORRELATION_USE_HIP)
  float *ptr = reinterpret_cast<float *>(0xDEADBEEF);
  hipError_t const err = hipMalloc(&ptr, sizeof(float) * 10);
  EXPECT_EQ(err, hipErrorNoDevice);
  EXPECT_EQ(ptr, nullptr);
#endif
}

TEST(GPUPortabilityTests, HostFallbackFreeReturnsSuccess) {
  hipError_t const err = hipFree(nullptr);
  EXPECT_EQ(err, hipSuccess);
}

TEST(GPUPortabilityTests, HostFallbackMemcpyAndSyncFailFast) {
#if !defined(CORRELATION_USE_CUDA) && !defined(CORRELATION_USE_HIP)
  float src = 1.0F;
  float dst = 0.0F;
  hipError_t const err_memcpy = hipMemcpy(&dst, &src, sizeof(float), hipMemcpyHostToDevice);
  EXPECT_EQ(err_memcpy, hipErrorNoDevice);

  hipError_t const err_sync = hipDeviceSynchronize();
  EXPECT_EQ(err_sync, hipErrorNoDevice);
#endif
}

TEST(GPUPortabilityTests, ErrorStringDescriptionsAreAccurate) {
  std::string const success_msg = hipGetErrorString(hipSuccess);
  EXPECT_FALSE(success_msg.empty());

#if !defined(CORRELATION_USE_CUDA) && !defined(CORRELATION_USE_HIP)
  std::string const no_dev_msg = hipGetErrorString(hipErrorNoDevice);
  EXPECT_NE(no_dev_msg.find("compiled without"), std::string::npos);
#endif
}

TEST(GPUPortabilityTests, HipCheckThrowsOnFailureAndPassesOnSuccess) {
  EXPECT_NO_THROW(hipCheck(hipSuccess));

#if !defined(CORRELATION_USE_CUDA) && !defined(CORRELATION_USE_HIP)
  EXPECT_THROW(
      {
        try {
          hipCheck(hipErrorNoDevice);
        } catch (const GPUError &ex) {
          EXPECT_NE(std::string(ex.what()).find("compiled without"), std::string::npos);
          throw;
        }
      },
      GPUError);
#endif
}

TEST(GPUPortabilityTests, DeviceBufferBehaviorAcrossArchitectures) {
  // Default constructor should always succeed and hold nullptr
  DeviceBuffer<double> empty_buf;
  EXPECT_EQ(empty_buf.get(), nullptr);
  EXPECT_EQ(empty_buf.size(), 0U);

#if !defined(CORRELATION_USE_CUDA) && !defined(CORRELATION_USE_HIP)
  // Non-empty allocation in host fallback must throw GPUError rather than silently returning
  // nullptr
  EXPECT_THROW((DeviceBuffer<double>(32)), GPUError);
#endif
}

} // namespace

} // namespace correlation::core::gpu::testing
