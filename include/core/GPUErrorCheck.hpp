/**
 * @file GPUErrorCheck.hpp
 * @brief Throwing error-check function for HIP/CUDA GPU API calls.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/GPUPortability.hpp"

#include <source_location>
#include <stdexcept>

namespace correlation::core::gpu {

/**
 * @brief Exception type for GPU runtime errors.
 */
class GPUError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

/**
 * @brief Validates a HIP/CUDA runtime call result and throws GPUError on failure.
 * @param[in] result The error code returned by the GPU runtime call.
 * @param[in] loc The source location where the check was invoked.
 * @throws GPUError If result is not equal to hipSuccess.
 */
inline void hipCheck(hipError_t result, std::source_location const loc = std::source_location::current()) {
#if defined(CORRELATION_USE_CUDA) || defined(CORRELATION_USE_HIP)
  if (result != hipSuccess) {
    throw GPUError(std::string("GPU error in ") + loc.file_name() + ":" + std::to_string(loc.line()) + " — " +
                   hipGetErrorString(result));
  }
#else
  (void)result;
  (void)loc;
#endif
}

} // namespace correlation::core::gpu
