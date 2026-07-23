/**
 * @file SIMDUtils.hpp
 * @brief Facade header for SIMD-accelerated kernels for distance, integration, and normalization.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/SIMDConfig.hpp"

namespace correlation::math {

/**
 * @brief Returns a string describing the compiled SIMD instruction set level.
 * @return String representation of the SIMD level ("AVX-512", "AVX2", or "Scalar").
 */
[[nodiscard]] inline const char *simd_level_string() noexcept {
#ifdef CORRELATION_SIMD_AVX512
  return "AVX-512";
#elif defined(CORRELATION_SIMD_AVX2)
  return "AVX2";
#else
  return "Scalar";
#endif
}

} // namespace correlation::math
