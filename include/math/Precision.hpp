/**
 * @file Precision.hpp
 * @brief Floating point precision type definitions for Correlation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

namespace correlation {

#if defined(CORRELATION_USE_SINGLE_PRECISION) && CORRELATION_USE_SINGLE_PRECISION
using real_t = float;
inline constexpr const char *REAL_TYPE_NAME = "float";
inline constexpr bool is_single_precision = true;
#else
using real_t = double;
inline constexpr const char *REAL_TYPE_NAME = "double";
inline constexpr bool is_single_precision = false;
#endif

/**
 * @struct KahanAccumulator
 * @brief Performs Kahan-Babuška-Neumaier compensated summation to prevent floating-point loss of precision.
 * @tparam T Floating point type (defaults to real_t).
 */
template <typename T = real_t> struct KahanAccumulator {
  T sum{0};
  T compensation{0};

  constexpr void add(T val) noexcept {
    T y_val = val - compensation;
    T t_val = sum + y_val;
    compensation = (t_val - sum) - y_val;
    sum = t_val;
  }

  constexpr void reset() noexcept {
    sum = static_cast<T>(0);
    compensation = static_cast<T>(0);
  }

  [[nodiscard]] constexpr T value() const noexcept { return sum; }

  constexpr operator T() const noexcept { return sum; }
};

} // namespace correlation
