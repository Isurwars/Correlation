/**
 * @file SIMDScattering.hpp
 * @brief SIMD-accelerated dot product, Debye sum, complex exponential, and Miller phase sum
 * kernels.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/simd/SIMDTypes.hpp"
#include "math/simd/detail/AVX2Kernels.hpp"
#include "math/simd/detail/AVX512Kernels.hpp" // IWYU pragma: export
#include "math/simd/detail/ScalarKernels.hpp" // IWYU pragma: export

#include <cmath>

namespace correlation::math {

/**
 * @brief Computes the dot product of two array buffers.
 * @tparam T Floating-point precision (float or double).
 */
template <typename T>
inline T simdDot(const T *CORRELATION_RESTRICT input_a, const T *CORRELATION_RESTRICT input_b,
                 std::size_t count) noexcept {
#ifdef CORRELATION_SIMD_AVX512
  return detail::avx512::simdDot(input_a, input_b, count);
#elifdef CORRELATION_SIMD_AVX2
  return detail::avx2::simdDot(input_a, input_b, count);
#else
  return detail::scalar::simdDot(input_a, input_b, count);
#endif
}

/**
 * @brief Computes dot products between a reference vector and a block of target vectors.
 * @tparam T Floating-point precision (float or double).
 */
template <typename T>
inline void dotBlock(T v1x, T v1y, T v1z, const T *CORRELATION_RESTRICT v2x,
                     const T *CORRELATION_RESTRICT v2y, const T *CORRELATION_RESTRICT v2z,
                     T *CORRELATION_RESTRICT out_dot, std::size_t count) noexcept {
#ifdef CORRELATION_SIMD_AVX512
  detail::avx512::dotBlock(v1x, v1y, v1z, v2x, v2y, v2z, out_dot, count);
#elifdef CORRELATION_SIMD_AVX2
  detail::avx2::dotBlock(v1x, v1y, v1z, v2x, v2y, v2z, out_dot, count);
#else
  detail::scalar::dotBlock(v1x, v1y, v1z, v2x, v2y, v2z, out_dot, count);
#endif
}

/**
 * @brief Computes the Debye scattering equation sum across distances.
 * @tparam T Floating-point precision (float or double).
 */
template <typename T>
inline T debyeSum(T q_magnitude, const T *CORRELATION_RESTRICT distances,
                  T *CORRELATION_RESTRICT scratch, std::size_t count) noexcept {
#ifdef CORRELATION_SIMD_AVX512
  return detail::avx512::debyeSum(q_magnitude, distances, scratch, count);
#elifdef CORRELATION_SIMD_AVX2
  return detail::avx2::debyeSum(q_magnitude, distances, scratch, count);
#else
  return detail::scalar::debyeSum(q_magnitude, distances, scratch, count);
#endif
}

/**
 * @brief Computes the sum of complex exponentials for a query q-vector.
 * @tparam T Coordinate scalar type (float or double).
 */
template <typename T>
inline void complexExpSum(T q_x, T q_y, T q_z, const T *CORRELATION_RESTRICT x_s,
                          const T *CORRELATION_RESTRICT y_s, const T *CORRELATION_RESTRICT z_s,
                          std::size_t count, ComplexExpSumResult<T> &result) noexcept {
  result.cos_sum = static_cast<T>(0.0);
  result.sin_sum = static_cast<T>(0.0);
  for (std::size_t idx = 0; idx < count; ++idx) {
    const T phase = q_x * x_s[idx] + q_y * y_s[idx] + q_z * z_s[idx];
    result.cos_sum += std::cos(phase);
    result.sin_sum += std::sin(phase);
  }
}

template <typename T>
inline void complexExpSum(const ComplexExpSumParams<T> &params,
                          ComplexExpSumResult<T> &result) noexcept {
  complexExpSum(params.q_x, params.q_y, params.q_z, params.x_s, params.y_s, params.z_s,
                params.count, result);
}

template <typename T = real_t>
inline ComplexExpSumResult<T> complexExpSum(T q_x, T q_y, T q_z, const T *CORRELATION_RESTRICT x_s,
                                            const T *CORRELATION_RESTRICT y_s,
                                            const T *CORRELATION_RESTRICT z_s,
                                            std::size_t count) noexcept {
  ComplexExpSumResult<T> result;
  complexExpSum(q_x, q_y, q_z, x_s, y_s, z_s, count, result);
  return result;
}

template <typename T = real_t>
inline ComplexExpSumResult<T> complexExpSum(const MillerExpSumParams<T> &params) noexcept {
  ComplexExpSumResult<T> result;
  complexExpSum(params, result);
  return result;
}

/**
 * @brief Computes the Miller phase sum across angle blocks.
 * @tparam T Floating-point precision (float or double).
 */
template <typename T>
inline void millerPhaseSum(const MillerPhaseSumParams<T> &params,
                           MillerPhaseSumResult<T> &result) noexcept {
#ifdef CORRELATION_SIMD_AVX512
  detail::avx512::millerPhaseSum(params, result);
#elifdef CORRELATION_SIMD_AVX2
  detail::avx2::millerPhaseSum(params, result);
#else
  detail::scalar::millerPhaseSum(params, result);
#endif
}

template <typename T = real_t>
inline MillerPhaseSumResult<T> millerPhaseSum(const MillerPhaseSumParams<T> &params) noexcept {
  MillerPhaseSumResult<T> result;
  millerPhaseSum(params, result);
  return result;
}

template <typename T>
inline void millerPhaseSum(const T *CORRELATION_RESTRICT cos1, const T *CORRELATION_RESTRICT sin1,
                           const T *CORRELATION_RESTRICT cos2, const T *CORRELATION_RESTRICT sin2,
                           const T *CORRELATION_RESTRICT cos3, const T *CORRELATION_RESTRICT sin3,
                           std::size_t count, MillerPhaseSumResult<T> &result) noexcept {
  millerPhaseSum(
      MillerPhaseSumParams<T>{
          .cos1 = cos1,
          .sin1 = sin1,
          .cos2 = cos2,
          .sin2 = sin2,
          .cos3 = cos3,
          .sin3 = sin3,
          .count = count,
      },
      result);
}

template <typename T = real_t>
inline MillerPhaseSumResult<T>
millerPhaseSum(const T *CORRELATION_RESTRICT cos1, const T *CORRELATION_RESTRICT sin1,
               const T *CORRELATION_RESTRICT cos2, const T *CORRELATION_RESTRICT sin2,
               const T *CORRELATION_RESTRICT cos3, const T *CORRELATION_RESTRICT sin3,
               std::size_t count) noexcept {
  MillerPhaseSumResult<T> result;
  millerPhaseSum(cos1, sin1, cos2, sin2, cos3, sin3, count, result);
  return result;
}

} // namespace correlation::math
