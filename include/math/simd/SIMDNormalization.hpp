/**
 * @file SIMDNormalization.hpp
 * @brief SIMD-accelerated integral, RDF normalization, and bin scaling routines.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/simd/SIMDTypes.hpp"
#include "math/simd/detail/AVX2Kernels.hpp"
#include "math/simd/detail/AVX512Kernels.hpp" // IWYU pragma: export
#include "math/simd/detail/ScalarKernels.hpp" // IWYU pragma: export

namespace correlation::math {

/**
 * @brief Computes a sinc-weighted integral over a range.
 * @tparam T Floating-point precision (float or double).
 */
template <typename T> inline T sincIntegral(const SincIntegralParams<T> &params) noexcept {
#ifdef CORRELATION_SIMD_AVX512
  return detail::avx512::sincIntegral(params);
#elifdef CORRELATION_SIMD_AVX2
  return detail::avx2::sincIntegral(params);
#else
  return detail::scalar::sincIntegral(params);
#endif
}

/**
 * @brief Computes a sinc-weighted integral for double precision (backward compatibility).
 */
inline double sincIntegral(double q_magnitude, const double *CORRELATION_RESTRICT integrand,
                           const double *CORRELATION_RESTRICT radial_bins,
                           double *CORRELATION_RESTRICT sinqr_scratch, std::size_t count) noexcept {
  return sincIntegral(SincIntegralParams<double>{
      .q_magnitude = q_magnitude,
      .integrand = integrand,
      .radial_bins = radial_bins,
      .sinqr_scratch = sinqr_scratch,
      .count = count,
  });
}

/**
 * @brief Computes a sinc-weighted integral for single precision (backward compatibility).
 */
inline float sincIntegral(float q_magnitude, const float *CORRELATION_RESTRICT integrand,
                          const float *CORRELATION_RESTRICT radial_bins,
                          float *CORRELATION_RESTRICT sinqr_scratch, std::size_t count) noexcept {
  return sincIntegral(SincIntegralParams<float>{
      .q_magnitude = q_magnitude,
      .integrand = integrand,
      .radial_bins = radial_bins,
      .sinqr_scratch = sinqr_scratch,
      .count = count,
  });
}

/**
 * @brief Normalizes Radial Distribution Function (RDF) bins.
 * @tparam T Floating-point precision (float or double).
 */
template <typename T>
inline void normalizeRdfBins(const RDFNormalizationParams<T> &params) noexcept {
#ifdef CORRELATION_SIMD_AVX512
  detail::avx512::normalizeRdfBins(params);
#elifdef CORRELATION_SIMD_AVX2
  detail::avx2::normalizeRdfBins(params);
#else
  detail::scalar::normalizeRdfBins(params);
#endif
}

/**
 * @brief Scales an array by a constant factor in-place using ScaleBinsParams container.
 * @tparam T Floating-point precision (float or double).
 */
template <typename T> inline void scaleBins(const ScaleBinsParams<T> &params) noexcept {
#ifdef CORRELATION_SIMD_AVX512
  detail::avx512::scaleBins(params);
#elifdef CORRELATION_SIMD_AVX2
  detail::avx2::scaleBins(params);
#else
  detail::scalar::scaleBins(params);
#endif
}

/**
 * @brief Scales an array by a constant factor in-place.
 * @tparam T Floating-point precision (float or double).
 */
template <typename T> inline void scaleBins(T *arr, T scale_factor, std::size_t count) noexcept {
  scaleBins(ScaleBinsParams<T>{
      .arr = arr,
      .scale_factor = scale_factor,
      .count = count,
  });
}

} // namespace correlation::math
