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
template <typename T> inline T sinc_integral(const SincIntegralParams<T> &params) noexcept {
#ifdef CORRELATION_SIMD_AVX512
  return detail::avx512::sinc_integral(params);
#elif defined(CORRELATION_SIMD_AVX2)
  return detail::avx2::sinc_integral(params);
#else
  return detail::scalar::sinc_integral(params);
#endif
}

/**
 * @brief Computes a sinc-weighted integral for double precision (backward compatibility).
 */
inline double sinc_integral(double q_magnitude, const double *CORRELATION_RESTRICT integrand,
                            const double *CORRELATION_RESTRICT radial_bins, double *CORRELATION_RESTRICT sinqr_scratch,
                            std::size_t count) noexcept {
  return sinc_integral(SincIntegralParams<double>{
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
inline float sinc_integral(float q_magnitude, const float *CORRELATION_RESTRICT integrand,
                           const float *CORRELATION_RESTRICT radial_bins, float *CORRELATION_RESTRICT sinqr_scratch,
                           std::size_t count) noexcept {
  return sinc_integral(SincIntegralParams<float>{
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
template <typename T> inline void normalize_rdf_bins(const RDFNormalizationParams<T> &params) noexcept {
#ifdef CORRELATION_SIMD_AVX512
  detail::avx512::normalize_rdf_bins(params);
#elif defined(CORRELATION_SIMD_AVX2)
  detail::avx2::normalize_rdf_bins(params);
#else
  detail::scalar::normalize_rdf_bins(params);
#endif
}

/**
 * @brief Scales an array by a constant factor in-place using ScaleBinsParams container.
 * @tparam T Floating-point precision (float or double).
 */
template <typename T> inline void scale_bins(const ScaleBinsParams<T> &params) noexcept {
#ifdef CORRELATION_SIMD_AVX512
  detail::avx512::scale_bins(params);
#elif defined(CORRELATION_SIMD_AVX2)
  detail::avx2::scale_bins(params);
#else
  detail::scalar::scale_bins(params);
#endif
}

/**
 * @brief Scales an array by a constant factor in-place.
 * @tparam T Floating-point precision (float or double).
 */
template <typename T> inline void scale_bins(T *arr, T scale_factor, std::size_t count) noexcept {
  scale_bins(ScaleBinsParams<T>{
      .arr = arr,
      .scale_factor = scale_factor,
      .count = count,
  });
}

} // namespace correlation::math
