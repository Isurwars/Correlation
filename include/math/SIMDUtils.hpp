/**
 * @file SIMDUtils.hpp
 * @brief SIMD-accelerated kernels for distance, integration, and normalization.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/Precision.hpp"
#include "math/SIMDConfig.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

namespace correlation::math {

// ---------------------------------------------------------------------------
// SoA block of atom positions
// ---------------------------------------------------------------------------
/**
 * @struct PositionBlockT
 * @brief Represents a Structure of Arrays (SoA) block of atom positions for SIMD processing.
 * @tparam T Floating point coordinate type (float or double).
 */
template <typename T = real_t> struct PositionBlockT {
  T *x{nullptr};        ///< Array of x-coordinates.
  T *y{nullptr};        ///< Array of y-coordinates.
  T *z{nullptr};        ///< Array of z-coordinates.
  std::size_t count{0}; ///< Number of atoms in this block.
};

/// Backward-compatible type alias for default real_t precision PositionBlock.
using PositionBlock = PositionBlockT<real_t>;

// ---------------------------------------------------------------------------
// Scalar helper
// ---------------------------------------------------------------------------
/**
 * @struct Point3T
 * @brief Lightweight POD for passing 3D coordinates to scalar SIMD helpers.
 * @tparam T Floating point coordinate type (float or double).
 */
template <typename T = real_t> struct Point3T {
  T x{static_cast<T>(0.0)}; ///< X-coordinate.
  T y{static_cast<T>(0.0)}; ///< Y-coordinate.
  T z{static_cast<T>(0.0)}; ///< Z-coordinate.
};

/// Backward-compatible type alias for default real_t precision Point3.
using Point3 = Point3T<real_t>;

/**
 * @brief Computes the squared distance between two 3D points.
 * @tparam T Floating point coordinate type.
 * @param[in] point_a The first point.
 * @param[in] point_b The second point.
 * @return The scalar squared distance.
 */
template <typename T> [[nodiscard]] inline T dist_sq_scalar(Point3T<T> point_a, Point3T<T> point_b) noexcept {
  const T dist_x = point_b.x - point_a.x;
  const T dist_y = point_b.y - point_a.y;
  const T dist_z = point_b.z - point_a.z;
  return dist_x * dist_x + dist_y * dist_y + dist_z * dist_z;
}

/**
 * @brief Non-templated overload of dist_sq_scalar for backward compatibility.
 * @param[in] point_a The first point.
 * @param[in] point_b The second point.
 * @return The scalar squared distance using real_t.
 */
[[nodiscard]] inline real_t dist_sq_scalar(Point3 point_a, Point3 point_b) noexcept {
  return dist_sq_scalar<real_t>(point_a, point_b);
}

// ---------------------------------------------------------------------------
// Core SIMD kernels
// ---------------------------------------------------------------------------

#ifdef CORRELATION_SIMD_AVX512

/**
 * @brief Computes the squared distances from a reference point to a single-precision block of positions (AVX-512).
 * @param[in] ref_x x-coordinate of reference point.
 * @param[in] ref_y y-coordinate of reference point.
 * @param[in] ref_z z-coordinate of reference point.
 * @param[in] block Structure of Arrays containing target positions.
 * @param[out] out_dsq Destination array for computed squared distances.
 */
inline void compute_dsq_block(float ref_x, float ref_y, float ref_z, const PositionBlockT<float> &block,
                              float *CORRELATION_RESTRICT out_dsq) noexcept {
  const __m512 va_x = _mm512_set1_ps(ref_x);
  const __m512 va_y = _mm512_set1_ps(ref_y);
  const __m512 va_z = _mm512_set1_ps(ref_z);

  std::size_t idx = 0;
  for (; idx + 16 <= block.count; idx += 16) {
    const __m512 dx = _mm512_sub_ps(_mm512_loadu_ps(block.x + idx), va_x);
    const __m512 dy = _mm512_sub_ps(_mm512_loadu_ps(block.y + idx), va_y);
    const __m512 dz = _mm512_sub_ps(_mm512_loadu_ps(block.z + idx), va_z);
    const __m512 dsq = _mm512_fmadd_ps(dx, dx, _mm512_fmadd_ps(dy, dy, _mm512_mul_ps(dz, dz)));
    _mm512_storeu_ps(out_dsq + idx, dsq);
  }
  if (idx < block.count) {
    const auto mask = static_cast<__mmask16>((1U << static_cast<std::uint32_t>(block.count - idx)) - 1U);
    const __m512 dx = _mm512_sub_ps(_mm512_maskz_loadu_ps(mask, block.x + idx), va_x);
    const __m512 dy = _mm512_sub_ps(_mm512_maskz_loadu_ps(mask, block.y + idx), va_y);
    const __m512 dz = _mm512_sub_ps(_mm512_maskz_loadu_ps(mask, block.z + idx), va_z);
    const __m512 dsq = _mm512_fmadd_ps(dx, dx, _mm512_fmadd_ps(dy, dy, _mm512_mul_ps(dz, dz)));
    _mm512_mask_storeu_ps(out_dsq + idx, mask, dsq);
  }
}

/**
 * @brief Computes the squared distances from a reference point to a double-precision block of positions (AVX-512).
 * @param[in] ref_x x-coordinate of reference point.
 * @param[in] ref_y y-coordinate of reference point.
 * @param[in] ref_z z-coordinate of reference point.
 * @param[in] block Structure of Arrays containing target positions.
 * @param[out] out_dsq Destination array for computed squared distances.
 */
inline void compute_dsq_block(double ref_x, double ref_y, double ref_z, const PositionBlockT<double> &block,
                              double *CORRELATION_RESTRICT out_dsq) noexcept {
  const __m512d va_x = _mm512_set1_pd(ref_x);
  const __m512d va_y = _mm512_set1_pd(ref_y);
  const __m512d va_z = _mm512_set1_pd(ref_z);

  std::size_t idx = 0;
  for (; idx + 8 <= block.count; idx += 8) {
    const __m512d dx = _mm512_sub_pd(_mm512_loadu_pd(block.x + idx), va_x);
    const __m512d dy = _mm512_sub_pd(_mm512_loadu_pd(block.y + idx), va_y);
    const __m512d dz = _mm512_sub_pd(_mm512_loadu_pd(block.z + idx), va_z);
    const __m512d dsq = _mm512_fmadd_pd(dx, dx, _mm512_fmadd_pd(dy, dy, _mm512_mul_pd(dz, dz)));
    _mm512_storeu_pd(out_dsq + idx, dsq);
  }
  if (idx < block.count) {
    const auto mask = static_cast<__mmask8>((1U << static_cast<std::uint32_t>(block.count - idx)) - 1U);
    const __m512d dx = _mm512_sub_pd(_mm512_maskz_loadu_pd(mask, block.x + idx), va_x);
    const __m512d dy = _mm512_sub_pd(_mm512_maskz_loadu_pd(mask, block.y + idx), va_y);
    const __m512d dz = _mm512_sub_pd(_mm512_maskz_loadu_pd(mask, block.z + idx), va_z);
    const __m512d dsq = _mm512_fmadd_pd(dx, dx, _mm512_fmadd_pd(dy, dy, _mm512_mul_pd(dz, dz)));
    _mm512_mask_storeu_pd(out_dsq + idx, mask, dsq);
  }
}

/**
 * @brief Computes the dot product of two double arrays (AVX-512 version).
 * @param[in] input_a First array.
 * @param[in] input_b Second array.
 * @param[in] count Number of elements.
 * @return The dot product value.
 */
inline double simd_dot(const double *CORRELATION_RESTRICT input_a, const double *CORRELATION_RESTRICT input_b,
                       std::size_t count) noexcept {
  __m512d vacc = _mm512_setzero_pd();
  std::size_t idx = 0;
  for (; idx + 8 <= count; idx += 8) {
    const __m512d va = _mm512_loadu_pd(input_a + idx);
    const __m512d vb = _mm512_loadu_pd(input_b + idx);
    vacc = _mm512_fmadd_pd(va, vb, vacc);
  }
  double acc = _mm512_reduce_add_pd(vacc);
  double carry = 0.0;
  for (; idx < count; ++idx) {
    const double val = input_a[idx] * input_b[idx];
    const double y_val = val - carry;
    const double temp_sum = acc + y_val;
    carry = (temp_sum - acc) - y_val;
    acc = temp_sum;
  }
  return acc;
}

/**
 * @brief Computes the dot product of two float arrays (AVX-512 version).
 * @param[in] input_a First array.
 * @param[in] input_b Second array.
 * @param[in] count Number of elements.
 * @return The dot product value.
 */
inline float simd_dot(const float *CORRELATION_RESTRICT input_a, const float *CORRELATION_RESTRICT input_b,
                      std::size_t count) noexcept {
  __m512 vacc = _mm512_setzero_ps();
  std::size_t idx = 0;
  for (; idx + 16 <= count; idx += 16) {
    const __m512 va = _mm512_loadu_ps(input_a + idx);
    const __m512 vb = _mm512_loadu_ps(input_b + idx);
    vacc = _mm512_fmadd_ps(va, vb, vacc);
  }
  float acc = _mm512_reduce_add_ps(vacc);
  for (; idx < count; ++idx) {
    acc += input_a[idx] * input_b[idx];
  }
  return acc;
}

/**
 * @brief Computes a sinc-weighted integral over a double range (AVX-512 version).
 * @param[in] q_magnitude The scattering vector magnitude.
 * @param[in] integrand Values to be integrated.
 * @param[in] radial_bins Radial distances for each point.
 * @param[out] sinqr_scratch Scratchpad memory for sinc calculations.
 * @param[in] count Number of data points.
 * @return The accumulated integral value.
 */
inline double sinc_integral(double q_magnitude, const double *CORRELATION_RESTRICT integrand,
                            const double *CORRELATION_RESTRICT radial_bins, double *CORRELATION_RESTRICT sinqr_scratch,
                            std::size_t count) noexcept {
  for (std::size_t idx = 0; idx < count; ++idx) {
    sinqr_scratch[idx] = std::sin(q_magnitude * radial_bins[idx]);
  }
  return simd_dot(integrand, sinqr_scratch, count);
}

/**
 * @brief Computes a sinc-weighted integral over a float range (AVX-512 version).
 * @param[in] q_magnitude The scattering vector magnitude.
 * @param[in] integrand Values to be integrated.
 * @param[in] radial_bins Radial distances for each point.
 * @param[out] sinqr_scratch Scratchpad memory for sinc calculations.
 * @param[in] count Number of data points.
 * @return The accumulated integral value.
 */
inline float sinc_integral(float q_magnitude, const float *CORRELATION_RESTRICT integrand,
                           const float *CORRELATION_RESTRICT radial_bins, float *CORRELATION_RESTRICT sinqr_scratch,
                           std::size_t count) noexcept {
  for (std::size_t idx = 0; idx < count; ++idx) {
    sinqr_scratch[idx] = std::sin(q_magnitude * radial_bins[idx]);
  }
  return simd_dot(integrand, sinqr_scratch, count);
}

/**
 * @brief Computes the Debye scattering equation sum for double precision (AVX-512 version).
 * @param[in] q_magnitude The scattering vector magnitude.
 * @param[in] distances Source array of distance values.
 * @param[out] scratch Scratchpad array for intermediate sinc evaluations.
 * @param[in] count Number of distances to process.
 * @return The computed Debye sum.
 */
inline double debye_sum(double q_magnitude, const double *CORRELATION_RESTRICT distances,
                        double *CORRELATION_RESTRICT scratch, std::size_t count) noexcept {
  if (q_magnitude < 1.0e-9) {
    return static_cast<double>(count);
  }

  for (std::size_t idx = 0; idx < count; ++idx) {
    const double val_x = q_magnitude * distances[idx];
    scratch[idx] = (val_x < 1.0e-4) ? (1.0 - (val_x * val_x) / 6.0) : (std::sin(val_x) / val_x);
  }

  __m512d vacc = _mm512_setzero_pd();
  std::size_t idx = 0;
  for (; idx + 8 <= count; idx += 8) {
    vacc = _mm512_add_pd(vacc, _mm512_loadu_pd(scratch + idx));
  }
  double acc = _mm512_reduce_add_pd(vacc);
  double carry = 0.0;
  for (; idx < count; ++idx) {
    const double val = scratch[idx];
    const double y_val = val - carry;
    const double temp_sum = acc + y_val;
    carry = (temp_sum - acc) - y_val;
    acc = temp_sum;
  }
  return acc;
}

/**
 * @brief Computes the Debye scattering equation sum for single precision (AVX-512 version).
 * @param[in] q_magnitude The scattering vector magnitude.
 * @param[in] distances Source array of distance values.
 * @param[out] scratch Scratchpad array for intermediate sinc evaluations.
 * @param[in] count Number of distances to process.
 * @return The computed Debye sum.
 */
inline float debye_sum(float q_magnitude, const float *CORRELATION_RESTRICT distances,
                       float *CORRELATION_RESTRICT scratch, std::size_t count) noexcept {
  if (q_magnitude < 1.0e-9F) {
    return static_cast<float>(count);
  }

  for (std::size_t idx = 0; idx < count; ++idx) {
    const float val_x = q_magnitude * distances[idx];
    scratch[idx] = (val_x < 1.0e-4F) ? (1.0F - (val_x * val_x) / 6.0F) : (std::sin(val_x) / val_x);
  }

  __m512 vacc = _mm512_setzero_ps();
  std::size_t idx = 0;
  for (; idx + 16 <= count; idx += 16) {
    vacc = _mm512_add_ps(vacc, _mm512_loadu_ps(scratch + idx));
  }
  float acc = _mm512_reduce_add_ps(vacc);
  for (; idx < count; ++idx) {
    acc += scratch[idx];
  }
  return acc;
}

/**
 * @brief Normalizes Radial Distribution Function (RDF) bins for double precision (AVX-512 version).
 * @param[in] hist_data Histogram array of raw bin counts.
 * @param[in] radial_bins Array of radial bin positions.
 * @param[in] g_norm Normalization factor for g(r).
 * @param[in] inv_Ni_dr Inverse scale factor times dr for species i.
 * @param[in] inv_Nj_dr Inverse scale factor times dr for species j.
 * @param[in] pi4_rho_j Density scaling factor (4 * pi * rho_j).
 * @param[out] g_out Output array for g(r).
 * @param[out] G_out Output array for G(r).
 * @param[out] J_out Output array for J(r).
 * @param[out] Jinv_out Output array for J^{-1}(r).
 * @param[in] count Total number of bins.
 */
inline void normalize_rdf_bins(const double *CORRELATION_RESTRICT hist_data,
                               const double *CORRELATION_RESTRICT radial_bins, double g_norm, double inv_Ni_dr,
                               double inv_Nj_dr, double pi4_rho_j, double *CORRELATION_RESTRICT g_out,
                               double *CORRELATION_RESTRICT G_out, double *CORRELATION_RESTRICT J_out,
                               double *CORRELATION_RESTRICT Jinv_out, std::size_t count) noexcept {
  const __m512d vg_norm = _mm512_set1_pd(g_norm);
  const __m512d v1 = _mm512_set1_pd(1.0);
  const __m512d vinNidr = _mm512_set1_pd(inv_Ni_dr);
  const __m512d vinNjdr = _mm512_set1_pd(inv_Nj_dr);
  const __m512d vpi4rho = _mm512_set1_pd(pi4_rho_j);

  std::size_t idx = 1;
  for (; idx + 8 <= count; idx += 8) {
    const __m512d vH = _mm512_loadu_pd(hist_data + idx);
    const __m512d vr = _mm512_loadu_pd(radial_bins + idx);
    const __m512d vr2 = _mm512_mul_pd(vr, vr);
    const __m512d vg = _mm512_div_pd(_mm512_mul_pd(vH, vg_norm), vr2);
    _mm512_storeu_pd(g_out + idx, vg);
    _mm512_storeu_pd(G_out + idx, _mm512_mul_pd(vpi4rho, _mm512_mul_pd(vr, _mm512_sub_pd(vg, v1))));
    _mm512_storeu_pd(J_out + idx, _mm512_mul_pd(vH, vinNidr));
    _mm512_storeu_pd(Jinv_out + idx, _mm512_mul_pd(vH, vinNjdr));
  }
  for (; idx < count; ++idx) {
    const double r_val = radial_bins[idx];
    if (r_val < 1.0e-9) {
      continue;
    }
    const double g_val = hist_data[idx] * g_norm / (r_val * r_val);
    g_out[idx] = g_val;
    G_out[idx] = pi4_rho_j * r_val * (g_val - 1.0);
    J_out[idx] = hist_data[idx] * inv_Ni_dr;
    Jinv_out[idx] = hist_data[idx] * inv_Nj_dr;
  }
}

/**
 * @brief Normalizes Radial Distribution Function (RDF) bins for single precision (AVX-512 version).
 * @param[in] hist_data Histogram array of raw bin counts.
 * @param[in] radial_bins Array of radial bin positions.
 * @param[in] g_norm Normalization factor for g(r).
 * @param[in] inv_Ni_dr Inverse scale factor times dr for species i.
 * @param[in] inv_Nj_dr Inverse scale factor times dr for species j.
 * @param[in] pi4_rho_j Density scaling factor (4 * pi * rho_j).
 * @param[out] g_out Output array for g(r).
 * @param[out] G_out Output array for G(r).
 * @param[out] J_out Output array for J(r).
 * @param[out] Jinv_out Output array for J^{-1}(r).
 * @param[in] count Total number of bins.
 */
inline void normalize_rdf_bins(const float *CORRELATION_RESTRICT hist_data,
                               const float *CORRELATION_RESTRICT radial_bins, float g_norm, float inv_Ni_dr,
                               float inv_Nj_dr, float pi4_rho_j, float *CORRELATION_RESTRICT g_out,
                               float *CORRELATION_RESTRICT G_out, float *CORRELATION_RESTRICT J_out,
                               float *CORRELATION_RESTRICT Jinv_out, std::size_t count) noexcept {
  const __m512 vg_norm = _mm512_set1_ps(g_norm);
  const __m512 v1 = _mm512_set1_ps(1.0F);
  const __m512 vinNidr = _mm512_set1_ps(inv_Ni_dr);
  const __m512 vinNjdr = _mm512_set1_ps(inv_Nj_dr);
  const __m512 vpi4rho = _mm512_set1_ps(pi4_rho_j);

  std::size_t idx = 1;
  for (; idx + 16 <= count; idx += 16) {
    const __m512 vH = _mm512_loadu_ps(hist_data + idx);
    const __m512 vr = _mm512_loadu_ps(radial_bins + idx);
    const __m512 vr2 = _mm512_mul_ps(vr, vr);
    const __m512 vg = _mm512_div_ps(_mm512_mul_ps(vH, vg_norm), vr2);
    _mm512_storeu_ps(g_out + idx, vg);
    _mm512_storeu_ps(G_out + idx, _mm512_mul_ps(vpi4rho, _mm512_mul_ps(vr, _mm512_sub_ps(vg, v1))));
    _mm512_storeu_ps(J_out + idx, _mm512_mul_ps(vH, vinNidr));
    _mm512_storeu_ps(Jinv_out + idx, _mm512_mul_ps(vH, vinNjdr));
  }
  for (; idx < count; ++idx) {
    const float r_val = radial_bins[idx];
    if (r_val < 1.0e-9F) {
      continue;
    }
    const float g_val = hist_data[idx] * g_norm / (r_val * r_val);
    g_out[idx] = g_val;
    G_out[idx] = pi4_rho_j * r_val * (g_val - 1.0F);
    J_out[idx] = hist_data[idx] * inv_Ni_dr;
    Jinv_out[idx] = hist_data[idx] * inv_Nj_dr;
  }
}

/**
 * @brief Scales a double array by a constant factor (AVX-512 version).
 * @param[in,out] arr Array to scale in-place.
 * @param[in] scale_factor Scalar multiplier.
 * @param[in] count Number of elements.
 */
inline void scale_bins(double *arr, double scale_factor, std::size_t count) noexcept {
  const __m512d vs = _mm512_set1_pd(scale_factor);
  std::size_t idx = 0;
  for (; idx + 8 <= count; idx += 8) {
    _mm512_storeu_pd(arr + idx, _mm512_mul_pd(_mm512_loadu_pd(arr + idx), vs));
  }
  for (; idx < count; ++idx) {
    arr[idx] *= scale_factor;
  }
}

/**
 * @brief Scales a float array by a constant factor (AVX-512 version).
 * @param[in,out] arr Array to scale in-place.
 * @param[in] scale_factor Scalar multiplier.
 * @param[in] count Number of elements.
 */
inline void scale_bins(float *arr, float scale_factor, std::size_t count) noexcept {
  const __m512 vs = _mm512_set1_ps(scale_factor);
  std::size_t idx = 0;
  for (; idx + 16 <= count; idx += 16) {
    _mm512_storeu_ps(arr + idx, _mm512_mul_ps(_mm512_loadu_ps(arr + idx), vs));
  }
  for (; idx < count; ++idx) {
    arr[idx] *= scale_factor;
  }
}

#elif defined(CORRELATION_SIMD_AVX2)

/**
 * @brief Computes the squared distances from a reference point to a single-precision block of positions (AVX2).
 * @param[in] ref_x x-coordinate of reference point.
 * @param[in] ref_y y-coordinate of reference point.
 * @param[in] ref_z z-coordinate of reference point.
 * @param[in] block Structure of Arrays containing target positions.
 * @param[out] out_dsq Destination array for computed squared distances.
 */
inline void compute_dsq_block(float ref_x, float ref_y, float ref_z, const PositionBlockT<float> &block,
                              float *CORRELATION_RESTRICT out_dsq) noexcept {
  const __m256 va_x = _mm256_set1_ps(ref_x);
  const __m256 va_y = _mm256_set1_ps(ref_y);
  const __m256 va_z = _mm256_set1_ps(ref_z);

  std::size_t idx = 0;
  for (; idx + 8 <= block.count; idx += 8) {
    const __m256 d_x = _mm256_sub_ps(_mm256_loadu_ps(block.x + idx), va_x);
    const __m256 d_y = _mm256_sub_ps(_mm256_loadu_ps(block.y + idx), va_y);
    const __m256 d_z = _mm256_sub_ps(_mm256_loadu_ps(block.z + idx), va_z);
#ifdef __FMA__
    const __m256 dsq = _mm256_fmadd_ps(d_x, d_x, _mm256_fmadd_ps(d_y, d_y, _mm256_mul_ps(d_z, d_z)));
#else
    const __m256 dsq =
        _mm256_add_ps(_mm256_mul_ps(d_x, d_x), _mm256_add_ps(_mm256_mul_ps(d_y, d_y), _mm256_mul_ps(d_z, d_z)));
#endif
    _mm256_storeu_ps(out_dsq + idx, dsq);
  }
  for (; idx < block.count; ++idx) {
    out_dsq[idx] = dist_sq_scalar<float>({.x = ref_x, .y = ref_y, .z = ref_z},
                                         {.x = block.x[idx], .y = block.y[idx], .z = block.z[idx]});
  }
}

/**
 * @brief Computes the squared distances from a reference point to a double-precision block of positions (AVX2).
 * @param[in] ref_x x-coordinate of reference point.
 * @param[in] ref_y y-coordinate of reference point.
 * @param[in] ref_z z-coordinate of reference point.
 * @param[in] block Structure of Arrays containing target positions.
 * @param[out] out_dsq Destination array for computed squared distances.
 */
inline void compute_dsq_block(double ref_x, double ref_y, double ref_z, const PositionBlockT<double> &block,
                              double *CORRELATION_RESTRICT out_dsq) noexcept {
  const __m256d va_x = _mm256_set1_pd(ref_x);
  const __m256d va_y = _mm256_set1_pd(ref_y);
  const __m256d va_z = _mm256_set1_pd(ref_z);

  std::size_t idx = 0;
  for (; idx + 4 <= block.count; idx += 4) {
    const __m256d d_x = _mm256_sub_pd(_mm256_loadu_pd(block.x + idx), va_x);
    const __m256d d_y = _mm256_sub_pd(_mm256_loadu_pd(block.y + idx), va_y);
    const __m256d d_z = _mm256_sub_pd(_mm256_loadu_pd(block.z + idx), va_z);
#ifdef __FMA__
    const __m256d dsq = _mm256_fmadd_pd(d_x, d_x, _mm256_fmadd_pd(d_y, d_y, _mm256_mul_pd(d_z, d_z)));
#else
    const __m256d dsq =
        _mm256_add_pd(_mm256_mul_pd(d_x, d_x), _mm256_add_pd(_mm256_mul_pd(d_y, d_y), _mm256_mul_pd(d_z, d_z)));
#endif
    _mm256_storeu_pd(out_dsq + idx, dsq);
  }
  for (; idx < block.count; ++idx) {
    out_dsq[idx] = dist_sq_scalar<double>({.x = ref_x, .y = ref_y, .z = ref_z},
                                          {.x = block.x[idx], .y = block.y[idx], .z = block.z[idx]});
  }
}

/**
 * @brief Computes the dot product of two double arrays (AVX2 version).
 * @param[in] input_a First array.
 * @param[in] input_b Second array.
 * @param[in] count Number of elements.
 * @return The dot product value.
 */
inline double simd_dot(const double *CORRELATION_RESTRICT input_a, const double *CORRELATION_RESTRICT input_b,
                       std::size_t count) noexcept {
  __m256d vacc = _mm256_setzero_pd();
  std::size_t idx = 0;
  for (; idx + 4 <= count; idx += 4) {
    const __m256d vec_a = _mm256_loadu_pd(input_a + idx);
    const __m256d vec_b = _mm256_loadu_pd(input_b + idx);
#ifdef __FMA__
    vacc = _mm256_fmadd_pd(vec_a, vec_b, vacc);
#else
    vacc = _mm256_add_pd(vacc, _mm256_mul_pd(vec_a, vec_b));
#endif
  }
  const __m128d low = _mm256_castpd256_pd128(vacc);
  const __m128d high = _mm256_extractf128_pd(vacc, 1);
  const __m128d sum2 = _mm_add_pd(low, high);
  const __m128d sum1 = _mm_hadd_pd(sum2, sum2);
  double acc = _mm_cvtsd_f64(sum1);
  double carry = 0.0;
  for (; idx < count; ++idx) {
    const double val = input_a[idx] * input_b[idx];
    const double y_val = val - carry;
    const double temp_sum = acc + y_val;
    carry = (temp_sum - acc) - y_val;
    acc = temp_sum;
  }
  return acc;
}

/**
 * @brief Computes the dot product of two float arrays (AVX2 version).
 * @param[in] input_a First array.
 * @param[in] input_b Second array.
 * @param[in] count Number of elements.
 * @return The dot product value.
 */
inline float simd_dot(const float *CORRELATION_RESTRICT input_a, const float *CORRELATION_RESTRICT input_b,
                      std::size_t count) noexcept {
  __m256 vacc = _mm256_setzero_ps();
  std::size_t idx = 0;
  for (; idx + 8 <= count; idx += 8) {
    const __m256 vec_a = _mm256_loadu_ps(input_a + idx);
    const __m256 vec_b = _mm256_loadu_ps(input_b + idx);
#ifdef __FMA__
    vacc = _mm256_fmadd_ps(vec_a, vec_b, vacc);
#else
    vacc = _mm256_add_ps(vacc, _mm256_mul_ps(vec_a, vec_b));
#endif
  }
  alignas(32) float float_buf[8];
  _mm256_storeu_ps(float_buf, vacc);
  float acc = float_buf[0] + float_buf[1] + float_buf[2] + float_buf[3] + float_buf[4] + float_buf[5] + float_buf[6] +
              float_buf[7];
  for (; idx < count; ++idx) {
    acc += input_a[idx] * input_b[idx];
  }
  return acc;
}

/**
 * @brief Computes a sinc-weighted integral over a double range (AVX2 version).
 * @param[in] q_magnitude The scattering vector magnitude.
 * @param[in] integrand Values to be integrated.
 * @param[in] radial_bins Radial distances for each point.
 * @param[out] sinqr_scratch Scratchpad memory for sinc calculations.
 * @param[in] count Number of data points.
 * @return The accumulated integral value.
 */
inline double sinc_integral(double q_magnitude, const double *CORRELATION_RESTRICT integrand,
                            const double *CORRELATION_RESTRICT radial_bins, double *CORRELATION_RESTRICT sinqr_scratch,
                            std::size_t count) noexcept {
  for (std::size_t idx = 0; idx < count; ++idx) {
    sinqr_scratch[idx] = std::sin(q_magnitude * radial_bins[idx]);
  }
  return simd_dot(integrand, sinqr_scratch, count);
}

/**
 * @brief Computes a sinc-weighted integral over a float range (AVX2 version).
 * @param[in] q_magnitude The scattering vector magnitude.
 * @param[in] integrand Values to be integrated.
 * @param[in] radial_bins Radial distances for each point.
 * @param[out] sinqr_scratch Scratchpad memory for sinc calculations.
 * @param[in] count Number of data points.
 * @return The accumulated integral value.
 */
inline float sinc_integral(float q_magnitude, const float *CORRELATION_RESTRICT integrand,
                           const float *CORRELATION_RESTRICT radial_bins, float *CORRELATION_RESTRICT sinqr_scratch,
                           std::size_t count) noexcept {
  for (std::size_t idx = 0; idx < count; ++idx) {
    sinqr_scratch[idx] = std::sin(q_magnitude * radial_bins[idx]);
  }
  return simd_dot(integrand, sinqr_scratch, count);
}

/**
 * @brief Computes the Debye scattering equation sum for double precision (AVX2 version).
 * @param[in] q_magnitude The scattering vector magnitude.
 * @param[in] distances Source array of distance values.
 * @param[out] scratch Scratchpad array for intermediate sinc evaluations.
 * @param[in] count Number of distances to process.
 * @return The computed Debye sum.
 */
inline double debye_sum(double q_magnitude, const double *CORRELATION_RESTRICT distances,
                        double *CORRELATION_RESTRICT scratch, std::size_t count) noexcept {
  if (q_magnitude < 1.0e-9) {
    return static_cast<double>(count);
  }

  for (std::size_t idx = 0; idx < count; ++idx) {
    const double val_x = q_magnitude * distances[idx];
    scratch[idx] = (val_x < 1.0e-4) ? (1.0 - (val_x * val_x) / 6.0) : (std::sin(val_x) / val_x);
  }

  __m256d vacc = _mm256_setzero_pd();
  std::size_t idx = 0;
  for (; idx + 4 <= count; idx += 4) {
    vacc = _mm256_add_pd(vacc, _mm256_loadu_pd(scratch + idx));
  }
  const __m128d lo_pd = _mm256_castpd256_pd128(vacc);
  const __m128d hi_pd = _mm256_extractf128_pd(vacc, 1);
  const __m128d sum2 = _mm_add_pd(lo_pd, hi_pd);
  const __m128d sum1 = _mm_hadd_pd(sum2, sum2);
  double acc = _mm_cvtsd_f64(sum1);
  double carry = 0.0;
  for (; idx < count; ++idx) {
    const double val = scratch[idx];
    const double y_val = val - carry;
    const double temp_sum = acc + y_val;
    carry = (temp_sum - acc) - y_val;
    acc = temp_sum;
  }
  return acc;
}

/**
 * @brief Computes the Debye scattering equation sum for single precision (AVX2 version).
 * @param[in] q_magnitude The scattering vector magnitude.
 * @param[in] distances Source array of distance values.
 * @param[out] scratch Scratchpad array for intermediate sinc evaluations.
 * @param[in] count Number of distances to process.
 * @return The computed Debye sum.
 */
inline float debye_sum(float q_magnitude, const float *CORRELATION_RESTRICT distances,
                       float *CORRELATION_RESTRICT scratch, std::size_t count) noexcept {
  if (q_magnitude < 1.0e-9F) {
    return static_cast<float>(count);
  }

  for (std::size_t idx = 0; idx < count; ++idx) {
    const float val_x = q_magnitude * distances[idx];
    scratch[idx] = (val_x < 1.0e-4F) ? (1.0F - (val_x * val_x) / 6.0F) : (std::sin(val_x) / val_x);
  }

  __m256 vacc = _mm256_setzero_ps();
  std::size_t idx = 0;
  for (; idx + 8 <= count; idx += 8) {
    vacc = _mm256_add_ps(vacc, _mm256_loadu_ps(scratch + idx));
  }
  alignas(32) float float_buf[8];
  _mm256_storeu_ps(float_buf, vacc);
  float acc = float_buf[0] + float_buf[1] + float_buf[2] + float_buf[3] + float_buf[4] + float_buf[5] + float_buf[6] +
              float_buf[7];
  for (; idx < count; ++idx) {
    acc += scratch[idx];
  }
  return acc;
}

/**
 * @brief Normalizes Radial Distribution Function (RDF) bins for double precision (AVX2 version).
 * @param[in] hist_data Histogram array of raw bin counts.
 * @param[in] radial_bins Array of radial bin positions.
 * @param[in] g_norm Normalization factor for g(r).
 * @param[in] inv_Ni_dr Inverse scale factor times dr for species i.
 * @param[in] inv_Nj_dr Inverse scale factor times dr for species j.
 * @param[in] pi4_rho_j Density scaling factor (4 * pi * rho_j).
 * @param[out] g_out Output array for g(r).
 * @param[out] G_out Output array for G(r).
 * @param[out] J_out Output array for J(r).
 * @param[out] Jinv_out Output array for J^{-1}(r).
 * @param[in] count Total number of bins.
 */
inline void normalize_rdf_bins(const double *CORRELATION_RESTRICT hist_data,
                               const double *CORRELATION_RESTRICT radial_bins, double g_norm, double inv_Ni_dr,
                               double inv_Nj_dr, double pi4_rho_j, double *CORRELATION_RESTRICT g_out,
                               double *CORRELATION_RESTRICT G_out, double *CORRELATION_RESTRICT J_out,
                               double *CORRELATION_RESTRICT Jinv_out, std::size_t count) noexcept {
  const __m256d vg_norm = _mm256_set1_pd(g_norm);
  const __m256d v1 = _mm256_set1_pd(1.0);
  const __m256d vinNidr = _mm256_set1_pd(inv_Ni_dr);
  const __m256d vinNjdr = _mm256_set1_pd(inv_Nj_dr);
  const __m256d vpi4rho = _mm256_set1_pd(pi4_rho_j);

  std::size_t idx = 1;
  for (; idx + 4 <= count; idx += 4) {
    const __m256d vH = _mm256_loadu_pd(hist_data + idx);
    const __m256d vr = _mm256_loadu_pd(radial_bins + idx);
    const __m256d vr2 = _mm256_mul_pd(vr, vr);
    const __m256d vg = _mm256_div_pd(_mm256_mul_pd(vH, vg_norm), vr2);
    _mm256_storeu_pd(g_out + idx, vg);
    _mm256_storeu_pd(G_out + idx, _mm256_mul_pd(vpi4rho, _mm256_mul_pd(vr, _mm256_sub_pd(vg, v1))));
    _mm256_storeu_pd(J_out + idx, _mm256_mul_pd(vH, vinNidr));
    _mm256_storeu_pd(Jinv_out + idx, _mm256_mul_pd(vH, vinNjdr));
  }
  for (; idx < count; ++idx) {
    const double r_val = radial_bins[idx];
    if (r_val < 1.0e-9) {
      continue;
    }
    const double g_val = hist_data[idx] * g_norm / (r_val * r_val);
    g_out[idx] = g_val;
    G_out[idx] = pi4_rho_j * r_val * (g_val - 1.0);
    J_out[idx] = hist_data[idx] * inv_Ni_dr;
    Jinv_out[idx] = hist_data[idx] * inv_Nj_dr;
  }
}

/**
 * @brief Normalizes Radial Distribution Function (RDF) bins for single precision (AVX2 version).
 * @param[in] hist_data Histogram array of raw bin counts.
 * @param[in] radial_bins Array of radial bin positions.
 * @param[in] g_norm Normalization factor for g(r).
 * @param[in] inv_Ni_dr Inverse scale factor times dr for species i.
 * @param[in] inv_Nj_dr Inverse scale factor times dr for species j.
 * @param[in] pi4_rho_j Density scaling factor (4 * pi * rho_j).
 * @param[out] g_out Output array for g(r).
 * @param[out] G_out Output array for G(r).
 * @param[out] J_out Output array for J(r).
 * @param[out] Jinv_out Output array for J^{-1}(r).
 * @param[in] count Total number of bins.
 */
inline void normalize_rdf_bins(const float *CORRELATION_RESTRICT hist_data,
                               const float *CORRELATION_RESTRICT radial_bins, float g_norm, float inv_Ni_dr,
                               float inv_Nj_dr, float pi4_rho_j, float *CORRELATION_RESTRICT g_out,
                               float *CORRELATION_RESTRICT G_out, float *CORRELATION_RESTRICT J_out,
                               float *CORRELATION_RESTRICT Jinv_out, std::size_t count) noexcept {
  const __m256 vg_norm = _mm256_set1_ps(g_norm);
  const __m256 v1 = _mm256_set1_ps(1.0F);
  const __m256 vinNidr = _mm256_set1_ps(inv_Ni_dr);
  const __m256 vinNjdr = _mm256_set1_ps(inv_Nj_dr);
  const __m256 vpi4rho = _mm256_set1_ps(pi4_rho_j);

  std::size_t idx = 1;
  for (; idx + 8 <= count; idx += 8) {
    const __m256 vH = _mm256_loadu_ps(hist_data + idx);
    const __m256 vr = _mm256_loadu_ps(radial_bins + idx);
    const __m256 vr2 = _mm256_mul_ps(vr, vr);
    const __m256 vg = _mm256_div_ps(_mm256_mul_ps(vH, vg_norm), vr2);
    _mm256_storeu_ps(g_out + idx, vg);
    _mm256_storeu_ps(G_out + idx, _mm256_mul_ps(vpi4rho, _mm256_mul_ps(vr, _mm256_sub_ps(vg, v1))));
    _mm256_storeu_ps(J_out + idx, _mm256_mul_ps(vH, vinNidr));
    _mm256_storeu_ps(Jinv_out + idx, _mm256_mul_ps(vH, vinNjdr));
  }
  for (; idx < count; ++idx) {
    const float r_val = radial_bins[idx];
    if (r_val < 1.0e-9F) {
      continue;
    }
    const float g_val = hist_data[idx] * g_norm / (r_val * r_val);
    g_out[idx] = g_val;
    G_out[idx] = pi4_rho_j * r_val * (g_val - 1.0F);
    J_out[idx] = hist_data[idx] * inv_Ni_dr;
    Jinv_out[idx] = hist_data[idx] * inv_Nj_dr;
  }
}

/**
 * @brief Scales a double array by a constant factor (AVX2 version).
 * @param[in,out] arr Array to scale in-place.
 * @param[in] scale_factor Scalar multiplier.
 * @param[in] count Number of elements.
 */
inline void scale_bins(double *arr, double scale_factor, std::size_t count) noexcept {
  const __m256d vs = _mm256_set1_pd(scale_factor);
  std::size_t idx = 0;
  for (; idx + 4 <= count; idx += 4) {
    _mm256_storeu_pd(arr + idx, _mm256_mul_pd(_mm256_loadu_pd(arr + idx), vs));
  }
  for (; idx < count; ++idx) {
    arr[idx] *= scale_factor;
  }
}

/**
 * @brief Scales a float array by a constant factor (AVX2 version).
 * @param[in,out] arr Array to scale in-place.
 * @param[in] scale_factor Scalar multiplier.
 * @param[in] count Number of elements.
 */
inline void scale_bins(float *arr, float scale_factor, std::size_t count) noexcept {
  const __m256 vs = _mm256_set1_ps(scale_factor);
  std::size_t idx = 0;
  for (; idx + 8 <= count; idx += 8) {
    _mm256_storeu_ps(arr + idx, _mm256_mul_ps(_mm256_loadu_ps(arr + idx), vs));
  }
  for (; idx < count; ++idx) {
    arr[idx] *= scale_factor;
  }
}

#else // Scalar fallback

/**
 * @brief Computes the squared distances from a reference point to a single-precision block of positions (Scalar
 * fallback).
 * @param[in] ref_x x-coordinate of reference point.
 * @param[in] ref_y y-coordinate of reference point.
 * @param[in] ref_z z-coordinate of reference point.
 * @param[in] block Structure of Arrays containing target positions.
 * @param[out] out_dsq Destination array for computed squared distances.
 */
inline void compute_dsq_block(float ref_x, float ref_y, float ref_z, const PositionBlockT<float> &block,
                              float *CORRELATION_RESTRICT out_dsq) noexcept {
  for (std::size_t idx = 0; idx < block.count; ++idx) {
    out_dsq[idx] = dist_sq_scalar<float>({ref_x, ref_y, ref_z}, {block.x[idx], block.y[idx], block.z[idx]});
  }
}

/**
 * @brief Computes the squared distances from a reference point to a double-precision block of positions (Scalar
 * fallback).
 * @param[in] ref_x x-coordinate of reference point.
 * @param[in] ref_y y-coordinate of reference point.
 * @param[in] ref_z z-coordinate of reference point.
 * @param[in] block Structure of Arrays containing target positions.
 * @param[out] out_dsq Destination array for computed squared distances.
 */
inline void compute_dsq_block(double ref_x, double ref_y, double ref_z, const PositionBlockT<double> &block,
                              double *CORRELATION_RESTRICT out_dsq) noexcept {
  for (std::size_t idx = 0; idx < block.count; ++idx) {
    out_dsq[idx] = dist_sq_scalar<double>({ref_x, ref_y, ref_z}, {block.x[idx], block.y[idx], block.z[idx]});
  }
}

/**
 * @brief Computes the dot product of two double arrays (Scalar fallback).
 * @param[in] input_a First array.
 * @param[in] input_b Second array.
 * @param[in] count Number of elements.
 * @return The dot product value.
 */
inline double simd_dot(const double *CORRELATION_RESTRICT input_a, const double *CORRELATION_RESTRICT input_b,
                       std::size_t count) noexcept {
  double acc = 0.0;
  double carry = 0.0;
  for (std::size_t idx = 0; idx < count; ++idx) {
    const double val = input_a[idx] * input_b[idx];
    const double y_val = val - carry;
    const double temp_sum = acc + y_val;
    carry = (temp_sum - acc) - y_val;
    acc = temp_sum;
  }
  return acc;
}

/**
 * @brief Computes the dot product of two float arrays (Scalar fallback).
 * @param[in] input_a First array.
 * @param[in] input_b Second array.
 * @param[in] count Number of elements.
 * @return The dot product value.
 */
inline float simd_dot(const float *CORRELATION_RESTRICT input_a, const float *CORRELATION_RESTRICT input_b,
                      std::size_t count) noexcept {
  float acc = 0.0F;
  for (std::size_t idx = 0; idx < count; ++idx) {
    acc += input_a[idx] * input_b[idx];
  }
  return acc;
}

/**
 * @brief Computes a sinc-weighted integral over a double range (Scalar fallback).
 * @param[in] q_magnitude The scattering vector magnitude.
 * @param[in] integrand Values to be integrated.
 * @param[in] radial_bins Radial distances for each point.
 * @param[out] sinqr_scratch Scratchpad memory for sinc calculations.
 * @param[in] count Number of data points.
 * @return The accumulated integral value.
 */
inline double sinc_integral(double q_magnitude, const double *CORRELATION_RESTRICT integrand,
                            const double *CORRELATION_RESTRICT radial_bins, double *CORRELATION_RESTRICT sinqr_scratch,
                            std::size_t count) noexcept {
  for (std::size_t idx = 0; idx < count; ++idx) {
    sinqr_scratch[idx] = std::sin(q_magnitude * radial_bins[idx]);
  }
  return simd_dot(integrand, sinqr_scratch, count);
}

/**
 * @brief Computes a sinc-weighted integral over a float range (Scalar fallback).
 * @param[in] q_magnitude The scattering vector magnitude.
 * @param[in] integrand Values to be integrated.
 * @param[in] radial_bins Radial distances for each point.
 * @param[out] sinqr_scratch Scratchpad memory for sinc calculations.
 * @param[in] count Number of data points.
 * @return The accumulated integral value.
 */
inline float sinc_integral(float q_magnitude, const float *CORRELATION_RESTRICT integrand,
                           const float *CORRELATION_RESTRICT radial_bins, float *CORRELATION_RESTRICT sinqr_scratch,
                           std::size_t count) noexcept {
  for (std::size_t idx = 0; idx < count; ++idx) {
    sinqr_scratch[idx] = std::sin(q_magnitude * radial_bins[idx]);
  }
  return simd_dot(integrand, sinqr_scratch, count);
}

/**
 * @brief Computes the Debye scattering equation sum for double precision (Scalar fallback).
 * @param[in] q_magnitude The scattering vector magnitude.
 * @param[in] distances Source array of distance values.
 * @param[out] scratch Scratchpad array for intermediate sinc evaluations.
 * @param[in] count Number of distances to process.
 * @return The computed Debye sum.
 */
inline double debye_sum(double q_magnitude, const double *CORRELATION_RESTRICT distances,
                        double *CORRELATION_RESTRICT /*scratch*/, std::size_t count) noexcept {
  if (q_magnitude < 1.0e-9) {
    return static_cast<double>(count);
  }
  double acc = 0.0;
  double carry = 0.0;
  for (std::size_t idx = 0; idx < count; ++idx) {
    const double val_x = q_magnitude * distances[idx];
    const double val = (val_x < 1.0e-4) ? (1.0 - (val_x * val_x) / 6.0) : (std::sin(val_x) / val_x);
    const double y_val = val - carry;
    const double temp_sum = acc + y_val;
    carry = (temp_sum - acc) - y_val;
    acc = temp_sum;
  }
  return acc;
}

/**
 * @brief Computes the Debye scattering equation sum for single precision (Scalar fallback).
 * @param[in] q_magnitude The scattering vector magnitude.
 * @param[in] distances Source array of distance values.
 * @param[out] scratch Scratchpad array for intermediate sinc evaluations.
 * @param[in] count Number of distances to process.
 * @return The computed Debye sum.
 */
inline float debye_sum(float q_magnitude, const float *CORRELATION_RESTRICT distances,
                       float *CORRELATION_RESTRICT /*scratch*/, std::size_t count) noexcept {
  if (q_magnitude < 1.0e-9F) {
    return static_cast<float>(count);
  }
  float acc = 0.0F;
  for (std::size_t idx = 0; idx < count; ++idx) {
    const float val_x = q_magnitude * distances[idx];
    const float val = (val_x < 1.0e-4F) ? (1.0F - (val_x * val_x) / 6.0F) : (std::sin(val_x) / val_x);
    acc += val;
  }
  return acc;
}

/**
 * @brief Normalizes Radial Distribution Function (RDF) bins for double precision (Scalar fallback).
 * @param[in] hist_data Histogram array of raw bin counts.
 * @param[in] radial_bins Array of radial bin positions.
 * @param[in] g_norm Normalization factor for g(r).
 * @param[in] inv_Ni_dr Inverse scale factor times dr for species i.
 * @param[in] inv_Nj_dr Inverse scale factor times dr for species j.
 * @param[in] pi4_rho_j Density scaling factor (4 * pi * rho_j).
 * @param[out] g_out Output array for g(r).
 * @param[out] G_out Output array for G(r).
 * @param[out] J_out Output array for J(r).
 * @param[out] Jinv_out Output array for J^{-1}(r).
 * @param[in] count Total number of bins.
 */
inline void normalize_rdf_bins(const double *CORRELATION_RESTRICT hist_data,
                               const double *CORRELATION_RESTRICT radial_bins, double g_norm, double inv_Ni_dr,
                               double inv_Nj_dr, double pi4_rho_j, double *CORRELATION_RESTRICT g_out,
                               double *CORRELATION_RESTRICT G_out, double *CORRELATION_RESTRICT J_out,
                               double *CORRELATION_RESTRICT Jinv_out, std::size_t count) noexcept {
  for (std::size_t idx = 1; idx < count; ++idx) {
    const double r_val = radial_bins[idx];
    if (r_val < 1.0e-9) {
      continue;
    }
    const double g_val = hist_data[idx] * g_norm / (r_val * r_val);
    g_out[idx] = g_val;
    G_out[idx] = pi4_rho_j * r_val * (g_val - 1.0);
    J_out[idx] = hist_data[idx] * inv_Ni_dr;
    Jinv_out[idx] = hist_data[idx] * inv_Nj_dr;
  }
}

/**
 * @brief Scales a double array by a constant factor (Scalar fallback).
 * @param[in,out] arr Array to scale in-place.
 * @param[in] scale_factor Scalar multiplier.
 * @param[in] count Number of elements.
 */
inline void scale_bins(double *arr, double scale_factor, std::size_t count) noexcept {
  for (std::size_t idx = 0; idx < count; ++idx) {
    arr[idx] *= scale_factor;
  }
}

#endif

/**
 * @brief Normalizes Radial Distribution Function (RDF) bins for generic scalar types.
 * @tparam T Floating point element type.
 * @param[in] hist_data Histogram array of raw bin counts.
 * @param[in] radial_bins Array of radial bin positions.
 * @param[in] g_norm Normalization factor for g(r).
 * @param[in] inv_Ni_dr Inverse scale factor times dr for species i.
 * @param[in] inv_Nj_dr Inverse scale factor times dr for species j.
 * @param[in] pi4_rho_j Density scaling factor (4 * pi * rho_j).
 * @param[out] g_out Output array for g(r).
 * @param[out] G_out Output array for G(r).
 * @param[out] J_out Output array for J(r).
 * @param[out] Jinv_out Output array for J^{-1}(r).
 * @param[in] count Total number of bins.
 */
template <typename T>
inline void normalize_rdf_bins(const T *CORRELATION_RESTRICT hist_data, const T *CORRELATION_RESTRICT radial_bins,
                               T g_norm, T inv_Ni_dr, T inv_Nj_dr, T pi4_rho_j, T *CORRELATION_RESTRICT g_out,
                               T *CORRELATION_RESTRICT G_out, T *CORRELATION_RESTRICT J_out,
                               T *CORRELATION_RESTRICT Jinv_out, std::size_t count) noexcept {
  for (std::size_t idx = 1; idx < count; ++idx) {
    const T r_val = radial_bins[idx];
    if (r_val < static_cast<T>(1.0e-9)) {
      continue;
    }
    const T g_val = hist_data[idx] * g_norm / (r_val * r_val);
    g_out[idx] = g_val;
    G_out[idx] = pi4_rho_j * r_val * (g_val - static_cast<T>(1.0));
    J_out[idx] = hist_data[idx] * inv_Ni_dr;
    Jinv_out[idx] = hist_data[idx] * inv_Nj_dr;
  }
}

/**
 * @brief Scales an array by a constant factor for generic scalar types.
 * @tparam T Floating point element type.
 * @param[in,out] arr Array to scale in-place.
 * @param[in] scale_factor Scalar multiplier.
 * @param[in] count Number of elements.
 */
template <typename T> inline void scale_bins(T *arr, T scale_factor, std::size_t count) noexcept {
  for (std::size_t idx = 0; idx < count; ++idx) {
    arr[idx] *= scale_factor;
  }
}

// ---------------------------------------------------------------------------
// dot_block
// ---------------------------------------------------------------------------

/**
 * @brief Computes dot products between a reference vector and a block of double vectors.
 * @param[in] v1x x-component of reference vector.
 * @param[in] v1y y-component of reference vector.
 * @param[in] v1z z-component of reference vector.
 * @param[in] v2x Array of x-components for target vectors.
 * @param[in] v2y Array of y-components for target vectors.
 * @param[in] v2z Array of z-components for target vectors.
 * @param[out] out_dot Output array for scalar dot products.
 * @param[in] count Number of elements.
 */
#ifdef CORRELATION_SIMD_AVX512
inline void dot_block(double v1x, double v1y, double v1z, const double *CORRELATION_RESTRICT v2x,
                      const double *CORRELATION_RESTRICT v2y, const double *CORRELATION_RESTRICT v2z,
                      double *CORRELATION_RESTRICT out_dot, std::size_t count) noexcept {
  const __m512d vv1x = _mm512_set1_pd(v1x);
  const __m512d vv1y = _mm512_set1_pd(v1y);
  const __m512d vv1z = _mm512_set1_pd(v1z);
  std::size_t idx = 0;
  for (; idx + 8 <= count; idx += 8) {
    const __m512d d_res = _mm512_fmadd_pd(
        vv1x, _mm512_loadu_pd(v2x + idx),
        _mm512_fmadd_pd(vv1y, _mm512_loadu_pd(v2y + idx), _mm512_mul_pd(vv1z, _mm512_loadu_pd(v2z + idx))));
    _mm512_storeu_pd(out_dot + idx, d_res);
  }
  for (; idx < count; ++idx) {
    out_dot[idx] = v1x * v2x[idx] + v1y * v2y[idx] + v1z * v2z[idx];
  }
}

/**
 * @brief Computes dot products between a reference vector and a block of float vectors.
 * @param[in] v1x x-component of reference vector.
 * @param[in] v1y y-component of reference vector.
 * @param[in] v1z z-component of reference vector.
 * @param[in] v2x Array of x-components for target vectors.
 * @param[in] v2y Array of y-components for target vectors.
 * @param[in] v2z Array of z-components for target vectors.
 * @param[out] out_dot Output array for scalar dot products.
 * @param[in] count Number of elements.
 */
inline void dot_block(float v1x, float v1y, float v1z, const float *CORRELATION_RESTRICT v2x,
                      const float *CORRELATION_RESTRICT v2y, const float *CORRELATION_RESTRICT v2z,
                      float *CORRELATION_RESTRICT out_dot, std::size_t count) noexcept {
  const __m512 vv1x = _mm512_set1_ps(v1x);
  const __m512 vv1y = _mm512_set1_ps(v1y);
  const __m512 vv1z = _mm512_set1_ps(v1z);
  std::size_t idx = 0;
  for (; idx + 16 <= count; idx += 16) {
    const __m512 d_res = _mm512_fmadd_ps(
        vv1x, _mm512_loadu_ps(v2x + idx),
        _mm512_fmadd_ps(vv1y, _mm512_loadu_ps(v2y + idx), _mm512_mul_ps(vv1z, _mm512_loadu_ps(v2z + idx))));
    _mm512_storeu_ps(out_dot + idx, d_res);
  }
  for (; idx < count; ++idx) {
    out_dot[idx] = v1x * v2x[idx] + v1y * v2y[idx] + v1z * v2z[idx];
  }
}
#elif defined(CORRELATION_SIMD_AVX2)
inline void dot_block(double v1x, double v1y, double v1z, const double *CORRELATION_RESTRICT v2x,
                      const double *CORRELATION_RESTRICT v2y, const double *CORRELATION_RESTRICT v2z,
                      double *CORRELATION_RESTRICT out_dot, std::size_t count) noexcept {
  const __m256d vv1x = _mm256_set1_pd(v1x);
  const __m256d vv1y = _mm256_set1_pd(v1y);
  const __m256d vv1z = _mm256_set1_pd(v1z);
  std::size_t idx = 0;
  for (; idx + 4 <= count; idx += 4) {
#ifdef __FMA__
    const __m256d d_res = _mm256_fmadd_pd(
        vv1x, _mm256_loadu_pd(v2x + idx),
        _mm256_fmadd_pd(vv1y, _mm256_loadu_pd(v2y + idx), _mm256_mul_pd(vv1z, _mm256_loadu_pd(v2z + idx))));
#else
    const __m256d d_res = _mm256_add_pd(_mm256_mul_pd(vv1x, _mm256_loadu_pd(v2x + idx)),
                                        _mm256_add_pd(_mm256_mul_pd(vv1y, _mm256_loadu_pd(v2y + idx)),
                                                      _mm256_mul_pd(vv1z, _mm256_loadu_pd(v2z + idx))));
#endif
    _mm256_storeu_pd(out_dot + idx, d_res);
  }
  for (; idx < count; ++idx) {
    out_dot[idx] = v1x * v2x[idx] + v1y * v2y[idx] + v1z * v2z[idx];
  }
}

inline void dot_block(float v1x, float v1y, float v1z, const float *CORRELATION_RESTRICT v2x,
                      const float *CORRELATION_RESTRICT v2y, const float *CORRELATION_RESTRICT v2z,
                      float *CORRELATION_RESTRICT out_dot, std::size_t count) noexcept {
  const __m256 vv1x = _mm256_set1_ps(v1x);
  const __m256 vv1y = _mm256_set1_ps(v1y);
  const __m256 vv1z = _mm256_set1_ps(v1z);
  std::size_t idx = 0;
  for (; idx + 8 <= count; idx += 8) {
#ifdef __FMA__
    const __m256 d_res = _mm256_fmadd_ps(
        vv1x, _mm256_loadu_ps(v2x + idx),
        _mm256_fmadd_ps(vv1y, _mm256_loadu_ps(v2y + idx), _mm256_mul_ps(vv1z, _mm256_loadu_ps(v2z + idx))));
#else
    const __m256 d_res = _mm256_add_ps(_mm256_mul_ps(vv1x, _mm256_loadu_ps(v2x + idx)),
                                       _mm256_add_ps(_mm256_mul_ps(vv1y, _mm256_loadu_ps(v2y + idx)),
                                                     _mm256_mul_ps(vv1z, _mm256_loadu_ps(v2z + idx))));
#endif
    _mm256_storeu_ps(out_dot + idx, d_res);
  }
  for (; idx < count; ++idx) {
    out_dot[idx] = v1x * v2x[idx] + v1y * v2y[idx] + v1z * v2z[idx];
  }
}
#else
inline void dot_block(double v1x, double v1y, double v1z, const double *CORRELATION_RESTRICT v2x,
                      const double *CORRELATION_RESTRICT v2y, const double *CORRELATION_RESTRICT v2z,
                      double *CORRELATION_RESTRICT out_dot, std::size_t count) noexcept {
  for (std::size_t idx = 0; idx < count; ++idx) {
    out_dot[idx] = v1x * v2x[idx] + v1y * v2y[idx] + v1z * v2z[idx];
  }
}
inline void dot_block(float v1x, float v1y, float v1z, const float *CORRELATION_RESTRICT v2x,
                      const float *CORRELATION_RESTRICT v2y, const float *CORRELATION_RESTRICT v2z,
                      float *CORRELATION_RESTRICT out_dot, std::size_t count) noexcept {
  for (std::size_t idx = 0; idx < count; ++idx) {
    out_dot[idx] = v1x * v2x[idx] + v1y * v2y[idx] + v1z * v2z[idx];
  }
}
#endif

// ---------------------------------------------------------------------------
// Utility
// ---------------------------------------------------------------------------

/**
 * @brief Populates vectors with standard atom block positions for SIMD.
 * @tparam AtomRange A range of atom objects (e.g., std::vector<Atom>).
 * @tparam T Coordinate scalar type (defaults to real_t).
 * @param[in] atoms The source atom collection.
 * @param[in] begin_idx Start index of the block.
 * @param[in] end_idx End index of the block.
 * @param[out] xs Output x-coordinate vector.
 * @param[out] ys Output y-coordinate vector.
 * @param[out] zs Output z-coordinate vector.
 * @return The number of atoms added to the block.
 */
template <typename AtomRange, typename T = real_t>
inline std::size_t fill_position_block(const AtomRange &atoms, std::size_t begin_idx, std::size_t end_idx,
                                       std::vector<T> &xs, std::vector<T> &ys, std::vector<T> &zs) noexcept {
  const std::size_t count = end_idx - begin_idx;
  xs.resize(count);
  ys.resize(count);
  zs.resize(count);
  for (std::size_t idx = 0; idx < count; ++idx) {
    const auto &pos = atoms[begin_idx + idx].position();
    xs[idx] = static_cast<T>(pos.x());
    ys[idx] = static_cast<T>(pos.y());
    zs[idx] = static_cast<T>(pos.z());
  }
  return count;
}

/**
 * @brief Computes the sum of complex exponentials for a given double query vector.
 * @param[in] qx x-component of q-vector.
 * @param[in] qy y-component of q-vector.
 * @param[in] qz z-component of q-vector.
 * @param[in] xs Array of atom x-positions.
 * @param[in] ys Array of atom y-positions.
 * @param[in] zs Array of atom z-positions.
 * @param[in] count Number of atoms.
 * @param[out] cos_sum Output for the real part (sum of cosines).
 * @param[out] sin_sum Output for the imaginary part (sum of sines).
 */
inline void complex_exp_sum(double qx, double qy, double qz, const double *CORRELATION_RESTRICT xs,
                            const double *CORRELATION_RESTRICT ys, const double *CORRELATION_RESTRICT zs,
                            std::size_t count, double &cos_sum, double &sin_sum) noexcept {
  cos_sum = 0.0;
  sin_sum = 0.0;
  for (std::size_t idx = 0; idx < count; ++idx) {
    const double phase = qx * xs[idx] + qy * ys[idx] + qz * zs[idx];
    cos_sum += std::cos(phase);
    sin_sum += std::sin(phase);
  }
}

/**
 * @brief Computes the sum of complex exponentials for a given float query vector.
 * @param[in] qx x-component of q-vector.
 * @param[in] qy y-component of q-vector.
 * @param[in] qz z-component of q-vector.
 * @param[in] xs Array of atom x-positions.
 * @param[in] ys Array of atom y-positions.
 * @param[in] zs Array of atom z-positions.
 * @param[in] count Number of atoms.
 * @param[out] cos_sum Output for the real part (sum of cosines).
 * @param[out] sin_sum Output for the imaginary part (sum of sines).
 */
inline void complex_exp_sum(float qx, float qy, float qz, const float *CORRELATION_RESTRICT xs,
                            const float *CORRELATION_RESTRICT ys, const float *CORRELATION_RESTRICT zs,
                            std::size_t count, float &cos_sum, float &sin_sum) noexcept {
  cos_sum = 0.0F;
  sin_sum = 0.0F;
  for (std::size_t idx = 0; idx < count; ++idx) {
    const float phase = qx * xs[idx] + qy * ys[idx] + qz * zs[idx];
    cos_sum += std::cos(phase);
    sin_sum += std::sin(phase);
  }
}

// ---------------------------------------------------------------------------
// miller_phase_sum
// ---------------------------------------------------------------------------

#ifdef CORRELATION_SIMD_AVX512

/**
 * @brief Computes the Miller phase sum across a block of double angles (AVX-512 version).
 * @param[in] cos1 Cosine of first angle block.
 * @param[in] sin1 Sine of first angle block.
 * @param[in] cos2 Cosine of second angle block.
 * @param[in] sin2 Sine of second angle block.
 * @param[in] cos3 Cosine of third angle block.
 * @param[in] sin3 Sine of third angle block.
 * @param[in] count Number of elements in the blocks.
 * @param[out] cos_sum Output for the cosine component of the sum.
 * @param[out] sin_sum Output for the sine component of the sum.
 */
inline void miller_phase_sum(const double *CORRELATION_RESTRICT cos1, const double *CORRELATION_RESTRICT sin1,
                             const double *CORRELATION_RESTRICT cos2, const double *CORRELATION_RESTRICT sin2,
                             const double *CORRELATION_RESTRICT cos3, const double *CORRELATION_RESTRICT sin3,
                             std::size_t count, double &cos_sum, double &sin_sum) noexcept {
  __m512d vc_sum = _mm512_setzero_pd();
  __m512d vs_sum = _mm512_setzero_pd();
  std::size_t idx = 0;
  for (; idx + 8 <= count; idx += 8) {
    const __m512d vc1 = _mm512_loadu_pd(cos1 + idx);
    const __m512d vs1 = _mm512_loadu_pd(sin1 + idx);
    const __m512d vc2 = _mm512_loadu_pd(cos2 + idx);
    const __m512d vs2 = _mm512_loadu_pd(sin2 + idx);
    const __m512d vc3 = _mm512_loadu_pd(cos3 + idx);
    const __m512d vs3 = _mm512_loadu_pd(sin3 + idx);
    const __m512d vc12 = _mm512_fmsub_pd(vc1, vc2, _mm512_mul_pd(vs1, vs2));
    const __m512d vs12 = _mm512_fmadd_pd(vs1, vc2, _mm512_mul_pd(vc1, vs2));
    const __m512d vc123 = _mm512_fmsub_pd(vc12, vc3, _mm512_mul_pd(vs12, vs3));
    const __m512d vs123 = _mm512_fmadd_pd(vs12, vc3, _mm512_mul_pd(vc12, vs3));
    vc_sum = _mm512_add_pd(vc_sum, vc123);
    vs_sum = _mm512_add_pd(vs_sum, vs123);
  }
  cos_sum = _mm512_reduce_add_pd(vc_sum);
  sin_sum = _mm512_reduce_add_pd(vs_sum);
  for (; idx < count; ++idx) {
    const double c12 = cos1[idx] * cos2[idx] - sin1[idx] * sin2[idx];
    const double s12 = sin1[idx] * cos2[idx] + cos1[idx] * sin2[idx];
    cos_sum += c12 * cos3[idx] - s12 * sin3[idx];
    sin_sum += s12 * cos3[idx] + c12 * sin3[idx];
  }
}

/**
 * @brief Computes the Miller phase sum across a block of float angles (AVX-512 version).
 * @param[in] cos1 Cosine of first angle block.
 * @param[in] sin1 Sine of first angle block.
 * @param[in] cos2 Cosine of second angle block.
 * @param[in] sin2 Sine of second angle block.
 * @param[in] cos3 Cosine of third angle block.
 * @param[in] sin3 Sine of third angle block.
 * @param[in] count Number of elements in the blocks.
 * @param[out] cos_sum Output for the cosine component of the sum.
 * @param[out] sin_sum Output for the sine component of the sum.
 */
inline void miller_phase_sum(const float *CORRELATION_RESTRICT cos1, const float *CORRELATION_RESTRICT sin1,
                             const float *CORRELATION_RESTRICT cos2, const float *CORRELATION_RESTRICT sin2,
                             const float *CORRELATION_RESTRICT cos3, const float *CORRELATION_RESTRICT sin3,
                             std::size_t count, float &cos_sum, float &sin_sum) noexcept {
  __m512 vc_sum = _mm512_setzero_ps();
  __m512 vs_sum = _mm512_setzero_ps();
  std::size_t idx = 0;
  for (; idx + 16 <= count; idx += 16) {
    const __m512 vc1 = _mm512_loadu_ps(cos1 + idx);
    const __m512 vs1 = _mm512_loadu_ps(sin1 + idx);
    const __m512 vc2 = _mm512_loadu_ps(cos2 + idx);
    const __m512 vs2 = _mm512_loadu_ps(sin2 + idx);
    const __m512 vc3 = _mm512_loadu_ps(cos3 + idx);
    const __m512 vs3 = _mm512_loadu_ps(sin3 + idx);
    const __m512 vc12 = _mm512_fmsub_ps(vc1, vc2, _mm512_mul_ps(vs1, vs2));
    const __m512 vs12 = _mm512_fmadd_ps(vs1, vc2, _mm512_mul_ps(vc1, vs2));
    const __m512 vc123 = _mm512_fmsub_ps(vc12, vc3, _mm512_mul_ps(vs12, vs3));
    const __m512 vs123 = _mm512_fmadd_ps(vs12, vc3, _mm512_mul_ps(vc12, vs3));
    vc_sum = _mm512_add_ps(vc_sum, vc123);
    vs_sum = _mm512_add_ps(vs_sum, vs123);
  }
  cos_sum = _mm512_reduce_add_ps(vc_sum);
  sin_sum = _mm512_reduce_add_ps(vs_sum);
  for (; idx < count; ++idx) {
    const float c12 = cos1[idx] * cos2[idx] - sin1[idx] * sin2[idx];
    const float s12 = sin1[idx] * cos2[idx] + cos1[idx] * sin2[idx];
    cos_sum += c12 * cos3[idx] - s12 * sin3[idx];
    sin_sum += s12 * cos3[idx] + c12 * sin3[idx];
  }
}

#elif defined(CORRELATION_SIMD_AVX2)

/**
 * @brief Computes the Miller phase sum across a block of double angles (AVX2 version).
 * @param[in] cos1 Cosine of first angle block.
 * @param[in] sin1 Sine of first angle block.
 * @param[in] cos2 Cosine of second angle block.
 * @param[in] sin2 Sine of second angle block.
 * @param[in] cos3 Cosine of third angle block.
 * @param[in] sin3 Sine of third angle block.
 * @param[in] count Number of elements in the blocks.
 * @param[out] cos_sum Output for the cosine component of the sum.
 * @param[out] sin_sum Output for the sine component of the sum.
 */
inline void miller_phase_sum(const double *CORRELATION_RESTRICT cos1, const double *CORRELATION_RESTRICT sin1,
                             const double *CORRELATION_RESTRICT cos2, const double *CORRELATION_RESTRICT sin2,
                             const double *CORRELATION_RESTRICT cos3, const double *CORRELATION_RESTRICT sin3,
                             std::size_t count, double &cos_sum, double &sin_sum) noexcept {
  __m256d vc_sum = _mm256_setzero_pd();
  __m256d vs_sum = _mm256_setzero_pd();
  std::size_t idx = 0;
  for (; idx + 4 <= count; idx += 4) {
    const __m256d vc1 = _mm256_loadu_pd(cos1 + idx);
    const __m256d vs1 = _mm256_loadu_pd(sin1 + idx);
    const __m256d vc2 = _mm256_loadu_pd(cos2 + idx);
    const __m256d vs2 = _mm256_loadu_pd(sin2 + idx);
    const __m256d vc3 = _mm256_loadu_pd(cos3 + idx);
    const __m256d vs3 = _mm256_loadu_pd(sin3 + idx);
#ifdef __FMA__
    const __m256d vc12 = _mm256_fmsub_pd(vc1, vc2, _mm256_mul_pd(vs1, vs2));
    const __m256d vs12 = _mm256_fmadd_pd(vs1, vc2, _mm256_mul_pd(vc1, vs2));
    const __m256d vc123 = _mm256_fmsub_pd(vc12, vc3, _mm256_mul_pd(vs12, vs3));
    const __m256d vs123 = _mm256_fmadd_pd(vs12, vc3, _mm256_mul_pd(vc12, vs3));
#else
    const __m256d vc12 = _mm256_sub_pd(_mm256_mul_pd(vc1, vc2), _mm256_mul_pd(vs1, vs2));
    const __m256d vs12 = _mm256_add_pd(_mm256_mul_pd(vs1, vc2), _mm256_mul_pd(vc1, vs2));
    const __m256d vc123 = _mm256_sub_pd(_mm256_mul_pd(vc12, vc3), _mm256_mul_pd(vs12, vs3));
    const __m256d vs123 = _mm256_add_pd(_mm256_mul_pd(vs12, vc3), _mm256_mul_pd(vc12, vs3));
#endif
    vc_sum = _mm256_add_pd(vc_sum, vc123);
    vs_sum = _mm256_add_pd(vs_sum, vs123);
  }
  const __m128d clo = _mm256_castpd256_pd128(vc_sum);
  const __m128d chi = _mm256_extractf128_pd(vc_sum, 1);
  const __m128d c_sum2 = _mm_add_pd(clo, chi);
  const __m128d c_sum1 = _mm_hadd_pd(c_sum2, c_sum2);
  cos_sum = _mm_cvtsd_f64(c_sum1);
  const __m128d slo = _mm256_castpd256_pd128(vs_sum);
  const __m128d shi = _mm256_extractf128_pd(vs_sum, 1);
  const __m128d s_sum2 = _mm_add_pd(slo, shi);
  const __m128d s_sum1 = _mm_hadd_pd(s_sum2, s_sum2);
  sin_sum = _mm_cvtsd_f64(s_sum1);
  for (; idx < count; ++idx) {
    const double c12 = cos1[idx] * cos2[idx] - sin1[idx] * sin2[idx];
    const double s12 = sin1[idx] * cos2[idx] + cos1[idx] * sin2[idx];
    cos_sum += c12 * cos3[idx] - s12 * sin3[idx];
    sin_sum += s12 * cos3[idx] + c12 * sin3[idx];
  }
}

/**
 * @brief Computes the Miller phase sum across a block of float angles (AVX2 version).
 * @param[in] cos1 Cosine of first angle block.
 * @param[in] sin1 Sine of first angle block.
 * @param[in] cos2 Cosine of second angle block.
 * @param[in] sin2 Sine of second angle block.
 * @param[in] cos3 Cosine of third angle block.
 * @param[in] sin3 Sine of third angle block.
 * @param[in] count Number of elements in the blocks.
 * @param[out] cos_sum Output for the cosine component of the sum.
 * @param[out] sin_sum Output for the sine component of the sum.
 */
inline void miller_phase_sum(const float *CORRELATION_RESTRICT cos1, const float *CORRELATION_RESTRICT sin1,
                             const float *CORRELATION_RESTRICT cos2, const float *CORRELATION_RESTRICT sin2,
                             const float *CORRELATION_RESTRICT cos3, const float *CORRELATION_RESTRICT sin3,
                             std::size_t count, float &cos_sum, float &sin_sum) noexcept {
  __m256 vc_sum = _mm256_setzero_ps();
  __m256 vs_sum = _mm256_setzero_ps();
  std::size_t idx = 0;
  for (; idx + 8 <= count; idx += 8) {
    const __m256 vc1 = _mm256_loadu_ps(cos1 + idx);
    const __m256 vs1 = _mm256_loadu_ps(sin1 + idx);
    const __m256 vc2 = _mm256_loadu_ps(cos2 + idx);
    const __m256 vs2 = _mm256_loadu_ps(sin2 + idx);
    const __m256 vc3 = _mm256_loadu_ps(cos3 + idx);
    const __m256 vs3 = _mm256_loadu_ps(sin3 + idx);
#ifdef __FMA__
    const __m256 vc12 = _mm256_fmsub_ps(vc1, vc2, _mm256_mul_ps(vs1, vs2));
    const __m256 vs12 = _mm256_fmadd_ps(vs1, vc2, _mm256_mul_ps(vc1, vs2));
    const __m256 vc123 = _mm256_fmsub_ps(vc12, vc3, _mm256_mul_ps(vs12, vs3));
    const __m256 vs123 = _mm256_fmadd_ps(vs12, vc3, _mm256_mul_ps(vc12, vs3));
#else
    const __m256 vc12 = _mm256_sub_ps(_mm256_mul_ps(vc1, vc2), _mm256_mul_ps(vs1, vs2));
    const __m256 vs12 = _mm256_add_ps(_mm256_mul_ps(vs1, vc2), _mm256_mul_ps(vc1, vs2));
    const __m256 vc123 = _mm256_sub_ps(_mm256_mul_ps(vc12, vc3), _mm256_mul_ps(vs12, vs3));
    const __m256 vs123 = _mm256_add_ps(_mm256_mul_ps(vs12, vc3), _mm256_mul_ps(vc12, vs3));
#endif
    vc_sum = _mm256_add_ps(vc_sum, vc123);
    vs_sum = _mm256_add_ps(vs_sum, vs123);
  }
  alignas(32) float c_buf[8];
  alignas(32) float s_buf[8];
  _mm256_storeu_ps(c_buf, vc_sum);
  _mm256_storeu_ps(s_buf, vs_sum);
  cos_sum = c_buf[0] + c_buf[1] + c_buf[2] + c_buf[3] + c_buf[4] + c_buf[5] + c_buf[6] + c_buf[7];
  sin_sum = s_buf[0] + s_buf[1] + s_buf[2] + s_buf[3] + s_buf[4] + s_buf[5] + s_buf[6] + s_buf[7];
  for (; idx < count; ++idx) {
    const float c12 = cos1[idx] * cos2[idx] - sin1[idx] * sin2[idx];
    const float s12 = sin1[idx] * cos2[idx] + cos1[idx] * sin2[idx];
    cos_sum += c12 * cos3[idx] - s12 * sin3[idx];
    sin_sum += s12 * cos3[idx] + c12 * sin3[idx];
  }
}

#else

/**
 * @brief Computes the Miller phase sum across a block of double angles (Scalar fallback).
 * @param[in] cos1 Cosine of first angle block.
 * @param[in] sin1 Sine of first angle block.
 * @param[in] cos2 Cosine of second angle block.
 * @param[in] sin2 Sine of second angle block.
 * @param[in] cos3 Cosine of third angle block.
 * @param[in] sin3 Sine of third angle block.
 * @param[in] count Number of elements in the blocks.
 * @param[out] cos_sum Output for the cosine component of the sum.
 * @param[out] sin_sum Output for the sine component of the sum.
 */
inline void miller_phase_sum(const double *CORRELATION_RESTRICT cos1, const double *CORRELATION_RESTRICT sin1,
                             const double *CORRELATION_RESTRICT cos2, const double *CORRELATION_RESTRICT sin2,
                             const double *CORRELATION_RESTRICT cos3, const double *CORRELATION_RESTRICT sin3,
                             std::size_t count, double &cos_sum, double &sin_sum) noexcept {
  cos_sum = 0.0;
  sin_sum = 0.0;
  for (std::size_t idx = 0; idx < count; ++idx) {
    const double c12 = cos1[idx] * cos2[idx] - sin1[idx] * sin2[idx];
    const double s12 = sin1[idx] * cos2[idx] + cos1[idx] * sin2[idx];
    cos_sum += c12 * cos3[idx] - s12 * sin3[idx];
    sin_sum += s12 * cos3[idx] + c12 * sin3[idx];
  }
}

/**
 * @brief Computes the Miller phase sum across a block of float angles (Scalar fallback).
 * @param[in] cos1 Cosine of first angle block.
 * @param[in] sin1 Sine of first angle block.
 * @param[in] cos2 Cosine of second angle block.
 * @param[in] sin2 Sine of second angle block.
 * @param[in] cos3 Cosine of third angle block.
 * @param[in] sin3 Sine of third angle block.
 * @param[in] count Number of elements in the blocks.
 * @param[out] cos_sum Output for the cosine component of the sum.
 * @param[out] sin_sum Output for the sine component of the sum.
 */
inline void miller_phase_sum(const float *CORRELATION_RESTRICT cos1, const float *CORRELATION_RESTRICT sin1,
                             const float *CORRELATION_RESTRICT cos2, const float *CORRELATION_RESTRICT sin2,
                             const float *CORRELATION_RESTRICT cos3, const float *CORRELATION_RESTRICT sin3,
                             std::size_t count, float &cos_sum, float &sin_sum) noexcept {
  cos_sum = 0.0F;
  sin_sum = 0.0F;
  for (std::size_t idx = 0; idx < count; ++idx) {
    const float c12 = cos1[idx] * cos2[idx] - sin1[idx] * sin2[idx];
    const float s12 = sin1[idx] * cos2[idx] + cos1[idx] * sin2[idx];
    cos_sum += c12 * cos3[idx] - s12 * sin3[idx];
    sin_sum += s12 * cos3[idx] + c12 * sin3[idx];
  }
}

#endif

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
