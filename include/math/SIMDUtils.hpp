/**
 * @file SIMDUtils.hpp
 * @brief SIMD-accelerated kernels for distance, integration, and normalization.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "math/SIMDConfig.hpp"
#include <cmath>
#include <cstddef>
#include <vector>

namespace correlation::math {

// ---------------------------------------------------------------------------
// SoA block of atom positions
// ---------------------------------------------------------------------------
/**
 * @struct PositionBlock
 * @brief Represents a Structure of Arrays (SoA) block of atom positions for SIMD processing.
 */
struct PositionBlock {
  double *x;
  double *y;
  double *z;
  std::size_t count;
};

// ---------------------------------------------------------------------------
// Scalar helper
// ---------------------------------------------------------------------------
/**
 * @brief Computes the squared distance between two 3D points.
 * 
 * @return The scalar squared distance.
 */
inline double dist_sq_scalar(double ax, double ay, double az, double bx,
                             double by, double bz) noexcept {
  double dx = bx - ax;
  double dy = by - ay;
  double dz = bz - az;
  return dx * dx + dy * dy + dz * dz;
}

// ---------------------------------------------------------------------------
// Core SIMD kernels
// ---------------------------------------------------------------------------

/**
 * @brief Computes the squared distances from a reference point to a block of positions.
 * 
 * @param ax x-coordinate of the reference point.
 * @param ay y-coordinate of the reference point.
 * @param az z-coordinate of the reference point.
 * @param block Structure of Arrays containing the target positions.
 * @param out_dsq Destination array for the computed squared distances.
 */
#if defined(CORRELATION_SIMD_AVX512)

inline void compute_dsq_block(double ax, double ay, double az,
                              const PositionBlock &block,
                              double *CORRELATION_RESTRICT out_dsq) noexcept {
  const __m512d va_x = _mm512_set1_pd(ax);
  const __m512d va_y = _mm512_set1_pd(ay);
  const __m512d va_z = _mm512_set1_pd(az);

  std::size_t k = 0;
  for (; k + 8 <= block.count; k += 8) {
    __m512d dx = _mm512_sub_pd(_mm512_loadu_pd(block.x + k), va_x);
    __m512d dy = _mm512_sub_pd(_mm512_loadu_pd(block.y + k), va_y);
    __m512d dz = _mm512_sub_pd(_mm512_loadu_pd(block.z + k), va_z);
    __m512d dsq =
        _mm512_fmadd_pd(dx, dx, _mm512_fmadd_pd(dy, dy, _mm512_mul_pd(dz, dz)));
    _mm512_storeu_pd(out_dsq + k, dsq);
  }
  if (k < block.count) {
    __mmask8 mask = static_cast<__mmask8>((1u << (block.count - k)) - 1u);
    __m512d dx = _mm512_sub_pd(_mm512_maskz_loadu_pd(mask, block.x + k), va_x);
    __m512d dy = _mm512_sub_pd(_mm512_maskz_loadu_pd(mask, block.y + k), va_y);
    __m512d dz = _mm512_sub_pd(_mm512_maskz_loadu_pd(mask, block.z + k), va_z);
    __m512d dsq =
        _mm512_fmadd_pd(dx, dx, _mm512_fmadd_pd(dy, dy, _mm512_mul_pd(dz, dz)));
    _mm512_mask_storeu_pd(out_dsq + k, mask, dsq);
  }
}

inline double sinc_integral(double Q,
                            const double *CORRELATION_RESTRICT integrand,
                            const double *CORRELATION_RESTRICT rbins,
                            double *CORRELATION_RESTRICT sinqr_scratch,
                            std::size_t count) noexcept {
  for (std::size_t j = 0; j < count; ++j) {
    sinqr_scratch[j] = std::sin(Q * rbins[j]);
  }
  __m512d vacc = _mm512_setzero_pd();
  std::size_t j = 0;
  for (; j + 8 <= count; j += 8) {
    __m512d vi = _mm512_loadu_pd(integrand + j);
    __m512d vs = _mm512_loadu_pd(sinqr_scratch + j);
    vacc = _mm512_fmadd_pd(vi, vs, vacc);
  }
  double acc = _mm512_reduce_add_pd(vacc);
  for (; j < count; ++j)
    acc += integrand[j] * sinqr_scratch[j];
  return acc;
}

#elif defined(CORRELATION_SIMD_AVX2)

inline void compute_dsq_block(double ax, double ay, double az,
                              const PositionBlock &block,
                              double *CORRELATION_RESTRICT out_dsq) noexcept {
  const __m256d va_x = _mm256_set1_pd(ax);
  const __m256d va_y = _mm256_set1_pd(ay);
  const __m256d va_z = _mm256_set1_pd(az);

  std::size_t k = 0;
  for (; k + 4 <= block.count; k += 4) {
    __m256d dx = _mm256_sub_pd(_mm256_loadu_pd(block.x + k), va_x);
    __m256d dy = _mm256_sub_pd(_mm256_loadu_pd(block.y + k), va_y);
    __m256d dz = _mm256_sub_pd(_mm256_loadu_pd(block.z + k), va_z);
#if defined(__FMA__)
    __m256d dsq =
        _mm256_fmadd_pd(dx, dx, _mm256_fmadd_pd(dy, dy, _mm256_mul_pd(dz, dz)));
#else
    __m256d dsq = _mm256_add_pd(
        _mm256_mul_pd(dx, dx),
        _mm256_add_pd(_mm256_mul_pd(dy, dy), _mm256_mul_pd(dz, dz)));
#endif
    _mm256_storeu_pd(out_dsq + k, dsq);
  }
  for (; k < block.count; ++k) {
    out_dsq[k] = dist_sq_scalar(ax, ay, az, block.x[k], block.y[k], block.z[k]);
  }
}

inline double sinc_integral(double Q,
                            const double *CORRELATION_RESTRICT integrand,
                            const double *CORRELATION_RESTRICT rbins,
                            double *CORRELATION_RESTRICT sinqr_scratch,
                            std::size_t count) noexcept {
  for (std::size_t j = 0; j < count; ++j) {
    sinqr_scratch[j] = std::sin(Q * rbins[j]);
  }
  __m256d vacc = _mm256_setzero_pd();
  std::size_t j = 0;
  for (; j + 4 <= count; j += 4) {
    __m256d vi = _mm256_loadu_pd(integrand + j);
    __m256d vs = _mm256_loadu_pd(sinqr_scratch + j);
#if defined(__FMA__)
    vacc = _mm256_fmadd_pd(vi, vs, vacc);
#else
    vacc = _mm256_add_pd(vacc, _mm256_mul_pd(vi, vs));
#endif
  }
  __m128d lo = _mm256_castpd256_pd128(vacc);
  __m128d hi = _mm256_extractf128_pd(vacc, 1);
  __m128d sum2 = _mm_add_pd(lo, hi);
  __m128d sum1 = _mm_hadd_pd(sum2, sum2);
  double acc = _mm_cvtsd_f64(sum1);
  for (; j < count; ++j)
    acc += integrand[j] * sinqr_scratch[j];
  return acc;
}

#else // Scalar fallback

inline void compute_dsq_block(double ax, double ay, double az,
                              const PositionBlock &block,
                              double *CORRELATION_RESTRICT out_dsq) noexcept {
  for (std::size_t k = 0; k < block.count; ++k) {
    out_dsq[k] = dist_sq_scalar(ax, ay, az, block.x[k], block.y[k], block.z[k]);
  }
}

inline double sinc_integral(double Q,
                            const double *CORRELATION_RESTRICT integrand,
                            const double *CORRELATION_RESTRICT rbins,
                            double *CORRELATION_RESTRICT /*sinqr_scratch*/,
                            std::size_t count) noexcept {
  double acc = 0.0;
  for (std::size_t j = 0; j < count; ++j) {
    acc += integrand[j] * std::sin(Q * rbins[j]);
  }
  return acc;
}

#endif

// ---------------------------------------------------------------------------
// Additional Kernels (Debye sum, normalization, etc.)
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// debye_sum kernel
// ---------------------------------------------------------------------------

/**
 * @brief Computes the Debye scattering equation sum.
 *
 * @param Q The scattering vector magnitude.
 * @param distances Source array of distance values.
 * @param scratch Scratchpad array for intermediate @f$ \text{sinc}(Qr) @f$ evaluations.
 * @param count Number of distances to process.
 * @return The computed Debye sum.
 */
#if defined(CORRELATION_SIMD_AVX512)

inline double debye_sum(double Q, const double *CORRELATION_RESTRICT distances,
                        double *CORRELATION_RESTRICT scratch,
                        std::size_t count) noexcept {
  if (Q < 1e-9)
    return static_cast<double>(count);

  // Pre-pass: compute sinc(Q*r) = sin(Q*r)/(Q*r) for every entry.
  // The Taylor guard (x < 1e-4 → 1 - x²/6) prevents NaN for r = 0
  // or near-zero distances before any vector arithmetic touches the data.
  for (std::size_t j = 0; j < count; ++j) {
    const double x = Q * distances[j];
    scratch[j] = (x < 1e-4) ? (1.0 - x * x / 6.0) : (std::sin(x) / x);
  }

  // SIMD body simply reduces scratch[] — no division, no NaN risk.
  __m512d vacc = _mm512_setzero_pd();
  std::size_t j = 0;
  for (; j + 8 <= count; j += 8)
    vacc = _mm512_add_pd(vacc, _mm512_loadu_pd(scratch + j));
  double acc = _mm512_reduce_add_pd(vacc);
  for (; j < count; ++j)
    acc += scratch[j];
  return acc;
}

#elif defined(CORRELATION_SIMD_AVX2)

inline double debye_sum(double Q, const double *CORRELATION_RESTRICT distances,
                        double *CORRELATION_RESTRICT scratch,
                        std::size_t count) noexcept {
  if (Q < 1e-9)
    return static_cast<double>(count);

  // Pre-pass: compute sinc(Q*r) = sin(Q*r)/(Q*r) for every entry.
  // The Taylor guard (x < 1e-4 → 1 - x²/6) prevents NaN for r = 0
  // or near-zero distances before any vector arithmetic touches the data.
  for (std::size_t j = 0; j < count; ++j) {
    const double x = Q * distances[j];
    scratch[j] = (x < 1e-4) ? (1.0 - x * x / 6.0) : (std::sin(x) / x);
  }

  // SIMD body simply reduces scratch[] — no division, no NaN risk.
  __m256d vacc = _mm256_setzero_pd();
  std::size_t j = 0;
  for (; j + 4 <= count; j += 4)
    vacc = _mm256_add_pd(vacc, _mm256_loadu_pd(scratch + j));
  __m128d lo = _mm256_castpd256_pd128(vacc);
  __m128d hi = _mm256_extractf128_pd(vacc, 1);
  __m128d sum2 = _mm_add_pd(lo, hi);
  __m128d sum1 = _mm_hadd_pd(sum2, sum2);
  double acc = _mm_cvtsd_f64(sum1);
  for (; j < count; ++j)
    acc += scratch[j];
  return acc;
}

#else

inline double debye_sum(double Q, const double *CORRELATION_RESTRICT distances,
                        double *CORRELATION_RESTRICT /*scratch*/,
                        std::size_t count) noexcept {
  if (Q < 1e-9)
    return static_cast<double>(count);
  double acc = 0.0;
  for (std::size_t j = 0; j < count; ++j) {
    const double x = Q * distances[j];
    if (x < 1e-4)
      acc += 1.0 - (x * x) / 6.0;
    else
      acc += std::sin(x) / x;
  }
  return acc;
}

#endif

// ---------------------------------------------------------------------------
// normalize_rdf_bins kernel
// ---------------------------------------------------------------------------

/**
 * @brief Normalizes Radial Distribution Function (RDF) bins.
 * 
 * @param H Histogram array of raw bin counts.
 * @param rbins Array of radial bin positions.
 * @param g_norm Normalization factor for g(r).
 * @param inv_Ni_dr Inverse scale factor times dr.
 * @param inv_Nj_dr Inverse scale factor times dr.
 * @param pi4_rho_j Density scaling factor (4 * pi * rho_j).
 * @param g_out Output array for g(r).
 * @param G_out Output array for G(r).
 * @param J_out Output array for J(r).
 * @param Jinv_out Output array for J^{-1}(r).
 * @param count Total number of bins.
 */
#if defined(CORRELATION_SIMD_AVX512)

inline void normalize_rdf_bins(const double *CORRELATION_RESTRICT H,
                               const double *CORRELATION_RESTRICT rbins,
                               double g_norm, double inv_Ni_dr,
                               double inv_Nj_dr, double pi4_rho_j,
                               double *g_out, double *G_out, double *J_out,
                               double *Jinv_out, std::size_t count) noexcept {
  const __m512d vg_norm = _mm512_set1_pd(g_norm);
  const __m512d v1 = _mm512_set1_pd(1.0);
  const __m512d vinNidr = _mm512_set1_pd(inv_Ni_dr);
  const __m512d vinNjdr = _mm512_set1_pd(inv_Nj_dr);
  const __m512d vpi4rho = _mm512_set1_pd(pi4_rho_j);

  std::size_t k = 1;
  for (; k + 8 <= count; k += 8) {
    __m512d vH = _mm512_loadu_pd(H + k);
    __m512d vr = _mm512_loadu_pd(rbins + k);
    __m512d vr2 = _mm512_mul_pd(vr, vr);
    __m512d vg = _mm512_div_pd(_mm512_mul_pd(vH, vg_norm), vr2);
    _mm512_storeu_pd(g_out + k, vg);
    _mm512_storeu_pd(
        G_out + k,
        _mm512_mul_pd(vpi4rho, _mm512_mul_pd(vr, _mm512_sub_pd(vg, v1))));
    _mm512_storeu_pd(J_out + k, _mm512_mul_pd(vH, vinNidr));
    _mm512_storeu_pd(Jinv_out + k, _mm512_mul_pd(vH, vinNjdr));
  }
  for (; k < count; ++k) {
    const double r = rbins[k];
    if (r < 1e-9)
      continue;
    const double g = H[k] * g_norm / (r * r);
    g_out[k] = g;
    G_out[k] = pi4_rho_j * r * (g - 1.0);
    J_out[k] = H[k] * inv_Ni_dr;
    Jinv_out[k] = H[k] * inv_Nj_dr;
  }
}

#elif defined(CORRELATION_SIMD_AVX2)

inline void normalize_rdf_bins(const double *CORRELATION_RESTRICT H,
                               const double *CORRELATION_RESTRICT rbins,
                               double g_norm, double inv_Ni_dr,
                               double inv_Nj_dr, double pi4_rho_j,
                               double *g_out, double *G_out, double *J_out,
                               double *Jinv_out, std::size_t count) noexcept {
  const __m256d vg_norm = _mm256_set1_pd(g_norm);
  const __m256d v1 = _mm256_set1_pd(1.0);
  const __m256d vinNidr = _mm256_set1_pd(inv_Ni_dr);
  const __m256d vinNjdr = _mm256_set1_pd(inv_Nj_dr);
  const __m256d vpi4rho = _mm256_set1_pd(pi4_rho_j);

  std::size_t k = 1;
  for (; k + 4 <= count; k += 4) {
    __m256d vH = _mm256_loadu_pd(H + k);
    __m256d vr = _mm256_loadu_pd(rbins + k);
    __m256d vr2 = _mm256_mul_pd(vr, vr);
    __m256d vg = _mm256_div_pd(_mm256_mul_pd(vH, vg_norm), vr2);
    _mm256_storeu_pd(g_out + k, vg);
    _mm256_storeu_pd(
        G_out + k,
        _mm256_mul_pd(vpi4rho, _mm256_mul_pd(vr, _mm256_sub_pd(vg, v1))));
    _mm256_storeu_pd(J_out + k, _mm256_mul_pd(vH, vinNidr));
    _mm256_storeu_pd(Jinv_out + k, _mm256_mul_pd(vH, vinNjdr));
  }
  for (; k < count; ++k) {
    const double r = rbins[k];
    if (r < 1e-9)
      continue;
    const double g = H[k] * g_norm / (r * r);
    g_out[k] = g;
    G_out[k] = pi4_rho_j * r * (g - 1.0);
    J_out[k] = H[k] * inv_Ni_dr;
    Jinv_out[k] = H[k] * inv_Nj_dr;
  }
}

#else

inline void normalize_rdf_bins(const double *CORRELATION_RESTRICT H,
                               const double *CORRELATION_RESTRICT rbins,
                               double g_norm, double inv_Ni_dr,
                               double inv_Nj_dr, double pi4_rho_j,
                               double *g_out, double *G_out, double *J_out,
                               double *Jinv_out, std::size_t count) noexcept {
  for (std::size_t k = 1; k < count; ++k) {
    const double r = rbins[k];
    if (r < 1e-9)
      continue;
    const double g = H[k] * g_norm / (r * r);
    g_out[k] = g;
    G_out[k] = pi4_rho_j * r * (g - 1.0);
    J_out[k] = H[k] * inv_Ni_dr;
    Jinv_out[k] = H[k] * inv_Nj_dr;
  }
}

#endif

// ---------------------------------------------------------------------------
// scale_bins
// ---------------------------------------------------------------------------

/**
 * @brief Scales an array by a constant scalar factor.
 * 
 * @param arr The array to scale in-place.
 * @param s The scalar multiplication factor.
 * @param count Number of elements to scale.
 */
#if defined(CORRELATION_SIMD_AVX512)
inline void scale_bins(double *arr, double s, std::size_t count) noexcept {
  const __m512d vs = _mm512_set1_pd(s);
  std::size_t k = 0;
  for (; k + 8 <= count; k += 8)
    _mm512_storeu_pd(arr + k, _mm512_mul_pd(_mm512_loadu_pd(arr + k), vs));
  for (; k < count; ++k)
    arr[k] *= s;
}
#elif defined(CORRELATION_SIMD_AVX2)
inline void scale_bins(double *arr, double s, std::size_t count) noexcept {
  const __m256d vs = _mm256_set1_pd(s);
  std::size_t k = 0;
  for (; k + 4 <= count; k += 4)
    _mm256_storeu_pd(arr + k, _mm256_mul_pd(_mm256_loadu_pd(arr + k), vs));
  for (; k < count; ++k)
    arr[k] *= s;
}
#else
inline void scale_bins(double *arr, double s, std::size_t count) noexcept {
  for (std::size_t k = 0; k < count; ++k)
    arr[k] *= s;
}
#endif

// ---------------------------------------------------------------------------
// dot_block
// ---------------------------------------------------------------------------

/**
 * @brief Computes dot products between a reference vector and a block of vectors.
 * 
 * @param v1x x-component of the reference vector.
 * @param v1y y-component of the reference vector.
 * @param v1z z-component of the reference vector.
 * @param v2x Array of x-components for the target vectors.
 * @param v2y Array of y-components for the target vectors.
 * @param v2z Array of z-components for the target vectors.
 * @param out Output array for the scalar dot products.
 * @param count Number of elements.
 */
#if defined(CORRELATION_SIMD_AVX512)
inline void dot_block(double v1x, double v1y, double v1z,
                      const double *CORRELATION_RESTRICT v2x,
                      const double *CORRELATION_RESTRICT v2y,
                      const double *CORRELATION_RESTRICT v2z,
                      double *CORRELATION_RESTRICT out,
                      std::size_t count) noexcept {
  const __m512d vv1x = _mm512_set1_pd(v1x);
  const __m512d vv1y = _mm512_set1_pd(v1y);
  const __m512d vv1z = _mm512_set1_pd(v1z);
  std::size_t k = 0;
  for (; k + 8 <= count; k += 8) {
    __m512d d = _mm512_fmadd_pd(
        vv1x, _mm512_loadu_pd(v2x + k),
        _mm512_fmadd_pd(vv1y, _mm512_loadu_pd(v2y + k),
                        _mm512_mul_pd(vv1z, _mm512_loadu_pd(v2z + k))));
    _mm512_storeu_pd(out + k, d);
  }
  for (; k < count; ++k)
    out[k] = v1x * v2x[k] + v1y * v2y[k] + v1z * v2z[k];
}
#elif defined(CORRELATION_SIMD_AVX2)
inline void dot_block(double v1x, double v1y, double v1z,
                      const double *CORRELATION_RESTRICT v2x,
                      const double *CORRELATION_RESTRICT v2y,
                      const double *CORRELATION_RESTRICT v2z,
                      double *CORRELATION_RESTRICT out,
                      std::size_t count) noexcept {
  const __m256d vv1x = _mm256_set1_pd(v1x);
  const __m256d vv1y = _mm256_set1_pd(v1y);
  const __m256d vv1z = _mm256_set1_pd(v1z);
  std::size_t k = 0;
  for (; k + 4 <= count; k += 4) {
#if defined(__FMA__)
    __m256d d = _mm256_fmadd_pd(
        vv1x, _mm256_loadu_pd(v2x + k),
        _mm256_fmadd_pd(vv1y, _mm256_loadu_pd(v2y + k),
                        _mm256_mul_pd(vv1z, _mm256_loadu_pd(v2z + k))));
#else
    __m256d d = _mm256_add_pd(
        _mm256_mul_pd(vv1x, _mm256_loadu_pd(v2x + k)),
        _mm256_add_pd(_mm256_mul_pd(vv1y, _mm256_loadu_pd(v2y + k)),
                      _mm256_mul_pd(vv1z, _mm256_loadu_pd(v2z + k))));
#endif
    _mm256_storeu_pd(out + k, d);
  }
  for (; k < count; ++k)
    out[k] = v1x * v2x[k] + v1y * v2y[k] + v1z * v2z[k];
}
#else
inline void dot_block(double v1x, double v1y, double v1z,
                      const double *CORRELATION_RESTRICT v2x,
                      const double *CORRELATION_RESTRICT v2y,
                      const double *CORRELATION_RESTRICT v2z,
                      double *CORRELATION_RESTRICT out,
                      std::size_t count) noexcept {
  for (std::size_t k = 0; k < count; ++k)
    out[k] = v1x * v2x[k] + v1y * v2y[k] + v1z * v2z[k];
}
#endif

// ---------------------------------------------------------------------------
// Utility
// ---------------------------------------------------------------------------
/**
 * @brief populates vectors with standard atom block positions for SIMD.
 */
template <typename AtomRange>
inline std::size_t
fill_position_block(const AtomRange &atoms, std::size_t begin, std::size_t end,
                    std::vector<double> &xs, std::vector<double> &ys,
                    std::vector<double> &zs) noexcept {
  const std::size_t count = end - begin;
  xs.resize(count);
  ys.resize(count);
  zs.resize(count);
  for (std::size_t i = 0; i < count; ++i) {
    const auto &pos = atoms[begin + i].position();
    xs[i] = pos.x();
    ys[i] = pos.y();
    zs[i] = pos.z();
  }
  return count;
}

// ---------------------------------------------------------------------------
// complex_exp_sum
// ---------------------------------------------------------------------------
/**
 * @brief Computes the sum of complex exponentials for a given query vector.
 */
inline void complex_exp_sum(double qx, double qy, double qz,
                            const double *CORRELATION_RESTRICT xs,
                            const double *CORRELATION_RESTRICT ys,
                            const double *CORRELATION_RESTRICT zs,
                            std::size_t count, double &cos_sum,
                            double &sin_sum) noexcept {
  cos_sum = 0.0;
  sin_sum = 0.0;
  for (std::size_t j = 0; j < count; ++j) {
    const double phase = qx * xs[j] + qy * ys[j] + qz * zs[j];
    cos_sum += std::cos(phase);
    sin_sum += std::sin(phase);
  }
}

// ---------------------------------------------------------------------------
// miller_phase_sum
// ---------------------------------------------------------------------------

/**
 * @brief Computes the Miller phase sum across a block of angles.
 */
#if defined(CORRELATION_SIMD_AVX512)

inline void miller_phase_sum(const double *CORRELATION_RESTRICT c1,
                             const double *CORRELATION_RESTRICT s1,
                             const double *CORRELATION_RESTRICT c2,
                             const double *CORRELATION_RESTRICT s2,
                             const double *CORRELATION_RESTRICT c3,
                             const double *CORRELATION_RESTRICT s3,
                             std::size_t count, double &c_sum,
                             double &s_sum) noexcept {
  __m512d vc_sum = _mm512_setzero_pd();
  __m512d vs_sum = _mm512_setzero_pd();
  std::size_t i = 0;
  for (; i + 8 <= count; i += 8) {
    __m512d vc1 = _mm512_loadu_pd(c1 + i);
    __m512d vs1 = _mm512_loadu_pd(s1 + i);
    __m512d vc2 = _mm512_loadu_pd(c2 + i);
    __m512d vs2 = _mm512_loadu_pd(s2 + i);
    __m512d vc3 = _mm512_loadu_pd(c3 + i);
    __m512d vs3 = _mm512_loadu_pd(s3 + i);
    __m512d vc12 = _mm512_fmsub_pd(vc1, vc2, _mm512_mul_pd(vs1, vs2));
    __m512d vs12 = _mm512_fmadd_pd(vs1, vc2, _mm512_mul_pd(vc1, vs2));
    __m512d vc123 = _mm512_fmsub_pd(vc12, vc3, _mm512_mul_pd(vs12, vs3));
    __m512d vs123 = _mm512_fmadd_pd(vs12, vc3, _mm512_mul_pd(vc12, vs3));
    vc_sum = _mm512_add_pd(vc_sum, vc123);
    vs_sum = _mm512_add_pd(vs_sum, vs123);
  }
  c_sum = _mm512_reduce_add_pd(vc_sum);
  s_sum = _mm512_reduce_add_pd(vs_sum);
  for (; i < count; ++i) {
    double c12 = c1[i] * c2[i] - s1[i] * s2[i];
    double s12 = s1[i] * c2[i] + c1[i] * s2[i];
    c_sum += c12 * c3[i] - s12 * s3[i];
    s_sum += s12 * c3[i] + c12 * s3[i];
  }
}

#elif defined(CORRELATION_SIMD_AVX2)

inline void miller_phase_sum(const double *CORRELATION_RESTRICT c1,
                             const double *CORRELATION_RESTRICT s1,
                             const double *CORRELATION_RESTRICT c2,
                             const double *CORRELATION_RESTRICT s2,
                             const double *CORRELATION_RESTRICT c3,
                             const double *CORRELATION_RESTRICT s3,
                             std::size_t count, double &c_sum,
                             double &s_sum) noexcept {
  __m256d vc_sum = _mm256_setzero_pd();
  __m256d vs_sum = _mm256_setzero_pd();
  std::size_t i = 0;
  for (; i + 4 <= count; i += 4) {
    __m256d vc1 = _mm256_loadu_pd(c1 + i);
    __m256d vs1 = _mm256_loadu_pd(s1 + i);
    __m256d vc2 = _mm256_loadu_pd(c2 + i);
    __m256d vs2 = _mm256_loadu_pd(s2 + i);
    __m256d vc3 = _mm256_loadu_pd(c3 + i);
    __m256d vs3 = _mm256_loadu_pd(s3 + i);
#if defined(__FMA__)
    __m256d vc12 = _mm256_fmsub_pd(vc1, vc2, _mm256_mul_pd(vs1, vs2));
    __m256d vs12 = _mm256_fmadd_pd(vs1, vc2, _mm256_mul_pd(vc1, vs2));
    __m256d vc123 = _mm256_fmsub_pd(vc12, vc3, _mm256_mul_pd(vs12, vs3));
    __m256d vs123 = _mm256_fmadd_pd(vs12, vc3, _mm256_mul_pd(vc12, vs3));
#else
    __m256d vc12 =
        _mm256_sub_pd(_mm256_mul_pd(vc1, vc2), _mm256_mul_pd(vs1, vs2));
    __m256d vs12 =
        _mm256_add_pd(_mm256_mul_pd(vs1, vc2), _mm256_mul_pd(vc1, vs2));
    __m256d vc123 =
        _mm256_sub_pd(_mm256_mul_pd(vc12, vc3), _mm256_mul_pd(vs12, vs3));
    __m256d vs123 =
        _mm256_add_pd(_mm256_mul_pd(vs12, vc3), _mm256_mul_pd(vc12, vs3));
#endif
    vc_sum = _mm256_add_pd(vc_sum, vc123);
    vs_sum = _mm256_add_pd(vs_sum, vs123);
  }
  __m128d clo = _mm256_castpd256_pd128(vc_sum);
  __m128d chi = _mm256_extractf128_pd(vc_sum, 1);
  __m128d c_sum2 = _mm_add_pd(clo, chi);
  __m128d c_sum1 = _mm_hadd_pd(c_sum2, c_sum2);
  c_sum = _mm_cvtsd_f64(c_sum1);
  __m128d slo = _mm256_castpd256_pd128(vs_sum);
  __m128d shi = _mm256_extractf128_pd(vs_sum, 1);
  __m128d s_sum2 = _mm_add_pd(slo, shi);
  __m128d s_sum1 = _mm_hadd_pd(s_sum2, s_sum2);
  s_sum = _mm_cvtsd_f64(s_sum1);
  for (; i < count; ++i) {
    double c12 = c1[i] * c2[i] - s1[i] * s2[i];
    double s12 = s1[i] * c2[i] + c1[i] * s2[i];
    c_sum += c12 * c3[i] - s12 * s3[i];
    s_sum += s12 * c3[i] + c12 * s3[i];
  }
}

#else

inline void miller_phase_sum(const double *CORRELATION_RESTRICT c1,
                             const double *CORRELATION_RESTRICT s1,
                             const double *CORRELATION_RESTRICT c2,
                             const double *CORRELATION_RESTRICT s2,
                             const double *CORRELATION_RESTRICT c3,
                             const double *CORRELATION_RESTRICT s3,
                             std::size_t count, double &c_sum,
                             double &s_sum) noexcept {
  c_sum = 0.0;
  s_sum = 0.0;
  for (std::size_t i = 0; i < count; ++i) {
    double c12 = c1[i] * c2[i] - s1[i] * s2[i];
    double s12 = s1[i] * c2[i] + c1[i] * s2[i];
    c_sum += c12 * c3[i] - s12 * s3[i];
    s_sum += s12 * c3[i] + c12 * s3[i];
  }
}

#endif

/**
 * @brief Returns a string describing the compiled SIMD instruction set level.
 * @return String representation of the SIMD level ("AVX-512", "AVX2", or "Scalar").
 */
inline const char *simd_level_string() noexcept {
#if defined(CORRELATION_SIMD_AVX512)
  return "AVX-512";
#elif defined(CORRELATION_SIMD_AVX2)
  return "AVX2";
#else
  return "Scalar";
#endif
}

} // namespace correlation::math
