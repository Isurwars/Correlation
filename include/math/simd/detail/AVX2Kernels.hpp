/**
 * @file AVX2Kernels.hpp
 * @brief Low-level AVX2 intrinsic kernel implementations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/simd/SIMDTypes.hpp"

#if defined(CORRELATION_SIMD_AVX2) && !defined(CORRELATION_SIMD_AVX512)

#include <immintrin.h>
#include <array>
#include <cmath>

namespace correlation::math::detail::avx2 {

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
  alignas(32) std::array<float, 8> float_buf{};
  _mm256_storeu_ps(float_buf.data(), vacc);
  float acc = float_buf[0] + float_buf[1] + float_buf[2] + float_buf[3] + float_buf[4] + float_buf[5] + float_buf[6] +
              float_buf[7];
  for (; idx < count; ++idx) {
    acc += input_a[idx] * input_b[idx];
  }
  return acc;
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
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

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
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

inline double sinc_integral(const SincIntegralParams<double> &params) noexcept {
  for (std::size_t idx = 0; idx < params.count; ++idx) {
    params.sinqr_scratch[idx] = std::sin(params.q_magnitude * params.radial_bins[idx]);
  }
  return simd_dot(params.integrand, params.sinqr_scratch, params.count);
}

inline float sinc_integral(const SincIntegralParams<float> &params) noexcept {
  for (std::size_t idx = 0; idx < params.count; ++idx) {
    params.sinqr_scratch[idx] = std::sin(params.q_magnitude * params.radial_bins[idx]);
  }
  return simd_dot(params.integrand, params.sinqr_scratch, params.count);
}

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
  alignas(32) std::array<float, 8> float_buf{};
  _mm256_storeu_ps(float_buf.data(), vacc);
  float acc = float_buf[0] + float_buf[1] + float_buf[2] + float_buf[3] + float_buf[4] + float_buf[5] + float_buf[6] +
              float_buf[7];
  for (; idx < count; ++idx) {
    acc += scratch[idx];
  }
  return acc;
}

inline void normalize_rdf_bins(const RDFNormalizationParams<double> &params) noexcept {
  if (params.count > 0) {
    params.g_out[0] = 0.0;
    params.G_out[0] = 0.0;
    params.J_out[0] = 0.0;
    params.Jinv_out[0] = 0.0;
  }
  const __m256d vg_norm = _mm256_set1_pd(params.g_norm);
  const __m256d v_1 = _mm256_set1_pd(1.0);
  const __m256d vinNidr = _mm256_set1_pd(params.inv_Ni_dr);
  const __m256d vinNjdr = _mm256_set1_pd(params.inv_Nj_dr);
  const __m256d vpi4rho = _mm256_set1_pd(params.pi4_rho_j);

  std::size_t idx = 1;
  for (; idx + 4 <= params.count; idx += 4) {
    const __m256d v_H = _mm256_loadu_pd(params.hist_data + idx);
    const __m256d v_r = _mm256_loadu_pd(params.radial_bins + idx);
    const __m256d vr2 = _mm256_mul_pd(v_r, v_r);
    const __m256d v_g = _mm256_div_pd(_mm256_mul_pd(v_H, vg_norm), vr2);
    _mm256_storeu_pd(params.g_out + idx, v_g);
    _mm256_storeu_pd(params.G_out + idx, _mm256_mul_pd(vpi4rho, _mm256_mul_pd(v_r, _mm256_sub_pd(v_g, v_1))));
    _mm256_storeu_pd(params.J_out + idx, _mm256_mul_pd(v_H, vinNidr));
    _mm256_storeu_pd(params.Jinv_out + idx, _mm256_mul_pd(v_H, vinNjdr));
  }
  for (; idx < params.count; ++idx) {
    const double r_val = params.radial_bins[idx];
    if (r_val < 1.0e-9) {
      params.g_out[idx] = 0.0;
      params.G_out[idx] = 0.0;
      params.J_out[idx] = 0.0;
      params.Jinv_out[idx] = 0.0;
      continue;
    }
    const double g_val = params.hist_data[idx] * params.g_norm / (r_val * r_val);
    params.g_out[idx] = g_val;
    params.G_out[idx] = params.pi4_rho_j * r_val * (g_val - 1.0);
    params.J_out[idx] = params.hist_data[idx] * params.inv_Ni_dr;
    params.Jinv_out[idx] = params.hist_data[idx] * params.inv_Nj_dr;
  }
}

inline void normalize_rdf_bins(const RDFNormalizationParams<float> &params) noexcept {
  if (params.count > 0) {
    params.g_out[0] = 0.0F;
    params.G_out[0] = 0.0F;
    params.J_out[0] = 0.0F;
    params.Jinv_out[0] = 0.0F;
  }
  const __m256 vg_norm = _mm256_set1_ps(params.g_norm);
  const __m256 v_1 = _mm256_set1_ps(1.0F);
  const __m256 vinNidr = _mm256_set1_ps(params.inv_Ni_dr);
  const __m256 vinNjdr = _mm256_set1_ps(params.inv_Nj_dr);
  const __m256 vpi4rho = _mm256_set1_ps(params.pi4_rho_j);

  std::size_t idx = 1;
  for (; idx + 8 <= params.count; idx += 8) {
    const __m256 v_H = _mm256_loadu_ps(params.hist_data + idx);
    const __m256 v_r = _mm256_loadu_ps(params.radial_bins + idx);
    const __m256 vr2 = _mm256_mul_ps(v_r, v_r);
    const __m256 v_g = _mm256_div_ps(_mm256_mul_ps(v_H, vg_norm), vr2);
    _mm256_storeu_ps(params.g_out + idx, v_g);
    _mm256_storeu_ps(params.G_out + idx, _mm256_mul_ps(vpi4rho, _mm256_mul_ps(v_r, _mm256_sub_ps(v_g, v_1))));
    _mm256_storeu_ps(params.J_out + idx, _mm256_mul_ps(v_H, vinNidr));
    _mm256_storeu_ps(params.Jinv_out + idx, _mm256_mul_ps(v_H, vinNjdr));
  }
  for (; idx < params.count; ++idx) {
    const float r_val = params.radial_bins[idx];
    if (r_val < 1.0e-9F) {
      params.g_out[idx] = 0.0F;
      params.G_out[idx] = 0.0F;
      params.J_out[idx] = 0.0F;
      params.Jinv_out[idx] = 0.0F;
      continue;
    }
    const float g_val = params.hist_data[idx] * params.g_norm / (r_val * r_val);
    params.g_out[idx] = g_val;
    params.G_out[idx] = params.pi4_rho_j * r_val * (g_val - 1.0F);
    params.J_out[idx] = params.hist_data[idx] * params.inv_Ni_dr;
    params.Jinv_out[idx] = params.hist_data[idx] * params.inv_Nj_dr;
  }
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
inline void scale_bins(double *arr, double scale_factor, std::size_t count) noexcept {
  const __m256d v_scale = _mm256_set1_pd(scale_factor);
  std::size_t idx = 0;
  for (; idx + 4 <= count; idx += 4) {
    _mm256_storeu_pd(arr + idx, _mm256_mul_pd(_mm256_loadu_pd(arr + idx), v_scale));
  }
  for (; idx < count; ++idx) {
    arr[idx] *= scale_factor;
  }
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
inline void scale_bins(float *arr, float scale_factor, std::size_t count) noexcept {
  const __m256 v_scale = _mm256_set1_ps(scale_factor);
  std::size_t idx = 0;
  for (; idx + 8 <= count; idx += 8) {
    _mm256_storeu_ps(arr + idx, _mm256_mul_ps(_mm256_loadu_ps(arr + idx), v_scale));
  }
  for (; idx < count; ++idx) {
    arr[idx] *= scale_factor;
  }
}

inline void miller_phase_sum(const MillerPhaseSumParams<double> &params, MillerPhaseSumResult<double> &result) noexcept {
  __m256d vc_sum = _mm256_setzero_pd();
  __m256d vs_sum = _mm256_setzero_pd();
  std::size_t idx = 0;
  const double *CORRELATION_RESTRICT cos1 = params.cos1;
  const double *CORRELATION_RESTRICT sin1 = params.sin1;
  const double *CORRELATION_RESTRICT cos2 = params.cos2;
  const double *CORRELATION_RESTRICT sin2 = params.sin2;
  const double *CORRELATION_RESTRICT cos3 = params.cos3;
  const double *CORRELATION_RESTRICT sin3 = params.sin3;
  const std::size_t count = params.count;

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
  result.cos_sum = _mm_cvtsd_f64(c_sum1);
  const __m128d slo = _mm256_castpd256_pd128(vs_sum);
  const __m128d shi = _mm256_extractf128_pd(vs_sum, 1);
  const __m128d s_sum2 = _mm_add_pd(slo, shi);
  const __m128d s_sum1 = _mm_hadd_pd(s_sum2, s_sum2);
  result.sin_sum = _mm_cvtsd_f64(s_sum1);
  for (; idx < count; ++idx) {
    const double c12 = cos1[idx] * cos2[idx] - sin1[idx] * sin2[idx];
    const double s12 = sin1[idx] * cos2[idx] + cos1[idx] * sin2[idx];
    result.cos_sum += c12 * cos3[idx] - s12 * sin3[idx];
    result.sin_sum += s12 * cos3[idx] + c12 * sin3[idx];
  }
}

inline void miller_phase_sum(const MillerPhaseSumParams<float> &params, MillerPhaseSumResult<float> &result) noexcept {
  __m256 vc_sum = _mm256_setzero_ps();
  __m256 vs_sum = _mm256_setzero_ps();
  std::size_t idx = 0;
  const float *CORRELATION_RESTRICT cos1 = params.cos1;
  const float *CORRELATION_RESTRICT sin1 = params.sin1;
  const float *CORRELATION_RESTRICT cos2 = params.cos2;
  const float *CORRELATION_RESTRICT sin2 = params.sin2;
  const float *CORRELATION_RESTRICT cos3 = params.cos3;
  const float *CORRELATION_RESTRICT sin3 = params.sin3;
  const std::size_t count = params.count;

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
  alignas(32) std::array<float, 8> c_buf{};
  alignas(32) std::array<float, 8> s_buf{};
  _mm256_storeu_ps(c_buf.data(), vc_sum);
  _mm256_storeu_ps(s_buf.data(), vs_sum);
  result.cos_sum = c_buf[0] + c_buf[1] + c_buf[2] + c_buf[3] + c_buf[4] + c_buf[5] + c_buf[6] + c_buf[7];
  result.sin_sum = s_buf[0] + s_buf[1] + s_buf[2] + s_buf[3] + s_buf[4] + s_buf[5] + s_buf[6] + s_buf[7];
  for (; idx < count; ++idx) {
    const float c12 = cos1[idx] * cos2[idx] - sin1[idx] * sin2[idx];
    const float s12 = sin1[idx] * cos2[idx] + cos1[idx] * sin2[idx];
    result.cos_sum += c12 * cos3[idx] - s12 * sin3[idx];
    result.sin_sum += s12 * cos3[idx] + c12 * sin3[idx];
  }
}

} // namespace correlation::math::detail::avx2

#endif
