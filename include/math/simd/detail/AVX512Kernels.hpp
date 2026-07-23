/**
 * @file AVX512Kernels.hpp
 * @brief Low-level AVX-512 intrinsic kernel implementations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#ifdef CORRELATION_SIMD_AVX512

#include <cmath>
#include <cstdint>
#include <immintrin.h>

namespace correlation::math::detail::avx512 {

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

inline void normalize_rdf_bins(const RDFNormalizationParams<double> &params) noexcept {
  if (params.count > 0) {
    params.g_out[0] = 0.0;
    params.G_out[0] = 0.0;
    params.J_out[0] = 0.0;
    params.Jinv_out[0] = 0.0;
  }
  const __m512d vg_norm = _mm512_set1_pd(params.g_norm);
  const __m512d v1 = _mm512_set1_pd(1.0);
  const __m512d vinNidr = _mm512_set1_pd(params.inv_Ni_dr);
  const __m512d vinNjdr = _mm512_set1_pd(params.inv_Nj_dr);
  const __m512d vpi4rho = _mm512_set1_pd(params.pi4_rho_j);

  std::size_t idx = 1;
  for (; idx + 8 <= params.count; idx += 8) {
    const __m512d vH = _mm512_loadu_pd(params.hist_data + idx);
    const __m512d vr = _mm512_loadu_pd(params.radial_bins + idx);
    const __m512d vr2 = _mm512_mul_pd(vr, vr);
    const __m512d vg = _mm512_div_pd(_mm512_mul_pd(vH, vg_norm), vr2);
    _mm512_storeu_pd(params.g_out + idx, vg);
    _mm512_storeu_pd(params.G_out + idx, _mm512_mul_pd(vpi4rho, _mm512_mul_pd(vr, _mm512_sub_pd(vg, v1))));
    _mm512_storeu_pd(params.J_out + idx, _mm512_mul_pd(vH, vinNidr));
    _mm512_storeu_pd(params.Jinv_out + idx, _mm512_mul_pd(vH, vinNjdr));
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
  const __m512 vg_norm = _mm512_set1_ps(params.g_norm);
  const __m512 v1 = _mm512_set1_ps(1.0F);
  const __m512 vinNidr = _mm512_set1_ps(params.inv_Ni_dr);
  const __m512 vinNjdr = _mm512_set1_ps(params.inv_Nj_dr);
  const __m512 vpi4rho = _mm512_set1_ps(params.pi4_rho_j);

  std::size_t idx = 1;
  for (; idx + 16 <= params.count; idx += 16) {
    const __m512 vH = _mm512_loadu_ps(params.hist_data + idx);
    const __m512 vr = _mm512_loadu_ps(params.radial_bins + idx);
    const __m512 vr2 = _mm512_mul_ps(vr, vr);
    const __m512 vg = _mm512_div_ps(_mm512_mul_ps(vH, vg_norm), vr2);
    _mm512_storeu_ps(params.g_out + idx, vg);
    _mm512_storeu_ps(params.G_out + idx, _mm512_mul_ps(vpi4rho, _mm512_mul_ps(vr, _mm512_sub_ps(vg, v1))));
    _mm512_storeu_ps(params.J_out + idx, _mm512_mul_ps(vH, vinNidr));
    _mm512_storeu_ps(params.Jinv_out + idx, _mm512_mul_ps(vH, vinNjdr));
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

inline void miller_phase_sum(const MillerPhaseSumParams<double> &params,
                             MillerPhaseSumResult<double> &result) noexcept {
  __m512d vc_sum = _mm512_setzero_pd();
  __m512d vs_sum = _mm512_setzero_pd();
  std::size_t idx = 0;
  const double *CORRELATION_RESTRICT cos1 = params.cos1;
  const double *CORRELATION_RESTRICT sin1 = params.sin1;
  const double *CORRELATION_RESTRICT cos2 = params.cos2;
  const double *CORRELATION_RESTRICT sin2 = params.sin2;
  const double *CORRELATION_RESTRICT cos3 = params.cos3;
  const double *CORRELATION_RESTRICT sin3 = params.sin3;
  const std::size_t count = params.count;

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
  result.cos_sum = _mm512_reduce_add_pd(vc_sum);
  result.sin_sum = _mm512_reduce_add_pd(vs_sum);
  for (; idx < count; ++idx) {
    const double c12 = cos1[idx] * cos2[idx] - sin1[idx] * sin2[idx];
    const double s12 = sin1[idx] * cos2[idx] + cos1[idx] * sin2[idx];
    result.cos_sum += c12 * cos3[idx] - s12 * sin3[idx];
    result.sin_sum += s12 * cos3[idx] + c12 * sin3[idx];
  }
}

inline void miller_phase_sum(const MillerPhaseSumParams<float> &params, MillerPhaseSumResult<float> &result) noexcept {
  __m512 vc_sum = _mm512_setzero_ps();
  __m512 vs_sum = _mm512_setzero_ps();
  std::size_t idx = 0;
  const float *CORRELATION_RESTRICT cos1 = params.cos1;
  const float *CORRELATION_RESTRICT sin1 = params.sin1;
  const float *CORRELATION_RESTRICT cos2 = params.cos2;
  const float *CORRELATION_RESTRICT sin2 = params.sin2;
  const float *CORRELATION_RESTRICT cos3 = params.cos3;
  const float *CORRELATION_RESTRICT sin3 = params.sin3;
  const std::size_t count = params.count;

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
  result.cos_sum = _mm512_reduce_add_ps(vc_sum);
  result.sin_sum = _mm512_reduce_add_ps(vs_sum);
  for (; idx < count; ++idx) {
    const float c12 = cos1[idx] * cos2[idx] - sin1[idx] * sin2[idx];
    const float s12 = sin1[idx] * cos2[idx] + cos1[idx] * sin2[idx];
    result.cos_sum += c12 * cos3[idx] - s12 * sin3[idx];
    result.sin_sum += s12 * cos3[idx] + c12 * sin3[idx];
  }
}

} // namespace correlation::math::detail::avx512

#endif
