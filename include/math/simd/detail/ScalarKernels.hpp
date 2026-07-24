/**
 * @file ScalarKernels.hpp
 * @brief Low-level Scalar fallback kernel implementations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/simd/SIMDTypes.hpp"

#include <cmath>

namespace correlation::math::detail::scalar {

inline void compute_dsq_block(float ref_x, float ref_y, float ref_z, const PositionBlockT<float> &block,
                              float *CORRELATION_RESTRICT out_dsq) noexcept {
  for (std::size_t idx = 0; idx < block.count; ++idx) {
    out_dsq[idx] = dist_sq_scalar<float>(
        {
            .x = ref_x,
            .y = ref_y,
            .z = ref_z,
        },
        {
            .x = block.x[idx],
            .y = block.y[idx],
            .z = block.z[idx],
        });
  }
}

inline void compute_dsq_block(double ref_x, double ref_y, double ref_z, const PositionBlockT<double> &block,
                              double *CORRELATION_RESTRICT out_dsq) noexcept {
  for (std::size_t idx = 0; idx < block.count; ++idx) {
    out_dsq[idx] = dist_sq_scalar<double>(
        {
            .x = ref_x,
            .y = ref_y,
            .z = ref_z,
        },
        {
            .x = block.x[idx],
            .y = block.y[idx],
            .z = block.z[idx],
        });
  }
}

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

inline float simd_dot(const float *CORRELATION_RESTRICT input_a, const float *CORRELATION_RESTRICT input_b,
                      std::size_t count) noexcept {
  float acc = 0.0F;
  for (std::size_t idx = 0; idx < count; ++idx) {
    acc += input_a[idx] * input_b[idx];
  }
  return acc;
}

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

  double acc = 0.0;
  double carry = 0.0;
  for (std::size_t idx = 0; idx < count; ++idx) {
    const double val_x = q_magnitude * distances[idx];
    const double val = (val_x < 1.0e-4) ? (1.0 - (val_x * val_x) / 6.0) : (std::sin(val_x) / val_x);
    scratch[idx] = val;
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

  float acc = 0.0F;
  for (std::size_t idx = 0; idx < count; ++idx) {
    const float val_x = q_magnitude * distances[idx];
    const float val = (val_x < 1.0e-4F) ? (1.0F - (val_x * val_x) / 6.0F) : (std::sin(val_x) / val_x);
    scratch[idx] = val;
    acc += val;
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
  for (std::size_t idx = 1; idx < params.count; ++idx) {
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
  for (std::size_t idx = 1; idx < params.count; ++idx) {
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

inline void scale_bins(const ScaleBinsParams<double> &params) noexcept {
  for (std::size_t idx = 0; idx < params.count; ++idx) {
    params.arr[idx] *= params.scale_factor;
  }
}

inline void scale_bins(const ScaleBinsParams<float> &params) noexcept {
  for (std::size_t idx = 0; idx < params.count; ++idx) {
    params.arr[idx] *= params.scale_factor;
  }
}

inline void miller_phase_sum(const MillerPhaseSumParams<double> &params,
                             MillerPhaseSumResult<double> &result) noexcept {
  result.cos_sum = 0.0;
  result.sin_sum = 0.0;
  for (std::size_t idx = 0; idx < params.count; ++idx) {
    const double c12 = params.cos1[idx] * params.cos2[idx] - params.sin1[idx] * params.sin2[idx];
    const double s12 = params.sin1[idx] * params.cos2[idx] + params.cos1[idx] * params.sin2[idx];
    result.cos_sum += c12 * params.cos3[idx] - s12 * params.sin3[idx];
    result.sin_sum += s12 * params.cos3[idx] + c12 * params.sin3[idx];
  }
}

inline void miller_phase_sum(const MillerPhaseSumParams<float> &params, MillerPhaseSumResult<float> &result) noexcept {
  result.cos_sum = 0.0F;
  result.sin_sum = 0.0F;
  for (std::size_t idx = 0; idx < params.count; ++idx) {
    const float c12 = params.cos1[idx] * params.cos2[idx] - params.sin1[idx] * params.sin2[idx];
    const float s12 = params.sin1[idx] * params.cos2[idx] + params.cos1[idx] * params.sin2[idx];
    result.cos_sum += c12 * params.cos3[idx] - s12 * params.sin3[idx];
    result.sin_sum += s12 * params.cos3[idx] + c12 * params.sin3[idx];
  }
}

} // namespace correlation::math::detail::scalar
