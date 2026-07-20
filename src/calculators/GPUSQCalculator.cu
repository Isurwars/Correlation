/**
 * @file GPUSQCalculator.cu
 * @brief CUDA implementation of GPU-accelerated S(Q) calculation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 *
 * This file is only compiled when BUILD_WITH_CUDA=ON.
 */

#include "calculators/CalculatorFactory.hpp"
#include "calculators/GPUSQCalculator.hpp"
#include "calculators/StructureFactorCalculator.hpp"
#include "math/Constants.hpp"

#include <cmath>
#include "core/GPUPortability.hpp"
#include <map>
#include <stdexcept>
#include <vector>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
// NOLINTNEXTLINE(cert-err58-cpp)
const bool registered = CalculatorFactory::registerTypeSafe<GPUSQCalculator>("GPUSQCalculator");

struct QVectorsData {
  std::vector<double> qx;
  std::vector<double> qy;
  std::vector<double> qz;
  std::vector<double> qmag;
};

struct DeviceAtoms {
  const double *__restrict__ x;
  const double *__restrict__ y;
  const double *__restrict__ z;
};

struct DeviceQVectors {
  const double *__restrict__ qx;
  const double *__restrict__ qy;
  const double *__restrict__ qz;
};

struct DeviceResults {
  double *__restrict__ rho_cos;
  double *__restrict__ rho_sin;
};

QVectorsData generateQVectors(const correlation::core::Cell &cell, double q_max) {
  const auto &inv = cell.inverseLatticeVectors();
  const double bx_x = correlation::math::two_pi * inv(0, 0);
  const double bx_y = correlation::math::two_pi * inv(0, 1);
  const double bx_z = correlation::math::two_pi * inv(0, 2);
  const double by_x = correlation::math::two_pi * inv(1, 0);
  const double by_y = correlation::math::two_pi * inv(1, 1);
  const double by_z = correlation::math::two_pi * inv(1, 2);
  const double bz_x = correlation::math::two_pi * inv(2, 0);
  const double bz_y = correlation::math::two_pi * inv(2, 1);
  const double bz_z = correlation::math::two_pi * inv(2, 2);

  const double b1_norm = std::sqrt(bx_x * bx_x + bx_y * bx_y + bx_z * bx_z);
  const double b2_norm = std::sqrt(by_x * by_x + by_y * by_y + by_z * by_z);
  const double b3_norm = std::sqrt(bz_x * bz_x + bz_y * bz_y + bz_z * bz_z);

  const int hmax = (b1_norm > 1e-10) ? static_cast<int>(std::ceil(q_max / b1_norm)) : 0;
  const int kmax = (b2_norm > 1e-10) ? static_cast<int>(std::ceil(q_max / b2_norm)) : 0;
  const int lmax = (b3_norm > 1e-10) ? static_cast<int>(std::ceil(q_max / b3_norm)) : 0;

  const double q_max_sq = q_max * q_max;
  QVectorsData q_data;
  for (int h_idx = -hmax; h_idx <= hmax; ++h_idx) {
    for (int k_idx = -kmax; k_idx <= kmax; ++k_idx) {
      for (int l_idx = -lmax; l_idx <= lmax; ++l_idx) {
        if (h_idx == 0 && k_idx == 0 && l_idx == 0) {
          continue;
        }
        double q_x = h_idx * bx_x + k_idx * by_x + l_idx * bz_x;
        double q_y = h_idx * bx_y + k_idx * by_y + l_idx * bz_y;
        double q_z = h_idx * bx_z + k_idx * by_z + l_idx * bz_z;
        double qmag_sq = q_x * q_x + q_y * q_y + q_z * q_z;
        if (qmag_sq <= q_max_sq) {
          q_data.qx.push_back(q_x);
          q_data.qy.push_back(q_y);
          q_data.qz.push_back(q_z);
          q_data.qmag.push_back(std::sqrt(qmag_sq));
        }
      }
    }
  }
  return q_data;
}

std::vector<double> averageBinnedSQ(const std::vector<double> &rho_cos, size_t num_q_bins,
                                    const std::vector<double> &rho_sin, double q_bin_width,
                                    const std::vector<double> &qmag, size_t num_atoms) {
  std::vector<double> total_sq(num_q_bins, 0.0);
  std::vector<size_t> total_count(num_q_bins, 0);
  const int num_q = static_cast<int>(qmag.size());

  for (int qi = 0; qi < num_q; ++qi) {
    auto bin = static_cast<size_t>(qmag[qi] / q_bin_width);
    if (bin >= num_q_bins) {
      continue;
    }
    double sq_val = (rho_cos[qi] * rho_cos[qi] + rho_sin[qi] * rho_sin[qi]) / static_cast<double>(num_atoms);
    total_sq[bin] += sq_val;
    total_count[bin] += 1;
  }

  for (size_t i = 0; i < num_q_bins; ++i) {
    if (total_count[i] > 0) {
      total_sq[i] /= static_cast<double>(total_count[i]);
    }
  }

  return total_sq;
}

// -------------------------------------------------------------------------
// CUDA kernel: compute rho_cos and rho_sin for a batch of q-vectors.
//
// Each thread handles one q-vector. For each q, it sums:
//   rho_cos[q] = sum_j cos(qx * xj + qy * yj + qz * zj)
//   rho_sin[q] = sum_j sin(qx * xj + qy * yj + qz * zj)
// -------------------------------------------------------------------------

__global__ void sq_kernel(DeviceAtoms atoms, int num_atoms, DeviceQVectors q_vecs, int num_q, DeviceResults results) {
  int q_i = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
  if (q_i >= num_q) {
    return;
  }

  double q_x = q_vecs.qx[q_i];
  double q_y = q_vecs.qy[q_i];
  double q_z = q_vecs.qz[q_i];

  double cos_sum = 0.0;
  double sin_sum = 0.0;

  for (int j = 0; j < num_atoms; ++j) {
    double phase = q_x * atoms.x[j] + q_y * atoms.y[j] + q_z * atoms.z[j];
    double sin_val = 0.0;
    double cos_val = 0.0;
    sincos(phase, &sin_val, &cos_val);
    cos_sum += cos_val;
    sin_sum += sin_val;
  }

  results.rho_cos[q_i] = cos_sum;
  results.rho_sin[q_i] = sin_sum;
}
} // namespace

// -------------------------------------------------------------------------
// Constructor — probes for a CUDA device.
// -------------------------------------------------------------------------
GPUSQCalculator::GPUSQCalculator() {
  int device_count = 0;
  hipError_t err = hipGetDeviceCount(&device_count);
  has_gpu_ = (err == hipSuccess && device_count > 0);
}

// -------------------------------------------------------------------------
// calculateFrame — GPU path or CPU fallback.
// -------------------------------------------------------------------------
void GPUSQCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                     const correlation::analysis::AnalysisSettings &settings) const {

  // ---------- CPU fallback ----------
  if (!has_gpu_) {
    StructureFactorCalculator cpu_calc;
    cpu_calc.calculateFrame(dists, settings);
    return;
  }

  // ---------- GPU path ----------
  if (settings.q_bin_width <= 0 || settings.q_max <= 0) {
    throw std::invalid_argument("Q-space parameters must be positive.");
  }

  const auto &cell = dists.cell();
  const auto &atoms = cell.atoms();
  const size_t num_atoms = atoms.size();
  if (num_atoms == 0) {
    return;
  }

  const double q_max = settings.q_max;
  const double q_bin_width = settings.q_bin_width;
  const auto num_q_bins = static_cast<size_t>(std::floor(q_max / q_bin_width));
  if (num_q_bins == 0) {
    throw std::invalid_argument("Q_max too small for Q_bin_width.");
  }

  // Generate q-vectors.
  QVectorsData q_data = generateQVectors(cell, q_max);
  if (q_data.qx.empty()) {
    return;
  }

  const int num_q = static_cast<int>(q_data.qx.size());

  // Prepare host position arrays.
  std::vector<double> h_x(num_atoms);
  std::vector<double> h_y(num_atoms);
  std::vector<double> h_z(num_atoms);
  for (size_t j = 0; j < num_atoms; ++j) {
    const auto &pos = atoms[j].position();
    h_x[j] = pos.x();
    h_y[j] = pos.y();
    h_z[j] = pos.z();
  }

  // Allocate device memory.
  double *d_x = nullptr;
  double *d_y = nullptr;
  double *d_z = nullptr;
  double *d_qx = nullptr;
  double *d_qy = nullptr;
  double *d_qz = nullptr;
  double *d_rho_cos = nullptr;
  double *d_rho_sin = nullptr;

  hipMalloc(&d_x, num_atoms * sizeof(double));
  hipMalloc(&d_y, num_atoms * sizeof(double));
  hipMalloc(&d_z, num_atoms * sizeof(double));
  hipMalloc(&d_qx, num_q * sizeof(double));
  hipMalloc(&d_qy, num_q * sizeof(double));
  hipMalloc(&d_qz, num_q * sizeof(double));
  hipMalloc(&d_rho_cos, num_q * sizeof(double));
  hipMalloc(&d_rho_sin, num_q * sizeof(double));

  hipMemcpy(d_x, h_x.data(), num_atoms * sizeof(double), hipMemcpyHostToDevice);
  hipMemcpy(d_y, h_y.data(), num_atoms * sizeof(double), hipMemcpyHostToDevice);
  hipMemcpy(d_z, h_z.data(), num_atoms * sizeof(double), hipMemcpyHostToDevice);
  hipMemcpy(d_qx, q_data.qx.data(), num_q * sizeof(double), hipMemcpyHostToDevice);
  hipMemcpy(d_qy, q_data.qy.data(), num_q * sizeof(double), hipMemcpyHostToDevice);
  hipMemcpy(d_qz, q_data.qz.data(), num_q * sizeof(double), hipMemcpyHostToDevice);

  // Launch kernel.
  const int block_size = 256;
  const int grid_size = (num_q + block_size - 1) / block_size;
  DeviceAtoms device_atoms{.x = d_x, .y = d_y, .z = d_z};
  DeviceQVectors device_q_vecs{.qx = d_qx, .qy = d_qy, .qz = d_qz};
  DeviceResults device_results{.rho_cos = d_rho_cos, .rho_sin = d_rho_sin};
  hipLaunchKernelGGL(sq_kernel, grid_size, block_size, 0, 0, device_atoms, static_cast<int>(num_atoms), device_q_vecs, num_q, device_results);
  hipDeviceSynchronize();

  // Copy results back.
  std::vector<double> rho_cos(num_q, 0.0);
  std::vector<double> rho_sin(num_q, 0.0);
  hipMemcpy(rho_cos.data(), d_rho_cos, num_q * sizeof(double), hipMemcpyDeviceToHost);
  hipMemcpy(rho_sin.data(), d_rho_sin, num_q * sizeof(double), hipMemcpyDeviceToHost);

  // Free device memory.
  hipFree(d_x);
  hipFree(d_y);
  hipFree(d_z);
  hipFree(d_qx);
  hipFree(d_qy);
  hipFree(d_qz);
  hipFree(d_rho_cos);
  hipFree(d_rho_sin);

  // Bin-average S(Q) = |rho(q)|^2 / N.
  std::vector<double> total_sq = averageBinnedSQ(rho_cos, num_q_bins, rho_sin, q_bin_width, q_data.qmag, num_atoms);

  // Build histogram.
  correlation::analysis::Histogram s_q_hist;
  s_q_hist.bins.resize(num_q_bins);
  s_q_hist.x_label = "Q";
  s_q_hist.title = "S(Q) — Structure Factor (GPU)";
  s_q_hist.y_label = "S(Q)";
  s_q_hist.x_unit = "Å⁻¹";
  s_q_hist.y_unit = "arbitrary units";
  s_q_hist.description = "Structure Factor S(Q) — GPU accelerated";
  s_q_hist.file_suffix = "_S_gpu";
  for (size_t q_idx = 0; q_idx < num_q_bins; ++q_idx) {
    s_q_hist.bins[q_idx] = static_cast<real_t>((static_cast<double>(q_idx) + 0.5) * q_bin_width);
  }
  s_q_hist.partials["Total"] = std::vector<real_t>(total_sq.begin(), total_sq.end());

  dists.addHistogram("S_q_gpu", std::move(s_q_hist));
}

} // namespace correlation::calculators
