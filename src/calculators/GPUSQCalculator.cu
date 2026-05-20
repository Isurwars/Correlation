/**
 * @file GPUSQCalculator.cu
 * @brief CUDA implementation of GPU-accelerated S(Q) calculation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 *
 * This file is only compiled when BUILD_WITH_CUDA=ON.
 */

#include "calculators/GPUSQCalculator.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "calculators/StructureFactorCalculator.hpp"
#include "math/Constants.hpp"

#include <cmath>
#include <cuda_runtime.h>
#include <map>
#include <stdexcept>
#include <vector>

namespace correlation::calculators {

namespace {
bool registered = CalculatorFactory::instance().registerCalculator(
    std::make_unique<GPUSQCalculator>());
} // namespace

// -------------------------------------------------------------------------
// CUDA kernel: compute rho_cos and rho_sin for a batch of q-vectors.
//
// Each thread handles one q-vector. For each q, it sums:
//   rho_cos[q] = sum_j cos(qx * xj + qy * yj + qz * zj)
//   rho_sin[q] = sum_j sin(qx * xj + qy * yj + qz * zj)
// -------------------------------------------------------------------------

__global__ void sq_kernel(const double *__restrict__ d_x,
                          const double *__restrict__ d_y,
                          const double *__restrict__ d_z, int N,
                          const double *__restrict__ d_qx,
                          const double *__restrict__ d_qy,
                          const double *__restrict__ d_qz, int num_q,
                          double *__restrict__ d_rho_cos,
                          double *__restrict__ d_rho_sin) {
  int qi = blockIdx.x * blockDim.x + threadIdx.x;
  if (qi >= num_q)
    return;

  double qx = d_qx[qi];
  double qy = d_qy[qi];
  double qz = d_qz[qi];

  double cos_sum = 0.0;
  double sin_sum = 0.0;

  for (int j = 0; j < N; ++j) {
    double phase = qx * d_x[j] + qy * d_y[j] + qz * d_z[j];
    double s, c;
    sincos(phase, &s, &c);
    cos_sum += c;
    sin_sum += s;
  }

  d_rho_cos[qi] = cos_sum;
  d_rho_sin[qi] = sin_sum;
}

// -------------------------------------------------------------------------
// Constructor — probes for a CUDA device.
// -------------------------------------------------------------------------
GPUSQCalculator::GPUSQCalculator() {
  int device_count = 0;
  cudaError_t err = cudaGetDeviceCount(&device_count);
  has_gpu_ = (err == cudaSuccess && device_count > 0);
}

// -------------------------------------------------------------------------
// calculateFrame — GPU path or CPU fallback.
// -------------------------------------------------------------------------
void GPUSQCalculator::calculateFrame(
    correlation::analysis::DistributionFunctions &df,
    const correlation::analysis::AnalysisSettings &settings) const {

  // ---------- CPU fallback ----------
  if (!has_gpu_) {
    StructureFactorCalculator cpu_calc;
    cpu_calc.calculateFrame(df, settings);
    return;
  }

  // ---------- GPU path ----------
  if (settings.q_bin_width <= 0 || settings.q_max <= 0) {
    throw std::invalid_argument("Q-space parameters must be positive.");
  }

  const auto &cell = df.cell();
  const auto &atoms = cell.atoms();
  const size_t N = atoms.size();
  if (N == 0)
    return;

  const double q_max = settings.q_max;
  const double q_bin_width = settings.q_bin_width;
  const size_t num_q_bins =
      static_cast<size_t>(std::floor(q_max / q_bin_width));
  if (num_q_bins == 0)
    throw std::invalid_argument("Q_max too small for Q_bin_width.");

  // Reciprocal lattice vectors.
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

  // Generate q-vectors.
  const double q_max_sq = q_max * q_max;
  std::vector<double> h_qx, h_qy, h_qz, h_qmag;
  for (int h = -hmax; h <= hmax; ++h) {
    for (int k = -kmax; k <= kmax; ++k) {
      for (int l = -lmax; l <= lmax; ++l) {
        if (h == 0 && k == 0 && l == 0)
          continue;
        double qx = h * bx_x + k * by_x + l * bz_x;
        double qy = h * bx_y + k * by_y + l * bz_y;
        double qz = h * bx_z + k * by_z + l * bz_z;
        double qmag_sq = qx * qx + qy * qy + qz * qz;
        if (qmag_sq <= q_max_sq) {
          h_qx.push_back(qx);
          h_qy.push_back(qy);
          h_qz.push_back(qz);
          h_qmag.push_back(std::sqrt(qmag_sq));
        }
      }
    }
  }
  if (h_qx.empty())
    return;

  const int num_q = static_cast<int>(h_qx.size());

  // Prepare host position arrays.
  std::vector<double> h_x(N), h_y(N), h_z(N);
  for (size_t j = 0; j < N; ++j) {
    const auto &p = atoms[j].position();
    h_x[j] = p.x();
    h_y[j] = p.y();
    h_z[j] = p.z();
  }

  // Allocate device memory.
  double *d_x, *d_y, *d_z;
  double *d_qx, *d_qy, *d_qz;
  double *d_rho_cos, *d_rho_sin;

  cudaMalloc(&d_x, N * sizeof(double));
  cudaMalloc(&d_y, N * sizeof(double));
  cudaMalloc(&d_z, N * sizeof(double));
  cudaMalloc(&d_qx, num_q * sizeof(double));
  cudaMalloc(&d_qy, num_q * sizeof(double));
  cudaMalloc(&d_qz, num_q * sizeof(double));
  cudaMalloc(&d_rho_cos, num_q * sizeof(double));
  cudaMalloc(&d_rho_sin, num_q * sizeof(double));

  cudaMemcpy(d_x, h_x.data(), N * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, h_y.data(), N * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_z, h_z.data(), N * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_qx, h_qx.data(), num_q * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_qy, h_qy.data(), num_q * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_qz, h_qz.data(), num_q * sizeof(double), cudaMemcpyHostToDevice);

  // Launch kernel.
  const int block_size = 256;
  const int grid_size = (num_q + block_size - 1) / block_size;
  sq_kernel<<<grid_size, block_size>>>(d_x, d_y, d_z, static_cast<int>(N),
                                       d_qx, d_qy, d_qz, num_q,
                                       d_rho_cos, d_rho_sin);
  cudaDeviceSynchronize();

  // Copy results back.
  std::vector<double> rho_cos(num_q), rho_sin(num_q);
  cudaMemcpy(rho_cos.data(), d_rho_cos, num_q * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(rho_sin.data(), d_rho_sin, num_q * sizeof(double),
             cudaMemcpyDeviceToHost);

  // Free device memory.
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_z);
  cudaFree(d_qx);
  cudaFree(d_qy);
  cudaFree(d_qz);
  cudaFree(d_rho_cos);
  cudaFree(d_rho_sin);

  // Bin-average S(Q) = |rho(q)|^2 / N.
  std::vector<double> total_sq(num_q_bins, 0.0);
  std::vector<size_t> total_count(num_q_bins, 0);

  for (int qi = 0; qi < num_q; ++qi) {
    size_t bin = static_cast<size_t>(h_qmag[qi] / q_bin_width);
    if (bin >= num_q_bins)
      continue;
    double sq_val =
        (rho_cos[qi] * rho_cos[qi] + rho_sin[qi] * rho_sin[qi]) /
        static_cast<double>(N);
    total_sq[bin] += sq_val;
    total_count[bin] += 1;
  }

  for (size_t i = 0; i < num_q_bins; ++i) {
    if (total_count[i] > 0)
      total_sq[i] /= static_cast<double>(total_count[i]);
  }

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
  for (size_t i = 0; i < num_q_bins; ++i)
    s_q_hist.bins[i] = (i + 0.5) * q_bin_width;
  s_q_hist.partials["Total"] = std::move(total_sq);

  df.addHistogram("S_q_gpu", std::move(s_q_hist));
}

} // namespace correlation::calculators
