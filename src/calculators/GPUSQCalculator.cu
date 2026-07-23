/**
 * @file GPUSQCalculator.cu
 * @brief CUDA/HIP implementation of GPU-accelerated S(Q) calculation supporting float and double precision.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 *
 * This file is only compiled when BUILD_WITH_CUDA=ON.
 */

#include "calculators/CalculatorFactory.hpp"
#include "calculators/GPUSQCalculator.hpp"
#include "calculators/StructureFactorCalculator.hpp"
#include "core/GPUPortability.hpp"
#include "math/Constants.hpp"
#include "math/Precision.hpp"

#include <cmath>
#include <map>
#include <stdexcept>
#include <vector>

namespace correlation::calculators {

namespace {
// Static registration of the calculator in the factory
const bool registered = CalculatorFactory::registerTypeSafe<GPUSQCalculator>("GPUSQCalculator");

template <typename T> struct QVectorsData {
  std::vector<T> qx;
  std::vector<T> qy;
  std::vector<T> qz;
  std::vector<T> qmag;
};

template <typename T> struct DeviceAtoms {
  const T *__restrict__ x;
  const T *__restrict__ y;
  const T *__restrict__ z;
};

template <typename T> struct DeviceQVectors {
  const T *__restrict__ qx;
  const T *__restrict__ qy;
  const T *__restrict__ qz;
};

template <typename T> struct DeviceResults {
  T *__restrict__ rho_cos;
  T *__restrict__ rho_sin;
};

template <typename T> __device__ __host__ inline void gpu_sincos(T val, T *sin_val, T *cos_val);

template <> __device__ __host__ inline void gpu_sincos<float>(float val, float *sin_val, float *cos_val) {
#if defined(__CUDA_ARCH__)
  __sincosf(val, sin_val, cos_val);
#elif defined(__HIP_DEVICE_COMPILE__)
  ::sincosf(val, sin_val, cos_val);
#elif defined(_GNU_SOURCE) || defined(__USE_GNU)
  ::sincosf(val, sin_val, cos_val);
#else
  *sin_val = std::sin(val);
  *cos_val = std::cos(val);
#endif
}

template <> __device__ __host__ inline void gpu_sincos<double>(double val, double *sin_val, double *cos_val) {
#if defined(__CUDA_ARCH__)
  __sincos(val, sin_val, cos_val);
#elif defined(__HIP_DEVICE_COMPILE__)
  ::sincos(val, sin_val, cos_val);
#elif defined(_GNU_SOURCE) || defined(__USE_GNU)
  ::sincos(val, sin_val, cos_val);
#else
  *sin_val = std::sin(val);
  *cos_val = std::cos(val);
#endif
}

template <typename T> QVectorsData<T> generateQVectors(const correlation::core::Cell &cell, T q_max) {
  const auto &inv = cell.inverseLatticeVectors();
  const T two_pi = static_cast<T>(correlation::math::two_pi);
  const T bx_x = two_pi * static_cast<T>(inv(0, 0));
  const T bx_y = two_pi * static_cast<T>(inv(0, 1));
  const T bx_z = two_pi * static_cast<T>(inv(0, 2));
  const T by_x = two_pi * static_cast<T>(inv(1, 0));
  const T by_y = two_pi * static_cast<T>(inv(1, 1));
  const T by_z = two_pi * static_cast<T>(inv(1, 2));
  const T bz_x = two_pi * static_cast<T>(inv(2, 0));
  const T bz_y = two_pi * static_cast<T>(inv(2, 1));
  const T bz_z = two_pi * static_cast<T>(inv(2, 2));

  const T b1_norm = std::sqrt(bx_x * bx_x + bx_y * bx_y + bx_z * bx_z);
  const T b2_norm = std::sqrt(by_x * by_x + by_y * by_y + by_z * by_z);
  const T b3_norm = std::sqrt(bz_x * bz_x + bz_y * bz_y + bz_z * bz_z);

  const int hmax = (b1_norm > static_cast<T>(1e-10)) ? static_cast<int>(std::ceil(q_max / b1_norm)) : 0;
  const int kmax = (b2_norm > static_cast<T>(1e-10)) ? static_cast<int>(std::ceil(q_max / b2_norm)) : 0;
  const int lmax = (b3_norm > static_cast<T>(1e-10)) ? static_cast<int>(std::ceil(q_max / b3_norm)) : 0;

  const T q_max_sq = q_max * q_max;
  QVectorsData<T> q_data;
  for (int h_idx = -hmax; h_idx <= hmax; ++h_idx) {
    for (int k_idx = -kmax; k_idx <= kmax; ++k_idx) {
      for (int l_idx = -lmax; l_idx <= lmax; ++l_idx) {
        if (h_idx == 0 && k_idx == 0 && l_idx == 0) {
          continue;
        }
        T q_x = static_cast<T>(h_idx) * bx_x + static_cast<T>(k_idx) * by_x + static_cast<T>(l_idx) * bz_x;
        T q_y = static_cast<T>(h_idx) * bx_y + static_cast<T>(k_idx) * by_y + static_cast<T>(l_idx) * bz_y;
        T q_z = static_cast<T>(h_idx) * bx_z + static_cast<T>(k_idx) * by_z + static_cast<T>(l_idx) * bz_z;
        T qmag_sq = q_x * q_x + q_y * q_y + q_z * q_z;
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

template <typename T>
std::vector<real_t> averageBinnedSQ(const std::vector<T> &rho_cos, size_t num_q_bins, const std::vector<T> &rho_sin,
                                    T q_bin_width, const std::vector<T> &qmag, size_t num_atoms) {
  std::vector<real_t> total_sq(num_q_bins, static_cast<real_t>(0.0));
  std::vector<size_t> total_count(num_q_bins, 0);
  const int num_q = static_cast<int>(qmag.size());

  for (int qi = 0; qi < num_q; ++qi) {
    auto bin = static_cast<size_t>(qmag[qi] / q_bin_width);
    if (bin >= num_q_bins) {
      continue;
    }
    T sq_val = (rho_cos[qi] * rho_cos[qi] + rho_sin[qi] * rho_sin[qi]) / static_cast<T>(num_atoms);
    total_sq[bin] += static_cast<real_t>(sq_val);
    total_count[bin] += 1;
  }

  for (size_t i = 0; i < num_q_bins; ++i) {
    if (total_count[i] > 0) {
      total_sq[i] /= static_cast<real_t>(total_count[i]);
    }
  }

  return total_sq;
}

// -------------------------------------------------------------------------
// CUDA kernel: compute rho_cos and rho_sin for a batch of q-vectors.
// -------------------------------------------------------------------------
template <typename T>
__global__ void sq_kernel(DeviceAtoms<T> atoms, int num_atoms, DeviceQVectors<T> q_vecs, int num_q,
                          DeviceResults<T> results) {
  int q_i = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
  if (q_i >= num_q) {
    return;
  }

  T q_x = q_vecs.qx[q_i];
  T q_y = q_vecs.qy[q_i];
  T q_z = q_vecs.qz[q_i];

  T cos_sum = static_cast<T>(0.0);
  T sin_sum = static_cast<T>(0.0);

  for (int j = 0; j < num_atoms; ++j) {
    T phase = q_x * atoms.x[j] + q_y * atoms.y[j] + q_z * atoms.z[j];
    T sin_val = static_cast<T>(0.0);
    T cos_val = static_cast<T>(0.0);
    gpu_sincos(phase, &sin_val, &cos_val);
    cos_sum += cos_val;
    sin_sum += sin_val;
  }

  results.rho_cos[q_i] = cos_sum;
  results.rho_sin[q_i] = sin_sum;
}
} // namespace

GPUSQCalculator::GPUSQCalculator() {
  int device_count = 0;
  hipError_t err = hipGetDeviceCount(&device_count);
  has_gpu_ = (err == hipSuccess && device_count > 0);
}

void GPUSQCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                     const correlation::analysis::AnalysisSettings &settings) const {

  if (!has_gpu_) {
    StructureFactorCalculator cpu_calc;
    cpu_calc.calculateFrame(dists, settings);
    return;
  }

  if (settings.q_bin_width <= 0 || settings.q_max <= 0) {
    throw std::invalid_argument("Q-space parameters must be positive.");
  }

  using T = real_t;

  const auto &cell = dists.cell();
  const auto &atoms = cell.atoms();
  const size_t num_atoms = atoms.size();
  if (num_atoms == 0) {
    return;
  }

  const T q_max = static_cast<T>(settings.q_max);
  const T q_bin_width = static_cast<T>(settings.q_bin_width);
  const auto num_q_bins = static_cast<size_t>(std::floor(q_max / q_bin_width));
  if (num_q_bins == 0) {
    throw std::invalid_argument("Q_max too small for Q_bin_width.");
  }

  QVectorsData<T> q_data = generateQVectors<T>(cell, q_max);
  const size_t num_q = q_data.qmag.size();
  if (num_q == 0) {
    return;
  }

  std::vector<T> h_x(num_atoms);
  std::vector<T> h_y(num_atoms);
  std::vector<T> h_z(num_atoms);

  for (size_t i = 0; i < num_atoms; ++i) {
    h_x[i] = static_cast<T>(atoms[i].position().x());
    h_y[i] = static_cast<T>(atoms[i].position().y());
    h_z[i] = static_cast<T>(atoms[i].position().z());
  }

  T *d_x = nullptr;
  T *d_y = nullptr;
  T *d_z = nullptr;
  T *d_qx = nullptr;
  T *d_qy = nullptr;
  T *d_qz = nullptr;
  T *d_rho_cos = nullptr;
  T *d_rho_sin = nullptr;

  hipMalloc(&d_x, num_atoms * sizeof(T));
  hipMalloc(&d_y, num_atoms * sizeof(T));
  hipMalloc(&d_z, num_atoms * sizeof(T));
  hipMalloc(&d_qx, num_q * sizeof(T));
  hipMalloc(&d_qy, num_q * sizeof(T));
  hipMalloc(&d_qz, num_q * sizeof(T));
  hipMalloc(&d_rho_cos, num_q * sizeof(T));
  hipMalloc(&d_rho_sin, num_q * sizeof(T));

  hipMemcpy(d_x, h_x.data(), num_atoms * sizeof(T), hipMemcpyHostToDevice);
  hipMemcpy(d_y, h_y.data(), num_atoms * sizeof(T), hipMemcpyHostToDevice);
  hipMemcpy(d_z, h_z.data(), num_atoms * sizeof(T), hipMemcpyHostToDevice);
  hipMemcpy(d_qx, q_data.qx.data(), num_q * sizeof(T), hipMemcpyHostToDevice);
  hipMemcpy(d_qy, q_data.qy.data(), num_q * sizeof(T), hipMemcpyHostToDevice);
  hipMemcpy(d_qz, q_data.qz.data(), num_q * sizeof(T), hipMemcpyHostToDevice);

  int block_size = 256;
  int grid_size = (static_cast<int>(num_q) + block_size - 1) / block_size;

  DeviceAtoms<T> dev_atoms{d_x, d_y, d_z};
  DeviceQVectors<T> dev_qvecs{d_qx, d_qy, d_qz};
  DeviceResults<T> dev_results{d_rho_cos, d_rho_sin};

  hipLaunchKernelGGL(sq_kernel<T>, grid_size, block_size, 0, 0, dev_atoms, static_cast<int>(num_atoms), dev_qvecs,
                     static_cast<int>(num_q), dev_results);
  hipDeviceSynchronize();

  std::vector<T> h_rho_cos(num_q);
  std::vector<T> h_rho_sin(num_q);

  hipMemcpy(h_rho_cos.data(), d_rho_cos, num_q * sizeof(T), hipMemcpyDeviceToHost);
  hipMemcpy(h_rho_sin.data(), d_rho_sin, num_q * sizeof(T), hipMemcpyDeviceToHost);

  hipFree(d_x);
  hipFree(d_y);
  hipFree(d_z);
  hipFree(d_qx);
  hipFree(d_qy);
  hipFree(d_qz);
  hipFree(d_rho_cos);
  hipFree(d_rho_sin);

  std::vector<real_t> s_q = averageBinnedSQ<T>(h_rho_cos, num_q_bins, h_rho_sin, q_bin_width, q_data.qmag, num_atoms);

  correlation::analysis::Histogram hist;
  hist.x_label = "Q";
  hist.title = "Structure Factor S(Q) (GPU)";
  hist.y_label = "S(Q)";
  hist.x_unit = "1/Å";
  hist.y_unit = "dimensionless";
  hist.description = "Static Structure Factor S(Q) computed on GPU";
  hist.file_suffix = "_sq_gpu";
  hist.bins.resize(num_q_bins);

  for (size_t i = 0; i < num_q_bins; ++i) {
    hist.bins[i] = (static_cast<real_t>(i) + static_cast<real_t>(0.5)) * q_bin_width;
  }

  hist.partials["Total"] = s_q;
  dists.addHistogram("S_Q_gpu", std::move(hist));
}

} // namespace correlation::calculators
