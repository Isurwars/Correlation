/**
 * @file GPUSteinhardtCalculator.cu
 * @brief CUDA/HIP implementation of GPU-accelerated Steinhardt parameter calculation supporting float and double.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/CalculatorFactory.hpp"
#include "calculators/GPUSteinhardtCalculator.hpp"
#include "calculators/SteinhardtCalculator.hpp"
#include "core/GPUPortability.hpp"
#include "math/Constants.hpp"
#include "math/Precision.hpp"

#include <array>
#include <cmath>
#include <vector>

namespace correlation::calculators {

namespace {

const bool registered = CalculatorFactory::registerTypeSafe<GPUSteinhardtCalculator>("GPUSteinhardtCalculator");

template <typename T> struct GPUPoint {
  T x;
  T y;
  T z;
};

template <typename T> struct SphericalHarmonicInput {
  T costheta;
  T phi;
};

template <typename T> struct SphericalHarmonicOutput {
  T *real_y;
  T *imag_y;
};

template <typename T> struct SteinhardtQ4Harmonics {
  std::array<T, 9> real{};
  std::array<T, 9> imag{};
};

template <typename T> struct SteinhardtQ6Harmonics {
  std::array<T, 13> real{};
  std::array<T, 13> imag{};
};

template <typename T> struct SteinhardtResult {
  T q4{static_cast<T>(0.0)};
  T q6{static_cast<T>(0.0)};
};

struct NeighborGraphPointers {
  const int *__restrict__ offsets;
  const int *__restrict__ indices;
};

template <typename T> struct SteinhardtOutputPointers {
  T *__restrict__ q4_out;
  T *__restrict__ q6_out;
};

// -------------------------------------------------------------------------
// Device helper: spherical harmonic Y_l^m for l=4 and l=6 (templated on T)
// -------------------------------------------------------------------------
template <typename T> __device__ void compute_y4m(SphericalHarmonicInput<T> input, SphericalHarmonicOutput<T> output) {
  T const costheta = input.costheta;
  T const phi = input.phi;
  T *real_y = output.real_y;
  T *imag_y = output.imag_y;

  T const sin2 = static_cast<T>(1.0) - costheta * costheta;
  T const sintheta = (sin2 > static_cast<T>(0.0)) ? sqrt(sin2) : static_cast<T>(0.0);

  T const p40 = static_cast<T>(0.125) * (static_cast<T>(35.0) * costheta * costheta * costheta * costheta -
                                         static_cast<T>(30.0) * costheta * costheta + static_cast<T>(3.0));
  T const p41 =
      static_cast<T>(-2.5) * costheta * (static_cast<T>(7.0) * costheta * costheta - static_cast<T>(3.0)) * sintheta;
  T const p42 = static_cast<T>(7.5) * (static_cast<T>(7.0) * costheta * costheta - static_cast<T>(1.0)) * sin2;
  T const p43 = static_cast<T>(-105.0) * costheta * sintheta * sin2;
  T const p44 = static_cast<T>(105.0) * sin2 * sin2;

  T const norm0 = static_cast<T>(0.1057855469369004);
  T const norm1 = static_cast<T>(0.0528927734684502);
  T const norm2 = static_cast<T>(0.0152687588147426);
  T const norm3 = static_cast<T>(0.0044077109282582);
  T const norm4 = static_cast<T>(0.0015583606622822);

  std::array<T, 5> const poly = {p40 * norm0, p41 * norm1, p42 * norm2, p43 * norm3, p44 * norm4};

  for (int m_idx = 0; m_idx <= 4; ++m_idx) {
    T cos_mphi = cos(static_cast<T>(m_idx) * phi);
    T sin_mphi = sin(static_cast<T>(m_idx) * phi);

    *(real_y + m_idx + 4) = *(poly.data() + m_idx) * cos_mphi;
    *(imag_y + m_idx + 4) = *(poly.data() + m_idx) * sin_mphi;

    if (m_idx > 0) {
      T const sign = (m_idx % 2 == 0) ? static_cast<T>(1.0) : static_cast<T>(-1.0);
      *(real_y - m_idx + 4) = sign * *(poly.data() + m_idx) * cos_mphi;
      *(imag_y - m_idx + 4) = -sign * *(poly.data() + m_idx) * sin_mphi;
    }
  }
}

template <typename T> __device__ void compute_y6m(SphericalHarmonicInput<T> input, SphericalHarmonicOutput<T> output) {
  T const costheta = input.costheta;
  T const phi = input.phi;
  T *real_y = output.real_y;
  T *imag_y = output.imag_y;

  T const sin2 = static_cast<T>(1.0) - costheta * costheta;
  T const sintheta = (sin2 > static_cast<T>(0.0)) ? sqrt(sin2) : static_cast<T>(0.0);

  T const cos2 = costheta * costheta;
  T const cos4 = cos2 * cos2;
  T const cos6 = cos4 * cos2;

  T const p60 = static_cast<T>(0.0625) * (static_cast<T>(231.0) * cos6 - static_cast<T>(315.0) * cos4 +
                                          static_cast<T>(105.0) * cos2 - static_cast<T>(5.0));
  T const p61 = static_cast<T>(-2.625) * costheta *
                (static_cast<T>(33.0) * cos4 - static_cast<T>(30.0) * cos2 + static_cast<T>(5.0)) * sintheta;
  T const p62 =
      static_cast<T>(13.125) * (static_cast<T>(33.0) * cos4 - static_cast<T>(18.0) * cos2 + static_cast<T>(1.0)) * sin2;
  T const p63 =
      static_cast<T>(-492.1875) * costheta * (static_cast<T>(11.0) * cos2 - static_cast<T>(3.0)) * sintheta * sin2;
  T const p64 = static_cast<T>(492.1875) * (static_cast<T>(11.0) * cos2 - static_cast<T>(1.0)) * sin2 * sin2;
  T const p65 = static_cast<T>(-10335.9375) * costheta * sintheta * sin2 * sin2;
  T const p66 = static_cast<T>(10335.9375) * sin2 * sin2 * sin2;

  T const norm0 = static_cast<T>(0.0906208035179047);
  T const norm1 = static_cast<T>(0.0453104017589523);
  T const norm2 = static_cast<T>(0.0101317187214731);
  T const norm3 = static_cast<T>(0.0022655200879476);
  T const norm4 = static_cast<T>(0.0006540003056029);
  T const norm5 = static_cast<T>(0.0002480962451000);
  T const norm6 = static_cast<T>(0.0001012848529241);

  std::array<T, 7> const poly = {p60 * norm0, p61 * norm1, p62 * norm2, p63 * norm3,
                                 p64 * norm4, p65 * norm5, p66 * norm6};

  for (int m_idx = 0; m_idx <= 6; ++m_idx) {
    T cos_mphi = cos(static_cast<T>(m_idx) * phi);
    T sin_mphi = sin(static_cast<T>(m_idx) * phi);

    *(real_y + m_idx + 6) = *(poly.data() + m_idx) * cos_mphi;
    *(imag_y + m_idx + 6) = *(poly.data() + m_idx) * sin_mphi;

    if (m_idx > 0) {
      T const sign = (m_idx % 2 == 0) ? static_cast<T>(1.0) : static_cast<T>(-1.0);
      *(real_y - m_idx + 6) = sign * *(poly.data() + m_idx) * cos_mphi;
      *(imag_y - m_idx + 6) = -sign * *(poly.data() + m_idx) * sin_mphi;
    }
  }
}

// -------------------------------------------------------------------------
// CUDA Kernel: compute Steinhardt parameters Q4 and Q6 per atom
// -------------------------------------------------------------------------
template <typename T>
__global__ void steinhardt_kernel(GPUPoint<T> const *__restrict__ atoms, NeighborGraphPointers graph, int num_atoms,
                                  SteinhardtOutputPointers<T> outputs) {

  int atom_idx = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
  if (atom_idx >= num_atoms) {
    return;
  }

  int const start = *(graph.offsets + atom_idx);
  int const end = *(graph.offsets + atom_idx + 1);
  int const num_neighbors = end - start;

  if (num_neighbors == 0) {
    *(outputs.q4_out + atom_idx) = static_cast<T>(0.0);
    *(outputs.q6_out + atom_idx) = static_cast<T>(0.0);
    return;
  }

  GPUPoint<T> const central_atom = *(atoms + atom_idx);

  std::array<T, 9> q4_real{};
  std::array<T, 9> q4_imag{};
  std::array<T, 13> q6_real{};
  std::array<T, 13> q6_imag{};

  for (int n_idx = start; n_idx < end; ++n_idx) {
    int const neighbor_idx = *(graph.indices + n_idx);
    GPUPoint<T> const neighbor_atom = *(atoms + neighbor_idx);

    T const delta_x = neighbor_atom.x - central_atom.x;
    T const delta_y = neighbor_atom.y - central_atom.y;
    T const delta_z = neighbor_atom.z - central_atom.z;

    T const r_sq = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
    T const r_val = sqrt(r_sq);

    if (r_val <= static_cast<T>(1e-12)) {
      continue;
    }

    T const costheta = delta_z / r_val;
    T const phi = atan2(delta_y, delta_x);

    std::array<T, 9> y4_r{};
    std::array<T, 9> y4_i{};
    std::array<T, 13> y6_r{};
    std::array<T, 13> y6_i{};

    compute_y4m<T>({.costheta = costheta, .phi = phi}, {.real_y = y4_r.data(), .imag_y = y4_i.data()});
    compute_y6m<T>({.costheta = costheta, .phi = phi}, {.real_y = y6_r.data(), .imag_y = y6_i.data()});

    for (int m_idx = 0; m_idx < 9; ++m_idx) {
      *(q4_real.data() + m_idx) += *(y4_r.data() + m_idx);
      *(q4_imag.data() + m_idx) += *(y4_i.data() + m_idx);
    }
    for (int m_idx = 0; m_idx < 13; ++m_idx) {
      *(q6_real.data() + m_idx) += *(y6_r.data() + m_idx);
      *(q6_imag.data() + m_idx) += *(y6_i.data() + m_idx);
    }
  }

  T const inv_n = static_cast<T>(1.0) / static_cast<T>(num_neighbors);
  T sum_sq4 = static_cast<T>(0.0);
  for (int m_idx = 0; m_idx < 9; ++m_idx) {
    *(q4_real.data() + m_idx) *= inv_n;
    *(q4_imag.data() + m_idx) *= inv_n;
    T const real4 = *(q4_real.data() + m_idx);
    T const imag4 = *(q4_imag.data() + m_idx);
    sum_sq4 += (real4 * real4 + imag4 * imag4);
  }

  T sum_sq6 = static_cast<T>(0.0);
  for (int m_idx = 0; m_idx < 13; ++m_idx) {
    *(q6_real.data() + m_idx) *= inv_n;
    *(q6_imag.data() + m_idx) *= inv_n;
    T const real6 = *(q6_real.data() + m_idx);
    T const imag6 = *(q6_imag.data() + m_idx);
    sum_sq6 += (real6 * real6 + imag6 * imag6);
  }

  T const pi_val = static_cast<T>(correlation::math::pi);
  T const norm_q4 = sqrt((static_cast<T>(4.0) * pi_val / static_cast<T>(9.0)) * sum_sq4);
  T const norm_q6 = sqrt((static_cast<T>(4.0) * pi_val / static_cast<T>(13.0)) * sum_sq6);

  *(outputs.q4_out + atom_idx) = norm_q4;
  *(outputs.q6_out + atom_idx) = norm_q6;
}

} // namespace

GPUSteinhardtCalculator::GPUSteinhardtCalculator() {
  int device_count = 0;
  hipError_t err = hipGetDeviceCount(&device_count);
  has_gpu_ = (err == hipSuccess && device_count > 0);
}

void GPUSteinhardtCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                             const correlation::analysis::AnalysisSettings &settings) const {

  if (!has_gpu_) {
    SteinhardtCalculator cpu_calc;
    cpu_calc.calculateFrame(dists, settings);
    return;
  }

  using T = real_t;

  if (dists.neighbors() == nullptr) {
    SteinhardtCalculator cpu_calc;
    cpu_calc.calculateFrame(dists, settings);
    return;
  }

  const auto &graph = dists.neighbors()->neighborGraph();
  const auto &cell = dists.cell();
  const auto &atoms = cell.atoms();
  size_t const num_atoms = atoms.size();

  if (num_atoms == 0) {
    return;
  }

  std::vector<GPUPoint<T>> h_atoms(num_atoms);
  for (size_t i = 0; i < num_atoms; ++i) {
    h_atoms[i] = GPUPoint<T>{.x = static_cast<T>(atoms[i].position().x()),
                             .y = static_cast<T>(atoms[i].position().y()),
                             .z = static_cast<T>(atoms[i].position().z())};
  }

  std::vector<int> h_offsets(num_atoms + 1, 0);
  std::vector<int> h_indices;

  for (size_t i = 0; i < num_atoms; ++i) {
    const auto &neighbors = graph.getNeighbors(i);
    h_offsets[i] = static_cast<int>(h_indices.size());
    for (const auto &neighbor : neighbors) {
      h_indices.push_back(static_cast<int>(neighbor.index));
    }
  }
  h_offsets[num_atoms] = static_cast<int>(h_indices.size());

  if (h_indices.empty()) {
    SteinhardtCalculator cpu_calc;
    cpu_calc.calculateFrame(dists, settings);
    return;
  }

  GPUPoint<T> *d_atoms = nullptr;
  int *d_offsets = nullptr;
  int *d_indices = nullptr;
  T *d_q4 = nullptr;
  T *d_q6 = nullptr;

  hipMalloc(&d_atoms, num_atoms * sizeof(GPUPoint<T>));
  hipMalloc(&d_offsets, (num_atoms + 1) * sizeof(int));
  hipMalloc(&d_indices, h_indices.size() * sizeof(int));
  hipMalloc(&d_q4, num_atoms * sizeof(T));
  hipMalloc(&d_q6, num_atoms * sizeof(T));

  hipMemcpy(d_atoms, h_atoms.data(), num_atoms * sizeof(GPUPoint<T>), hipMemcpyHostToDevice);
  hipMemcpy(d_offsets, h_offsets.data(), (num_atoms + 1) * sizeof(int), hipMemcpyHostToDevice);
  hipMemcpy(d_indices, h_indices.data(), h_indices.size() * sizeof(int), hipMemcpyHostToDevice);

  int const block_size = 256;
  int const grid_size = (static_cast<int>(num_atoms) + block_size - 1) / block_size;
  hipLaunchKernelGGL(steinhardt_kernel<T>, grid_size, block_size, 0, 0, d_atoms,
                     NeighborGraphPointers{d_offsets, d_indices}, static_cast<int>(num_atoms),
                     SteinhardtOutputPointers<T>{d_q4, d_q6});
  hipDeviceSynchronize();

  std::vector<T> h_q4(num_atoms, static_cast<T>(0.0));
  std::vector<T> h_q6(num_atoms, static_cast<T>(0.0));
  hipMemcpy(h_q4.data(), d_q4, num_atoms * sizeof(T), hipMemcpyDeviceToHost);
  hipMemcpy(h_q6.data(), d_q6, num_atoms * sizeof(T), hipMemcpyDeviceToHost);

  hipFree(d_atoms);
  hipFree(d_offsets);
  hipFree(d_indices);
  hipFree(d_q4);
  hipFree(d_q6);

  size_t const bins_Q = 100;
  T const Q_max = static_cast<T>(1.0);
  T const d_q = Q_max / static_cast<T>(bins_Q);

  correlation::analysis::Histogram hist_Q4;
  hist_Q4.x_label = "Q4";
  hist_Q4.title = "Steinhardt Q4 Interface Parameter (GPU)";
  hist_Q4.y_label = "Probability";
  hist_Q4.x_unit = "arbitrary units";
  hist_Q4.y_unit = "counts";
  hist_Q4.description = "Steinhardt Q4 Bond Orientational Order Parameter (GPU)";
  hist_Q4.file_suffix = "_Q4_gpu";
  hist_Q4.bins.resize(bins_Q);
  std::vector<real_t> q4_bins(bins_Q, static_cast<real_t>(0.0));

  correlation::analysis::Histogram hist_Q6;
  hist_Q6.x_label = "Q6";
  hist_Q6.title = "Steinhardt Q6 Interface Parameter (GPU)";
  hist_Q6.y_label = "Probability";
  hist_Q6.x_unit = "arbitrary units";
  hist_Q6.y_unit = "counts";
  hist_Q6.description = "Steinhardt Q6 Bond Orientational Order Parameter (GPU)";
  hist_Q6.file_suffix = "_Q6_gpu";
  hist_Q6.bins.resize(bins_Q);
  std::vector<real_t> q6_bins(bins_Q, static_cast<real_t>(0.0));

  for (size_t bin_idx = 0; bin_idx < bins_Q; ++bin_idx) {
    *(hist_Q4.bins.data() + bin_idx) = static_cast<real_t>((static_cast<T>(bin_idx) + static_cast<T>(0.5)) * d_q);
    *(hist_Q6.bins.data() + bin_idx) = static_cast<real_t>((static_cast<T>(bin_idx) + static_cast<T>(0.5)) * d_q);
  }

  for (size_t i = 0; i < num_atoms; ++i) {
    T const q4_val = *(h_q4.data() + i);
    T const q6_val = *(h_q6.data() + i);
    if (q4_val >= static_cast<T>(0.0) && q4_val < Q_max) {
      auto bin = static_cast<size_t>(q4_val / d_q);
      if (bin < bins_Q) {
        *(q4_bins.data() + bin) += static_cast<real_t>(1.0);
      }
    }
    if (q6_val >= static_cast<T>(0.0) && q6_val < Q_max) {
      auto bin = static_cast<size_t>(q6_val / d_q);
      if (bin < bins_Q) {
        *(q6_bins.data() + bin) += static_cast<real_t>(1.0);
      }
    }
  }

  hist_Q4.partials["Total"] = q4_bins;
  hist_Q6.partials["Total"] = q6_bins;

  dists.addHistogram("Q4_gpu", std::move(hist_Q4));
  dists.addHistogram("Q6_gpu", std::move(hist_Q6));
}

} // namespace correlation::calculators
