/**
 * @file GPUXRDCalculator.cu
 * @brief CUDA/HIP implementation of GPU-accelerated X-Ray Diffraction (XRD) using direct Debye scattering.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "calculators/CalculatorFactory.hpp"
#include "calculators/GPUXRDCalculator.hpp"
#include "calculators/XRDCalculator.hpp"
#include "core/DeviceBuffer.hpp"
#include "core/GPUPortability.hpp"
#include "math/Constants.hpp"
#include "math/Precision.hpp"
#include "physics/PhysicalData.hpp"

#include <cmath>
#include <map>
#include <stdexcept>
#include <vector>

namespace correlation::calculators {

namespace {

const bool registered = CalculatorFactory::registerTypeSafe<GPUXRDCalculator>("GPUXRDCalculator");

bool has_gpu_device() {
  int device_count = 0;
  hipError_t const err = hipGetDeviceCount(&device_count);
  return (err == hipSuccess && device_count > 0);
}

template <typename T> struct CromerMannCoeffs {
  std::array<T, 4> a{};
  std::array<T, 4> b{};
  T c{};
};

template <typename T> struct DeviceAtomData {
  const T *CORRELATION_RESTRICT x;
  const T *CORRELATION_RESTRICT y;
  const T *CORRELATION_RESTRICT z;
  const int *CORRELATION_RESTRICT type_idx;
};

template <typename T> struct XRDGridConfig {
  T lambda;
  T theta_min;
  T theta_max;
  T bin_width;
  int num_bins;
};

template <typename T> CORRELATION_DEVICE T eval_form_factor(const CromerMannCoeffs<T> &coeffs, T q_val) {
  T const s_val = q_val / static_cast<T>(correlation::math::four_pi);
  T const s_sq = s_val * s_val;
  T form_factor = coeffs.c;
  for (int k = 0; k < 4; ++k) {
    form_factor += coeffs.a[k] * std::exp(-coeffs.b[k] * s_sq);
  }
  return form_factor;
}

template <typename T>
CORRELATION_GLOBAL void
debye_xrd_kernel(DeviceAtomData<T> atoms, int num_atoms, const CromerMannCoeffs<T> *CORRELATION_RESTRICT element_coeffs,
                 int num_elements, XRDGridConfig<T> grid_cfg, T *CORRELATION_RESTRICT out_intensity) {
  int const bin_idx = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
  if (bin_idx >= grid_cfg.num_bins) {
    return;
  }

  T const two_theta = grid_cfg.theta_min + static_cast<T>(bin_idx) * grid_cfg.bin_width;
  T const theta_rad = (two_theta / static_cast<T>(2.0)) * static_cast<T>(correlation::math::deg_to_rad);
  T const q_val = static_cast<T>(correlation::math::four_pi) * std::sin(theta_rad) / grid_cfg.lambda;

  if (q_val < static_cast<T>(1e-6)) {
    out_intensity[bin_idx] = static_cast<T>(0.0);
    return;
  }

  std::array<T, 16> f_values{};
  for (int elem_idx = 0; elem_idx < num_elements && elem_idx < 16; ++elem_idx) {
    f_values[static_cast<size_t>(elem_idx)] = eval_form_factor(element_coeffs[elem_idx], q_val);
  }

  T intensity = static_cast<T>(0.0);

  for (int i = 0; i < num_atoms; ++i) {
    int const type_i = atoms.type_idx[i];
    T const f_i = f_values[static_cast<size_t>(type_i)];
    intensity += f_i * f_i;

    T const x_i = atoms.x[i];
    T const y_i = atoms.y[i];
    T const z_i = atoms.z[i];

    for (int j = i + 1; j < num_atoms; ++j) {
      int const type_j = atoms.type_idx[j];
      T const f_j = f_values[static_cast<size_t>(type_j)];

      T const delta_x = x_i - atoms.x[j];
      T const delta_y = y_i - atoms.y[j];
      T const delta_z = z_i - atoms.z[j];
      T const r_ij = std::sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);

      T const q_dist = q_val * r_ij;
      T const sinc_val = (q_dist > static_cast<T>(1e-6)) ? (std::sin(q_dist) / q_dist) : static_cast<T>(1.0);

      intensity += static_cast<T>(2.0) * f_i * f_j * sinc_val;
    }
  }

  if (num_atoms > 0) {
    out_intensity[bin_idx] = intensity / static_cast<T>(num_atoms);
  } else {
    out_intensity[bin_idx] = static_cast<T>(0.0);
  }
}

template <typename T> std::vector<CromerMannCoeffs<T>> build_element_coeffs(const correlation::core::Cell &cell) {
  const auto &elements = cell.elements();
  size_t const num_elements = elements.size();
  std::vector<CromerMannCoeffs<T>> coeffs(num_elements);

  for (size_t i = 0; i < num_elements; ++i) {
    const auto &form_factors = correlation::physics::getAtomicFormFactors(elements[i].symbol);
    coeffs[i].a[0] = static_cast<T>(form_factors[0]);
    coeffs[i].b[0] = static_cast<T>(form_factors[1]);
    coeffs[i].a[1] = static_cast<T>(form_factors[2]);
    coeffs[i].b[1] = static_cast<T>(form_factors[3]);
    coeffs[i].a[2] = static_cast<T>(form_factors[4]);
    coeffs[i].b[2] = static_cast<T>(form_factors[5]);
    coeffs[i].a[3] = static_cast<T>(form_factors[6]);
    coeffs[i].b[3] = static_cast<T>(form_factors[7]);
    coeffs[i].c = static_cast<T>(form_factors[8]);
  }
  return coeffs;
}

template <typename T> struct HostAtomBuffer {
  std::vector<T> x;
  std::vector<T> y;
  std::vector<T> z;
  std::vector<int> type_idx;
};

template <typename T> HostAtomBuffer<T> extract_host_atoms(const correlation::core::Cell &cell) {
  const auto &atoms = cell.atoms();
  const size_t num_atoms = atoms.size();
  const auto &elements = cell.elements();

  std::map<std::string, int> element_map;
  for (size_t i = 0; i < elements.size(); ++i) {
    element_map[elements[i].symbol] = static_cast<int>(i);
  }

  HostAtomBuffer<T> h_atoms{
      .x = std::vector<T>(num_atoms),
      .y = std::vector<T>(num_atoms),
      .z = std::vector<T>(num_atoms),
      .type_idx = std::vector<int>(num_atoms),
  };

  for (size_t i = 0; i < num_atoms; ++i) {
    const auto &pos = atoms[i].position();
    h_atoms.x[i] = static_cast<T>(pos.x());
    h_atoms.y[i] = static_cast<T>(pos.y());
    h_atoms.z[i] = static_cast<T>(pos.z());
    h_atoms.type_idx[i] = element_map.at(atoms[i].element().symbol);
  }
  return h_atoms;
}

template <typename T>
std::vector<T> execute_gpu_debye(const correlation::core::Cell &cell, const XRDGridConfig<T> &grid_cfg) {
  size_t const num_atoms = cell.atoms().size();
  size_t const num_elements = cell.elements().size();
  auto const num_bins = static_cast<size_t>(grid_cfg.num_bins);

  HostAtomBuffer<T> const h_atoms = extract_host_atoms<T>(cell);
  std::vector<CromerMannCoeffs<T>> const h_coeffs = build_element_coeffs<T>(cell);

  using correlation::core::gpu::DeviceBuffer;
  DeviceBuffer<T> d_x(num_atoms);
  DeviceBuffer<T> d_y(num_atoms);
  DeviceBuffer<T> d_z(num_atoms);
  DeviceBuffer<int> d_type_idx(num_atoms);
  DeviceBuffer<CromerMannCoeffs<T>> d_coeffs(num_elements);
  DeviceBuffer<T> const d_intensity(num_bins);

  d_x.copyFromHost(h_atoms.x.data(), num_atoms);
  d_y.copyFromHost(h_atoms.y.data(), num_atoms);
  d_z.copyFromHost(h_atoms.z.data(), num_atoms);
  d_type_idx.copyFromHost(h_atoms.type_idx.data(), num_atoms);
  d_coeffs.copyFromHost(h_coeffs.data(), num_elements);

  int const block_size = 256;
  int const grid_size = (static_cast<int>(num_bins) + block_size - 1) / block_size;

  DeviceAtomData<T> const dev_atoms{
      .x = d_x.get(),
      .y = d_y.get(),
      .z = d_z.get(),
      .type_idx = d_type_idx.get(),
  };

  hipLaunchKernelGGL(debye_xrd_kernel<T>, grid_size, block_size, 0, 0, dev_atoms, static_cast<int>(num_atoms),
                     d_coeffs.get(), static_cast<int>(num_elements), grid_cfg, d_intensity.get());
  correlation::core::gpu::hipCheck(hipDeviceSynchronize());

  std::vector<T> h_intensity(num_bins);
  d_intensity.copyToHost(h_intensity.data(), num_bins);
  return h_intensity;
}

correlation::analysis::Histogram assemble_histogram(const std::vector<real_t> &intensities,
                                                    const GPUXRDParams &params) {
  size_t const num_bins = intensities.size();
  correlation::analysis::Histogram hist;
  hist.x_label = "2θ";
  hist.title = "XRD Pattern (GPU Debye)";
  hist.y_label = "Intensity";
  hist.x_unit = "°";
  hist.y_unit = "arbitrary units";
  hist.description = "X-Ray Diffraction pattern computed on GPU using direct Debye scattering";
  hist.file_suffix = "_XRD_gpu";
  hist.bins.resize(num_bins);

  for (size_t i = 0; i < num_bins; ++i) {
    hist.bins[i] = params.theta_min + static_cast<real_t>(i) * params.bin_width;
  }
  hist.partials["Total"] = intensities;
  return hist;
}

} // namespace

namespace gpu {

correlation::analysis::Histogram compute_xrd_gpu(const correlation::core::Cell &cell, const GPUXRDParams &params) {
  if (params.bin_width <= static_cast<real_t>(0.0)) {
    throw std::invalid_argument("Angular resolution bin_width must be strictly positive.");
  }
  if (params.theta_max <= params.theta_min) {
    throw std::invalid_argument("theta_max must be greater than theta_min.");
  }
  if (params.lambda <= static_cast<real_t>(0.0)) {
    throw std::invalid_argument("Wavelength lambda must be strictly positive.");
  }

  size_t const num_bins = static_cast<size_t>((params.theta_max - params.theta_min) / params.bin_width) + 1;

  if (cell.atomCount() == 0) {
    std::vector<real_t> const zero_intensity(num_bins, static_cast<real_t>(0.0));
    return assemble_histogram(zero_intensity, params);
  }

  if (!has_gpu_device()) {
    // Fallback to CPU calculation via S(Q) / XRDCalculator
    XRDCalculator const cpu_calc;
    correlation::analysis::DistributionFunctions dists(cell);
    correlation::analysis::AnalysisSettings const settings{};
    cpu_calc.calculateFrame(dists, settings);
    if (dists.getAllHistograms().contains("XRD")) {
      return dists.getHistogram("XRD");
    }
    std::vector<real_t> const fallback_intensity(num_bins, static_cast<real_t>(0.0));
    return assemble_histogram(fallback_intensity, params);
  }

  XRDGridConfig<real_t> const grid_cfg{
      .lambda = params.lambda,
      .theta_min = params.theta_min,
      .theta_max = params.theta_max,
      .bin_width = params.bin_width,
      .num_bins = static_cast<int>(num_bins),
  };

  std::vector<real_t> const intensities = execute_gpu_debye<real_t>(cell, grid_cfg);
  return assemble_histogram(intensities, params);
}

} // namespace gpu

GPUXRDCalculator::GPUXRDCalculator() : has_gpu_(has_gpu_device()) {}

void GPUXRDCalculator::calculateFrame(correlation::analysis::DistributionFunctions &dists,
                                      const correlation::analysis::AnalysisSettings & /*settings*/) const {
  GPUXRDParams const params{
      .lambda = static_cast<real_t>(1.5406),
      .theta_min = static_cast<real_t>(10.0),
      .theta_max = static_cast<real_t>(80.0),
      .bin_width = static_cast<real_t>(0.1),
  };

  correlation::analysis::Histogram hist = gpu::compute_xrd_gpu(dists.cell(), params);
  dists.addHistogram("XRD_gpu", std::move(hist));
}

} // namespace correlation::calculators
