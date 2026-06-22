/**
 * @file Smoothing.hpp
 * @brief Kernel-based smoothing algorithms for analysis data.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/Constants.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace correlation::math {

/**
 * @enum KernelType
 * @brief Selection of kernel types for smoothing.
 */
enum class KernelType : std::uint8_t { Gaussian, Bump, Triweight, Epanechnikov, Cosine, Biweight };

/**
 * @brief Parameters for kernel generation.
 */
struct KernelGenerationParams {
  size_t size;
  double bin_width;
  double sigma;
  KernelType type;
};
namespace detail {

inline void fillGaussian(std::vector<double> &kernel, const KernelGenerationParams &params, double center) {
  const double prefactor = 1.0 / (std::sqrt(correlation::math::two_pi) * params.sigma);
  const double exp_coeff = -1.0 / (2.0 * params.sigma * params.sigma);
  for (size_t i = 0; i < params.size; ++i) {
    const double distance = (static_cast<double>(i) - center) * params.bin_width;
    kernel[i] = prefactor * std::exp(exp_coeff * distance * distance);
  }
}

inline void fillTriweight(std::vector<double> &kernel, const KernelGenerationParams &params, double center) {
  constexpr double factor = 35.0 / 32.0;
  for (size_t i = 0; i < params.size; ++i) {
    const double norm_dist = ((static_cast<double>(i) - center) * params.bin_width) / params.sigma;
    if (std::abs(norm_dist) <= 1.0) {
      const double norm_dist_sq = norm_dist * norm_dist;
      kernel[i] = factor * std::pow(1.0 - norm_dist_sq, 3);
    } else {
      kernel[i] = 0.0;
    }
  }
}

inline void fillBump(std::vector<double> &kernel, const KernelGenerationParams &params, double center) {
  for (size_t i = 0; i < params.size; ++i) {
    const double norm_dist = ((static_cast<double>(i) - center) * params.bin_width) / params.sigma;
    if (std::abs(norm_dist) < 1.0) {
      const double norm_dist_sq = norm_dist * norm_dist;
      kernel[i] = std::exp(-1.0 / (1.0 - norm_dist_sq));
    } else {
      kernel[i] = 0.0;
    }
  }
}

inline void fillEpanechnikov(std::vector<double> &kernel, const KernelGenerationParams &params, double center) {
  constexpr double factor = 3.0 / 4.0;
  for (size_t i = 0; i < params.size; ++i) {
    const double norm_dist = ((static_cast<double>(i) - center) * params.bin_width) / params.sigma;
    if (std::abs(norm_dist) <= 1.0) {
      kernel[i] = factor * (1.0 - norm_dist * norm_dist);
    } else {
      kernel[i] = 0.0;
    }
  }
}

inline void fillCosine(std::vector<double> &kernel, const KernelGenerationParams &params, double center) {
  const double cos_factor = correlation::math::pi / 4.0;
  for (size_t i = 0; i < params.size; ++i) {
    const double norm_dist = ((static_cast<double>(i) - center) * params.bin_width) / params.sigma;
    if (std::abs(norm_dist) <= 1.0) {
      kernel[i] = cos_factor * std::cos(correlation::math::pi * norm_dist / 2.0);
    } else {
      kernel[i] = 0.0;
    }
  }
}

inline void fillBiweight(std::vector<double> &kernel, const KernelGenerationParams &params, double center) {
  constexpr double factor = 15.0 / 16.0;
  for (size_t i = 0; i < params.size; ++i) {
    const double norm_dist = ((static_cast<double>(i) - center) * params.bin_width) / params.sigma;
    if (std::abs(norm_dist) <= 1.0) {
      const double norm_dist_sq = norm_dist * norm_dist;
      kernel[i] = factor * (1.0 - norm_dist_sq) * (1.0 - norm_dist_sq);
    } else {
      kernel[i] = 0.0;
    }
  }
}

} // namespace detail

/**
 * @brief Generates a smoothing kernel.
 * @param params The parameters for kernel generation.
 * @return A vector containing the normalized kernel.
 */
[[nodiscard]] inline std::vector<double> generateKernel(const KernelGenerationParams &params) {
  const size_t size = params.size;
  std::vector<double> kernel(size);
  const double center = static_cast<double>(size - 1) / 2.0;

  switch (params.type) {
  case KernelType::Gaussian:
    detail::fillGaussian(kernel, params, center);
    break;
  case KernelType::Triweight:
    detail::fillTriweight(kernel, params, center);
    break;
  case KernelType::Bump:
    detail::fillBump(kernel, params, center);
    break;
  case KernelType::Epanechnikov:
    detail::fillEpanechnikov(kernel, params, center);
    break;
  case KernelType::Cosine:
    detail::fillCosine(kernel, params, center);
    break;
  case KernelType::Biweight:
    detail::fillBiweight(kernel, params, center);
    break;
  default:
    throw std::invalid_argument("Unsupported kernel type for smoothing.");
  }

  // Normalize the discrete kernel to ensure the sum of its elements is 1.
  const double sum = std::accumulate(kernel.begin(), kernel.end(), 0.0);
  if (sum > 1e-9) {
    for (double &val : kernel) {
      val /= sum;
    }
  }
  return kernel;
}

/**
 * @brief Applies kernel smoothing to a data series.
 *
 * This function convolves the chosen kernel with `y_values` to produce a smoothed
 * version. The bin width `bin_width` fully characterises the sampling grid; the full
 * grid vector is not required and is not accepted to keep the interface lean.
 *
 * @param bin_width    Uniform bin width of the data (e.g., Δr in Å).
 * @param y_values     The dependent variable values to be smoothed.
 * @param sigma Kernel bandwidth (standard deviation for Gaussian), in the
 *              same physical units as `bin_width`.
 * @param type  Kernel type to use.
 * @return A vector containing the smoothed data (same length as `y_values`).
 */
inline std::vector<double> KernelSmoothing(double bin_width, const std::vector<double> &y_values, double sigma,
                                           KernelType type) {
  if (y_values.empty()) {
    return {};
  }
  if (bin_width <= 0) {
    throw std::invalid_argument("'bin_width' must be a positive bin width for smoothing.");
  }
  if (sigma <= 0.0) {
    throw std::invalid_argument("Smoothing sigma must be positive.");
  }

  // Calculate the kernel radius in bins by dividing the physical sigma by bin_width.
  // We use 4.0 * sigma to cover approximately 99.99% of the Gaussian area.
  const size_t num_points = y_values.size();
  auto kernel_radius_bins = static_cast<size_t>(4.0 * sigma / bin_width);

  // Clamp the radius to prevent the kernel from exceeding half the data size.
  const size_t max_kernel_radius = num_points / 2;
  size_t kernel_radius = std::min(kernel_radius_bins, max_kernel_radius);

  // Ensure a minimum kernel size (3 bins: -1, 0, +1) if sigma is very small.
  if (kernel_radius == 0) {
    kernel_radius = 1;
  }

  const size_t kernel_size = 2 * kernel_radius + 1;

  // Pass the original physical sigma (in distance units) to the kernel
  // generator.
  const auto kernel = generateKernel({kernel_size, bin_width, sigma, type});
  std::vector<double> smoothed(num_points, 0.0);

  // Phase 1 – left boundary (bins 0 .. kernel_radius-1): needs index clamping.
  for (size_t i = 0; i < std::min(kernel_radius, num_points); ++i) {
    for (size_t j = 0; j < kernel_size; ++j) {
      const long long idx =
          static_cast<long long>(i) + static_cast<long long>(j) - static_cast<long long>(kernel_radius);
      smoothed[i] += y_values[static_cast<size_t>(std::max(0LL, idx))] * kernel[j];
    }
  }

  // Phase 2 – interior (bins kernel_radius .. num_points-kernel_radius-1):
  // index  i - kernel_radius + j  is always in [0, num_points-1],  no branch needed.
  // The clean inner loop (sum += src[j] * kernel[j]) is auto-vectorised by
  // GCC/Clang at -O2 using whatever SIMD width the TU was compiled with.
  if (num_points > 2 * kernel_radius) {
    for (size_t i = kernel_radius; i < num_points - kernel_radius; ++i) {
      const double *src = y_values.data() + (i - kernel_radius);
      double sum = 0.0;
      for (size_t j = 0; j < kernel_size; ++j) {
        sum += src[j] * kernel[j];
      }
      smoothed[i] = sum;
    }
  }

  // Phase 3 – right boundary (bins num_points-kernel_radius .. num_points-1): needs clamping.
  const size_t right_start = (num_points > kernel_radius) ? num_points - kernel_radius : num_points;
  for (size_t i = right_start; i < num_points; ++i) {
    for (size_t j = 0; j < kernel_size; ++j) {
      const long long idx =
          static_cast<long long>(i) + static_cast<long long>(j) - static_cast<long long>(kernel_radius);
      const long long clamped = std::min(idx, static_cast<long long>(num_points) - 1);
      smoothed[i] += y_values[static_cast<size_t>(clamped)] * kernel[j];
    }
  }

  return smoothed;
}

/**
 * @brief Compatibility overload: derives `bin_width` from the grid and delegates.
 *
 * Prefer the `(bin_width, y_values, sigma, type)` overload when the bin width is already
 * known to avoid recomputing it.
 *
 * @param r_values The independent variable grid (e.g., radial distances).
 * @param y_values The dependent variable values to be smoothed.
 * @param sigma Kernel bandwidth.
 * @param type Kernel type to use.
 * @return A vector containing the smoothed data.
 */
inline std::vector<double> KernelSmoothing(const std::vector<double> &r_values, const std::vector<double> &y_values,
                                           double sigma, KernelType type) {
  if (r_values.size() != y_values.size() || r_values.size() < 2) {
    return {};
  }
  const double bin_width = r_values[1] - r_values[0];
  return KernelSmoothing(bin_width, y_values, sigma, type);
}

} // namespace correlation::math
