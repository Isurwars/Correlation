/**
 * @file Smoothing.hpp
 * @brief Kernel-based smoothing algorithms for analysis data.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/Constants.hpp"
#include "math/Precision.hpp"
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
  size_t size;      ///< Discrete kernel array length (odd integer).
  real_t bin_width; ///< Bin width of the underlying data.
  real_t sigma;     ///< Physical kernel standard deviation / bandwidth.
  KernelType type;  ///< Mathematical kernel function type.
};

/**
 * @brief Options bundling for data series kernel smoothing.
 */
struct KernelSmoothingParams {
  real_t bin_width = 0.0;                 ///< Uniform bin width of the data (e.g., Δr in Å).
  real_t sigma = 0.1;                     ///< Kernel bandwidth standard deviation.
  KernelType type = KernelType::Gaussian; ///< Kernel function selection.
};
namespace detail {

inline void fillGaussian(std::vector<real_t> &kernel, const KernelGenerationParams &params, real_t center) {
  const real_t prefactor = static_cast<real_t>(1.0) / static_cast<real_t>(std::sqrt(two_pi) * params.sigma);
  const real_t exp_coeff = static_cast<real_t>(-1.0) / (static_cast<real_t>(2.0) * params.sigma * params.sigma);
  for (size_t i = 0; i < params.size; ++i) {
    const real_t distance = (static_cast<real_t>(i) - center) * params.bin_width;
    kernel[i] = prefactor * std::exp(exp_coeff * distance * distance);
  }
}

inline void fillTriweight(std::vector<real_t> &kernel, const KernelGenerationParams &params, real_t center) {
  constexpr real_t factor = static_cast<real_t>(35.0) / static_cast<real_t>(32.0);
  for (size_t i = 0; i < params.size; ++i) {
    const real_t norm_dist = ((static_cast<real_t>(i) - center) * params.bin_width) / params.sigma;
    if (std::abs(norm_dist) <= 1.0) {
      const real_t norm_dist_sq = norm_dist * norm_dist;
      kernel[i] = factor * std::pow(static_cast<real_t>(1.0) - norm_dist_sq, static_cast<real_t>(3.0));
    } else {
      kernel[i] = 0.0;
    }
  }
}

inline void fillBump(std::vector<real_t> &kernel, const KernelGenerationParams &params, real_t center) {
  for (size_t i = 0; i < params.size; ++i) {
    const real_t norm_dist = ((static_cast<real_t>(i) - center) * params.bin_width) / params.sigma;
    if (std::abs(norm_dist) < 1.0) {
      const real_t norm_dist_sq = norm_dist * norm_dist;
      kernel[i] = std::exp(-static_cast<real_t>(1.0) / (static_cast<real_t>(1.0) - norm_dist_sq));
    } else {
      kernel[i] = 0.0;
    }
  }
}

inline void fillEpanechnikov(std::vector<real_t> &kernel, const KernelGenerationParams &params, real_t center) {
  const real_t factor = static_cast<real_t>(3.0) / static_cast<real_t>(4.0);
  for (size_t i = 0; i < params.size; ++i) {
    const real_t norm_dist = ((static_cast<real_t>(i) - center) * params.bin_width) / params.sigma;
    if (std::abs(norm_dist) <= 1.0) {
      kernel[i] = factor * (static_cast<real_t>(1.0) - norm_dist * norm_dist);
    } else {
      kernel[i] = 0.0;
    }
  }
}

inline void fillCosine(std::vector<real_t> &kernel, const KernelGenerationParams &params, real_t center) {
  const real_t cos_factor = pi / static_cast<real_t>(4.0);
  for (size_t i = 0; i < params.size; ++i) {
    const real_t norm_dist = ((static_cast<real_t>(i) - center) * params.bin_width) / params.sigma;
    if (std::abs(norm_dist) <= 1.0) {
      kernel[i] = static_cast<real_t>(cos_factor * std::cos(pi * norm_dist / static_cast<real_t>(2.0)));
    } else {
      kernel[i] = 0.0;
    }
  }
}

inline void fillBiweight(std::vector<real_t> &kernel, const KernelGenerationParams &params, real_t center) {
  const real_t factor = static_cast<real_t>(15.0) / static_cast<real_t>(16.0);
  for (size_t i = 0; i < params.size; ++i) {
    const real_t norm_dist = ((static_cast<real_t>(i) - center) * params.bin_width) / params.sigma;
    if (std::abs(norm_dist) <= 1.0) {
      const real_t norm_dist_sq = norm_dist * norm_dist;
      kernel[i] = factor * (static_cast<real_t>(1.0) - norm_dist_sq) * (static_cast<real_t>(1.0) - norm_dist_sq);
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
[[nodiscard]] inline std::vector<real_t> generateKernel(const KernelGenerationParams &params) {
  const size_t size = params.size;
  std::vector<real_t> kernel(size);
  const real_t center = static_cast<real_t>(size - 1) / static_cast<real_t>(2.0);

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
  const real_t sum = std::accumulate(kernel.begin(), kernel.end(), static_cast<real_t>(0.0));
  if (sum > 1e-9) {
    for (real_t &val : kernel) {
      val /= sum;
    }
  }
  return kernel;
}

/**
 * @brief Applies kernel smoothing to a data series.
 *
 * This function convolves the chosen kernel with `y_values` to produce a smoothed
 * version. The bin width in `params` characterises the sampling grid.
 *
 * @param y_values The dependent variable values to be smoothed.
 * @param params   Kernel smoothing configuration parameters.
 * @return A vector containing the smoothed data (same length as `y_values`).
 */
inline std::vector<real_t> KernelSmoothing(const std::vector<real_t> &y_values, KernelSmoothingParams params) {
  real_t const bin_width = params.bin_width;
  real_t const sigma = params.sigma;
  KernelType const type = params.type;

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
  std::vector<real_t> smoothed(num_points, 0.0);

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
      const real_t *src = y_values.data() + (i - kernel_radius);
      real_t sum = 0.0;
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
 * @brief Grid overload: derives `bin_width` from the independent variable grid.
 *
 * @param r_values The independent variable grid (e.g., radial distances).
 * @param y_values The dependent variable values to be smoothed.
 * @param params   Kernel smoothing configuration parameters.
 * @return A vector containing the smoothed data.
 */
inline std::vector<real_t> KernelSmoothing(const std::vector<real_t> &r_values, const std::vector<real_t> &y_values,
                                           KernelSmoothingParams params) {
  if (r_values.size() != y_values.size() || r_values.size() < 2) {
    return {};
  }
  params.bin_width = r_values[1] - r_values[0];
  return KernelSmoothing(y_values, params);
}

} // namespace correlation::math
