/**
 * @file Smoothing.hpp
 * @brief Kernel-based smoothing algorithms for analysis data.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @license SPDX-License-Identifier: MIT
 */

#pragma once

#include "math/Constants.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace correlation::math {

/**
 * @enum KernelType
 * @brief Selection of kernel types for smoothing.
 */
enum class KernelType { Gaussian, Bump, Triweight };

/**
 * @brief Generates a smoothing kernel.
 * @param size The size of the kernel in bins (must be odd).
 * @param dx The bin width of the data (e.g., delta r).
 * @param sigma The bandwidth (standard deviation) of the kernel.
 * @param type The type of kernel to generate.
 * @return A vector containing the normalized kernel.
 */

inline std::vector<double> generateKernel(size_t size, double dx, double sigma,
                                          KernelType type) {
  std::vector<double> kernel(size);
  // Center of the discrete kernel in bin index units.
  const double center = static_cast<double>(size - 1) / 2.0;

  if (type == KernelType::Gaussian) {
    // 1 / (sigma * sqrt(2*pi))
    const double a = 1.0 / (std::sqrt(correlation::math::two_pi) * sigma);
    const double b = -1.0 / (2.0 * sigma * sigma);
    for (size_t i = 0; i < size; ++i) {
      // Calculate x in distance units (r units)
      const double x = (static_cast<double>(i) - center) * dx;
      kernel[i] = a * std::exp(b * x * x);
    }
  } else if (type == KernelType::Triweight) {
    constexpr double factor = 35.0 / 32.0;
    for (size_t i = 0; i < size; ++i) {
      // Calculate u (distance normalized by sigma)
      const double u = ((static_cast<double>(i) - center) * dx) / sigma;
      if (std::abs(u) <= 1.0) {
        const double u2 = u * u;
        kernel[i] = factor * std::pow(1 - u2, 3);
      } else {
        kernel[i] = 0.0;
      }
    }
  } else if (type == KernelType::Bump) {
    for (size_t i = 0; i < size; ++i) {
      // Calculate u (distance normalized by sigma)
      const double u = ((static_cast<double>(i) - center) * dx) / sigma;
      if (std::abs(u) < 1.0) {
        const double u2 = u * u;
        kernel[i] = std::exp(-1.0 / (1.0 - u2));
      } else {
        kernel[i] = 0.0;
      }
    }
  } else {
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
 * This function convolves the chosen kernel with `y` to produce a smoothed
 * version. The bin width `dx` fully characterises the sampling grid; the full
 * grid vector is not required and is not accepted to keep the interface lean.
 *
 * @param dx    Uniform bin width of the data (e.g., Δr in Å).
 * @param y     The dependent variable values to be smoothed.
 * @param sigma Kernel bandwidth (standard deviation for Gaussian), in the
 *              same physical units as `dx`.
 * @param type  Kernel type to use.
 * @return A vector containing the smoothed data (same length as `y`).
 */
inline std::vector<double> KernelSmoothing(double dx,
                                           const std::vector<double> &y,
                                           double sigma, KernelType type) {
  if (y.empty()) {
    return {};
  }
  if (dx <= 0) {
    throw std::invalid_argument(
        "'dx' must be a positive bin width for smoothing.");
  }
  if (sigma <= 0.0) {
    throw std::invalid_argument("Smoothing sigma must be positive.");
  }

  // Calculate the kernel radius in bins by dividing the physical sigma by dx.
  // We use 4.0 * sigma to cover approximately 99.99% of the Gaussian area.
  const size_t n = y.size();
  size_t kernel_radius_bins = static_cast<size_t>(4.0 * sigma / dx);

  // Clamp the radius to prevent the kernel from exceeding half the data size.
  const size_t max_kernel_radius = n / 2;
  size_t kernel_radius = std::min(kernel_radius_bins, max_kernel_radius);

  // Ensure a minimum kernel size (3 bins: -1, 0, +1) if sigma is very small.
  if (kernel_radius == 0) {
    kernel_radius = 1;
  }

  const size_t kernel_size = 2 * kernel_radius + 1;

  // Pass the original physical sigma (in distance units) to the kernel
  // generator.
  const auto kernel = generateKernel(kernel_size, dx, sigma, type);
  std::vector<double> smoothed(n, 0.0);

  // Phase 1 – left boundary (bins 0 .. kernel_radius-1): needs index clamping.
  for (size_t i = 0; i < std::min(kernel_radius, n); ++i) {
    for (size_t j = 0; j < kernel_size; ++j) {
      const long long idx = static_cast<long long>(i) +
                            static_cast<long long>(j) -
                            static_cast<long long>(kernel_radius);
      smoothed[i] += y[static_cast<size_t>(std::max(0LL, idx))] * kernel[j];
    }
  }

  // Phase 2 – interior (bins kernel_radius .. n-kernel_radius-1):
  // index  i - kernel_radius + j  is always in [0, n-1],  no branch needed.
  // The clean inner loop (sum += src[j] * kernel[j]) is auto-vectorised by
  // GCC/Clang at -O2 using whatever SIMD width the TU was compiled with.
  if (n > 2 * kernel_radius) {
    for (size_t i = kernel_radius; i < n - kernel_radius; ++i) {
      const double *src = y.data() + (i - kernel_radius);
      double sum = 0.0;
      for (size_t j = 0; j < kernel_size; ++j)
        sum += src[j] * kernel[j];
      smoothed[i] = sum;
    }
  }

  // Phase 3 – right boundary (bins n-kernel_radius .. n-1): needs clamping.
  const size_t right_start = (n > kernel_radius) ? n - kernel_radius : n;
  for (size_t i = right_start; i < n; ++i) {
    for (size_t j = 0; j < kernel_size; ++j) {
      const long long idx = static_cast<long long>(i) +
                            static_cast<long long>(j) -
                            static_cast<long long>(kernel_radius);
      const long long clamped = std::min(idx, static_cast<long long>(n) - 1);
      smoothed[i] += y[static_cast<size_t>(clamped)] * kernel[j];
    }
  }

  return smoothed;
}

/**
 * @brief Compatibility overload: derives `dx` from the grid and delegates.
 *
 * Prefer the `(dx, y, sigma, type)` overload when the bin width is already
 * known to avoid recomputing it.
 */
inline std::vector<double> KernelSmoothing(const std::vector<double> &r,
                                           const std::vector<double> &y,
                                           double sigma, KernelType type) {
  if (r.size() != y.size() || r.size() < 2)
    return {};
  const double dx = r[1] - r[0];
  return KernelSmoothing(dx, y, sigma, type);
}

} // namespace correlation::math
