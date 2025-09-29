// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "../include/PhysicalData.hpp"

// Enum class for type-safe selection of kernel types.
// This prevents passing invalid integer values.
enum class KernelType { Gaussian, Bump, Triweight };

/**
 * @brief Applies kernel smoothing to a data series.
 *
 * This function convolves a given kernel with the input data `y` to produce a
 * smoothed version.
 *
 * @param r The independent variable values (x-axis).
 * @param y The dependent variable values (y-axis) to be smoothed.
 * @param sigma The bandwidth (standard deviation for Gaussian) of the kernel,
 * in the same units as the 'r' bins (e.g., Å or degrees).
 * @param type The type of kernel to use for smoothing.
 * @return A vector containing the smoothed data.
 */

inline std::vector<double> generateKernel(size_t size, double dx, double sigma,
                                          KernelType type) {
  std::vector<double> kernel(size);
  // Center of the discrete kernel in bin index units.
  const double center = static_cast<double>(size - 1) / 2.0;

  if (type == KernelType::Gaussian) {
    // 1 / (sigma * sqrt(2*pi))
    const double a = 1.0 / (std::sqrt(2 * constants::pi) * sigma);
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

inline std::vector<double> KernelSmoothing(const std::vector<double> &r,
                                           const std::vector<double> &y,
                                           double sigma, KernelType type) {
  if (r.size() != y.size() || r.empty()) {
    return {}; // Return empty vector for invalid input
  }
  if (sigma <= 0.0) {
    throw std::invalid_argument("Smoothing sigma must be positive.");
  }

  const size_t n = r.size();
  const double dx = (n > 1) ? (r[1] - r[0]) : 0;
  if (dx <= 0) {
    throw std::invalid_argument(
        "Input 'r' must be uniformly increasing for smoothing.");
  }

  // Calculate the kernel radius in bins by dividing the physical sigma by dx.
  // We use 4.0 * sigma to cover approximately 99.99% of the Gaussian area.
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

  // Apply the convolution.
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < kernel_size; ++j) {
      const long long data_idx = static_cast<long long>(i) +
                                 static_cast<long long>(j) -
                                 static_cast<long long>(kernel_radius);

      // Handle boundary conditions: clamp the index to [0, n-1].
      // This ensures safe access and uses the boundary values for convolution
      // with kernel points outside the data range.
      const long long clamped_idx = std::clamp<long long>(data_idx, 0, n - 1);

      smoothed[i] += y[clamped_idx] * kernel[j];
    }
  }

  return smoothed;
}
