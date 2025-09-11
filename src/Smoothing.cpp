// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "../include/Smoothing.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "../include/Constants.hpp"

namespace { // Anonymous namespace for internal linkage

// A helper function to generate a normalized kernel vector centered at zero.
std::vector<double> generateKernel(size_t size, double dx, double sigma,
                                   KernelType type) {
  std::vector<double> kernel(size);
  const double center = static_cast<double>(size - 1) / 2.0;

  if (type == KernelType::Gaussian) {
    const double a = 1.0 / (std::sqrt(2 * constants::pi) * sigma);
    const double b = -1.0 / (2.0 * sigma * sigma);
    for (size_t i = 0; i < size; ++i) {
      const double x = (static_cast<double>(i) - center) * dx;
      kernel[i] = a * std::exp(b * x * x);
    }
  } else if (type == KernelType::Triweight) {
    constexpr double factor = 35.0 / 32.0;
    for (size_t i = 0; i < size; ++i) {
      const double u = ((static_cast<double>(i) - center) * dx) / sigma;
      if (std::abs(u) <= 1.0) {
        const double u2 = u * u;
        kernel[i] = factor * std::pow(1 - u2, 3);
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
} // namespace

std::vector<double> KernelSmoothing(const std::vector<double> &r,
                                    const std::vector<double> &y, double sigma,
                                    KernelType type) {
  if (r.size() != y.size() || r.empty()) {
    return {}; // Return empty vector for invalid input
  }

  const size_t n = r.size();
  const double dx = (n > 1) ? (r[1] - r[0]) : 0;
  if (dx <= 0) {
    throw std::invalid_argument(
        "Input 'r' must be uniformly increasing for smoothing.");
  }

  const size_t kernel_radius =
      std::min(n / 2, static_cast<size_t>(4.0 * sigma / dx));
  const size_t kernel_size = 2 * kernel_radius + 1;

  const auto kernel = generateKernel(kernel_size, dx, sigma, type);
  std::vector<double> smoothed(n, 0.0);

  // Apply the convolution.
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < kernel_size; ++j) {
      const long long data_idx = static_cast<long long>(i) +
                                 static_cast<long long>(j) -
                                 static_cast<long long>(kernel_radius);

      // Handle boundary conditions (points near the start/end of the data)
      const long long clamped_idx = std::clamp<long long>(data_idx, 0, n - 1);

      smoothed[i] += y[clamped_idx] * kernel[j];
    }
  }

  return smoothed;
}
