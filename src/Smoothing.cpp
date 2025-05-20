// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "../include/Smoothing.hpp"

#include <cmath> // For exp and pow

#include "../include/Constants.hpp"
#include "../include/Templates.hpp"

std::vector<double> GaussianKernel(const std::vector<double> &r_vals,
          double r0, double sigma) {
  // Normalized Gaussian Kernel
  double a = 1 / (std::sqrt(2 * constants::pi) * sigma);
  double b = -1 / (2 * sigma * sigma);
  std::vector<double> aux(r_vals.size());
  for (std::size_t i = 0; i < r_vals.size(); ++i) {
    aux[i] = a * std::exp(b * std::pow(r_vals[i] - r0, 2));
  }
  return aux;
} // GaussianKernel

std::vector<double> BumpKernel(const std::vector<double> &r_vals,
            double r0, double radius) {
  // Bump Function
  std::vector<double> aux(r_vals.size(), 0);
  for (std::size_t i = 0; i < r_vals.size(); ++i) {
    double x = (r_vals[i] - r0) / radius;
    if (std::abs(x) < 1) {
      aux[i] = std::exp(-1 / (1 - std::pow(x, 2)));
    } else {
      aux[i] = 0.0;
    }
  }
  return aux;
}

// Triweight Kernel function
std::vector<double> TriweightKernel(const std::vector<double> &r_vals,
           double r0, double radius) {
  std::vector<double> aux(r_vals.size(), 0);
  constexpr double factor = 35.0 / 32.0;

  for (std::size_t i = 0; i < r_vals.size(); ++i) {
    double u = (r_vals[i] - r0) / radius;
    if (std::abs(u) <= 1) {
      double u2 = u * u;
      aux[i] = factor * std::pow(1 - u2, 3);
    } else {
      aux[i] = 0.0;
    }
  }
  return aux;
}

std::vector<double> KernelSmoothing(const std::vector<double> &r,
            const std::vector<double> &y, double sigma,
            int _kernel_) {
  // Kernel Smoothing
  int n_ = r.size();
  std::vector<double> kernel_at_pos(n_, 0);
  std::vector<double> smoothed(n_, 0);
  for (int i = 0; i < n_; i++) {
    if (_kernel_ == 1) {
      kernel_at_pos = GaussianKernel(r, r[i], sigma);
      kernel_at_pos = NormalizeVector(kernel_at_pos);
    } else if (_kernel_ == 2) {
      kernel_at_pos = BumpKernel(r, r[i], sigma);
      kernel_at_pos = NormalizeVector(kernel_at_pos);
    } else if (_kernel_ == 3) {
      kernel_at_pos = TriweightKernel(r, r[i], sigma);
      kernel_at_pos = NormalizeVector(kernel_at_pos);
    }
    for (int j = 0; j < n_; j++) {
      smoothed[i] += kernel_at_pos[j] * y[j];
    }
  }
  return smoothed;
} // Kernel Smoothing
