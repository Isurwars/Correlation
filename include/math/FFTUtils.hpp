// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "math/Constants.hpp"
#include <algorithm>
#include <complex>
#include <stdexcept>
#include <vector>

namespace correlation::math::fft {

/**
 * @brief Simple Radix-2 Cooley-Tukey FFT implementation with bit-reversal
 * optimization.
 */
inline void computeFFT(std::vector<std::complex<double>> &a, bool invert) {
  size_t n = a.size();
  if (n == 0)
    return;

  if ((n & (n - 1)) != 0) {
    throw std::invalid_argument("FFT length must be a power of 2");
  }

  // Bit-reversal permutation (Optimized)
  for (size_t i = 1, j = 0; i < n; i++) {
    size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) {
      j ^= bit;
    }
    j ^= bit;
    if (i < j) {
      std::swap(a[i], a[j]);
    }
  }

  // Cooley-Tukey with precomputed twiddle basics
  for (size_t len = 2; len <= n; len <<= 1) {
    double angle = 2.0 * constants::pi / len * (invert ? -1 : 1);
    std::complex<double> wlen(std::cos(angle), std::sin(angle));
    for (size_t i = 0; i < n; i += len) {
      std::complex<double> w(1.0, 0.0);
      for (size_t j = 0; j < len / 2; j++) {
        std::complex<double> u = a[i + j];
        std::complex<double> v = a[i + j + len / 2] * w;
        a[i + j] = u + v;
        a[i + j + len / 2] = u - v;
        w *= wlen;
      }
    }
  }

  if (invert) {
    for (auto &x : a) {
      x /= static_cast<double>(n);
    }
  }
}

/**
 * @brief Computes the autocorrelation of a 1D real sequence using FFT.
 */
inline std::vector<double> autocorrelate(const std::vector<double> &x) {
  size_t n = x.size();
  if (n == 0)
    return {};

  size_t len = 1;
  while (len < 2 * n)
    len <<= 1;

  std::vector<std::complex<double>> cx(len, {0.0, 0.0});
  for (size_t i = 0; i < n; ++i)
    cx[i] = {x[i], 0.0};

  computeFFT(cx, false);

  for (size_t i = 0; i < len; ++i) {
    double mag_sq = cx[i].real() * cx[i].real() + cx[i].imag() * cx[i].imag();
    cx[i] = {mag_sq, 0.0};
  }

  computeFFT(cx, true);

  std::vector<double> result(n);
  for (size_t i = 0; i < n; ++i) {
    result[i] = cx[i].real();
  }
  return result;
}

} // namespace correlation::math::fft
