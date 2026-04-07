/**
 * @file FFTUtils.hpp
 * @brief Fast Fourier Transform (FFT) and autocorrelation utilities.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "math/Constants.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

namespace correlation::math {

/**
 * @brief Simple Radix-2 Cooley-Tukey FFT implementation with bit-reversal
 * optimization.
 * 
 * @param a The input/output vector of complex numbers. Modifies in-place.
 * @param invert If true, performs an inverse FFT and scales the result by 1/N.
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
    double angle = correlation::math::two_pi / len * (invert ? -1 : 1);
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
 *
 * @param x         Input signal.
 * @param workspace Reusable scratch buffer. Resized automatically when needed.
 *                  Pass the same buffer across repeated calls (e.g. per-atom
 *                  loops) to avoid repeated heap allocation.
 * @return Autocorrelation array of length n.
 */
inline std::vector<double>
autocorrelate(const std::vector<double> &x,
              std::vector<std::complex<double>> &workspace) {
  const size_t n = x.size();
  if (n == 0)
    return {};

  size_t len = 1;
  while (len < 2 * n)
    len <<= 1;

  // Reuse caller-supplied workspace; only reallocates when length grows.
  workspace.assign(len, {0.0, 0.0});
  for (size_t i = 0; i < n; ++i)
    workspace[i] = {x[i], 0.0};

  computeFFT(workspace, false);

  // In-place |X[k]|² — no intermediate vector needed.
  for (size_t i = 0; i < len; ++i) {
    const double re = workspace[i].real();
    const double im = workspace[i].imag();
    workspace[i] = {re * re + im * im, 0.0};
  }

  computeFFT(workspace, true);

  std::vector<double> result(n);
  for (size_t i = 0; i < n; ++i)
    result[i] = workspace[i].real();
  return result;
}

/**
 * @brief Convenience overload — allocates its own workspace.
 *
 * Prefer the two-argument overload in tight loops to recycle the allocation.
 */
inline std::vector<double> autocorrelate(const std::vector<double> &x) {
  std::vector<std::complex<double>> workspace;
  return autocorrelate(x, workspace);
}

} // namespace correlation::math
