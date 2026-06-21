/**
 * @file FFTUtils.hpp
 * @brief Fast Fourier Transform (FFT) and autocorrelation utilities.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/Constants.hpp"
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

namespace correlation::math {

/**
 * @brief Simple Radix-2 Cooley-Tukey FFT implementation with bit-reversal
 * optimization.
 *
 * @param data The input/output vector of complex numbers. Modifies in-place.
 * @param invert If true, performs an inverse FFT and scales the result by 1/N.
 */
inline void computeFFT(std::vector<std::complex<double>> &data, bool invert) {
  size_t size = data.size();
  if (size == 0) {
    return;
  }

  if ((size & (size - 1)) != 0) {
    throw std::invalid_argument("FFT length must be a power of 2");
  }

  // Bit-reversal permutation (Optimized)
  for (size_t i = 1, j = 0; i < size; i++) {
    size_t bit = size >> 1;
    for (; (j & bit) != 0; bit >>= 1) {
      j ^= bit;
    }
    j ^= bit;
    if (i < j) {
      std::swap(data[i], data[j]);
    }
  }

  // Cooley-Tukey with precomputed twiddle basics
  for (size_t len = 2; len <= size; len <<= 1) {
    double angle = correlation::math::two_pi / static_cast<double>(len) * (invert ? -1 : 1);
    std::complex<double> wlen(std::cos(angle), std::sin(angle));
    for (size_t i = 0; i < size; i += len) {
      std::complex<double> twiddle(1.0, 0.0);
      for (size_t j = 0; j < len / 2; j++) {
        std::complex<double> even_val = data[i + j];
        std::complex<double> odd_val = data[i + j + len / 2] * twiddle;
        data[i + j] = even_val + odd_val;
        data[i + j + len / 2] = even_val - odd_val;
        twiddle *= wlen;
      }
    }
  }

  if (invert) {
    for (auto &val : data) {
      val /= static_cast<double>(size);
    }
  }
}

/**
 * @brief Computes the autocorrelation of a 1D real sequence using FFT.
 *
 * @param signal    Input signal.
 * @param workspace Reusable scratch buffer. Resized automatically when needed.
 *                  Pass the same buffer across repeated calls (e.g. per-atom
 *                  loops) to avoid repeated heap allocation.
 * @return Autocorrelation array of length n.
 */
inline std::vector<double> autocorrelate(const std::vector<double> &signal,
                                         std::vector<std::complex<double>> &workspace) {
  const size_t size = signal.size();
  if (size == 0) {
    return {};
  }

  size_t len = 1;
  while (len < 2 * size) {
    len <<= 1;
  }

  // Reuse caller-supplied workspace; only reallocates when length grows.
  workspace.assign(len, {0.0, 0.0});
  for (size_t i = 0; i < size; ++i) {
    workspace[i] = {signal[i], 0.0};
  }

  computeFFT(workspace, false);

  // In-place |X[k]|² — no intermediate vector needed.
  for (size_t i = 0; i < len; ++i) {
    const double real_part = workspace[i].real();
    const double imag_part = workspace[i].imag();
    workspace[i] = {real_part * real_part + imag_part * imag_part, 0.0};
  }

  computeFFT(workspace, true);

  std::vector<double> result(size);
  for (size_t i = 0; i < size; ++i) {
    result[i] = workspace[i].real();
  }
  return result;
}

/**
 * @brief Convenience overload — allocates its own workspace.
 *
 * Prefer the two-argument overload in tight loops to recycle the allocation.
 *
 * @param signal Input signal vector.
 * @return The autocorrelation vector.
 */
inline std::vector<double> autocorrelate(const std::vector<double> &signal) {
  std::vector<std::complex<double>> workspace;
  return autocorrelate(signal, workspace);
}

} // namespace correlation::math
