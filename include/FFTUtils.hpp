// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace FFTUtils {

/**
 * @brief Simple Radix-2 Cooley-Tukey FFT implementation.
 * @param a The complex array to evaluate. Mutated in-place.
 * @param invert If true, computes the Inverse FFT.
 */
inline void computeFFT(std::vector<std::complex<double>>& a, bool invert) {
    size_t n = a.size();
    if (n == 0) return;
    
    // Ensure n is a power of 2
    if ((n & (n - 1)) != 0) {
        throw std::invalid_argument("FFT length must be a power of 2");
    }

    // Bit-reversal permutation
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

    // Cooley-Tukey
    const double pi = std::acos(-1.0);
    for (size_t len = 2; len <= n; len <<= 1) {
        double angle = 2.0 * pi / len * (invert ? -1 : 1);
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
        for (std::complex<double>& x : a) {
            x /= n;
        }
    }
}

/**
 * @brief Computes the autocorrelation of a 1D real sequence using FFT.
 * Returns the correlation for positive lags from 0 to x.size()-1.
 */
inline std::vector<double> autocorrelate(const std::vector<double>& x) {
    size_t n = x.size();
    if (n == 0) return {};
    
    // Pad to power of 2 >= 2*n - 1 to prevent circular correlation wrap-around
    size_t len = 1;
    while (len < 2 * n) len <<= 1;
    
    std::vector<std::complex<double>> cx(len, {0.0, 0.0});
    for(size_t i = 0; i < n; ++i) cx[i] = {x[i], 0.0};
    
    computeFFT(cx, false);
    
    for(size_t i = 0; i < len; ++i) {
        double mag = std::abs(cx[i]);
        cx[i] = {mag * mag, 0.0};
    }
    
    computeFFT(cx, true);
    
    std::vector<double> result(n);
    for(size_t i = 0; i < n; ++i) {
        result[i] = cx[i].real();
    }
    return result;
}

} // namespace FFTUtils
