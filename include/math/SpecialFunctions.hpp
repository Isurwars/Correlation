// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "math/SpecialFunctions.hpp"
#include "math/Constants.hpp"
#include <cmath>

namespace correlation::math::special {

/**
 * @brief Inline helper for factorials up to 20.
 * Beyond 20, the values exceed the capacity of a 64-bit integer,
 * so we return doubles for consistency.
 */
inline double factorial(int n) {
  static const double fact[] = {1.0,
                                1.0,
                                2.0,
                                6.0,
                                24.0,
                                120.0,
                                720.0,
                                5040.0,
                                40320.0,
                                362880.0,
                                3628800.0,
                                39916800.0,
                                479001600.0,
                                6227020800.0,
                                87178291200.0,
                                1307674368000.0,
                                20922789888000.0,
                                355687428096000.0,
                                6402373705728000.0,
                                121645100408832000.0,
                                2432902008176640000.0};
  if (n < 0)
    return 0.0;
  if (n <= 20)
    return fact[n];
  return std::tgamma(n + 1.0);
}

/**
 * @brief Computes the spherical associated Legendre polynomial
 * corresponding to std::sph_legendre from C++17.
 *
 * Specifically, computes Y_l^m(theta, 0) without the Condon-Shortley phase.
 */
inline double sph_legendre(int l, int m, double theta) {
  if (m < 0 || m > l) {
    return 0.0;
  }

  double x = std::cos(theta);

  // Compute P_m^m(x)
  double p_mm = 1.0;
  if (m > 0) {
    double somx2 = std::sqrt((1.0 - x) * (1.0 + x));
    double fact = 1.0;
    for (int i = 1; i <= m; ++i) {
      p_mm *= -fact * somx2;
      fact += 2.0;
    }
  }

  double p_lm = p_mm;

  if (l != m) {
    // Compute P_{m+1}^m(x)
    double p_mp1m = x * (2 * m + 1) * p_mm;
    p_lm = p_mp1m;

    if (l > m + 1) {
      // Compute P_l^m(x) using recurrence relation
      double p_ll_2 = p_mm;   // P_{l-2}^m
      double p_ll_1 = p_mp1m; // P_{l-1}^m

      for (int ll = m + 2; ll <= l; ++ll) {
        p_lm = (x * (2 * ll - 1) * p_ll_1 - (ll + m - 1) * p_ll_2) / (ll - m);
        p_ll_2 = p_ll_1;
        p_ll_1 = p_lm;
      }
    }
  }

  // Normalization factor
  double norm = std::sqrt((2.0 * l + 1.0) / (4.0 * constants::pi) *
                          factorial(l - m) / factorial(l + m));

  // Cancel Condon-Shortley phase to match std::sph_legendre
  if (m % 2 != 0) {
    p_lm = -p_lm;
  }

  return p_lm * norm;
}

/**
 * @brief Vectorized version of sph_legendre.
 * Computes the polynomial for a range of angles.
 */
inline void sph_legendre_batch(int l, int m, const double *theta,
                               double *results, size_t count) {
  for (size_t i = 0; i < count; ++i) {
    results[i] = sph_legendre(l, m, theta[i]);
  }
  // TODO: Add SIMD optimization for specific platforms if needed.
}

} // namespace correlation::math::special
