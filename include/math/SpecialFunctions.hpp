/**
 * @file SpecialFunctions.hpp
 * @brief Implementation of specialized mathematical functions (factorials, Legendre polynomials).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @license SPDX-License-Identifier: MIT
 */

#pragma once

#include "math/SIMDConfig.hpp"
#include <cmath>
#include <numbers>

namespace correlation::math {

/**
 * @brief Inline helper for factorials up to 20.
 * Beyond 20, the values exceed the capacity of a 64-bit integer,
 * so we return doubles for consistency.
 * 
 * @param n the number to compute factorial for.
 * @return The factorial of n.
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
 * 
 * @param l The degree of the polynomial.
 * @param m The order of the polynomial.
 * @param theta The colatitudinal angle in radians.
 * @return The evaluated spherical Legendre polynomial.
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
  double norm = std::sqrt((2.0 * l + 1.0) / (4.0 * std::numbers::pi) *
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
 * 
 * @param l The degree of the polynomial.
 * @param m The order of the polynomial.
 * @param theta Array of colatitudinal angles.
 * @param results Output array where computed polynomials will be stored.
 * @param count Number of angles to process.
 */
inline void sph_legendre_batch(int l, int m, const double *CORRELATION_RESTRICT theta,
                               double *CORRELATION_RESTRICT results, size_t count) {
  if (m < 0 || m > l) {
    for (size_t i = 0; i < count; ++i) {
      results[i] = 0.0;
    }
    return;
  }

  // Precompute normalization factor and Condon-Shortley phase
  double norm = std::sqrt((2.0 * l + 1.0) / (4.0 * std::numbers::pi) *
                          factorial(l - m) / factorial(l + m));
  if (m % 2 != 0) {
    norm = -norm;
  }

#if defined(CORRELATION_SIMD_AVX512)
  size_t i = 0;
  const __m512d vnorm = _mm512_set1_pd(norm);
  for (; i + 8 <= count; i += 8) {
    double x_arr[8], somx2_arr[8];
    for (int k = 0; k < 8; ++k) {
      x_arr[k] = std::cos(theta[i + k]);
      somx2_arr[k] = std::sqrt((1.0 - x_arr[k]) * (1.0 + x_arr[k]));
    }
    __m512d vx = _mm512_loadu_pd(x_arr);
    __m512d vsomx2 = _mm512_loadu_pd(somx2_arr);

    __m512d vp_mm = _mm512_set1_pd(1.0);
    if (m > 0) {
      double fact = 1.0;
      for (int k = 1; k <= m; ++k) {
        vp_mm = _mm512_mul_pd(vp_mm, _mm512_mul_pd(_mm512_set1_pd(-fact), vsomx2));
        fact += 2.0;
      }
    }

    __m512d vp_lm = vp_mm;
    if (l != m) {
      __m512d vp_mp1m = _mm512_mul_pd(_mm512_mul_pd(vx, _mm512_set1_pd(2.0 * m + 1.0)), vp_mm);
      vp_lm = vp_mp1m;

      if (l > m + 1) {
        __m512d vp_ll_2 = vp_mm;
        __m512d vp_ll_1 = vp_mp1m;

        for (int ll = m + 2; ll <= l; ++ll) {
          __m512d term1 = _mm512_mul_pd(_mm512_mul_pd(vx, _mm512_set1_pd(2.0 * ll - 1.0)), vp_ll_1);
          __m512d term2 = _mm512_mul_pd(_mm512_set1_pd(ll + m - 1.0), vp_ll_2);
          __m512d diff = _mm512_sub_pd(term1, term2);
          vp_lm = _mm512_mul_pd(diff, _mm512_set1_pd(1.0 / (ll - m)));
          vp_ll_2 = vp_ll_1;
          vp_ll_1 = vp_lm;
        }
      }
    }
    
    __m512d vres = _mm512_mul_pd(vp_lm, vnorm);
    _mm512_storeu_pd(results + i, vres);
  }
  for (; i < count; ++i) {
    results[i] = sph_legendre(l, m, theta[i]);
  }
#elif defined(CORRELATION_SIMD_AVX2)
  size_t i = 0;
  const __m256d vnorm = _mm256_set1_pd(norm);
  for (; i + 4 <= count; i += 4) {
    double x_arr[4], somx2_arr[4];
    for (int k = 0; k < 4; ++k) {
      x_arr[k] = std::cos(theta[i + k]);
      somx2_arr[k] = std::sqrt((1.0 - x_arr[k]) * (1.0 + x_arr[k]));
    }
    __m256d vx = _mm256_loadu_pd(x_arr);
    __m256d vsomx2 = _mm256_loadu_pd(somx2_arr);

    __m256d vp_mm = _mm256_set1_pd(1.0);
    if (m > 0) {
      double fact = 1.0;
      for (int k = 1; k <= m; ++k) {
        vp_mm = _mm256_mul_pd(vp_mm, _mm256_mul_pd(_mm256_set1_pd(-fact), vsomx2));
        fact += 2.0;
      }
    }

    __m256d vp_lm = vp_mm;
    if (l != m) {
      __m256d vp_mp1m = _mm256_mul_pd(_mm256_mul_pd(vx, _mm256_set1_pd(2.0 * m + 1.0)), vp_mm);
      vp_lm = vp_mp1m;

      if (l > m + 1) {
        __m256d vp_ll_2 = vp_mm;
        __m256d vp_ll_1 = vp_mp1m;

        for (int ll = m + 2; ll <= l; ++ll) {
          __m256d term1 = _mm256_mul_pd(_mm256_mul_pd(vx, _mm256_set1_pd(2.0 * ll - 1.0)), vp_ll_1);
          __m256d term2 = _mm256_mul_pd(_mm256_set1_pd(ll + m - 1.0), vp_ll_2);
          __m256d diff = _mm256_sub_pd(term1, term2);
          vp_lm = _mm256_mul_pd(diff, _mm256_set1_pd(1.0 / (ll - m)));
          vp_ll_2 = vp_ll_1;
          vp_ll_1 = vp_lm;
        }
      }
    }

    __m256d vres = _mm256_mul_pd(vp_lm, vnorm);
    _mm256_storeu_pd(results + i, vres);
  }
  for (; i < count; ++i) {
    results[i] = sph_legendre(l, m, theta[i]);
  }
#else
  for (size_t i = 0; i < count; ++i) {
    results[i] = sph_legendre(l, m, theta[i]);
  }
#endif
}

} // namespace correlation::math
