/**
 * @file SpecialFunctions.hpp
 * @brief Implementation of specialized mathematical functions (factorials, Legendre polynomials).
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/Constants.hpp"
#include "math/SIMDConfig.hpp"
#include <array>
#include <cmath>

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
  static const std::array<double, 21> fact = {1.0,
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
  if (n < 0) {
    return 0.0;
  }
  if (n <= 20) {
    return fact.at(static_cast<size_t>(n));
  }
  return std::tgamma(n + 1.0);
}

/**
 * @brief Parameters for Legendre polynomial evaluation (degree and order).
 */
struct LegendreParams {
  int degree;
  int order;
};

/**
 * @brief Computes the spherical associated Legendre polynomial
 * corresponding to std::sph_legendre from C++17.
 *
 * Specifically, computes Y_l^m(theta, 0) without the Condon-Shortley phase.
 *
 * @param params Degree and order parameters.
 * @param theta The colatitudinal angle in radians.
 * @return The evaluated spherical Legendre polynomial.
 */
inline double sph_legendre(LegendreParams params, double theta) {
  const int degree = params.degree;
  const int order = params.order;
  if (order < 0 || order > degree) {
    return 0.0;
  }

  double cos_theta = std::cos(theta);

  // Compute P_m^m(x)
  double p_mm = 1.0;
  if (order > 0) {
    double somx2 = std::sqrt((1.0 - cos_theta) * (1.0 + cos_theta));
    double fact = 1.0;
    for (int i = 1; i <= order; ++i) {
      p_mm *= -fact * somx2;
      fact += 2.0;
    }
  }

  double p_lm = p_mm;

  if (degree != order) {
    // Compute P_{m+1}^m(x)
    double p_mp1m = cos_theta * (2 * order + 1) * p_mm;
    p_lm = p_mp1m;

    if (degree > order + 1) {
      // Compute P_l^m(x) using recurrence relation
      double p_ll_2 = p_mm;   // P_{l-2}^m
      double p_ll_1 = p_mp1m; // P_{l-1}^m

      for (int ll = order + 2; ll <= degree; ++ll) {
        p_lm = (cos_theta * (2 * ll - 1) * p_ll_1 - (ll + order - 1) * p_ll_2) / (ll - order);
        p_ll_2 = p_ll_1;
        p_ll_1 = p_lm;
      }
    }
  }

  // Normalization factor
  double norm = std::sqrt((2.0 * degree + 1.0) / (4.0 * correlation::math::pi) * factorial(degree - order) /
                          factorial(degree + order));

  // Cancel Condon-Shortley phase to match std::sph_legendre
  if (order % 2 != 0) {
    p_lm = -p_lm;
  }

  return p_lm * norm;
}

/**
 * @brief Vectorized version of sph_legendre.
 * Computes the polynomial for a range of angles.
 *
 * @param params Degree and order parameters.
 * @param theta Array of colatitudinal angles.
 * @param results Output array where computed polynomials will be stored.
 * @param count Number of angles to process.
 */
inline void sph_legendre_batch(LegendreParams params, const double *CORRELATION_RESTRICT theta,
                               double *CORRELATION_RESTRICT results, size_t count) {
  const int degree = params.degree;
  const int order = params.order;
  if (order < 0 || order > degree) {
    for (size_t idx = 0; idx < count; ++idx) {
      results[idx] = 0.0;
    }
    return;
  }

  // Precompute normalization factor and Condon-Shortley phase
  double norm = std::sqrt((2.0 * degree + 1.0) / (4.0 * correlation::math::pi) * factorial(degree - order) /
                          factorial(degree + order));
  if (order % 2 != 0) {
    norm = -norm;
  }

#if defined(CORRELATION_SIMD_AVX512)
  size_t idx = 0;
  const __m512d vnorm = _mm512_set1_pd(norm);
  for (; idx + 8 <= count; idx += 8) {
    std::array<double, 8> x_arr{};
    std::array<double, 8> somx2_arr{};
    for (size_t k = 0; k < 8; ++k) {
      x_arr.at(k) = std::cos(theta[idx + k]);
      somx2_arr.at(k) = std::sqrt((1.0 - x_arr.at(k)) * (1.0 + x_arr.at(k)));
    }
    __m512d v_x = _mm512_loadu_pd(x_arr.data());
    __m512d v_somx2 = _mm512_loadu_pd(somx2_arr.data());

    __m512d vp_mm = _mm512_set1_pd(1.0);
    if (order > 0) {
      double fact = 1.0;
      for (int k = 1; k <= order; ++k) {
        vp_mm = _mm512_mul_pd(vp_mm, _mm512_mul_pd(_mm512_set1_pd(-fact), v_somx2));
        fact += 2.0;
      }
    }

    __m512d vp_lm = vp_mm;
    if (degree != order) {
      __m512d vp_mp1m = _mm512_mul_pd(_mm512_mul_pd(v_x, _mm512_set1_pd(2.0 * order + 1.0)), vp_mm);
      vp_lm = vp_mp1m;

      if (degree > order + 1) {
        __m512d vp_ll_2 = vp_mm;
        __m512d vp_ll_1 = vp_mp1m;

        for (int ll = order + 2; ll <= degree; ++ll) {
          __m512d term1 = _mm512_mul_pd(_mm512_mul_pd(v_x, _mm512_set1_pd(2.0 * ll - 1.0)), vp_ll_1);
          __m512d term2 = _mm512_mul_pd(_mm512_set1_pd(ll + order - 1.0), vp_ll_2);
          __m512d diff = _mm512_sub_pd(term1, term2);
          vp_lm = _mm512_mul_pd(diff, _mm512_set1_pd(1.0 / (ll - order)));
          vp_ll_2 = vp_ll_1;
          vp_ll_1 = vp_lm;
        }
      }
    }

    __m512d vres = _mm512_mul_pd(vp_lm, vnorm);
    _mm512_storeu_pd(results + idx, vres);
  }
  for (; idx < count; ++idx) {
    results[idx] = sph_legendre({.degree = degree, .order = order}, theta[idx]);
  }
#elif defined(CORRELATION_SIMD_AVX2)
  size_t idx = 0;
  const __m256d vnorm = _mm256_set1_pd(norm);
  for (; idx + 4 <= count; idx += 4) {
    std::array<double, 4> x_arr{};
    std::array<double, 4> somx2_arr{};
    for (size_t k = 0; k < 4; ++k) {
      x_arr.at(k) = std::cos(theta[idx + k]);
      somx2_arr.at(k) = std::sqrt((1.0 - x_arr.at(k)) * (1.0 + x_arr.at(k)));
    }
    __m256d v_x = _mm256_loadu_pd(x_arr.data());
    __m256d v_somx2 = _mm256_loadu_pd(somx2_arr.data());

    __m256d vp_mm = _mm256_set1_pd(1.0);
    if (order > 0) {
      double fact = 1.0;
      for (int k = 1; k <= order; ++k) {
        vp_mm = _mm256_mul_pd(vp_mm, _mm256_mul_pd(_mm256_set1_pd(-fact), v_somx2));
        fact += 2.0;
      }
    }

    __m256d vp_lm = vp_mm;
    if (degree != order) {
      __m256d vp_mp1m = _mm256_mul_pd(_mm256_mul_pd(v_x, _mm256_set1_pd(2.0 * order + 1.0)), vp_mm);
      vp_lm = vp_mp1m;

      if (degree > order + 1) {
        __m256d vp_ll_2 = vp_mm;
        __m256d vp_ll_1 = vp_mp1m;

        for (int ll = order + 2; ll <= degree; ++ll) {
          __m256d term1 = _mm256_mul_pd(_mm256_mul_pd(v_x, _mm256_set1_pd(2.0 * ll - 1.0)), vp_ll_1);
          __m256d term2 = _mm256_mul_pd(_mm256_set1_pd(ll + order - 1.0), vp_ll_2);
          __m256d diff = _mm256_sub_pd(term1, term2);
          vp_lm = _mm256_mul_pd(diff, _mm256_set1_pd(1.0 / (ll - order)));
          vp_ll_2 = vp_ll_1;
          vp_ll_1 = vp_lm;
        }
      }
    }

    __m256d vres = _mm256_mul_pd(vp_lm, vnorm);
    _mm256_storeu_pd(results + idx, vres);
  }
  for (; idx < count; ++idx) {
    results[idx] = sph_legendre({.degree = degree, .order = order}, theta[idx]);
  }
#else
  for (size_t idx = 0; idx < count; ++idx) {
    results[idx] = sph_legendre({.degree = degree, .order = order}, theta[idx]);
  }
#endif
}

} // namespace correlation::math
