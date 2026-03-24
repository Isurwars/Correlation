// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

// ---------------------------------------------------------------------------
// SIMD level detection
// ---------------------------------------------------------------------------
#if defined(__AVX512F__)
#define CORRELATION_SIMD_AVX512
#elif defined(__AVX2__)
#define CORRELATION_SIMD_AVX2
#endif

#if defined(CORRELATION_SIMD_AVX512)
#include <immintrin.h>
#define CORRELATION_SIMD_WIDTH 8 // 8 doubles per AVX-512 register
#elif defined(CORRELATION_SIMD_AVX2)
#include <immintrin.h>
#define CORRELATION_SIMD_WIDTH 4 // 4 doubles per AVX2 register
#else
#define CORRELATION_SIMD_WIDTH 1 // scalar fallback
#endif

// ---------------------------------------------------------------------------
// Portable restrict keyword
// ---------------------------------------------------------------------------
#if defined(_MSC_VER)
#define CORRELATION_RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
#define CORRELATION_RESTRICT __restrict__
#else
#define CORRELATION_RESTRICT
#endif

#include <cmath>
#include <cstddef>
#include <vector>

namespace simd_utils {

// ---------------------------------------------------------------------------
// SoA block of atom positions
// Caller fills these arrays before invoking the SIMD kernel.
// ---------------------------------------------------------------------------
struct PositionBlock {
  // Aligned arrays for SoA layout. Padding to the next multiple of
  // CORRELATION_SIMD_WIDTH is the caller's responsibility.
  double *x;
  double *y;
  double *z;
  std::size_t count; // actual number of atoms (without padding)
};

// ---------------------------------------------------------------------------
// Scalar helper: squared distance between atom A (single) and atom B[k]
// ---------------------------------------------------------------------------
inline double dist_sq_scalar(double ax, double ay, double az, double bx,
                             double by, double bz) noexcept {
  double dx = bx - ax;
  double dy = by - ay;
  double dz = bz - az;
  return dx * dx + dy * dy + dz * dz;
}

// ---------------------------------------------------------------------------
// Core SIMD kernel: compute squared distances between a single atom A
// and a contiguous block of atoms B stored in SoA layout.
//
// out_dsq[k] = ||A - B[k]||^2  for k in [0, block.count)
//
// The caller must ensure out_dsq has at least block.count elements.
// Masked processing handles tails so no out-of-bounds reads occur.
// ---------------------------------------------------------------------------

#if defined(CORRELATION_SIMD_AVX512)

inline void compute_dsq_block(double ax, double ay, double az,
                              const PositionBlock &block,
                              double *CORRELATION_RESTRICT out_dsq) noexcept {
  const __m512d va_x = _mm512_set1_pd(ax);
  const __m512d va_y = _mm512_set1_pd(ay);
  const __m512d va_z = _mm512_set1_pd(az);

  std::size_t k = 0;
  // Process 8 atoms at a time
  for (; k + 8 <= block.count; k += 8) {
    __m512d dx = _mm512_sub_pd(_mm512_loadu_pd(block.x + k), va_x);
    __m512d dy = _mm512_sub_pd(_mm512_loadu_pd(block.y + k), va_y);
    __m512d dz = _mm512_sub_pd(_mm512_loadu_pd(block.z + k), va_z);
    __m512d dsq =
        _mm512_fmadd_pd(dx, dx, _mm512_fmadd_pd(dy, dy, _mm512_mul_pd(dz, dz)));
    _mm512_storeu_pd(out_dsq + k, dsq);
  }
  // Masked tail: remaining atoms (0-7)
  if (k < block.count) {
    // Build mask for the remaining lanes
    __mmask8 mask = static_cast<__mmask8>((1u << (block.count - k)) - 1u);
    __m512d dx = _mm512_sub_pd(_mm512_maskz_loadu_pd(mask, block.x + k), va_x);
    __m512d dy = _mm512_sub_pd(_mm512_maskz_loadu_pd(mask, block.y + k), va_y);
    __m512d dz = _mm512_sub_pd(_mm512_maskz_loadu_pd(mask, block.z + k), va_z);
    __m512d dsq =
        _mm512_fmadd_pd(dx, dx, _mm512_fmadd_pd(dy, dy, _mm512_mul_pd(dz, dz)));
    _mm512_mask_storeu_pd(out_dsq + k, mask, dsq);
  }
}

#elif defined(CORRELATION_SIMD_AVX2)

inline void compute_dsq_block(double ax, double ay, double az,
                              const PositionBlock &block,
                              double *CORRELATION_RESTRICT out_dsq) noexcept {
  const __m256d va_x = _mm256_set1_pd(ax);
  const __m256d va_y = _mm256_set1_pd(ay);
  const __m256d va_z = _mm256_set1_pd(az);

  std::size_t k = 0;
  // Process 4 atoms at a time
  for (; k + 4 <= block.count; k += 4) {
    __m256d dx = _mm256_sub_pd(_mm256_loadu_pd(block.x + k), va_x);
    __m256d dy = _mm256_sub_pd(_mm256_loadu_pd(block.y + k), va_y);
    __m256d dz = _mm256_sub_pd(_mm256_loadu_pd(block.z + k), va_z);
    // Use FMA if available (AVX2 CPUs typically have FMA3 as well)
#if defined(__FMA__)
    __m256d dsq =
        _mm256_fmadd_pd(dx, dx, _mm256_fmadd_pd(dy, dy, _mm256_mul_pd(dz, dz)));
#else
    __m256d dsq = _mm256_add_pd(
        _mm256_mul_pd(dx, dx),
        _mm256_add_pd(_mm256_mul_pd(dy, dy), _mm256_mul_pd(dz, dz)));
#endif
    _mm256_storeu_pd(out_dsq + k, dsq);
  }
  // Scalar tail
  for (; k < block.count; ++k) {
    out_dsq[k] = dist_sq_scalar(ax, ay, az, block.x[k], block.y[k], block.z[k]);
  }
}

#else // Scalar fallback

inline void compute_dsq_block(double ax, double ay, double az,
                              const PositionBlock &block,
                              double *CORRELATION_RESTRICT out_dsq) noexcept {
  for (std::size_t k = 0; k < block.count; ++k) {
    out_dsq[k] = dist_sq_scalar(ax, ay, az, block.x[k], block.y[k], block.z[k]);
  }
}

#endif

// ---------------------------------------------------------------------------
// sinc_integral kernel
//
// Computes: sum_{j=0}^{count-1} integrand[j] * sin(Q * rbins[j])
//
// Strategy (portable, no SVML/sleef dependency):
//   Phase 1 – fill a caller-provided scratch buffer sinQr[j] = sin(Q*rbins[j])
//             using a tight scalar loop.  Modern compilers (GCC/Clang -O2 +
//             -ffast-math) auto-vectorize this with the platform sin library.
//   Phase 2 – AVX2/AVX-512 FMA dot-product: acc += integrand[j] * sinQr[j].
//
// The caller must supply `sinqr_scratch` with capacity >= count.
// ---------------------------------------------------------------------------

#if defined(CORRELATION_SIMD_AVX512)

inline double sinc_integral(double Q,
                            const double *CORRELATION_RESTRICT integrand,
                            const double *CORRELATION_RESTRICT rbins,
                            double *CORRELATION_RESTRICT sinqr_scratch,
                            std::size_t count) noexcept {
  // Phase 1: scalar sin (auto-vectorized by the compiler)
  for (std::size_t j = 0; j < count; ++j) {
    sinqr_scratch[j] = std::sin(Q * rbins[j]);
  }

  // Phase 2: AVX-512 FMA dot product
  __m512d vacc = _mm512_setzero_pd();
  std::size_t j = 0;
  for (; j + 8 <= count; j += 8) {
    __m512d vi = _mm512_loadu_pd(integrand + j);
    __m512d vs = _mm512_loadu_pd(sinqr_scratch + j);
    vacc = _mm512_fmadd_pd(vi, vs, vacc);
  }
  double acc = _mm512_reduce_add_pd(vacc);
  // Scalar tail
  for (; j < count; ++j) {
    acc += integrand[j] * sinqr_scratch[j];
  }
  return acc;
}

#elif defined(CORRELATION_SIMD_AVX2)

inline double sinc_integral(double Q,
                            const double *CORRELATION_RESTRICT integrand,
                            const double *CORRELATION_RESTRICT rbins,
                            double *CORRELATION_RESTRICT sinqr_scratch,
                            std::size_t count) noexcept {
  // Phase 1: scalar sin (auto-vectorized by the compiler)
  for (std::size_t j = 0; j < count; ++j) {
    sinqr_scratch[j] = std::sin(Q * rbins[j]);
  }

  // Phase 2: AVX2 FMA dot product (4 doubles per register)
  __m256d vacc = _mm256_setzero_pd();
  std::size_t j = 0;
  for (; j + 4 <= count; j += 4) {
    __m256d vi = _mm256_loadu_pd(integrand + j);
    __m256d vs = _mm256_loadu_pd(sinqr_scratch + j);
#if defined(__FMA__)
    vacc = _mm256_fmadd_pd(vi, vs, vacc);
#else
    vacc = _mm256_add_pd(vacc, _mm256_mul_pd(vi, vs));
#endif
  }
  // Horizontal reduction of the 4-lane AVX2 accumulator
  __m128d lo = _mm256_castpd256_pd128(vacc);
  __m128d hi = _mm256_extractf128_pd(vacc, 1);
  __m128d sum2 = _mm_add_pd(lo, hi);
  __m128d sum1 = _mm_hadd_pd(sum2, sum2);
  double acc = _mm_cvtsd_f64(sum1);
  // Scalar tail
  for (; j < count; ++j) {
    acc += integrand[j] * sinqr_scratch[j];
  }
  return acc;
}

#else // Scalar fallback

inline double sinc_integral(double Q,
                            const double *CORRELATION_RESTRICT integrand,
                            const double *CORRELATION_RESTRICT rbins,
                            double *CORRELATION_RESTRICT /*sinqr_scratch*/,
                            std::size_t count) noexcept {
  double acc = 0.0;
  for (std::size_t j = 0; j < count; ++j) {
    acc += integrand[j] * std::sin(Q * rbins[j]);
  }
  return acc;
}

#endif

// ---------------------------------------------------------------------------
// debye_sum kernel
//
// Computes: sum_{j=0}^{count-1} sin(Q * distances[j]) / (Q * distances[j])
//
// Strategy:
//   Phase 1 - fill scratch buffer with sin(Q*distances[j]) using auto-vectorized loop.
//   Phase 2 - AVX2/AVX-512 division and sum.
// ---------------------------------------------------------------------------

#if defined(CORRELATION_SIMD_AVX512)

inline double debye_sum(double Q,
                        const double *CORRELATION_RESTRICT distances,
                        double *CORRELATION_RESTRICT scratch,
                        std::size_t count) noexcept {
  const double Q2 = Q * Q;
  if (Q < 1e-9) {
      return static_cast<double>(count);
  }
  const double invQ = 1.0 / Q;
  for (std::size_t j = 0; j < count; ++j) {
    scratch[j] = std::sin(Q * distances[j]);
  }

  __m512d vacc = _mm512_setzero_pd();
  const __m512d vinvQ = _mm512_set1_pd(invQ);
  std::size_t j = 0;
  for (; j + 8 <= count; j += 8) {
    __m512d vr = _mm512_loadu_pd(distances + j);
    __m512d vs = _mm512_loadu_pd(scratch + j);
    // sinc(x) = sin(x)/x. vs is sin(Q*r), vr is r. 
    // vs * invQ / vr = sin(Q*r) / (Q*r)
    __m512d term = _mm512_div_pd(_mm512_mul_pd(vs, vinvQ), vr);
    vacc = _mm512_add_pd(vacc, term);
  }
  double acc = _mm512_reduce_add_pd(vacc);
  for (; j < count; ++j) {
    const double x = Q * distances[j];
    if (x < 1e-4) {
      acc += 1.0 - (x * x) / 6.0;
    } else {
      acc += std::sin(x) / x;
    }
  }
  return acc;
}

#elif defined(CORRELATION_SIMD_AVX2)

inline double debye_sum(double Q,
                        const double *CORRELATION_RESTRICT distances,
                        double *CORRELATION_RESTRICT scratch,
                        std::size_t count) noexcept {
  if (Q < 1e-9) {
      return static_cast<double>(count);
  }
  const double invQ = 1.0 / Q;
  for (std::size_t j = 0; j < count; ++j) {
    scratch[j] = std::sin(Q * distances[j]);
  }

  __m256d vacc = _mm256_setzero_pd();
  const __m256d vinvQ = _mm256_set1_pd(invQ);
  std::size_t j = 0;
  for (; j + 4 <= count; j += 4) {
    __m256d vr = _mm256_loadu_pd(distances + j);
    __m256d vs = _mm256_loadu_pd(scratch + j);
    __m256d term = _mm256_div_pd(_mm256_mul_pd(vs, vinvQ), vr);
    vacc = _mm256_add_pd(vacc, term);
  }
  __m128d lo = _mm256_castpd256_pd128(vacc);
  __m128d hi = _mm256_extractf128_pd(vacc, 1);
  __m128d sum2 = _mm_add_pd(lo, hi);
  __m128d sum1 = _mm_hadd_pd(sum2, sum2);
  double acc = _mm_cvtsd_f64(sum1);
  for (; j < count; ++j) {
    const double x = Q * distances[j];
    if (x < 1e-4) {
      acc += 1.0 - (x * x) / 6.0;
    } else {
      acc += std::sin(x) / x;
    }
  }
  return acc;
}

#else // Scalar fallback

inline double debye_sum(double Q,
                        const double *CORRELATION_RESTRICT distances,
                        double *CORRELATION_RESTRICT /*scratch*/,
                        std::size_t count) noexcept {
  if (Q < 1e-9) {
      return static_cast<double>(count);
  }
  double acc = 0.0;
  for (std::size_t j = 0; j < count; ++j) {
    const double x = Q * distances[j];
    if (x < 1e-4) {
      acc += 1.0 - (x * x) / 6.0;
    } else {
      acc += std::sin(x) / x;
    }
  }
  return acc;
}

#endif

// ---------------------------------------------------------------------------
// normalize_rdf_bins kernel
//
// Computes for k in [1, count)  (bin 0 is left as-is / zero):
//   g[k]    = H[k] * g_norm / (r[k]^2)
//   G[k]    = pi4_rho_j * r[k] * (g[k] - 1)
//   J[k]    = H[k] * inv_Ni_dr
//   Jinv[k] = H[k] * inv_Nj_dr
//
// NOTE: J_out and Jinv_out may alias (self-pair case i==j) — no restrict.
// H and rbins are read-only and must not alias any output.
// ---------------------------------------------------------------------------

#if defined(CORRELATION_SIMD_AVX512)

inline void normalize_rdf_bins(const double *CORRELATION_RESTRICT H,
                               const double *CORRELATION_RESTRICT rbins,
                               double g_norm, double inv_Ni_dr,
                               double inv_Nj_dr, double pi4_rho_j,
                               double *g_out, double *G_out, double *J_out,
                               double *Jinv_out, std::size_t count) noexcept {
  const __m512d vg_norm = _mm512_set1_pd(g_norm);
  const __m512d v1 = _mm512_set1_pd(1.0);
  const __m512d vinNidr = _mm512_set1_pd(inv_Ni_dr);
  const __m512d vinNjdr = _mm512_set1_pd(inv_Nj_dr);
  const __m512d vpi4rho = _mm512_set1_pd(pi4_rho_j);

  std::size_t k = 1; // bin 0 stays 0
  for (; k + 8 <= count; k += 8) {
    __m512d vH = _mm512_loadu_pd(H + k);
    __m512d vr = _mm512_loadu_pd(rbins + k);
    __m512d vr2 = _mm512_mul_pd(vr, vr);
    __m512d vg = _mm512_div_pd(_mm512_mul_pd(vH, vg_norm), vr2);
    _mm512_storeu_pd(g_out + k, vg);
    _mm512_storeu_pd(
        G_out + k,
        _mm512_mul_pd(vpi4rho, _mm512_mul_pd(vr, _mm512_sub_pd(vg, v1))));
    _mm512_storeu_pd(J_out + k, _mm512_mul_pd(vH, vinNidr));
    _mm512_storeu_pd(Jinv_out + k, _mm512_mul_pd(vH, vinNjdr));
  }
  for (; k < count; ++k) {
    const double r = rbins[k];
    if (r < 1e-9)
      continue;
    const double g = H[k] * g_norm / (r * r);
    g_out[k] = g;
    G_out[k] = pi4_rho_j * r * (g - 1.0);
    J_out[k] = H[k] * inv_Ni_dr;
    Jinv_out[k] = H[k] * inv_Nj_dr;
  }
}

#elif defined(CORRELATION_SIMD_AVX2)

inline void normalize_rdf_bins(const double *CORRELATION_RESTRICT H,
                               const double *CORRELATION_RESTRICT rbins,
                               double g_norm, double inv_Ni_dr,
                               double inv_Nj_dr, double pi4_rho_j,
                               double *g_out, double *G_out, double *J_out,
                               double *Jinv_out, std::size_t count) noexcept {
  const __m256d vg_norm = _mm256_set1_pd(g_norm);
  const __m256d v1 = _mm256_set1_pd(1.0);
  const __m256d vinNidr = _mm256_set1_pd(inv_Ni_dr);
  const __m256d vinNjdr = _mm256_set1_pd(inv_Nj_dr);
  const __m256d vpi4rho = _mm256_set1_pd(pi4_rho_j);

  std::size_t k = 1; // bin 0 stays 0
  for (; k + 4 <= count; k += 4) {
    __m256d vH = _mm256_loadu_pd(H + k);
    __m256d vr = _mm256_loadu_pd(rbins + k);
    __m256d vr2 = _mm256_mul_pd(vr, vr);
    __m256d vg = _mm256_div_pd(_mm256_mul_pd(vH, vg_norm), vr2);
    _mm256_storeu_pd(g_out + k, vg);
    _mm256_storeu_pd(
        G_out + k,
        _mm256_mul_pd(vpi4rho, _mm256_mul_pd(vr, _mm256_sub_pd(vg, v1))));
    _mm256_storeu_pd(J_out + k, _mm256_mul_pd(vH, vinNidr));
    _mm256_storeu_pd(Jinv_out + k, _mm256_mul_pd(vH, vinNjdr));
  }
  // Scalar tail
  for (; k < count; ++k) {
    const double r = rbins[k];
    if (r < 1e-9)
      continue;
    const double g = H[k] * g_norm / (r * r);
    g_out[k] = g;
    G_out[k] = pi4_rho_j * r * (g - 1.0);
    J_out[k] = H[k] * inv_Ni_dr;
    Jinv_out[k] = H[k] * inv_Nj_dr;
  }
}

#else // Scalar fallback

inline void normalize_rdf_bins(const double *CORRELATION_RESTRICT H,
                               const double *CORRELATION_RESTRICT rbins,
                               double g_norm, double inv_Ni_dr,
                               double inv_Nj_dr, double pi4_rho_j,
                               double *g_out, double *G_out, double *J_out,
                               double *Jinv_out, std::size_t count) noexcept {
  for (std::size_t k = 1; k < count; ++k) {
    const double r = rbins[k];
    if (r < 1e-9)
      continue;
    const double g = H[k] * g_norm / (r * r);
    g_out[k] = g;
    G_out[k] = pi4_rho_j * r * (g - 1.0);
    J_out[k] = H[k] * inv_Ni_dr;
    Jinv_out[k] = H[k] * inv_Nj_dr;
  }
}

#endif

// ---------------------------------------------------------------------------
// scale_bins: multiply every element of arr[0..count) by scalar `s`
// Used for the self-pair H(r) x2 correction.
// Auto-vectorized by the compiler; explicit AVX2 for clarity and alignment
// with the rest of this header.
// ---------------------------------------------------------------------------

#if defined(CORRELATION_SIMD_AVX512)

inline void scale_bins(double *arr, double s, std::size_t count) noexcept {
  const __m512d vs = _mm512_set1_pd(s);
  std::size_t k = 0;
  for (; k + 8 <= count; k += 8)
    _mm512_storeu_pd(arr + k, _mm512_mul_pd(_mm512_loadu_pd(arr + k), vs));
  for (; k < count; ++k)
    arr[k] *= s;
}

#elif defined(CORRELATION_SIMD_AVX2)

inline void scale_bins(double *arr, double s, std::size_t count) noexcept {
  const __m256d vs = _mm256_set1_pd(s);
  std::size_t k = 0;
  for (; k + 4 <= count; k += 4)
    _mm256_storeu_pd(arr + k, _mm256_mul_pd(_mm256_loadu_pd(arr + k), vs));
  for (; k < count; ++k)
    arr[k] *= s;
}

#else

inline void scale_bins(double *arr, double s, std::size_t count) noexcept {
  for (std::size_t k = 0; k < count; ++k)
    arr[k] *= s;
}

#endif

// ---------------------------------------------------------------------------
// dot_block kernel
//
// Computes: out[k] = v1x*v2x[k] + v1y*v2y[k] + v1z*v2z[k]  for k in [0, count)
//
// Used by AngleCalculator to batch all cos(θ) numerators for a fixed
// neighbor j against the remaining neighbors [j+1, CN) in one SIMD pass.
// Caller still applies scalar acos per element after dividing by distances.
// ---------------------------------------------------------------------------

#if defined(CORRELATION_SIMD_AVX512)

inline void dot_block(double v1x, double v1y, double v1z,
                      const double *CORRELATION_RESTRICT v2x,
                      const double *CORRELATION_RESTRICT v2y,
                      const double *CORRELATION_RESTRICT v2z,
                      double *CORRELATION_RESTRICT out,
                      std::size_t count) noexcept {
  const __m512d vv1x = _mm512_set1_pd(v1x);
  const __m512d vv1y = _mm512_set1_pd(v1y);
  const __m512d vv1z = _mm512_set1_pd(v1z);
  std::size_t k = 0;
  for (; k + 8 <= count; k += 8) {
    __m512d d = _mm512_fmadd_pd(vv1x, _mm512_loadu_pd(v2x + k),
                _mm512_fmadd_pd(vv1y, _mm512_loadu_pd(v2y + k),
                _mm512_mul_pd (vv1z, _mm512_loadu_pd(v2z + k))));
    _mm512_storeu_pd(out + k, d);
  }
  for (; k < count; ++k)
    out[k] = v1x * v2x[k] + v1y * v2y[k] + v1z * v2z[k];
}

#elif defined(CORRELATION_SIMD_AVX2)

inline void dot_block(double v1x, double v1y, double v1z,
                      const double *CORRELATION_RESTRICT v2x,
                      const double *CORRELATION_RESTRICT v2y,
                      const double *CORRELATION_RESTRICT v2z,
                      double *CORRELATION_RESTRICT out,
                      std::size_t count) noexcept {
  const __m256d vv1x = _mm256_set1_pd(v1x);
  const __m256d vv1y = _mm256_set1_pd(v1y);
  const __m256d vv1z = _mm256_set1_pd(v1z);
  std::size_t k = 0;
  for (; k + 4 <= count; k += 4) {
#if defined(__FMA__)
    __m256d d = _mm256_fmadd_pd(vv1x, _mm256_loadu_pd(v2x + k),
                _mm256_fmadd_pd(vv1y, _mm256_loadu_pd(v2y + k),
                _mm256_mul_pd (vv1z, _mm256_loadu_pd(v2z + k))));
#else
    __m256d d = _mm256_add_pd(
        _mm256_mul_pd(vv1x, _mm256_loadu_pd(v2x + k)),
        _mm256_add_pd(_mm256_mul_pd(vv1y, _mm256_loadu_pd(v2y + k)),
                      _mm256_mul_pd(vv1z, _mm256_loadu_pd(v2z + k))));
#endif
    _mm256_storeu_pd(out + k, d);
  }
  for (; k < count; ++k)
    out[k] = v1x * v2x[k] + v1y * v2y[k] + v1z * v2z[k];
}

#else // Scalar fallback

inline void dot_block(double v1x, double v1y, double v1z,
                      const double *CORRELATION_RESTRICT v2x,
                      const double *CORRELATION_RESTRICT v2y,
                      const double *CORRELATION_RESTRICT v2z,
                      double *CORRELATION_RESTRICT out,
                      std::size_t count) noexcept {
  for (std::size_t k = 0; k < count; ++k)
    out[k] = v1x * v2x[k] + v1y * v2y[k] + v1z * v2z[k];
}

#endif

// ---------------------------------------------------------------------------
// Utility: build a SoA position block from a range of atoms.
// Writes into caller-provided aligned storage. Returns the count.
// ---------------------------------------------------------------------------
template <typename AtomRange>
inline std::size_t
fill_position_block(const AtomRange &atoms, std::size_t begin, std::size_t end,
                    std::vector<double> &xs, std::vector<double> &ys,
                    std::vector<double> &zs) noexcept {
  const std::size_t count = end - begin;
  xs.resize(count);
  ys.resize(count);
  zs.resize(count);
  for (std::size_t i = 0; i < count; ++i) {
    const auto &pos = atoms[begin + i].position();
    xs[i] = pos.x();
    ys[i] = pos.y();
    zs[i] = pos.z();
  }
  return count;
}

// ---------------------------------------------------------------------------
// complex_exp_sum kernel
//
// Computes the real and imaginary parts of: rho(q) = sum_j exp(i q.r_j)
//   cos_sum = sum_j cos(qx*x_j + qy*y_j + qz*z_j)
//   sin_sum = sum_j sin(qx*x_j + qy*y_j + qz*z_j)
//
// Then |rho|^2 = cos_sum^2 + sin_sum^2.
// ---------------------------------------------------------------------------
inline void complex_exp_sum(double qx, double qy, double qz,
                            const double *CORRELATION_RESTRICT xs,
                            const double *CORRELATION_RESTRICT ys,
                            const double *CORRELATION_RESTRICT zs,
                            std::size_t count,
                            double &cos_sum, double &sin_sum) noexcept {
  cos_sum = 0.0;
  sin_sum = 0.0;
  for (std::size_t j = 0; j < count; ++j) {
    const double phase = qx * xs[j] + qy * ys[j] + qz * zs[j];
    cos_sum += std::cos(phase);
    sin_sum += std::sin(phase);
  }
}

// ---------------------------------------------------------------------------
// Report which SIMD level was compiled in (for logging/debugging)
// ---------------------------------------------------------------------------
inline const char *simd_level_string() noexcept {
#if defined(CORRELATION_SIMD_AVX512)
  return "AVX-512";
#elif defined(CORRELATION_SIMD_AVX2)
  return "AVX2";
#else
  return "Scalar";
#endif
}

} // namespace simd_utils
