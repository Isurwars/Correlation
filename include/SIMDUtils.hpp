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
  #define CORRELATION_SIMD_WIDTH 8   // 8 doubles per AVX-512 register
#elif defined(CORRELATION_SIMD_AVX2)
  #include <immintrin.h>
  #define CORRELATION_SIMD_WIDTH 4   // 4 doubles per AVX2 register
#else
  #define CORRELATION_SIMD_WIDTH 1   // scalar fallback
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
inline double dist_sq_scalar(double ax, double ay, double az,
                              double bx, double by, double bz) noexcept {
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
    __m512d dsq = _mm512_fmadd_pd(dx, dx,
                  _mm512_fmadd_pd(dy, dy,
                  _mm512_mul_pd(dz, dz)));
    _mm512_storeu_pd(out_dsq + k, dsq);
  }
  // Masked tail: remaining atoms (0-7)
  if (k < block.count) {
    // Build mask for the remaining lanes
    __mmask8 mask = static_cast<__mmask8>((1u << (block.count - k)) - 1u);
    __m512d dx = _mm512_sub_pd(_mm512_maskz_loadu_pd(mask, block.x + k), va_x);
    __m512d dy = _mm512_sub_pd(_mm512_maskz_loadu_pd(mask, block.y + k), va_y);
    __m512d dz = _mm512_sub_pd(_mm512_maskz_loadu_pd(mask, block.z + k), va_z);
    __m512d dsq = _mm512_fmadd_pd(dx, dx,
                  _mm512_fmadd_pd(dy, dy,
                  _mm512_mul_pd(dz, dz)));
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
    __m256d dsq = _mm256_fmadd_pd(dx, dx,
                  _mm256_fmadd_pd(dy, dy,
                  _mm256_mul_pd(dz, dz)));
#else
    __m256d dsq = _mm256_add_pd(_mm256_mul_pd(dx, dx),
                  _mm256_add_pd(_mm256_mul_pd(dy, dy),
                                _mm256_mul_pd(dz, dz)));
#endif
    _mm256_storeu_pd(out_dsq + k, dsq);
  }
  // Scalar tail
  for (; k < block.count; ++k) {
    out_dsq[k] = dist_sq_scalar(ax, ay, az,
                                 block.x[k], block.y[k], block.z[k]);
  }
}

#else // Scalar fallback

inline void compute_dsq_block(double ax, double ay, double az,
                               const PositionBlock &block,
                               double *CORRELATION_RESTRICT out_dsq) noexcept {
  for (std::size_t k = 0; k < block.count; ++k) {
    out_dsq[k] = dist_sq_scalar(ax, ay, az,
                                 block.x[k], block.y[k], block.z[k]);
  }
}

#endif

// ---------------------------------------------------------------------------
// Utility: build a SoA position block from a range of atoms.
// Writes into caller-provided aligned storage. Returns the count.
// ---------------------------------------------------------------------------
template <typename AtomRange>
inline std::size_t fill_position_block(const AtomRange &atoms,
                                       std::size_t begin, std::size_t end,
                                       std::vector<double> &xs,
                                       std::vector<double> &ys,
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
