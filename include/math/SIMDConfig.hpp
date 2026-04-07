/**
 * @file SIMDConfig.hpp
 * @brief Configuration macros for SIMD instruction sets and compiler specifics.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @license SPDX-License-Identifier: MIT
 */

#pragma once

/**
 * @file SIMDConfig.hpp
 * @brief Configuration macros for SIMD instruction sets and compiler specifics.
 */

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

// ---------------------------------------------------------------------------
// Alignment macro
// ---------------------------------------------------------------------------
#if defined(_MSC_VER)
#define CORRELATION_ALIGN(n) __declspec(align(n))
#else
#define CORRELATION_ALIGN(n) __attribute__((aligned(n)))
#endif
