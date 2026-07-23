/**
 * @file SIMDConfig.hpp
 * @brief Configuration macros for SIMD instruction sets and compiler specifics.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/CompilerPortability.hpp" // IWYU pragma: export

// ---------------------------------------------------------------------------
// SIMD level detection
// ---------------------------------------------------------------------------
#if defined(__AVX512F__)
/** @brief Defined if AVX-512 instruction set is supported and enabled. */
#define CORRELATION_SIMD_AVX512
#elif defined(__AVX2__)
/** @brief Defined if AVX2 instruction set is supported and enabled. */
#define CORRELATION_SIMD_AVX2
#endif

#if defined(CORRELATION_SIMD_AVX512)
#include <immintrin.h>
/** @brief Number of 64-bit float elements per SIMD register. */
#define CORRELATION_SIMD_WIDTH 8
#elif defined(CORRELATION_SIMD_AVX2)
#include <immintrin.h>
/** @brief Number of 64-bit float elements per SIMD register. */
#define CORRELATION_SIMD_WIDTH 4
#else
/** @brief Number of 64-bit float elements per SIMD register. */
#define CORRELATION_SIMD_WIDTH 1
#endif
