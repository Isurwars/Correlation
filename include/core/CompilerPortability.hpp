/**
 * @file CompilerPortability.hpp
 * @brief Unified compiler extensions, alignment, restrict, and GPU execution qualifiers.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

// ---------------------------------------------------------------------------
// GPU / Compiler Execution Qualifiers
// ---------------------------------------------------------------------------
#ifndef CORRELATION_DEVICE
#if defined(__CUDACC__) || defined(__HIPCC__)
#define CORRELATION_DEVICE __device__
#define CORRELATION_HOST __host__
#define CORRELATION_GLOBAL __global__
#define CORRELATION_FORCEINLINE __forceinline__
#else
#define CORRELATION_DEVICE
#define CORRELATION_HOST
#define CORRELATION_GLOBAL
#define CORRELATION_FORCEINLINE inline
#endif
#endif

// ---------------------------------------------------------------------------
// Portable restrict keyword macro
// ---------------------------------------------------------------------------
#ifndef CORRELATION_RESTRICT
#if defined(__CUDACC__) || defined(__HIPCC__)
#define CORRELATION_RESTRICT __restrict__
#elif defined(_MSC_VER)
#define CORRELATION_RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
#define CORRELATION_RESTRICT __restrict__
#else
#define CORRELATION_RESTRICT
#endif
#endif

// ---------------------------------------------------------------------------
// Memory Alignment Macro
// ---------------------------------------------------------------------------
#ifndef CORRELATION_ALIGN
#if defined(_MSC_VER)
#define CORRELATION_ALIGN(n) __declspec(align(n))
#else
#define CORRELATION_ALIGN(n) __attribute__((aligned(n)))
#endif
#endif
