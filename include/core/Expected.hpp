/**
 * @file Expected.hpp
 * @brief Polyfill and type alias for std::expected and std::unexpected.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#if defined(__cpp_lib_expected) && __cpp_lib_expected >= 202211L
#include <expected>
namespace correlation {
/**
 * @brief Represents an expected value or an unexpected error.
 * @tparam T The value type.
 * @tparam E The error type.
 */
template <typename T, typename E> using expected = std::expected<T, E>;

/**
 * @brief Wraps an unexpected error value.
 * @tparam E The error type.
 */
template <typename E> using unexpected = std::unexpected<E>;

using unexpect_t = std::unexpect_t;
inline constexpr std::unexpect_t unexpect{};
} // namespace correlation
#else
#include "core/detail/tl_expected.hpp"
namespace correlation {
/**
 * @brief Represents an expected value or an unexpected error (portable fallback).
 * @tparam T The value type.
 * @tparam E The error type.
 */
template <typename T, typename E> using expected = tl::expected<T, E>;

/**
 * @brief Wraps an unexpected error value (portable fallback).
 * @tparam E The error type.
 */
template <typename E> using unexpected = tl::unexpected<E>;

using unexpect_t = tl::unexpect_t;
inline constexpr tl::unexpect_t unexpect{};
} // namespace correlation
#endif
