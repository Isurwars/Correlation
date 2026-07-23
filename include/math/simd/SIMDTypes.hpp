/**
 * @file SIMDTypes.hpp
 * @brief Common parameter containers, POD types, and scalar helpers for SIMD operations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/Precision.hpp"
#include "math/SIMDConfig.hpp"

#include <cstddef>
#include <vector>

namespace correlation::math {

// ---------------------------------------------------------------------------
// SoA block of atom positions
// ---------------------------------------------------------------------------
/**
 * @struct PositionBlockT
 * @brief Represents a Structure of Arrays (SoA) block of atom positions for SIMD processing.
 * @tparam T Floating point coordinate type (float or double).
 */
template <typename T = real_t> struct PositionBlockT {
  T *x{nullptr};        ///< Array of x-coordinates.
  T *y{nullptr};        ///< Array of y-coordinates.
  T *z{nullptr};        ///< Array of z-coordinates.
  std::size_t count{0}; ///< Number of atoms in this block.
};

/// Backward-compatible type alias for default real_t precision PositionBlock.
using PositionBlock = PositionBlockT<real_t>;

// ---------------------------------------------------------------------------
// Scalar helper
// ---------------------------------------------------------------------------
/**
 * @struct Point3T
 * @brief Lightweight POD for passing 3D coordinates to scalar SIMD helpers.
 * @tparam T Floating point coordinate type (float or double).
 */
template <typename T = real_t> struct Point3T {
  T x{static_cast<T>(0.0)}; ///< X-coordinate.
  T y{static_cast<T>(0.0)}; ///< Y-coordinate.
  T z{static_cast<T>(0.0)}; ///< Z-coordinate.
};

/// Backward-compatible type alias for default real_t precision Point3.
using Point3 = Point3T<real_t>;

/**
 * @brief Parameter container for sinc_integral kernel invocations.
 * @tparam T Floating-point precision (float or double).
 */
template <typename T> struct SincIntegralParams {
  T q_magnitude{static_cast<T>(0.0)};
  const T *CORRELATION_RESTRICT integrand{nullptr};
  const T *CORRELATION_RESTRICT radial_bins{nullptr};
  T *CORRELATION_RESTRICT sinqr_scratch{nullptr};
  std::size_t count{0};
};

/**
 * @brief Parameter container for normalize_rdf_bins kernel invocations.
 * @tparam T Floating-point precision (float or double).
 */
template <typename T> struct RDFNormalizationParams {
  const T *CORRELATION_RESTRICT hist_data{nullptr};
  const T *CORRELATION_RESTRICT radial_bins{nullptr};
  T g_norm{static_cast<T>(0.0)};
  T inv_Ni_dr{static_cast<T>(0.0)};
  T inv_Nj_dr{static_cast<T>(0.0)};
  T pi4_rho_j{static_cast<T>(0.0)};
  T *CORRELATION_RESTRICT g_out{nullptr};
  T *CORRELATION_RESTRICT G_out{nullptr};
  T *CORRELATION_RESTRICT J_out{nullptr};
  T *CORRELATION_RESTRICT Jinv_out{nullptr};
  std::size_t count{0};
};

/**
 * @brief Parameter container for scale_bins kernel invocations.
 * @tparam T Floating-point precision (float or double).
 */
template <typename T> struct ScaleBinsParams {
  T *arr{nullptr};
  T scale_factor{static_cast<T>(1.0)};
  std::size_t count{0};
};

/**
 * @brief Parameter container for miller_phase_sum kernel invocations.
 * @tparam T Floating-point precision (float or double).
 */
template <typename T> struct MillerPhaseSumParams {
  const T *CORRELATION_RESTRICT cos1{nullptr};
  const T *CORRELATION_RESTRICT sin1{nullptr};
  const T *CORRELATION_RESTRICT cos2{nullptr};
  const T *CORRELATION_RESTRICT sin2{nullptr};
  const T *CORRELATION_RESTRICT cos3{nullptr};
  const T *CORRELATION_RESTRICT sin3{nullptr};
  std::size_t count{0};
};

/**
 * @brief Parameter container for miller_exp_sum / complex_exp_sum kernel invocations.
 * @tparam T Coordinate scalar type (defaults to real_t).
 */
template <typename T = real_t> struct MillerExpSumParams {
  T q_x{static_cast<T>(0.0)};
  T q_y{static_cast<T>(0.0)};
  T q_z{static_cast<T>(0.0)};
  const T *CORRELATION_RESTRICT x_s{nullptr};
  const T *CORRELATION_RESTRICT y_s{nullptr};
  const T *CORRELATION_RESTRICT z_s{nullptr};
  std::size_t count{0};
};

/// Alias for MillerExpSumParams used with complex_exp_sum kernel invocations.
template <typename T = real_t> using ComplexExpSumParams = MillerExpSumParams<T>;

/**
 * @brief Output result container holding cosine and sine phase sums.
 * @tparam T Floating-point precision (defaults to real_t).
 */
template <typename T = real_t> struct ComplexExpSumResult {
  T cos_sum{static_cast<T>(0.0)};
  T sin_sum{static_cast<T>(0.0)};
};

/// Alias for PhaseSumResult used with phase sum kernel invocations.
template <typename T = real_t> using PhaseSumResult = ComplexExpSumResult<T>;

/// Alias for MillerPhaseSumResult used with miller_phase_sum kernel invocations.
template <typename T = real_t> using MillerPhaseSumResult = ComplexExpSumResult<T>;

/**
 * @brief Parameter container for fill_position_block kernel invocations.
 * @tparam AtomRange A range of atom objects.
 * @tparam T Coordinate scalar type (defaults to real_t).
 */
template <typename AtomRange, typename T = real_t> struct FillPositionBlockParams {
  const AtomRange *atoms{nullptr};
  std::size_t begin_idx{0};
  std::size_t end_idx{0};
  std::vector<T> *x_s{nullptr};
  std::vector<T> *y_s{nullptr};
  std::vector<T> *z_s{nullptr};
};

/**
 * @brief Computes the squared distance between two 3D points.
 * @tparam T Floating point coordinate type.
 * @param[in] point_a The first point.
 * @param[in] point_b The second point.
 * @return The scalar squared distance.
 */
template <typename T> [[nodiscard]] inline T dist_sq_scalar(Point3T<T> point_a, Point3T<T> point_b) noexcept {
  const T dist_x = point_b.x - point_a.x;
  const T dist_y = point_b.y - point_a.y;
  const T dist_z = point_b.z - point_a.z;
  return dist_x * dist_x + dist_y * dist_y + dist_z * dist_z;
}

/**
 * @brief Non-templated overload of dist_sq_scalar for backward compatibility.
 * @param[in] point_a The first point.
 * @param[in] point_b The second point.
 * @return The scalar squared distance using real_t.
 */
[[nodiscard]] inline real_t dist_sq_scalar(Point3 point_a, Point3 point_b) noexcept {
  return dist_sq_scalar<real_t>(point_a, point_b);
}

} // namespace correlation::math
