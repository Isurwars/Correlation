/**
 * @file PeakIntegrationService.hpp
 * @brief Numerical integration, peak feature extraction, and coordination number (CN) computation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/Precision.hpp"

#include <optional>
#include <span>
#include <string>

namespace correlation::analysis {

/**
 * @struct PeakIntegralResult
 * @brief Quantitative metrics resulting from peak detection and integration.
 */
struct PeakIntegralResult {
  real_t r_peak{0.0};                        ///< Coordinate of peak maximum (Å or unit).
  real_t r_min{0.0};                         ///< Lower integration boundary (Å or unit).
  real_t r_max{0.0};                         ///< Upper integration boundary (Å or unit).
  real_t peak_height{0.0};                   ///< Amplitude/intensity at peak maximum.
  real_t area{0.0};                          ///< Standard trapezoidal integral: ∫ y(x) dx.
  std::optional<real_t> coordination_number; ///< 4πρ ∫ r² g(r) dr for radial distributions.
  real_t fwhm{0.0};                          ///< Full width at half maximum.
  real_t centroid{0.0};                      ///< Weighted center of mass: ∫ x y(x) dx / ∫ y(x) dx.
};

/**
 * @class PeakIntegrationService
 * @brief High-performance numerical integration and automated peak feature detection.
 */
class PeakIntegrationService {
public:
  /**
   * @brief Numerically integrates a distribution across a specified domain range.
   *
   * @param[in] x_bins Domain bin centers.
   * @param[in] y_vals Function values.
   * @param[in] r_min Lower bound of integration range.
   * @param[in] r_max Upper bound of integration range.
   * @param[in] number_density Optional macroscopic atomic number density ρ₀ (atoms/Å³).
   * @param[in] is_radial True if distribution is radial g(r) or H(r).
   * @return Computed integral result metrics.
   */
  [[nodiscard]] static PeakIntegralResult
  integrateRange(std::span<const real_t> x_bins, std::span<const real_t> y_vals, real_t r_min,
                 real_t r_max, std::optional<real_t> number_density = std::nullopt,
                 bool is_radial = false);

  /**
   * @brief Detects the nearest peak maximum to a clicked coordinate and integrates its shell.
   *
   * @param[in] x_bins Domain bin centers.
   * @param[in] y_vals Function values.
   * @param[in] target_x Selected coordinate to search around.
   * @param[in] number_density Optional macroscopic atomic number density ρ₀ (atoms/Å³).
   * @param[in] is_radial True if distribution is radial g(r) or H(r).
   * @return Detected peak result, or std::nullopt if no valid peak exists.
   */
  [[nodiscard]] static std::optional<PeakIntegralResult>
  autoDetectPeakNear(std::span<const real_t> x_bins, std::span<const real_t> y_vals,
                     real_t target_x, std::optional<real_t> number_density = std::nullopt,
                     bool is_radial = false);

  /**
   * @brief Calculates Full Width at Half Maximum (FWHM) around a identified peak index.
   *
   * @param[in] x_bins Domain coordinates.
   * @param[in] y_vals Amplitudes.
   * @param[in] peak_idx Index of peak maximum.
   * @param[in] min_idx Lower search boundary index.
   * @param[in] max_idx Upper search boundary index.
   * @return Computed FWHM in coordinate units.
   */
  [[nodiscard]] static real_t calculateFWHM(std::span<const real_t> x_bins,
                                            std::span<const real_t> y_vals, std::size_t peak_idx,
                                            std::size_t min_idx, std::size_t max_idx) noexcept;
};

} // namespace correlation::analysis
