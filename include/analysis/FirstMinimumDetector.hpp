/**
 * @file FirstMinimumDetector.hpp
 * @brief Algorithmic detection of first coordination shell peak, onset, and minimum in distribution
 * functions.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "math/Precision.hpp"

#include <expected>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace correlation::analysis {

/**
 * @enum class DetectionError
 * @brief Error categories encountered during first minimum or peak detection.
 */
enum class DetectionError {
  InsufficientData, ///< Bins or intensities array has fewer than required points.
  MismatchedSizes,  ///< Bins and intensities spans have different lengths.
  NoPeakFound,      ///< No primary coordination peak detected above noise threshold.
  NoMinimumFound,   ///< No inflection point or trough detected following primary peak.
  InvalidRange,     ///< Minimum distance exceeds maximum distance or window is unphysical.
  PartialNotFound   ///< Requested partial key was not found in Histogram container.
};

/**
 * @struct CutoffResult
 * @brief Calculated radial cutoff boundaries and key feature locations.
 */
struct CutoffResult {
  real_t r_min{0.0};       ///< Lower cutoff bound (onset of first coordination shell) in Å.
  real_t r_max{0.0};       ///< Upper cutoff bound (first minimum following primary peak) in Å.
  real_t r_peak{0.0};      ///< Radial position of primary peak maximum in Å.
  real_t peak_height{0.0}; ///< Intensity at primary peak maximum.
  real_t min_height{0.0};  ///< Intensity at detected first minimum.
};

/**
 * @struct DetectionOptions
 * @brief Tuning parameters for first minimum and cutoff detection.
 */
struct DetectionOptions {
  real_t min_search_radius{0.5}; ///< Minimum radius (Å) to start searching for the primary peak.
  real_t max_search_radius{8.0}; ///< Maximum radius (Å) to search for the primary peak.
  real_t search_window{1.8}; ///< Radial search window (Å) after peak to locate the first minimum.
  real_t onset_threshold_frac{0.05}; ///< Fraction of peak height defining peak onset (r_min).
  real_t core_floor_radius{0.5};     ///< Hard-core lower bound floor for r_min in Å.
  bool apply_smoothing{true}; ///< Whether to apply 3-point filter before derivative analysis.
};

/**
 * @class FirstMinimumDetector
 * @brief Evaluates radial distribution functions to locate coordination shell boundaries.
 */
class FirstMinimumDetector {
public:
  /**
   * @brief Detects first coordination shell cutoffs from raw or smoothed spans.
   * @param[in] r_bins Radial grid coordinates in Ångströms.
   * @param[in] g_r Radial distribution intensities.
   * @param[in] opts Detection parameters and tuning options.
   * @return Detected CutoffResult or DetectionError.
   */
  [[nodiscard]] static std::expected<CutoffResult, DetectionError>
  detectCutoffs(std::span<const real_t> r_bins, std::span<const real_t> g_r,
                const DetectionOptions &opts = {});

  /**
   * @brief Convenience overload extracting partial curve from Histogram container.
   * @param[in] hist Distribution function histogram container.
   * @param[in] partial_key Partial key identifier (e.g. "Si-O", "Total").
   * @param[in] opts Detection parameters and tuning options.
   * @return Detected CutoffResult or DetectionError.
   */
  [[nodiscard]] static std::expected<CutoffResult, DetectionError>
  detectCutoffs(const Histogram &hist, std::string_view partial_key,
                const DetectionOptions &opts = {});

  /**
   * @brief Converts a DetectionError to a human-readable diagnostic message.
   * @param[in] error Target error code.
   * @return Descriptive error string.
   */
  [[nodiscard]] static std::string_view formatError(DetectionError error) noexcept;
};

} // namespace correlation::analysis
