/**
 * @file ValidationRuleService.hpp
 * @brief Stateless validation and parsing service for user inputs and simulation parameters.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include <expected>
#include <string>
#include <string_view>

namespace correlation::app {

/**
 * @class ValidationRuleService
 * @brief Stateless utility class providing functional parsing and validation rules.
 */
class ValidationRuleService {
public:
  /**
   * @brief Parses and validates a positive floating point value (> 0).
   * @param[in] str Input text to parse.
   * @return Parsed float value on success, or error message on failure.
   */
  [[nodiscard]] static std::expected<float, std::string> parsePositiveFloat(std::string_view str);

  /**
   * @brief Parses and validates a non-negative floating point value (>= 0).
   * @param[in] str Input text to parse.
   * @return Parsed float value on success, or error message on failure.
   */
  [[nodiscard]] static std::expected<float, std::string>
  parseNonNegativeFloat(std::string_view str);

  /**
   * @brief Parses and validates a positive integer value (> 0).
   * @param[in] str Input text to parse.
   * @return Parsed integer value on success, or error message on failure.
   */
  [[nodiscard]] static std::expected<int, std::string> parsePositiveInt(std::string_view str);

  /**
   * @brief Parses starting frame input ("Start", "End", or 1-based numeric).
   * @param[in] frame_s Input text string.
   * @param[in] total_frames Total number of available frames in trajectory.
   * @return 0-based frame index on success, or error message on failure.
   */
  [[nodiscard]] static std::expected<int, std::string> parseMinFrame(std::string_view frame_s,
                                                                     int total_frames);

  /**
   * @brief Parses ending frame input ("End", "Start", or 1-based numeric).
   * @param[in] frame_s Input text string.
   * @param[in] total_frames Total number of available frames in trajectory.
   * @return 0-based frame index (-1 for end) on success, or error message on failure.
   */
  [[nodiscard]] static std::expected<int, std::string> parseMaxFrame(std::string_view frame_s,
                                                                     int total_frames);

  /**
   * @brief Parses frame stride input (positive integer >= 1).
   * @param[in] stride_s Input text string.
   * @return Stride integer on success, or error message on failure.
   */
  [[nodiscard]] static std::expected<int, std::string> parseFrameStride(std::string_view stride_s);

  /**
   * @brief Validates that a histogram bin width does not exceed its maximum cutoff.
   * @param[in] bin_val Evaluated bin width.
   * @param[in] max_val Evaluated maximum cutoff.
   * @param[in] max_label Name of the maximum property (e.g. "r_max", "q_max").
   * @return Void on success, or error message on failure.
   */
  [[nodiscard]] static std::expected<void, std::string>
  validateBinWithinMax(float bin_val, float max_val, std::string_view max_label);

  /**
   * @brief Validates that an angular bin width does not exceed 180 degrees.
   * @param[in] angle_val Evaluated angle in degrees.
   * @return Void on success, or error message on failure.
   */
  [[nodiscard]] static std::expected<void, std::string> validateAngleDegrees(float angle_val);

  /**
   * @brief Validates X-ray diffraction 2-theta minimum and maximum angular bounds.
   * @param[in] theta_min Evaluated minimum 2-theta angle.
   * @param[in] theta_max Evaluated maximum 2-theta angle.
   * @return Void on success, or error message on failure.
   */
  [[nodiscard]] static std::expected<void, std::string> validateXrdTheta(float theta_min,
                                                                         float theta_max);
};

} // namespace correlation::app
