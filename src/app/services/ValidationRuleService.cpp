/**
 * @file ValidationRuleService.cpp
 * @brief Implementation of stateless validation and parsing service.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/services/ValidationRuleService.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <format>
#include <string>

namespace correlation::app {

namespace {

[[nodiscard]] std::string toLowerStr(std::string_view str) {
  std::string data(str);
  std::ranges::transform(data, data.begin(),
                         [](unsigned char chr) { return static_cast<char>(std::tolower(chr)); });
  return data;
}

} // namespace

std::expected<float, std::string> ValidationRuleService::parsePositiveFloat(std::string_view str) {
  if (str.empty()) {
    return std::unexpected("Must be a positive number");
  }
  try {
    size_t idx = 0;
    const std::string text(str);
    const float val = std::stof(text, &idx);
    if (idx < text.size() || val <= 0.0F || std::isnan(val) || std::isinf(val)) {
      return std::unexpected("Must be a positive number");
    }
    return val;
  } catch (const std::exception &) {
    return std::unexpected("Must be a positive number");
  }
}

std::expected<float, std::string>
ValidationRuleService::parseNonNegativeFloat(std::string_view str) {
  if (str.empty()) {
    return std::unexpected("Must be non-negative");
  }
  try {
    size_t idx = 0;
    const std::string text(str);
    const float val = std::stof(text, &idx);
    if (idx < text.size() || val < 0.0F || std::isnan(val) || std::isinf(val)) {
      return std::unexpected("Must be non-negative");
    }
    return val;
  } catch (const std::exception &) {
    return std::unexpected("Must be non-negative");
  }
}

std::expected<int, std::string> ValidationRuleService::parsePositiveInt(std::string_view str) {
  if (str.empty()) {
    return std::unexpected("Must be a positive integer");
  }
  try {
    size_t idx = 0;
    const std::string text(str);
    const int val = std::stoi(text, &idx);
    if (idx < text.size() || val <= 0) {
      return std::unexpected("Must be a positive integer");
    }
    return val;
  } catch (const std::exception &) {
    return std::unexpected("Must be a positive integer");
  }
}

std::expected<int, std::string> ValidationRuleService::parseMinFrame(std::string_view frame_s,
                                                                     int total_frames) {
  const std::string lower = toLowerStr(frame_s);
  if (lower == "start" || frame_s.empty()) {
    return 0;
  }
  if (lower == "end") {
    return total_frames > 0 ? total_frames - 1 : 0;
  }

  const auto parsed = parsePositiveInt(frame_s);
  if (!parsed) {
    return std::unexpected("Must be positive integer, 'Start', or 'End'");
  }
  if (total_frames > 0 && *parsed > total_frames) {
    return std::unexpected(std::format("Must be ≤ total frames ({})", total_frames));
  }
  return *parsed - 1;
}

std::expected<int, std::string> ValidationRuleService::parseMaxFrame(std::string_view frame_s,
                                                                     int total_frames) {
  const std::string lower = toLowerStr(frame_s);
  if (lower == "end" || frame_s.empty()) {
    return total_frames > 0 ? total_frames - 1 : -1;
  }
  if (lower == "start") {
    return 0;
  }

  const auto parsed = parsePositiveInt(frame_s);
  if (!parsed) {
    return std::unexpected("Must be positive integer or 'End'");
  }
  if (total_frames > 0 && *parsed > total_frames) {
    return std::unexpected(std::format("Must be ≤ total frames ({})", total_frames));
  }
  return *parsed - 1;
}

std::expected<int, std::string> ValidationRuleService::parseFrameStride(std::string_view stride_s) {
  if (stride_s.empty() || stride_s == "1") {
    return 1;
  }
  const auto parsed = parsePositiveInt(stride_s);
  if (!parsed || *parsed <= 0) {
    return std::unexpected("Must be a positive integer (≥ 1)");
  }
  return *parsed;
}

std::expected<void, std::string>
ValidationRuleService::validateBinWithinMax(float bin_val, float max_val,
                                            std::string_view max_label) {
  if (max_val > 0.0F && bin_val > max_val) {
    return std::unexpected(std::format("Must be ≤ {}", max_label));
  }
  return {};
}

std::expected<void, std::string> ValidationRuleService::validateAngleDegrees(float angle_val) {
  if (angle_val > 180.0F) {
    return std::unexpected("Must be ≤ 180°");
  }
  return {};
}

std::expected<void, std::string> ValidationRuleService::validateXrdTheta(float theta_min,
                                                                         float theta_max) {
  if (theta_min >= 180.0F) {
    return std::unexpected("Must be < 180°");
  }
  if (theta_max > 180.0F) {
    return std::unexpected("Must be ≤ 180°");
  }
  if (theta_min >= theta_max) {
    return std::unexpected("Must be > Min 2θ");
  }
  return {};
}

} // namespace correlation::app
