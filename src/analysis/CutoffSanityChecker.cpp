/**
 * @file CutoffSanityChecker.cpp
 * @brief Implementation of CutoffSanityChecker.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "analysis/CutoffSanityChecker.hpp"
#include "physics/PhysicalData.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <optional>

namespace correlation::analysis {
namespace {

constexpr real_t MIN_PLAUSIBLE_CUTOFF = 0.8;
constexpr real_t MAX_PLAUSIBLE_CUTOFF = 6.0;
constexpr real_t PBC_SAFETY_MARGIN = 0.48;

[[nodiscard]] real_t safeCovalentRadius(std::string_view symbol) {
  try {
    return physics::getCovalentRadius(symbol);
  } catch (const std::out_of_range &) {
    return static_cast<real_t>(1.5);
  }
}

[[nodiscard]] auto getShortestBoxWidth(const core::Cell *cell) -> std::optional<real_t> {
  if (cell == nullptr || cell->volume() <= 0.0) {
    return std::nullopt;
  }
  const auto widths = cell->perpendicularWidths();
  const real_t shortest = std::min({widths[0], widths[1], widths[2]});
  if (shortest <= 0.0) {
    return std::nullopt;
  }
  return shortest;
}

[[nodiscard]] auto formatPairName(const std::vector<core::Element> &elements, size_t idx_a,
                                  size_t idx_b) -> std::string {
  if (idx_a < elements.size() && idx_b < elements.size()) {
    return std::format("{}-{}", elements[idx_a].symbol, elements[idx_b].symbol);
  }
  return std::format("Pair({},{})", idx_a, idx_b);
}

[[nodiscard]] auto checkEmptyOrIncomplete(const BondCutoffMatrix &matrix,
                                          const std::vector<core::Element> &elements)
    -> std::optional<SanityReport> {
  const size_t num_elements = elements.size();
  if (num_elements == 0) {
    return SanityReport{
        .severity = SanitySeverity::Warning,
        .title = "No Elements",
        .message = "No elements loaded. Cutoff validation is inactive.",
        .has_fix = false,
        .fix_action = SanityFixAction::None,
        .fix_label = "",
        .problematic_pair = "",
        .suggested_value = 0.0,
    };
  }
  if (matrix.size() < num_elements) {
    return SanityReport{
        .severity = SanitySeverity::Error,
        .title = "Cutoffs Not Initialized",
        .message = "Bond cutoff matrix is empty. Neighbor and angular analyses are disabled.",
        .has_fix = true,
        .fix_action = SanityFixAction::ResetToCovalent,
        .fix_label = "Initialize",
        .problematic_pair = "",
        .suggested_value = 0.0,
    };
  }
  return std::nullopt;
}

[[nodiscard]] auto checkInvertedBounds(const BondCutoffMatrix &matrix,
                                       const std::vector<core::Element> &elements)
    -> std::optional<SanityReport> {
  const size_t num_elements = elements.size();
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      const real_t min_d = std::sqrt(matrix[i][j].min_sq);
      const real_t max_d = std::sqrt(matrix[i][j].max_sq);
      if (min_d > 0.0 && max_d > 0.0 && min_d >= max_d) {
        const std::string pair_str = formatPairName(elements, i, j);
        return SanityReport{
            .severity = SanitySeverity::Error,
            .title = "Inverted Cutoff Bounds",
            .message =
                std::format("Cutoff range for {} is inverted (r_min {:.2f} Å >= r_max {:.2f} Å).",
                            pair_str, min_d, max_d),
            .has_fix = true,
            .fix_action = SanityFixAction::ResetToCovalent,
            .fix_label = "Reset",
            .problematic_pair = pair_str,
            .suggested_value = max_d * static_cast<real_t>(0.6),
        };
      }
    }
  }
  return std::nullopt;
}

[[nodiscard]] auto checkBoxHalfViolations(const BondCutoffMatrix &matrix,
                                          const std::vector<core::Element> &elements,
                                          real_t shortest_box_width)
    -> std::optional<SanityReport> {
  const real_t half_box = shortest_box_width * static_cast<real_t>(0.5);
  const size_t num_elements = elements.size();
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      const real_t max_d = std::sqrt(matrix[i][j].max_sq);
      if (max_d > half_box) {
        const std::string pair_str = formatPairName(elements, i, j);
        const real_t safe_target = shortest_box_width * PBC_SAFETY_MARGIN;
        return SanityReport{
            .severity = SanitySeverity::Error,
            .title = "PBC Violation (r_max > L/2)",
            .message = std::format(
                "Cutoff for {} ({:.2f} Å) exceeds half-box length (L/2 = {:.2f} Å). Unphysical "
                "self-interactions will occur.",
                pair_str, max_d, half_box),
            .has_fix = true,
            .fix_action = SanityFixAction::ClampToBoxHalf,
            .fix_label = "Clamp to L/2",
            .problematic_pair = pair_str,
            .suggested_value = safe_target,
        };
      }
    }
  }
  return std::nullopt;
}

[[nodiscard]] auto checkMissingCutoffs(const BondCutoffMatrix &matrix,
                                       const std::vector<core::Element> &elements)
    -> std::optional<SanityReport> {
  const size_t num_elements = elements.size();
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      const real_t max_d = std::sqrt(matrix[i][j].max_sq);
      if (max_d <= 0.05) {
        const std::string pair_str = formatPairName(elements, i, j);
        return SanityReport{
            .severity = SanitySeverity::Warning,
            .title = "Unset Cutoff Distance",
            .message =
                std::format("Cutoff for {} is unset or zero (r_max = {:.2f} Å).", pair_str, max_d),
            .has_fix = true,
            .fix_action = SanityFixAction::AutoDetectFirstMin,
            .fix_label = "Auto-Detect",
            .problematic_pair = pair_str,
            .suggested_value = 2.5,
        };
      }
    }
  }
  return std::nullopt;
}

[[nodiscard]] auto checkPlausibilityLimits(const BondCutoffMatrix &matrix,
                                           const std::vector<core::Element> &elements)
    -> std::optional<SanityReport> {
  const size_t num_elements = elements.size();
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      const real_t max_d = std::sqrt(matrix[i][j].max_sq);
      if (max_d > MAX_PLAUSIBLE_CUTOFF) {
        const std::string pair_str = formatPairName(elements, i, j);
        return SanityReport{
            .severity = SanitySeverity::Warning,
            .title = "Unusually Large Cutoff",
            .message = std::format(
                "Cutoff for {} ({:.2f} Å) exceeds typical first coordination shell bounds (> {:.1f} "
                "Å).",
                pair_str, max_d, MAX_PLAUSIBLE_CUTOFF),
            .has_fix = true,
            .fix_action = SanityFixAction::AutoDetectFirstMin,
            .fix_label = "Auto-Detect",
            .problematic_pair = pair_str,
            .suggested_value = 3.2,
        };
      }
      if (max_d < MIN_PLAUSIBLE_CUTOFF) {
        const std::string pair_str = formatPairName(elements, i, j);
        return SanityReport{
            .severity = SanitySeverity::Warning,
            .title = "Unusually Short Cutoff",
            .message = std::format(
                "Cutoff for {} ({:.2f} Å) is shorter than physical covalent bonds (< {:.1f} Å).",
                pair_str, max_d, MIN_PLAUSIBLE_CUTOFF),
            .has_fix = true,
            .fix_action = SanityFixAction::ResetToCovalent,
            .fix_label = "Reset",
            .problematic_pair = pair_str,
            .suggested_value = 1.5,
        };
      }
    }
  }
  return std::nullopt;
}

auto applyClampToBoxHalf(BondCutoffMatrix &matrix, real_t shortest_box_width) -> bool {
  const real_t clamp_limit = shortest_box_width * PBC_SAFETY_MARGIN;
  const real_t clamp_limit_sq = clamp_limit * clamp_limit;
  bool modified = false;

  for (auto &row : matrix) {
    for (auto &range : row) {
      if (range.max_sq > clamp_limit_sq) {
        range.max_sq = clamp_limit_sq;
        if (range.min_sq >= range.max_sq) {
          range.min_sq = range.max_sq * static_cast<real_t>(0.36); // 0.6^2
        }
        modified = true;
      }
    }
  }
  return modified;
}

auto applyResetToCovalent(BondCutoffMatrix &matrix, const std::vector<core::Element> &elements)
    -> bool {
  const size_t num_elements = elements.size();
  if (num_elements == 0) {
    return false;
  }
  if (matrix.size() != num_elements) {
    matrix.resize(num_elements);
  }

  for (size_t i = 0; i < num_elements; ++i) {
    const real_t r_cov_i = safeCovalentRadius(elements[i].symbol);
    for (size_t j = i; j < num_elements; ++j) {
      const real_t r_cov_j = safeCovalentRadius(elements[j].symbol);
      const real_t base_dist = r_cov_i + r_cov_j;
      const real_t min_dist = base_dist * static_cast<real_t>(0.6);
      const real_t max_dist = base_dist * static_cast<real_t>(1.2);

      const BondCutoffRange range{
          .min_sq = min_dist * min_dist,
          .max_sq = max_dist * max_dist,
      };
      matrix[i][j] = range;
      matrix[j][i] = range;
    }
  }
  return true;
}

} // namespace

auto CutoffSanityChecker::checkCutoffs(const BondCutoffMatrix &matrix,
                                       const std::vector<core::Element> &elements,
                                       const core::Cell *cell) -> SanityReport {
  if (auto empty_rep = checkEmptyOrIncomplete(matrix, elements); empty_rep.has_value()) {
    return std::move(*empty_rep);
  }

  if (auto inv_rep = checkInvertedBounds(matrix, elements); inv_rep.has_value()) {
    return std::move(*inv_rep);
  }

  if (const auto box_width = getShortestBoxWidth(cell); box_width.has_value()) {
    if (auto pbc_rep = checkBoxHalfViolations(matrix, elements, *box_width); pbc_rep.has_value()) {
      return std::move(*pbc_rep);
    }
  }

  if (auto missing_rep = checkMissingCutoffs(matrix, elements); missing_rep.has_value()) {
    return std::move(*missing_rep);
  }

  if (auto plaus_rep = checkPlausibilityLimits(matrix, elements); plaus_rep.has_value()) {
    return std::move(*plaus_rep);
  }

  return SanityReport{
      .severity = SanitySeverity::Optimal,
      .title = "Optimal Cutoffs",
      .message =
          "All bond cutoffs satisfy periodic boundary conditions and physical coordination bounds.",
      .has_fix = false,
      .fix_action = SanityFixAction::None,
      .fix_label = "",
      .problematic_pair = "",
      .suggested_value = 0.0,
  };
}

auto CutoffSanityChecker::applyFix(BondCutoffMatrix &matrix,
                                   const std::vector<core::Element> &elements,
                                   const core::Cell *cell, SanityFixAction action) -> bool {
  if (action == SanityFixAction::ClampToBoxHalf) {
    if (const auto box_width = getShortestBoxWidth(cell); box_width.has_value()) {
      return applyClampToBoxHalf(matrix, *box_width);
    }
    return false;
  }
  if (action == SanityFixAction::ResetToCovalent) {
    return applyResetToCovalent(matrix, elements);
  }
  return false;
}

} // namespace correlation::analysis
