/**
 * @file CutoffSanityChecker.hpp
 * @brief Real-time sanity checker and auto-repair engine for atomic bond cutoffs.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/AnalysisTypes.hpp"
#include "core/Atom.hpp"
#include "core/Cell.hpp"
#include "math/Precision.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace correlation::analysis {

/**
 * @enum SanitySeverity
 * @brief Severity level for scientific cutoff sanity diagnostics.
 */
enum class SanitySeverity : std::uint8_t {
  Optimal = 0, ///< All cutoffs satisfy physical coordination and PBC requirements.
  Warning = 1, ///< Sub-optimal, missing, or potentially problematic cutoff bounds.
  Error = 2,   ///< Unphysical inversion or hard minimum-image convention violation.
};

/**
 * @enum SanityFixAction
 * @brief Automated one-click remediation actions for cutoff sanity warnings.
 */
enum class SanityFixAction : std::uint8_t {
  None = 0,               ///< No action required or available.
  ClampToBoxHalf = 1,     ///< Clamp maximum cutoffs to 48% of the shortest perpendicular box width.
  ResetToCovalent = 2,    ///< Reset problematic bounds to standard covalent radius sums.
  AutoDetectFirstMin = 3, ///< Trigger auto-detection of the first RDF minimum.
};

/**
 * @struct SanityReport
 * @brief Diagnostic findings produced by the cutoff sanity checker.
 */
struct SanityReport {
  SanitySeverity severity = SanitySeverity::Optimal; ///< Diagnostic severity code.
  std::string title;    ///< Concise header string (e.g. "Optimal Cutoffs").
  std::string message;  ///< Detailed explanation of findings.
  bool has_fix = false; ///< Whether an automated fix action can be performed.
  SanityFixAction fix_action = SanityFixAction::None; ///< Associated automated remediation code.
  std::string fix_label;        ///< Button label for UI (e.g. "Clamp to L/2").
  std::string problematic_pair; ///< Element pair label triggering the warning (e.g. "Si-O").
  real_t suggested_value = 0.0; ///< Numeric cutoff distance suggested for repair.
};

/**
 * @class CutoffSanityChecker
 * @brief Validates bond cutoff matrices against simulation cell periodic boundary conditions
 * and fundamental physical coordination chemistry constraints.
 */
class CutoffSanityChecker {
public:
  CutoffSanityChecker() = delete;

  /**
   * @brief Evaluates the sanity of an entire bond cutoff matrix against periodic box geometry
   * and physical constraints.
   * @param[in] matrix The bond cutoff matrix to evaluate.
   * @param[in] elements List of elements present in the system.
   * @param[in] cell Pointer to the simulation cell (optional, used for PBC checks).
   * @return A SanityReport containing severity, title, diagnostic message, and fix action.
   */
  [[nodiscard]] static auto checkCutoffs(const BondCutoffMatrix &matrix,
                                         const std::vector<core::Element> &elements,
                                         const core::Cell *cell = nullptr) -> SanityReport;

  /**
   * @brief Automatically applies a recommended repair action to the bond cutoff matrix.
   * @param[in,out] matrix The bond cutoff matrix to modify.
   * @param[in] elements List of elements present in the system.
   * @param[in] cell Pointer to the simulation cell.
   * @param[in] action The remediation action code to execute.
   * @return True if the matrix was successfully modified, false otherwise.
   */
  static auto applyFix(BondCutoffMatrix &matrix, const std::vector<core::Element> &elements,
                       const core::Cell *cell, SanityFixAction action) -> bool;
};

} // namespace correlation::analysis
