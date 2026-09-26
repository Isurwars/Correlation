/**
 * @file HBondCalculator.hpp
 * @brief Hydrogen Bond Analysis calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */
#pragma once

#include "BaseCalculator.hpp"

namespace correlation::calculators {

/**
 * @class HBondCalculator
 * @brief Computes hydrogen bond statistics based on geometric criteria.
 */
class HBondCalculator : public BaseCalculator {
public:
  [[nodiscard]] std::string_view getName() const override { return "Hydrogen Bond"; }
  [[nodiscard]] std::string_view getShortName() const override { return "HBond"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Structural"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "Computes hydrogen bond statistics (counts and distribution).";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Performs Hydrogen Bond analysis.
   * @param[in] cell The simulation cell.
   * @param[in] neighbors The structural analyzer.
   * @return A histogram of H-bond counts.
   */
  [[nodiscard]] static correlation::analysis::Histogram
  calculate(const correlation::core::Cell &cell,
            const correlation::analysis::StructureAnalyzer *neighbors);
};

} // namespace correlation::calculators
