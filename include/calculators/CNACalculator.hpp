/**
 * @file CNACalculator.hpp
 * @brief Common Neighbor Analysis (CNA) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */
#pragma once

#include "BaseCalculator.hpp"

namespace correlation::calculators {

/**
 * @class CNACalculator
 * @brief Performs Common Neighbor Analysis (CNA) to classify local structures.
 */
class CNACalculator : public BaseCalculator {
public:
  [[nodiscard]] std::string getName() const override { return "CNA"; }
  [[nodiscard]] std::string getShortName() const override { return "CNA"; }
  [[nodiscard]] std::string getGroup() const override { return "Structural"; }
  [[nodiscard]] std::string getDescription() const override {
    return "Common Neighbor Analysis for local environment classification.";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Performs CNA on a single cell.
   * @param cell The simulation cell.
   * @param neighbors The structural analyzer with neighbor info.
   * @return A histogram of CNA indices.
   */
  static correlation::analysis::Histogram calculate(const correlation::core::Cell &cell,
                                                    const correlation::analysis::StructureAnalyzer *neighbors);
};

} // namespace correlation::calculators
