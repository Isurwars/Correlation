/**
 * @file CNCalculator.hpp
 * @brief Coordination number (CN) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"

namespace correlation::calculators {

/**
 * @class CNCalculator
 * @brief Computes the Coordination Number (CN) distribution.
 */
class CNCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "CN"; }
  std::string getShortName() const override { return "CN"; }
  std::string getGroup() const override { return "Structural"; }
  std::string getDescription() const override {
    return "Computes the Coordination Number (CN).";
  }

  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &df,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  static correlation::analysis::Histogram calculate(const correlation::core::Cell &cell,
                             const correlation::analysis::StructureAnalyzer *neighbors);
};

} // namespace correlation::calculators
