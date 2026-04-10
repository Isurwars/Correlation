/**
 * @file DADCalculator.hpp
 * @brief Dihedral angle distribution (DAD) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"

namespace correlation::calculators {

/**
 * @class DADCalculator
 * @brief Computes the Dihedral Angle Distribution (DAD).
 */
class DADCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "DAD"; }
  std::string getShortName() const override { return "DAD"; }
  std::string getGroup() const override { return "Angular"; }
  std::string getDescription() const override {
    return "Computes the Dihedral Angle Distribution (DAD).";
  }

  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &df,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  static correlation::analysis::Histogram calculate(const correlation::core::Cell &cell,
                             const correlation::analysis::StructureAnalyzer *neighbors,
                             double bin_width);
};

} // namespace correlation::calculators
