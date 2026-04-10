/**
 * @file PADCalculator.hpp
 * @brief Plane angle distribution (PAD) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"

namespace correlation::calculators {

/**
 * @class PADCalculator
 * @brief Computes the Bond Angle Distribution (BAD/PAD).
 */
class PADCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "PAD"; }
  std::string getShortName() const override { return "PAD"; }
  std::string getGroup() const override { return "Angular"; }
  std::string getDescription() const override {
    return "Computes the Plane-Angle Distribution (PAD).";
  }

  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(
      correlation::analysis::DistributionFunctions &df,
      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief High-performance computation of the Plane Angle Distribution (PAD).
   * 
   * @param cell The periodic cell.
   * @param neighbors Structural analyzer containing the neighbor graph.
   * @param bin_width Angular resolution (radians).
   * @return A histogram representing the PAD distribution.
   */
  static correlation::analysis::Histogram
  calculate(const correlation::core::Cell &cell,
            const correlation::analysis::StructureAnalyzer *neighbors,
            double bin_width);
};

} // namespace correlation::calculators
