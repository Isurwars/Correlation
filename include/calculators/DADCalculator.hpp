/**
 * @file DADCalculator.hpp
 * @brief Dihedral angle distribution (DAD) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
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
  [[nodiscard]] std::string getName() const override { return "DAD"; }
  [[nodiscard]] std::string getShortName() const override { return "DAD"; }
  [[nodiscard]] std::string getGroup() const override { return "Angular"; }
  [[nodiscard]] std::string getDescription() const override {
    return "Computes the Dihedral Angle Distribution (DAD).";
  }

  [[nodiscard]] bool isFrameCalculator() const override { return true; }
  [[nodiscard]] bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(correlation::analysis::DistributionFunctions &dists,
                      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief High-performance computation of the Dihedral Angle Distribution (DAD).
   *
   * @param cell The periodic cell.
   * @param neighbors Structural analyzer containing the neighbor graph.
   * @param bin_width Angular resolution (radians).
   * @return A histogram representing the DAD distribution.
   */
  static correlation::analysis::Histogram calculate(const correlation::core::Cell &cell,
                                                    const correlation::analysis::StructureAnalyzer *neighbors,
                                                    double bin_width);
};

} // namespace correlation::calculators
