// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "BaseCalculator.hpp"
#include "DistributionFunctions.hpp"

/**
 * @class DADCalculator
 * @brief Computes the Dihedral Angle Distribution (DAD).
 */
class DADCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "DAD"; }
  std::string getGroup() const override { return "Angular"; }
  std::string getDescription() const override {
    return "Computes the Dihedral Angle Distribution.";
  }

  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(DistributionFunctions &df,
                      const AnalysisSettings &settings) const override;

  static Histogram calculate(const Cell &cell,
                             const StructureAnalyzer *neighbors,
                             double bin_width);
};
