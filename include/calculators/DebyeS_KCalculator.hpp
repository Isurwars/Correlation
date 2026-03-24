// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "BaseCalculator.hpp"
#include "DistributionFunctions.hpp"

class DebyeS_KCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "Debye_SK"; }
  std::string getGroup() const override { return "Radial"; }
  std::string getDescription() const override {
    return "Computes the Structure Factor S(k) using the Debye direct summation formula (ideal for non-periodic/cluster systems).";
  }

  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(DistributionFunctions &df,
                      const AnalysisSettings &settings) const override;
};
