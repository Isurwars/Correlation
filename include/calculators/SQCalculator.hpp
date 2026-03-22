// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "BaseCalculator.hpp"
#include "DistributionFunctions.hpp"
#include <map>
#include <string>

/**
 * @class SQCalculator
 * @brief Computes the Structure Factor S(Q) from g(r).
 */
class SQCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "SQ"; }
  std::string getGroup() const override { return "Radial"; }
  std::string getDescription() const override {
    return "Computes the Structure Factor S(Q) from the pair distribution "
           "function g(r).";
  }

  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(DistributionFunctions &df,
                      const AnalysisSettings &settings) const override;

  static Histogram
  calculate(const Histogram &g_r_hist, const Cell &cell,
            const std::map<std::string, double> &ashcroft_weights, double q_max,
            double q_bin_width, double r_integration_max);
};
