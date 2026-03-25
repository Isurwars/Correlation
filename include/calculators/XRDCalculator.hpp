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
 * @class XRDCalculator
 * @brief Computes the X-Ray Diffraction (XRD) pattern from g(r).
 */
class XRDCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "XRD"; }
  std::string getShortName() const override { return "XRD"; }
  std::string getGroup() const override { return "Radial"; }
  std::string getDescription() const override {
    return "Computes the X-Ray Diffraction Pattern (XRD).";
  }

  bool isFrameCalculator() const override { return true; }
  bool isTrajectoryCalculator() const override { return false; }

  void calculateFrame(DistributionFunctions &df,
                      const AnalysisSettings &settings) const override;

  static Histogram
  calculate(const Histogram &g_r_hist, const Cell &cell,
            const std::map<std::string, double> &ashcroft_weights,
            double lambda, double theta_min, double theta_max,
            double bin_width);
};
