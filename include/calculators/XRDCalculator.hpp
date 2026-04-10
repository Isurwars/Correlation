/**
 * @file XRDCalculator.hpp
 * @brief X-ray diffraction (XRD) pattern calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"

#include <map>
#include <string>

namespace correlation::calculators {

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

  void calculateFrame(
      correlation::analysis::DistributionFunctions &df,
      const correlation::analysis::AnalysisSettings &settings) const override;

  static correlation::analysis::Histogram
  calculate(const correlation::analysis::Histogram &g_r_hist,
            const correlation::core::Cell &cell,
            const std::map<std::string, double> &ashcroft_weights,
            double lambda, double theta_min, double theta_max,
            double bin_width);
};

} // namespace correlation::calculators
