/**
 * @file VDOSCalculator.hpp
 * @brief Vibrational density of states (VDOS) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "BaseCalculator.hpp"
#include "DistributionFunctions.hpp"

/**
 * @class VDOSCalculator
 * @brief Computes the Vibrational Density of States (VDOS) from the VACF.
 */
class VDOSCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "vDoS"; }
  std::string getShortName() const override { return "vDoS"; }
  std::string getGroup() const override { return "Dynamic"; }
  std::string getDescription() const override {
    return "Computes the Vibrational Density of States (vDoS).";
  }

  bool isFrameCalculator() const override { return false; }
  bool isTrajectoryCalculator() const override { return true; }

  void calculateTrajectory(DistributionFunctions &df, const Trajectory &traj,
                           const AnalysisSettings &settings) const override;

  static Histogram calculate(const Histogram &vacf_hist);
};
