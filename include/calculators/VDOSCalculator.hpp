/**
 * @file VDOSCalculator.hpp
 * @brief Vibrational density of states (VDOS) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "BaseCalculator.hpp"
#include "analysis/DistributionFunctions.hpp"

namespace correlation::calculators {

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

  void calculateTrajectory(
      correlation::analysis::DistributionFunctions &df,
      const correlation::core::Trajectory &traj,
      const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Computes the Vibrational Density of States (VDOS).
   * @param vacf_hist Input Velocity Autocorrelation Function (VACF) histogram.
   * @return A histogram representing intensity vs frequency (THz/cm^-1).
   */
  static correlation::analysis::Histogram
  calculate(const correlation::analysis::Histogram &vacf_hist);
};

} // namespace correlation::calculators
