/**
 * @file VACFCalculator.hpp
 * @brief Velocity autocorrelation function (VACF) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "BaseCalculator.hpp"
#include "DistributionFunctions.hpp"
#include "core/Trajectory.hpp"

#include <map>
#include <string>

/**
 * @class VACFCalculator
 * @brief Computes the Velocity Autocorrelation Function (VACF).
 */
class VACFCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "VACF"; }
  std::string getShortName() const override { return "VACF"; }
  std::string getGroup() const override { return "Dynamic"; }
  std::string getDescription() const override {
    return "Computes the Velocity Autocorrelation Function (VACF).";
  }

  bool isFrameCalculator() const override { return false; }
  bool isTrajectoryCalculator() const override { return true; }

  void calculateTrajectory(DistributionFunctions &df,
                           const correlation::core::Trajectory &traj,
                           const AnalysisSettings &settings) const override;

  static std::map<std::string, Histogram>
  calculate(const correlation::core::Trajectory &traj,
            int max_correlation_frames, size_t start_frame, size_t end_frame);
};
