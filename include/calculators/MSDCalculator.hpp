/**
 * @file MSDCalculator.hpp
 * @brief Mean Squared Displacement (MSD) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "BaseCalculator.hpp"
#include "DistributionFunctions.hpp"
#include "Trajectory.hpp"
#include <map>
#include <string>

/**
 * @class MSDCalculator
 * @brief Computes the Mean Squared Displacement (MSD) and the diffusion
 *        coefficient via the Einstein relation.
 *
 * MSD(t) = <|r_i(t0 + t) - r_i(t0)|^2>
 *
 * Atomic positions are unwrapped using the minimum image convention so that
 * periodic boundary condition crossings are handled correctly.
 * The diffusion coefficient can be extracted from the long-time slope:
 *   D = MSD(t) / (6 * t)   (3D isotropic)
 */
class MSDCalculator : public BaseCalculator {
public:
  std::string getName() const override { return "MSD"; }
  std::string getShortName() const override { return "MSD"; }
  std::string getGroup() const override { return "Dynamic"; }
  std::string getDescription() const override {
    return "Computes the Mean Squared Displacement (MSD) and diffusion "
           "coefficient via the Einstein relation.";
  }

  bool isFrameCalculator() const override { return false; }
  bool isTrajectoryCalculator() const override { return true; }

  void calculateTrajectory(DistributionFunctions &df, const Trajectory &traj,
                           const AnalysisSettings &settings) const override;

  /**
   * @brief Core MSD calculation returning named histograms.
   *
   * Returns:
   *  - "MSD"   : raw MSD (Å²) vs time (fs)
   *  - "D_eff" : running diffusion coefficient D(t) = MSD(t) / (6t)  (Å²/fs)
   *
   * @param traj               The trajectory.
   * @param max_correlation_frames  Maximum lag in frames (-1 = half trajectory).
   * @param start_frame        First frame to use.
   * @param end_frame          One-past-last frame to use.
   */
  static std::map<std::string, Histogram>
  calculate(const Trajectory &traj, int max_correlation_frames,
            size_t start_frame, size_t end_frame);
};
