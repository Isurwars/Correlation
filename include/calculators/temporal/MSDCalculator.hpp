/**
 * @file MSDCalculator.hpp
 * @brief Mean Squared Displacement (MSD) calculator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "calculators/contracts/ITemporalCalculator.hpp"
#include "core/Trajectory.hpp"

#include <map>
#include <string>

namespace correlation::calculators {

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
class MSDCalculator : public ITemporalCalculator {
public:
  [[nodiscard]] std::string_view getName() const override { return "MSD"; }
  [[nodiscard]] std::string_view getShortName() const override { return "MSD"; }
  [[nodiscard]] std::string_view getGroup() const override { return "Dynamic"; }
  [[nodiscard]] std::string_view getDescription() const override {
    return "Computes the Mean Squared Displacement (MSD) and diffusion "
           "coefficient via the Einstein relation.";
  }

  /**
   * @brief Dispatches the calculation for an entire trajectory.
   *
   * @param[in,out] dists Distribution functions container to append results to.
   * @param[in] traj Atomic trajectory container.
   * @param[in] settings Current analysis configuration settings.
   */
  void calculateTrajectory(correlation::analysis::DistributionFunctions &dists,
                           const correlation::core::Trajectory &traj,
                           const correlation::analysis::AnalysisSettings &settings) const override;

  /**
   * @brief Core MSD calculation returning named histograms.
   *
   * @param[in] traj               The trajectory.
   * @param[in] max_correlation_frames  Maximum lag in frames (-1 = half trajectory).
   * @param[in] start_frame        First frame to use.
   * @param[in] end_frame          One-past-last frame to use.
   * @return A map of histograms:
   *  - "MSD"   : raw MSD (Angstrom^2) vs time (fs)
   *  - "D_eff" : running diffusion coefficient D(t) = MSD(t) / (6t) (Angstrom^2/fs)
   */
  static std::map<std::string, correlation::analysis::Histogram>
  calculate(const correlation::core::Trajectory &traj,
            correlation::analysis::MaxFrames max_correlation_frames,
            correlation::analysis::StartFrame start_frame,
            correlation::analysis::EndFrame end_frame);
};

} // namespace correlation::calculators
