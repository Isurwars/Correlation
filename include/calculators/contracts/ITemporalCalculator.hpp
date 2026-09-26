/**
 * @file ITemporalCalculator.hpp
 * @brief Contract for multi-frame time-series trajectory calculators.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "calculators/BaseCalculator.hpp"

namespace correlation::calculators {

/**
 * @brief Abstract contract for calculators operating across multi-frame trajectories over time.
 * Enforces isFrameCalculator() = false, isTrajectoryCalculator() = true, and provides pure virtual calculateTrajectory.
 */
class ITemporalCalculator : public BaseCalculator {
public:
  ~ITemporalCalculator() override = default;

  [[nodiscard]] bool isFrameCalculator() const noexcept override { return false; }
  [[nodiscard]] bool isTrajectoryCalculator() const noexcept override { return true; }

  /**
   * @brief Calculate multi-frame properties across the full trajectory.
   * @param[in,out] dists The DistributionFunctions container to store trajectory-wide results.
   * @param[in] traj The complete trajectory data.
   * @param[in] settings The analysis settings and parameters.
   */
  void calculateTrajectory(correlation::analysis::DistributionFunctions &dists,
                           const correlation::core::Trajectory &traj,
                           const correlation::analysis::AnalysisSettings &settings) const override = 0;
};

} // namespace correlation::calculators
