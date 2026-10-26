/**
 * @file ITrajectoryCalculator.hpp
 * @brief Abstract interface for multi-frame time-series trajectory calculators.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "calculators/contracts/ICalculator.hpp"

namespace correlation::core {
class Trajectory;
} // namespace correlation::core

namespace correlation::analysis {
class DistributionFunctions;
struct AnalysisSettings;
} // namespace correlation::analysis

namespace correlation::calculators {

/**
 * @brief Segregated contract for calculators operating on multi-frame time-series trajectories.
 */
class ITrajectoryCalculator : public virtual ICalculator {
public:
  ~ITrajectoryCalculator() override = default;

  /**
   * @brief Calculates trajectory-wide dynamic or temporal properties.
   * @param[in,out] dists DistributionFunctions container to accumulate trajectory-wide results.
   * @param[in] traj Source trajectory data spanning multiple time steps.
   * @param[in] settings Configuration parameters for analysis execution.
   */
  virtual void calculateTrajectory(correlation::analysis::DistributionFunctions &dists,
                                   const correlation::core::Trajectory &traj,
                                   const correlation::analysis::AnalysisSettings &settings) const = 0;
};

} // namespace correlation::calculators
