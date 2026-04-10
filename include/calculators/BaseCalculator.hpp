/**
 * @file BaseCalculator.hpp
 * @brief Abstract base class for all analysis calculators.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "core/Trajectory.hpp"

#include <string>

namespace correlation::calculators {

/**
 * @brief Base class for all analysis calculators.
 */
class BaseCalculator {
public:
  virtual ~BaseCalculator() = default;

  /**
   * @brief Returns the unique identifier/name of the calculator.
   * @return The full name of the calculator (e.g. "Radial Distribution Function").
   */
  virtual std::string getName() const = 0;

  /**
   * @brief Returns a short, UI-friendly name of the calculator (e.g. "g_r",
   * "S_q").
   * @return The short name/abbreviation of the calculator.
   */
  virtual std::string getShortName() const = 0;

  /**
   * @brief Returns the UI group this calculator belongs to (e.g., "Radial",
   * "Angular", "Dynamic", "Rings").
   * @return The group name for UI categorization.
   */
  virtual std::string getGroup() const = 0;

  /**
   * @brief Returns a brief description of the calculator's purpose.
   * @return A human-readable description string.
   */
  virtual std::string getDescription() const = 0;

  /**
   * @brief Check if this calculator runs per-frame (e.g. RDF, PAD).
   * @return True if it calculates properties for individual snapshots.
   */
  virtual bool isFrameCalculator() const = 0;

  /**
   * @brief Check if this calculator runs on the whole trajectory (e.g. VACF).
   * @return True if it calculates time-dependent or multi-frame properties.
   */
  virtual bool isTrajectoryCalculator() const = 0;

  /**
   * @brief Calculate per-frame properties.
   *        Called concurrently on different DistributionFunctions objects.
   *        Must be thread-safe (const).
   * @param df The DistributionFunctions object to store the results.
   * @param settings The analysis settings and parameters.
   */
  virtual void calculateFrame(
      correlation::analysis::DistributionFunctions &df,
      const correlation::analysis::AnalysisSettings &settings) const {}

  /**
   * @brief Calculate multi-frame properties or post-processing.
   *        Called once after all frames are accumulated.
   * @param df The DistributionFunctions container to store trajectory-wide results.
   * @param traj The complete trajectory data.
   * @param settings The analysis settings and parameters.
   */
  virtual void calculateTrajectory(
      correlation::analysis::DistributionFunctions &df,
      const correlation::core::Trajectory &traj,
      const correlation::analysis::AnalysisSettings &settings) const {}
};

} // namespace correlation::calculators
