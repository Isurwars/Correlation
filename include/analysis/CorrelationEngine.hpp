/**
 * @file CorrelationEngine.hpp
 * @brief Analysis pipeline orchestration engine for molecular simulations.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "core/Trajectory.hpp"

#include <expected>
#include <functional>
#include <memory>
#include <string>

namespace correlation::analysis {

/**
 * @brief Configuration parameters for running the correlation analysis engine.
 */
struct CorrelationEngineConfig {
  AnalysisSettings settings;
  BondCutoffMatrix bond_cutoffs;
  int min_frame{0};
  int max_frame{-1};
  real_t time_step{1.0};
};

/**
 * @class CorrelationEngine
 * @brief Standalone analysis pipeline orchestration for molecular simulations.
 *
 * Encapsulates the execution pipeline (frame-by-frame structural analysis,
 * trajectory-wide dynamics, and post-processing) decoupled from UI and storage concerns.
 */
class CorrelationEngine {
public:
  /**
   * @brief Validates analysis configuration parameters.
   * @param config The configuration to validate.
   * @return Error message, or empty string if configuration is valid.
   */
  [[nodiscard]] static std::string validateConfig(const CorrelationEngineConfig &config);

  /**
   * @brief Executes the complete analysis pipeline on a trajectory.
   * @param trajectory The trajectory to analyze.
   * @param config Configuration parameters for the run.
   * @param progress_callback Optional progress callback function.
   * @return Unique pointer to computed DistributionFunctions, or error message string.
   */
  [[nodiscard]] static std::expected<std::unique_ptr<DistributionFunctions>, std::string>
  runAnalysis(correlation::core::Trajectory &trajectory, const CorrelationEngineConfig &config,
              std::function<void(float, const std::string &)> progress_callback = nullptr);

  /**
   * @brief Computes trajectory-wide dynamic properties (MSD, VACF, relaxation time, Deborah number).
   * @param distribution_functions Target DistributionFunctions container to update.
   */
  static void calculateDynamicProperties(DistributionFunctions &distribution_functions);

  /**
   * @brief Dispatches active trajectory calculators on the trajectory.
   * @param distribution_functions Target DistributionFunctions container.
   * @param trajectory The trajectory to calculate on.
   * @param settings Active analysis settings.
   */
  static void runTrajectoryCalculators(DistributionFunctions &distribution_functions,
                                       correlation::core::Trajectory &trajectory, const AnalysisSettings &settings);
};

} // namespace correlation::analysis
