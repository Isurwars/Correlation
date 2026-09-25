/**
 * @file BondCutoffController.hpp
 * @brief Controller managing bond cutoff matrix operations and UI synchronization.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "AppWindow.h"
#include "analysis/DistributionFunctions.hpp"
#include "app/BondCutoffService.hpp"
#include "app/TrajectoryLoader.hpp"
#include "app/core/AppOptions.hpp"

namespace correlation::app {

class AppBackend;

/**
 * @enum FactorBound
 * @brief Specifies which cutoff bound (Min or Max) to adjust by a covalent factor.
 */
enum class FactorBound { Min, Max };

/**
 * @class BondCutoffController
 * @brief Controller that manages element-pair bond cutoff matrices and synchronizes with the UI.
 */
class BondCutoffController {
public:
  /**
   * @brief Constructs a BondCutoffController with injected services.
   * @param[in,out] window Target UI window.
   * @param[in,out] loader Trajectory loader providing cell and trajectory.
   * @param[in,out] options Program options holding bond cutoffs.
   */
  BondCutoffController(AppWindow &window, TrajectoryLoader &loader, ProgramOptions &options);

  /**
   * @brief Constructs a BondCutoffController with injected services (transitional).
   * @param[in,out] window Target UI window.
   * @param[in,out] loader Trajectory loader providing cell and trajectory.
   * @param[in,out] cutoff_service Cutoff calculation service.
   * @param[in,out] options Program options holding bond cutoffs.
   */
  BondCutoffController(AppWindow &window, TrajectoryLoader &loader,
                       BondCutoffService &cutoff_service, ProgramOptions &options);

  /**
   * @brief Constructs a BondCutoffController with the given UI window and backend (transitional).
   * @param[in,out] window Target UI window.
   * @param[in,out] backend Application backend holding atomic cell and cutoffs.
   */
  BondCutoffController(AppWindow &window, AppBackend &backend);

  /**
   * @brief Populates the UI model with default recommended bond cutoffs based on atomic elements.
   */
  void setBondCutoffs();

  /**
   * @brief Parses the current UI model into a BondCutoffMatrix.
   * @return Active bond cutoff matrix.
   */
  [[nodiscard]] correlation::analysis::BondCutoffMatrix getBondCutoffs() const;

  /**
   * @brief Scales all current bond cutoffs by a multiplicative factor.
   * @param[in] scale_factor Scaling factor (> 0).
   */
  void applyScaledCutoffs(float scale_factor);

  /**
   * @brief Sets uniform maximum cutoff across all element pairs.
   * @param[in] max_cutoff Maximum cutoff distance in Angstroms (> 0).
   */
  void setUniformCutoff(float max_cutoff);

  /**
   * @brief Adjusts min or max bond cutoffs using covalent radii sum scaled by factor.
   * @param[in] factor Covalent scale factor (> 0).
   * @param[in] bound Target bound (Min or Max).
   */
  void applyCovalentFactor(float factor, FactorBound bound);

  /**
   * @brief Applies a uniform global cutoff to all element pairs.
   * @param[in] global_cutoff Cutoff distance in Angstroms (> 0).
   */
  void applyGlobalCutoff(float global_cutoff);

private:
  AppWindow &window_;
  TrajectoryLoader &loader_;
  ProgramOptions &options_;
};

} // namespace correlation::app
