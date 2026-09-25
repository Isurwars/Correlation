/**
 * @file OptionsResetService.hpp
 * @brief Stateless service restoring default analysis options in UI cards.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "AppWindow.h"

namespace correlation::app {

class AppBackend;
class TrajectoryLoader;

/**
 * @class OptionsResetService
 * @brief Stateless utility class providing default options restoration for UI analysis cards.
 */
class OptionsResetService {
public:
  /**
   * @brief Restores default radial distribution (RDF) options in the UI.
   * @param[in,out] window Target UI window.
   */
  static void resetRDF(AppWindow &window);

  /**
   * @brief Restores default bond and dihedral angle options in the UI.
   * @param[in,out] window Target UI window.
   */
  static void resetAngle(AppWindow &window);

  /**
   * @brief Restores default structure factor (SQ) options in the UI.
   * @param[in,out] window Target UI window.
   */
  static void resetSQ(AppWindow &window);

  /**
   * @brief Restores default powder XRD options in the UI.
   * @param[in,out] window Target UI window.
   */
  static void resetXRD(AppWindow &window);

  /**
   * @brief Restores default topological rings options in the UI.
   * @param[in,out] window Target UI window.
   */
  static void resetRings(AppWindow &window);

  /**
   * @brief Restores default smoothing options in the UI.
   * @param[in,out] window Target UI window.
   */
  static void resetSmoothing(AppWindow &window);

  /**
   * @brief Restores default advanced order parameter options in the UI.
   * @param[in,out] window Target UI window.
   */
  static void resetAdvanced(AppWindow &window);

  /**
   * @brief Restores default trajectory analysis options in the UI.
   * @param[in,out] window Target UI window.
   * @param[in] loader Source TrajectoryLoader containing frame count and trajectory data.
   */
  static void resetTrajectory(AppWindow &window, const TrajectoryLoader &loader);

  /**
   * @brief Restores default trajectory analysis options in the UI (transitional).
   * @param[in,out] window Target UI window.
   * @param[in] backend Source backend containing frame count and trajectory data.
   */
  static void resetTrajectory(AppWindow &window, const AppBackend &backend);

  /**
   * @brief Restores default publication export settings in the UI.
   * @param[in,out] window Target UI window.
   */
  static void resetExportSettings(AppWindow &window);

  /**
   * @brief Restores default material type in the UI.
   * @param[in,out] window Target UI window.
   */
  static void resetMaterialType(AppWindow &window);
};

} // namespace correlation::app
