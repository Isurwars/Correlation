/**
 * @file AppController.hpp
 * @brief High-level application controller orchestrating analysis workflows.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "AppWindow.h"
#include "app/core/AppOptions.hpp"
#include "app/core/SettingsManager.hpp"
#include "app/services/AnalysisDispatcher.hpp"
#include "app/services/BondCutoffService.hpp"
#include "app/services/TrajectoryLoader.hpp"
#include "app/viewmodel/BondCutoffController.hpp"
#include <nfd.h>

#include <memory>

class AppControllerTests;
class AppWindow;

namespace correlation::app {

class AnalysisRunner;
class FileIOHandler;
class InputValidator;
class PlotController;
class PresetController;

/**
 * @class AppController
 * @brief Controller class for the application orchestrating services and UI.
 *
 * This class handles the interaction between the User Interface (AppWindow)
 * and domain services (TrajectoryLoader, AnalysisDispatcher, BondCutoffService).
 * It manages event handling, threading for analysis, and data synchronization.
 */
class AppController {
public:
  /** @name Constructors & Methods */
  ///@{

  /**
   * @brief Constructs the AppController with injected domain services.
   * @param window Reference to the main application window.
   * @param loader Reference to the trajectory loading service.
   * @param dispatcher Reference to the analysis dispatching service.
   * @param options Reference to application configuration options.
   */
  AppController(AppWindow &window, TrajectoryLoader &loader, AnalysisDispatcher &dispatcher,
                ProgramOptions &options);

  /**
   * @brief Constructs the AppController with injected domain services (transitional).
   * @param window Reference to the main application window.
   * @param loader Reference to the trajectory loading service.
   * @param dispatcher Reference to the analysis dispatching service.
   * @param cutoff_service Reference to the bond cutoff service.
   * @param options Reference to application configuration options.
   */
  AppController(AppWindow &window, TrajectoryLoader &loader, AnalysisDispatcher &dispatcher,
                BondCutoffService &cutoff_service, ProgramOptions &options);

  /**
   * @brief Destructor. Ensures analysis threads are joined before destruction.
   */
  ~AppController();

  /**
   * @brief Loads application layout settings and applies window/column geometry.
   */
  void loadSettings();

  /**
   * @brief Saves current layout settings and window geometry to disk.
   */
  void saveSettings() const;

  /**
   * @brief Retrieves the current user-selected options from the UI.
   * @return The populated ProgramOptions struct.
   */
  ProgramOptions handleOptionsfromUI();

  /**
   * @brief Populates the UI with options from the backend.
   */
  void handleOptionstoUI();

  /**
   * @brief Updates the UI with the recommended bond cutoffs from the backend.
   */
  void setBondCutoffs();

  /**
   * @brief Parses and retrieves the user-modified bond cutoff matrix from the UI.
   * @return A BondCutoffMatrix representing both min and max squared cutoffs.
   */
  correlation::analysis::BondCutoffMatrix getBondCutoffs();

  /**
   * @brief Restores default radial distribution (RDF) options in the UI.
   */
  void handleResetRDFOptions();

  /**
   * @brief Restores default bond and dihedral angle options in the UI.
   */
  void handleResetAngleOptions();

  /**
   * @brief Restores default structure factor (SQ) options in the UI.
   */
  void handleResetSQOptions();

  /**
   * @brief Restores default powder XRD options in the UI.
   */
  void handleResetXRDOptions();

  /**
   * @brief Updates XRD wavelength based on radiation source preset.
   * @param preset_idx Selected preset index.
   */
  void handleXRDPresetChanged(int preset_idx);

  /**
   * @brief Scales recommended covalent cutoffs and updates the UI matrix.
   * @param scale_factor Multiplier applied to covalent distances.
   */
  void handleApplyScaledCutoffs(float scale_factor);

  /**
   * @brief Sets uniform cutoffs across all pairs and updates the UI matrix.
   * @param max_cutoff Maximum cutoff distance in Å.
   */
  void handleSetUniformCutoff(float max_cutoff);

  /**
   * @brief Applies min factor multiplier to covalent radii sums for minimum bond cutoffs.
   * @param min_factor Multiplier applied to sum of covalent radii for min distance.
   */
  void handleApplyMinFactor(float min_factor);

  /**
   * @brief Applies max factor multiplier to covalent radii sums for maximum bond cutoffs.
   * @param max_factor Multiplier applied to sum of covalent radii for max distance.
   */
  void handleApplyMaxFactor(float max_factor);

  /**
   * @brief Applies global uniform cutoff across all atom pairs.
   * @param global_cutoff Maximum cutoff distance in Å.
   */
  void handleApplyGlobalCutoff(float global_cutoff);

  /**
   * @brief Restores default topological rings options in the UI.
   */
  void handleResetRingsOptions();

  /**
   * @brief Restores default smoothing options in the UI.
   */
  void handleResetSmoothingOptions();

  /**
   * @brief Restores default advanced order parameter options in the UI.
   */
  void handleResetAdvancedOptions();

  /**
   * @brief Restores default trajectory analysis options in the UI.
   */
  void handleResetTrajectoryOptions();

  /**
   * @brief Restores default publication export settings in the UI.
   */
  void handleResetExportSettings();

  /**
   * @brief Restores default calculator selections in the UI.
   */
  void handleResetAnalysesSelection();

  /**
   * @brief Restores default material type in the UI.
   */
  void handleResetMaterialType();

  /**
   * @brief Clears all comparison curves from the plot overlay.
   */
  void handleClearComparisonCurves();

  /**
   * @brief Populates the UI calculator groups from CalculatorFactory.
   */
  void populateCalculatorGroups();

  /**
   * @brief Updates active group flags in the UI based on active calculators.
   */
  void updateActiveGroupFlags();

  /**
   * @brief Returns the InputValidator instance.
   */
  InputValidator *getInputValidator() { return input_validator_.get(); }

  /**
   * @brief Returns the PlotController instance.
   */
  PlotController *getPlotController() { return plot_controller_.get(); }

  /**
   * @brief Returns the BondCutoffController instance.
   */
  BondCutoffController &getBondCutoffController() { return bond_cutoff_controller_; }

  /**
   * @brief Returns the const BondCutoffController instance.
   */
  const BondCutoffController &getBondCutoffController() const { return bond_cutoff_controller_; }

  ///@}

private:
  friend class ::AppControllerTests;
  AppWindow &window_;
  TrajectoryLoader &loader_;
  AnalysisDispatcher &dispatcher_;
  ProgramOptions &options_;

  BondCutoffController bond_cutoff_controller_;
  std::unique_ptr<AnalysisRunner> analysis_runner_;
  std::unique_ptr<FileIOHandler> file_io_handler_;
  std::unique_ptr<InputValidator> input_validator_;
  std::unique_ptr<PlotController> plot_controller_;
  std::unique_ptr<PresetController> preset_controller_;
  AppSettings settings_;
};
} // namespace correlation::app
