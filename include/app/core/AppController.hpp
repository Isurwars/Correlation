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

#include <functional>
#include <memory>

class AppControllerTests;
class AppWindow;

namespace correlation::app {

class AnalysisRunner;
class FileIOHandler;
class InputValidator;
class PlotController;
class PresetController;
class Structure3DController;

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
   * @brief Loads a structure or trajectory file asynchronously.
   * @param filepath Path to the structure or trajectory file.
   * @param on_complete Optional callback invoked on event loop when loading finishes.
   */
  void loadTrajectory(const std::string &filepath, std::function<void(bool)> on_complete = nullptr);

  /**
   * @brief Executes configured analyses on the active trajectory.
   */
  void runAnalysis();

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
   * @brief Populates the export analyses selection list with computed histograms.
   */
  void populateExportAnalyses();

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

  /**
   * @brief Returns the Structure3DController instance.
   */
  [[nodiscard]] Structure3DController *getStructure3DController() noexcept {
    return structure_3d_controller_.get();
  }

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
  std::unique_ptr<Structure3DController> structure_3d_controller_;
  AppSettings settings_;

  void registerActionCallbacks();
  void registerPlotCallbacks();
  void registerStructure3DCallbacks();
  void registerPresetAndOptionCallbacks();
  void registerExportCallbacks();
};
} // namespace correlation::app
