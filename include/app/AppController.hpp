/**
 * @file AppController.hpp
 * @brief High-level application controller orchestrating analysis workflows.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "AppWindow.h"
#include "app/AppBackend.hpp"
#include <nfd.h>

#include <memory>

class AppControllerTests;

namespace correlation::app {

class AnalysisRunner;
class FileIOHandler;
class InputValidator;
class PlotController;
class PresetController;

/**
 * @class AppController
 * @brief Controller class for the application.
 *
 * This class handles the interaction between the User Interface (AppWindow)
 * and the logic Backend (AppBackend). It manages event handling, threading for
 * analysis, and data synchronization between UI and Backend.
 */
class AppController {
public:
  /** @name Constructors & Methods */
  ///@{

  /**
   * @brief Constructs the AppController.
   * @param window Reference to the main application window.
   * @param backend Reference to the application backend.
   */
  AppController(AppWindow &window, AppBackend &backend);

  /**
   * @brief Destructor. Ensures analysis threads are joined before destruction.
   */
  ~AppController();

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
   * @brief Parses and retrieves the user-modified bond cutoffs from the UI.
   * @return A nested vector representing the bond cutoff matrix.
   */
  std::vector<std::vector<real_t>> getBondCutoffs();
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

  ///@}

private:
  friend class ::AppControllerTests;
  AppWindow &window_;   ///< Reference to the managed UI window.
  AppBackend &backend_; ///< Reference to the logic backend.

  std::unique_ptr<AnalysisRunner> analysis_runner_;
  std::unique_ptr<FileIOHandler> file_io_handler_;
  std::unique_ptr<InputValidator> input_validator_;
  std::unique_ptr<PlotController> plot_controller_;
  std::unique_ptr<PresetController> preset_controller_;
};
} // namespace correlation::app
