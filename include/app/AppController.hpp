/**
 * @file AppController.hpp
 * @brief High-level application controller orchestrating analysis workflows.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "AppWindow.h"
#include "app/AppBackend.hpp"
#include "app/PresetManager.hpp"
#include "plotters/SvgPlotter.hpp"
#include <nfd.h>

#include <memory>
#include <thread>

class AppControllerTests;

namespace correlation::app {

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
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//

  /**
   * @brief Constructs the AppController.
   * @param ui Reference to the main application window.
   * @param backend Reference to the application backend.
   */
  AppController(AppWindow &ui, AppBackend &backend);

  /**
   * @brief Destructor. Ensures analysis threads are joined before destruction.
   */
  ~AppController();

private:
  friend class ::AppControllerTests;
  AppWindow &ui_;       ///< Reference to the managed UI window.
  AppBackend &backend_; ///< Reference to the logic backend.



  std::thread analysis_thread_; ///< Handle for the background analysis computation.
  std::thread load_thread_;     ///< Handle for the background file loading process.

  std::vector<std::string> available_plot_keys_; ///< Map of UI indices to plot data keys.

  struct PinnedRun {
    std::string label;
    std::map<std::string, correlation::analysis::Histogram> histograms;
  };
  std::vector<PinnedRun> pinned_runs_; ///< Pinned runs for comparison plots.
  std::vector<Preset> presets_;        ///< List of loaded parameter presets.

  float last_mouse_x_ = -1.0f;
  float last_mouse_y_ = -1.0f;
  bool mouse_hover_ = false;
  float last_plot_width_ = 0.0f;
  float last_plot_height_ = 0.0f;

  /**
   * @brief Helper to update the UI progress bar and status text safely.
   * @param p Progress value between 0.0 and 1.0.
   * @param msg Status message to display in the UI.
   */
  void updateProgress(float p, const std::string &msg);

  /**
   * @brief Propagates the current backend options to the UI elements.
   * @param ui Reference to the main application window.
   */
  void handleOptionstoUI(AppWindow &ui);

  /**
   * @brief Updates active group flags in the UI based on active calculators.
   * @param ui Reference to the main application window.
   */
  void updateActiveGroupFlags(AppWindow &ui);

  /**
   * @brief Retrieves the current user-selected options from the UI.
   * @param ui Reference to the main application window.
   * @return The populated ProgramOptions struct.
   */
  ProgramOptions handleOptionsfromUI(AppWindow &ui);

  /**
   * @brief Updates the UI with the recommended bond cutoffs from the backend.
   * @param ui Reference to the main application window.
   */
  void setBondCutoffs(AppWindow &ui);

  /**
   * @brief Constructs a SvgPlotter PlotConfig based on current UI settings.
   */
  correlation::plotters::PlotConfig buildPlotConfigFromUI();

  /**
   * @brief Parses and retrieves the user-modified bond cutoffs from the UI.
   * @param ui Reference to the main application window.
   * @return A nested vector representing the bond cutoff matrix.
   */
  std::vector<std::vector<double>> getBondCutoffs(AppWindow &ui);

  /**
   * @brief Populates the UI calculator groups from CalculatorFactory.
   * @param ui Reference to the main application window.
   */
  void populateCalculatorGroups(AppWindow &ui);

  /**
   * @brief Validates all numeric input fields and pushes error states to the UI.
   * @return true if all inputs are valid, false otherwise.
   */
  bool validateInputs();

  /**
   * @brief Generates and updates the equivalent CLI command in the UI.
   */
  void updateCliCommand();

  /**
   * @brief Copies the equivalent CLI command to the clipboard.
   */
  void handleCopyCliCommand();

  //-------------------------------------------------------------------------//
  //-------------------------------- Methods --------------------------------//
  //-------------------------------------------------------------------------//

  /**
   * @brief Handles the "Run Analysis" signal from the UI.
   * Starts the analysis in a separate thread.
   */
  void handleRunAnalysis();

  /**
   * @brief Handles the "Write Files" signal from the UI.
   * Opens a save file dialog to select the output location.
   */
  void handleWriteFiles();

  /**
   * @brief Handles the "Browse File" signal from the UI.
   * Opens an open file dialog to select the input trajectory/structure.
   */
  void handleBrowseFile();



  /**
   * @brief Populates the UI plot dropdown with the available histogram names
   *        from the last completed analysis.
   */
  void populatePlotList();

  /**
   * @brief Handles the "Select Plot" signal from the UI.
   *        Generates an SVG image for the requested histogram and pushes it
   *        to the UI via `ui_.set_preview_plot()`.
   * @param index The index of the selected plot in the dropdown menu.
   */
  void handleSelectPlot(int index);

  /**
   * @brief Handles mouse movements over the preview plot area.
   */
  void handleMouseMove(float mx, float my, bool hover, float w, float h);

  /**
   * @brief Handles the "Save Plot" signal from the UI.
   * Opens a save file dialog to export the currently selected plot (SVG or PDF).
   */
  void handleSavePlot();

  /**
   * @brief Handles pinning the current run for comparison.
   */
  void handlePinRun();

  /**
   * @brief Handles clearing all pinned runs.
   */
  void handleClearPinnedRuns();

  /**
   * @brief Preset management handlers.
   */
  void handleLoadPreset(int index);
  void handleSavePreset(const std::string &name);
  void handleDeletePreset(int index);
  void refreshPresetList();
};
} // namespace correlation::app
