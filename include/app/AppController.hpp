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
#include "app/PresetManager.hpp"
#include "plotters/SvgPlotter.hpp"
#include <nfd.h>

#include <atomic>
#include <chrono>
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
   * @param window Reference to the main application window.
   * @param backend Reference to the application backend.
   */
  AppController(AppWindow &window, AppBackend &backend);

  /**
   * @brief Destructor. Ensures analysis threads are joined before destruction.
   */
  ~AppController();

private:
  friend class ::AppControllerTests;
  AppWindow &window_;       ///< Reference to the managed UI window.
  AppBackend &backend_; ///< Reference to the logic backend.

  std::thread analysis_thread_; ///< Handle for the background analysis computation.
  std::thread load_thread_;     ///< Handle for the background file loading process.
  std::thread render_thread_;   ///< Handle for the background SVG rendering process.

  struct RenderTaskData {
    correlation::analysis::Histogram active_hist;
    std::vector<std::pair<std::string, correlation::analysis::Histogram>> comparison_hists;
    correlation::plotters::PlotConfig config;
    correlation::plotters::HoverInfo hover;
    std::map<std::string, double> ashcroft_weights;
  };

  std::atomic<bool> is_rendering_{false};   ///< Whether a background render is currently active.
  std::atomic<bool> render_pending_{false}; ///< Whether another render request is pending.

  std::vector<std::string> available_plot_keys_; ///< Map of UI indices to plot data keys.

  struct PinnedRun {
    std::string label;
    std::map<std::string, correlation::analysis::Histogram> histograms;
  };
  std::vector<PinnedRun> pinned_runs_; ///< Pinned runs for comparison plots.
  std::vector<Preset> presets_;        ///< List of loaded parameter presets.

  float last_mouse_x_ = -1.0F;
  float last_mouse_y_ = -1.0F;
  bool mouse_hover_ = false;
  float last_plot_width_ = 0.0F;
  float last_plot_height_ = 0.0F;

  // Cache and throttling variables for mouse hover / select plot
  std::chrono::steady_clock::time_point last_replot_time_;
  slint::Timer hover_timer_;
  slint::Timer update_timer_;
  bool update_scheduled_ = false;
  int pending_plot_index_ = -1;
  bool needs_redraw_ = false;

  int last_rendered_index_ = -1;
  correlation::plotters::PlotConfig last_config_;
  correlation::plotters::HoverInfo last_hover_;
  std::size_t last_pinned_runs_count_ = 0;

  std::shared_ptr<std::string> current_svg_;

  /**
   * @brief Helper to update the UI progress bar and status text safely.
   * @param progress Progress value between 0.0 and 1.0.
   * @param msg Status message to display in the UI.
   */
  void updateProgress(float progress, const std::string &msg);

  /**
   * @brief Propagates the current backend options to the UI elements.
   */
  void handleOptionstoUI();

  /**
   * @brief Updates active group flags in the UI based on active calculators.
   */
  void updateActiveGroupFlags();

  /**
   * @brief Retrieves the current user-selected options from the UI.
   * @return The populated ProgramOptions struct.
   */
  ProgramOptions handleOptionsfromUI();

  /**
   * @brief Updates the UI with the recommended bond cutoffs from the backend.
   */
  void setBondCutoffs();

  /**
   * @brief Constructs a SvgPlotter PlotConfig based on current UI settings.
   */
  correlation::plotters::PlotConfig buildPlotConfigFromUI();

  /**
   * @brief Parses and retrieves the user-modified bond cutoffs from the UI.
   * @return A nested vector representing the bond cutoff matrix.
   */
  std::vector<std::vector<double>> getBondCutoffs();

  /**
   * @brief Populates the UI calculator groups from CalculatorFactory.
   */
  void populateCalculatorGroups();

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
   * @brief Handles mouse movements over the preview plot area.
   */
  void handleMouseMove(float mouse_x, float mouse_y, bool hover, float width, float height);

  /**
   * @brief Requests a plot update, throttling redraws to every 200ms unless immediate is true.
   * @param index The index of the selected plot.
   * @param immediate If true, bypasses the throttle and redraws immediately.
   */
  void requestPlotUpdate(int index, bool immediate = false);

  /**
   * @brief Checks if the given plot parameters match the previously rendered plot.
   */
  bool isPlotCacheHit(int index, const correlation::plotters::PlotConfig &config,
                      const correlation::plotters::HoverInfo &hover) const;

  /**
   * @brief Executes the SVG rendering in a background thread and updates the UI.
   */
  void executePlotRender(RenderTaskData data);

  /**
   * @brief Handles the "Save Plot" signal from the UI.
   * Opens a save file dialog to export the currently selected plot (SVG or PDF).
   */
  void handleSavePlot();

  /**
   * @brief Executes the saving logic for a plot to reduce handleSavePlot complexity.
   */
  void executeSavePlot(const std::string &filepath, const correlation::analysis::Histogram *hist,
                       const std::string &name);

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
  void handleMaterialTypeChanged(int type);
};
} // namespace correlation::app
