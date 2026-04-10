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
#include "external/PortableFileDialogs.hpp"

#include <memory>
#include <thread>

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
  AppWindow &ui_;         ///< Reference to the managed UI window.
  AppBackend &backend_;   ///< Reference to the logic backend.

  std::unique_ptr<pfd::open_file> current_file_dialog_; ///< Active file open dialog handler.
  std::unique_ptr<pfd::save_file> current_save_dialog_; ///< Active file save dialog handler.

  std::thread analysis_thread_; ///< Handle for the background analysis computation.
  std::thread load_thread_;     ///< Handle for the background file loading process.

  std::vector<std::string> available_plot_keys_; ///< Map of UI indices to plot data keys.

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
   * @brief Periodic check for file dialog status.
   * Since file dialogs might be non-blocking or running in a separate
   * process/thread, this method checks if they have returned a result.
   */
  void handleCheckFileDialogStatus();

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
};
} // namespace correlation::app
