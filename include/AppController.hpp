// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <memory>
#include <thread>

#include "AppBackend.hpp"
#include "AppWindow.h"
#include "PortableFileDialogs.hpp"

/**
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
  // Pointers to the UI and backend objects
  AppWindow &ui_;
  AppBackend &backend_;

  // Handle to the file dialog, moved here from main.cpp
  std::unique_ptr<pfd::open_file> current_file_dialog_;
  std::unique_ptr<pfd::save_file> current_save_dialog_;

  // Analysis thread
  std::thread analysis_thread_;

  // Progress helper
  void updateProgress(float p);

  // functions to handle options
  void handleOptionstoUI(AppWindow &ui);
  ProgramOptions handleOptionsfromUI(AppWindow &ui);

  // functios to handle Bond_Cutoffs_Sq_
  void setBondCutoffs(AppWindow &ui);
  std::vector<std::vector<double>> getBondCutoffs(AppWindow &ui);

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
};
