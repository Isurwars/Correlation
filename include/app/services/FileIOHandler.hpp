/**
 * @file FileIOHandler.hpp
 * @brief Handles file operations and native file dialogs.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "app/AnalysisDispatcher.hpp"
#include "app/AppOptions.hpp"
#include "app/TrajectoryLoader.hpp"
#include <atomic>
#include <string>
#include <thread>

class AppWindow;

namespace correlation::app {

class AppBackend;
class AppController;

/**
 * @class FileIOHandler
 * @brief Manages open/save file dialogs and background file loading.
 */
class FileIOHandler {
public:
  /**
   * @brief Constructs the FileIOHandler with injected services.
   * @param[in,out] window Reference to the UI window.
   * @param[in,out] loader Reference to the trajectory loader.
   * @param[in,out] dispatcher Reference to the analysis dispatcher.
   * @param[in,out] options Reference to program options.
   * @param[in,out] controller Reference to the main AppController.
   */
  FileIOHandler(::AppWindow &window, TrajectoryLoader &loader, AnalysisDispatcher &dispatcher,
                ProgramOptions &options, AppController &controller);

  /**
   * @brief Constructs the FileIOHandler (transitional).
   * @param[in,out] window Reference to the UI window.
   * @param[in,out] backend Reference to the application backend.
   * @param[in,out] controller Reference to the main AppController.
   */
  FileIOHandler(::AppWindow &window, AppBackend &backend, AppController &controller);

  /**
   * @brief Destructor. Ensures background threads are joined.
   */
  ~FileIOHandler();

  /**
   * @brief Displays file open dialog asynchronously and triggers background loading of selected
   * structure/trajectory file.
   */
  void handleBrowseFile();

  /**
   * @brief Reloads the currently loaded structure/trajectory file from disk asynchronously.
   */
  void handleReloadFile();

  /**
   * @brief Displays file save dialog asynchronously and exports active structural data.
   */
  void handleWriteFiles();

private:
  /**
   * @brief Initiates asynchronous trajectory loading on the background worker thread.
   * @param[in] filepath Absolute path to structure or trajectory file.
   */
  void startLoadingTrajectory(const std::string &filepath);

  /**
   * @brief Updates simulation cell diagnostics (volume, dimensions, density) on the UI.
   * @param[in] cell Pointer to the simulation cell, or nullptr if unavailable.
   */
  void updateBoxDiagnostics(const core::Cell *cell);

  /**
   * @brief Updates loaded file metadata (basename, frame counts, atom counts) on the UI.
   * @param[in] filepath Path to the loaded structure or trajectory file.
   */
  void updateFileMetadata(const std::string &filepath);

  /**
   * @brief Executes export of active structural datasets based on the selected file path extension.
   * @param[in] filepath Absolute base path for writing files.
   */
  void executeWriteFiles(const std::string &filepath);

  ::AppWindow &window_;
  TrajectoryLoader &loader_;
  AnalysisDispatcher &dispatcher_;
  ProgramOptions &options_;
  AppController &controller_;

  std::thread dialog_thread_; ///< Background worker thread for native file dialogs
  std::thread load_thread_;   ///< Background thread for loading files without blocking UI
  std::atomic<bool> dialog_active_{
      false}; ///< Concurrency guard preventing duplicate dialog launches
};

} // namespace correlation::app
