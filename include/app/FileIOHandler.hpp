/**
 * @file FileIOHandler.hpp
 * @brief Handles file operations and native file dialogs.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "app/AppBackend.hpp"
#include <atomic>
#include <string>
#include <thread>

class AppWindow;

namespace correlation::app {

class AppController;

/**
 * @class FileIOHandler
 * @brief Manages open/save file dialogs and background file loading.
 */
class FileIOHandler {
public:
  /**
   * @brief Constructs the FileIOHandler.
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
   * @brief Displays file open dialog asynchronously and triggers background loading of selected structure/trajectory
   * file.
   */
  void handleBrowseFile();

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
   * @brief Executes export of active structural datasets based on the selected file path extension.
   * @param[in] filepath Absolute base path for writing files.
   */
  void executeWriteFiles(const std::string &filepath);

  ::AppWindow &window_;
  AppBackend &backend_;
  AppController &controller_;

  std::thread dialog_thread_;              ///< Background worker thread for native file dialogs
  std::thread load_thread_;                ///< Background thread for loading files without blocking UI
  std::atomic<bool> dialog_active_{false}; ///< Concurrency guard preventing duplicate dialog launches
};

} // namespace correlation::app
