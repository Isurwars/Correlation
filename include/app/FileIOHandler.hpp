/**
 * @file FileIOHandler.hpp
 * @brief Handles file operations and native file dialogs.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "app/AppBackend.hpp"
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
   * @brief Destructor. Ensures background load thread is joined.
   */
  ~FileIOHandler();

  /**
   * @brief Displays file open dialog and triggers background loading of selected structure/trajectory file.
   */
  void handleBrowseFile();

  /**
   * @brief Displays file save dialog and exports active structural data.
   */
  void handleWriteFiles();

private:
  ::AppWindow &window_;
  AppBackend &backend_;
  AppController &controller_;

  std::thread load_thread_; ///< Background thread for loading files without blocking UI
};

} // namespace correlation::app
