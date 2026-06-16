/**
 * @file FileIOHandler.hpp
 * @brief Handles file operations and native file dialogs.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "AppWindow.h"
#include "app/AppBackend.hpp"
#include <thread>

namespace correlation::app {

class AppController;

/**
 * @class FileIOHandler
 * @brief Manages open/save file dialogs and background file loading.
 */
class FileIOHandler {
public:
  FileIOHandler(AppWindow &window, AppBackend &backend, AppController &controller);
  ~FileIOHandler();

  void handleBrowseFile();
  void handleWriteFiles();

private:
  AppWindow &window_;
  AppBackend &backend_;
  AppController &controller_;

  std::thread load_thread_; ///< Background thread for loading files without blocking UI
};

} // namespace correlation::app
