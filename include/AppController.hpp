// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <memory>

#include "AppBackend.hpp"
#include "PortableFileDialogs.hpp"
#include "app_window.h"

class AppController {
public:
  // Constructor to set up the UI and backend references
  AppController(AppWindow &ui, AppBackend &backend);

private:
  // Pointers to the UI and backend objects
  AppWindow &ui_;
  AppBackend &backend_;

  // Handle to the file dialog, moved here from main.cpp
  std::unique_ptr<pfd::open_file> current_file_dialog_;

  // Member functions for the UI callbacks
  void handleRunAnalysis();
  void handleBrowseFile();
  void handleCheckFileDialogStatus();
};
