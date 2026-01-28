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

class AppController {
public:
  // Constructor to set up the UI and backend references
  AppController(AppWindow &ui, AppBackend &backend);
  ~AppController();

private:
  // Pointers to the UI and backend objects
  AppWindow &ui_;
  AppBackend &backend_;

  // Handle to the file dialog, moved here from main.cpp
  std::unique_ptr<pfd::open_file> current_file_dialog_;
  
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

  // Member functions for the UI callbacks
  void handleRunAnalysis();
  void handleWriteFiles();
  void handleBrowseFile();
  void handleCheckFileDialogStatus();
};
