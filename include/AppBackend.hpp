// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <memory>

#include "Cell.hpp"
#include "DistributionFunctions.hpp"
#include "ProgramOptions.hpp"

// This class holds the application state and logic.

class AppBackend {
public:
  // Constructor to set up the UI and backend references
  AppBackend();

  // Member functions
  void load_file(const std::string &path);
  void run_analysis(const ProgramOptions &options);

private:
  // Member functions
  void analysis_thread_func();

  // Pointers to the Cell and DF
  std::unique_ptr<Cell> cell_;
  std::unique_ptr<DistributionFunctions> df_;

  ProgramOptions options_;
};
