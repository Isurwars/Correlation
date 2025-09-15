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
// The Slint UI will interact with this object.
class AppBackend {
public:
  AppBackend();
  void load_file(const std::string &path);
  void run_analysis(const ProgramOptions &options);

private:
  void analysis_thread_func();

  std::unique_ptr<Cell> cell_;
  ProgramOptions options_;
  std::unique_ptr<DistributionFunctions> df_;
};
