#ifndef APP_BACKEND_HPP_
#define APP_BACKEND_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

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
  void run_analysis();

private:
  void analysis_thread_func();

  std::unique_ptr<Cell> cell_;
  ProgramOptions options_;
  std::unique_ptr<DistributionFunctions> df_;
};

#endif // APP_BACKEND_HPP_
