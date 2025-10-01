// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <memory>

#include "Cell.hpp"
#include "DistributionFunctions.hpp"

// Encapsulates all command-line configurable options for the application.
struct ProgramOptions {
  std::string input_file;
  std::string output_file_base;
  bool normalize = true;
  bool smoothing = false;
  double r_max = 20.0;
  double r_bin_width = 0.02;
  double q_max = 20.0;
  double q_bin_width = 0.02;
  double r_int_max = 10.0;
  double angle_max = 180.0;
  double angle_bin_width = 1.0;
  double bond_factor = 1.2;
  double smoothing_sigma = 0.1;
  KernelType smoothing_kernel = KernelType::Gaussian;
};

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
