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
  bool smoothing = false;
  double r_max = 20.0;
  double r_bin_width = 0.02;
  double q_max = 20.0;
  double q_bin_width = 0.02;
  double r_int_max = 10.0;
  double angle_bin_width = 1.0;
  double smoothing_sigma = 0.1;
  KernelType smoothing_kernel = KernelType::Gaussian;
  std::vector<std::vector<double>> bond_cutoffs_sq_;
};

class AppBackend {
public:
  // Constructor to set up the UI and backend references
  AppBackend();

  // Accesors
  void setOptions(ProgramOptions opt) { options_ = opt; }
  ProgramOptions options() { return options_; }
  const std::unique_ptr<Cell> &cell() const { return cell_; }

  // Member functions
  std::string load_file(const std::string &path);
  void run_analysis();
  std::map<std::string, int> getAtomCounts() const;
  std::vector<std::vector<double>> getRecommendedBondCutoffs() const;
  double getBondCutoff(int, int);
  void setBondCutoffs(const std::vector<std::vector<double>> &cutoffs);

private:
  // Member functions
  void analysis_thread_func();

  // Pointers to the Cell and DF
  std::unique_ptr<Cell> cell_;
  std::unique_ptr<DistributionFunctions> df_;

  ProgramOptions options_;
};
