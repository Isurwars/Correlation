// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <map>
#include <memory>
#include <functional>

#include "Smoothing.hpp"
#include "Trajectory.hpp"
#include "TrajectoryAnalyzer.hpp"
#include "DistributionFunctions.hpp"

// Encapsulates all command-line configurable options for the application.
struct ProgramOptions {
  std::string input_file;
  std::string output_file_base;
  bool smoothing = true;
  bool use_hdf5 = true;
  bool use_csv = false;
  double r_max = 20.0;
  double r_bin_width = 0.02;
  double q_max = 20.0;
  double q_bin_width = 0.02;
  double r_int_max = 10.0;
  double angle_bin_width = 1.0;
  double smoothing_sigma = 0.1;
  KernelType smoothing_kernel = KernelType::Gaussian;
  int min_frame = 0;
  int max_frame = -1;
  std::vector<std::vector<double>> bond_cutoffs_sq_;
};

class AppBackend {
public:
  // Constructor to set up the UI and backend references
  AppBackend();

  // Accesors
  void setOptions(const ProgramOptions &opt) { options_ = opt; }
  ProgramOptions options() { return options_; }
  const Cell *cell() const {
    if (trajectory_ && !trajectory_->getFrames().empty()) {
      return &trajectory_->getFrames()[0];
    }
    return nullptr;
  }

  // Member functions
  std::string load_file(const std::string &path);
  void run_analysis();
  void write_files();
  std::map<std::string, int> getAtomCounts() const;
  int getFrameCount() const;
  int getTotalAtomCount() const;
  size_t getRemovedFrameCount() const;
  std::vector<std::vector<double>> getRecommendedBondCutoffs() const;
  double getBondCutoff(int, int);
  void setBondCutoffs(const std::vector<std::vector<double>> &cutoffs);

  // Callbacks
  void setProgressCallback(std::function<void(float)> cb) { progress_callback_ = cb; }

private:
  // Member functions
  void analysis_thread_func();

  // Pointers to the Trajectory and DF
  std::unique_ptr<Trajectory> trajectory_;
  std::unique_ptr<TrajectoryAnalyzer> trajectory_analyzer_;
  std::unique_ptr<DistributionFunctions> df_;

  ProgramOptions options_;
  std::function<void(float)> progress_callback_;
};
