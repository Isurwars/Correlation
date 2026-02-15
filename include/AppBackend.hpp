// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <functional>
#include <map>
#include <memory>

#include "DistributionFunctions.hpp"
#include "Smoothing.hpp"
#include "Trajectory.hpp"
#include "TrajectoryAnalyzer.hpp"

/**
 * @brief Encapsulates all command-line configurable options for the application.
 */
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
  double time_step = 1.0;
  std::vector<std::vector<double>> bond_cutoffs_sq_;
};

/**
 * @brief The main backend class for the application.
 *
 * This class orchestrates the loading of files, setting of options, running of
 * analyses, and writing of results. It acts as the bridge between the UI/Controller
 * and the core data/analysis logic.
 */
class AppBackend {
public:
  //-------------------------------------------------------------------------//
  //----------------------------- Constructors ------------------------------//
  //-------------------------------------------------------------------------//

  /**
   * @brief Default constructor.
   */
  AppBackend();

  //-------------------------------------------------------------------------//
  //------------------------------- Accessors -------------------------------//
  //-------------------------------------------------------------------------//

  /**
   * @brief Sets the program options.
   * @param opt The options to set.
   */
  void setOptions(const ProgramOptions &opt) { options_ = opt; }

  /**
   * @brief Gets the current program options.
   * @return Copy of the current options.
   */
  ProgramOptions options() { return options_; }

  /**
   * @brief Gets a pointer to the current cell (first frame of trajectory).
   * @return Pointer to the Cell, or nullptr if no trajectory loaded.
   */
  const Cell *cell() const {
    if (trajectory_ && !trajectory_->getFrames().empty()) {
      return &trajectory_->getFrames()[0];
    }
    return nullptr;
  }

  //-------------------------------------------------------------------------//
  //-------------------------------- Methods --------------------------------//
  //-------------------------------------------------------------------------//

  /**
   * @brief Loads a file from a given path.
   * @param path Absolute path to the file.
   * @return A status message indicating success or details about the loaded file.
   */
  std::string load_file(const std::string &path);

  /**
   * @brief Runs the analysis based on the current options and loaded trajectory.
   *
   * This method performs the following:
   * 1. Sets up bond cutoffs and time steps.
   * 2. Initializes the TrajectoryAnalyzer.
   * 3. Computes the mean Distribution Functions (RDF, ADF, etc.).
   * 4. Calculates VACF and VDOS if applicable.
   */
  void run_analysis();

  /**
   * @brief Writes the analysis results to files (CSV, HDF5) as specified in options.
   */
  void write_files();

  /**
   * @brief Gets the atom counts for the current structure.
   * @return A map of Element Symbol -> Count.
   */
  std::map<std::string, int> getAtomCounts() const;

  /**
   * @brief Gets the total number of frames in the trajectory.
   * @return Number of frames.
   */
  int getFrameCount() const;

  /**
   * @brief Gets the total number of atoms in the first frame.
   * @return Number of atoms.
   */
  int getTotalAtomCount() const;

  /**
   * @brief Gets the count of frames removed/skipped during loading/processing.
   * @return Number of removed frames.
   */
  size_t getRemovedFrameCount() const;

  /**
   * @brief Gets the time step of the trajectory.
   * @return Time step in femtoseconds.
   */
  double getTimeStep() const;

  /**
   * @brief Calculates recommended bond cutoffs based on the first pair density minimum.
   * @return A matrix of cutoffs where entry [i][j] is the cutoff for pair i-j.
   */
  std::vector<std::vector<double>> getRecommendedBondCutoffs() const;

  /**
   * @brief Gets the bond cutoff for a specific pair of element types.
   * @param type1 Index of the first element type.
   * @param type2 Index of the second element type.
   * @return The cutoff distance.
   */
  double getBondCutoff(int type1, int type2);

  /**
   * @brief Sets the bond cutoffs to be used in analysis.
   * @param cutoffs Matrix of cutoffs.
   */
  void setBondCutoffs(const std::vector<std::vector<double>> &cutoffs);

  // Callbacks
  void setProgressCallback(std::function<void(float)> cb) {
    progress_callback_ = cb;
  }

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
