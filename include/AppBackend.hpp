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
 * @brief Default values and messages for the application.
 */
struct AppDefaults {
  static constexpr double R_MAX = 20.0;
  static constexpr double R_BIN_WIDTH = 0.02;
  static constexpr double Q_MAX = 20.0;
  static constexpr double Q_BIN_WIDTH = 0.02;
  static constexpr double R_INT_MAX = 10.0;
  static constexpr double ANGLE_BIN_WIDTH = 1.0;
  static constexpr double SMOOTHING_SIGMA = 0.1;
  static constexpr decltype(KernelType::Gaussian) SMOOTHING_KERNEL =
      KernelType::Gaussian;
  static constexpr double TIME_STEP = 1.0;

  static constexpr const char *MSG_RUNNING_ANALYSIS = "Running Analysis...";
  static constexpr const char *MSG_ANALYSIS_ENDED = "Analysis ended.";
  static constexpr const char *MSG_SELECTING_OUTPUT =
      "Selecting output file...";
  static constexpr const char *MSG_FILE_SELECTION_CANCELLED =
      "File selection cancelled.";
  static constexpr const char *MSG_ERROR_LOADING = "Error loading file: ";
  static constexpr const char *MSG_FILES_WRITTEN = "Files Written.";
  static constexpr const char *MSG_SAVE_CANCELLED = "Save cancelled.";
};

/**
 * @brief Encapsulates all configurable options for the application.
 */
struct ProgramOptions {
  std::string input_file;
  std::string output_file_base;
  bool smoothing = true;
  bool use_hdf5 = true;
  bool use_csv = false;
  double r_max = AppDefaults::R_MAX;
  double r_bin_width = AppDefaults::R_BIN_WIDTH;
  double q_max = AppDefaults::Q_MAX;
  double q_bin_width = AppDefaults::Q_BIN_WIDTH;
  double r_int_max = AppDefaults::R_INT_MAX;
  double angle_bin_width = AppDefaults::ANGLE_BIN_WIDTH;
  double smoothing_sigma = AppDefaults::SMOOTHING_SIGMA;
  KernelType smoothing_kernel = AppDefaults::SMOOTHING_KERNEL;
  int min_frame = 0;
  int max_frame = -1;
  double time_step = AppDefaults::TIME_STEP;
  std::vector<std::vector<double>> bond_cutoffs_sq_;
};

/**
 * @brief The main backend class for the application.
 *
 * This class orchestrates the loading of files, setting of options, running of
 * analyses, and writing of results. It acts as the bridge between the
 * UI/Controller and the core data/analysis logic.
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
  [[nodiscard]] ProgramOptions options() const { return options_; }

  /**
   * @brief Gets a pointer to the current cell (first frame of trajectory).
   * @return Pointer to the Cell, or nullptr if no trajectory loaded.
   */
  [[nodiscard]] const Cell *cell() const {
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
   * @return A status message indicating success or details about the loaded
   * file.
   */
  std::string load_file(const std::string &path);

  /**
   * @brief Runs the analysis based on the current options and loaded
   * trajectory.
   *
   * This method performs the following:
   * 1. Sets up bond cutoffs and time steps.
   * 2. Initializes the TrajectoryAnalyzer.
   * 3. Computes the mean Distribution Functions (RDF, ADF, etc.).
   * 4. Calculates VACF and VDOS if applicable.
   */
  void run_analysis();

  /**
   * @brief Writes the analysis results to files (CSV, HDF5) as specified in
   * options.
   */
  void write_files();

  /**
   * @brief Gets the atom counts for the current structure.
   * @return A map of Element Symbol -> Count.
   */
  [[nodiscard]] std::map<std::string, int> getAtomCounts() const;

  /**
   * @brief Gets the total number of frames in the trajectory.
   * @return Number of frames.
   */
  [[nodiscard]] int getFrameCount() const;

  /**
   * @brief Gets the total number of atoms in the first frame.
   * @return Number of atoms.
   */
  [[nodiscard]] int getTotalAtomCount() const;

  /**
   * @brief Gets the count of frames removed/skipped during loading/processing.
   * @return Number of removed frames.
   */
  [[nodiscard]] size_t getRemovedFrameCount() const;

  /**
   * @brief Gets the time step of the trajectory.
   * @return Time step in femtoseconds.
   */
  [[nodiscard]] double getTimeStep() const;

  /**
   * @brief Calculates recommended bond cutoffs based on the first pair density
   * minimum.
   * @return A matrix of cutoffs where entry [i][j] is the cutoff for pair i-j.
   */
  [[nodiscard]] std::vector<std::vector<double>>
  getRecommendedBondCutoffs() const;

  /**
   * @brief Gets the bond cutoff for a specific pair of element types.
   * @param type1 Index of the first element type.
   * @param type2 Index of the second element type.
   * @return The cutoff distance.
   */
  [[nodiscard]] double getBondCutoff(int type1, int type2);

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
