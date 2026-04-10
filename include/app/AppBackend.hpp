/**
 * @file AppBackend.hpp
 * @brief Application backend interface between the UI and analysis engine.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "analysis/TrajectoryAnalyzer.hpp"
#include "core/Trajectory.hpp"
#include "math/Smoothing.hpp"

#include <functional>
#include <map>
#include <memory>
namespace correlation::app {

/**
 * @brief Default values and messages for the application.
 */
struct AppDefaults {
  static constexpr double R_MAX = 20.0;           ///< Default max radius for RDF (Angstrom).
  static constexpr double R_BIN_WIDTH = 0.02;     ///< Default bin width for RDF (Angstrom).
  static constexpr double Q_MAX = 20.0;           ///< Default max q for S(Q) (Angstrom^-1).
  static constexpr double Q_BIN_WIDTH = 0.02;     ///< Default bin width for S(Q) (Angstrom^-1).
  static constexpr double R_INT_MAX = 10.0;       ///< Default max radius for integration (Angstrom).
  static constexpr double ANGLE_BIN_WIDTH = 1.0;  ///< Default bin width for ADF (Degrees).
  static constexpr double SMOOTHING_SIGMA = 0.1;  ///< Default Gaussian smoothing sigma.

  /** @brief Default smoothing kernel. */
  static constexpr decltype(correlation::math::KernelType::Gaussian)
      SMOOTHING_KERNEL = correlation::math::KernelType::Gaussian;

  static constexpr double TIME_STEP = 1.0; ///< Default time step (fs).

  // --- Status Messages ---
  static constexpr const char *MSG_RUNNING_ANALYSIS = "Running Analysis...";      ///< Status: Computation in progress.
  static constexpr const char *MSG_ANALYSIS_ENDED = "Analysis ended.";              ///< Status: Successfully completed.
  static constexpr const char *MSG_SELECTING_OUTPUT = "Selecting output file...";   ///< UI: File picker open.
  static constexpr const char *MSG_FILE_SELECTION_CANCELLED = "File selection cancelled."; ///< UI: User closed picker.
  static constexpr const char *MSG_ERROR_LOADING = "Error loading file: ";          ///< Error: IO or parsing failure.
  static constexpr const char *MSG_FILES_WRITTEN = "Files Written.";                ///< Success: Data exported.
  static constexpr const char *MSG_SAVE_CANCELLED = "Save cancelled.";              ///< UI: User aborted save.
  static constexpr const char *MSG_ANALYSIS_ABORTED = "Analysis aborted: No trajectory loaded."; ///< Error: Missing data.
  static constexpr const char *MSG_ERROR_ANALYSIS = "Error during analysis: ";       ///< Error: Computation failure.
  static constexpr const char *MSG_ERROR_WRITING = "Error during file writing: ";    ///< Error: Export failure.
  static constexpr const char *MSG_NO_DATA_TO_WRITE = "No analysis data to write.";  ///< Error: Empty results.
};

/**
 * @brief Encapsulates all configurable options for the application.
 */
struct ProgramOptions {
  std::string input_file;       ///< Path to the input trajectory file.
  std::string output_file_base; ///< Base path/name for output files.
  bool smoothing = true;        ///< Whether to apply Gaussian smoothing to results.
  bool use_hdf5 = true;         ///< Enable HDF5 output format.
  bool use_csv = false;         ///< Enable CSV output format.
  bool use_parquet = false;     ///< Enable Parquet output format.
  double r_max = AppDefaults::R_MAX;             ///< Max distance for RDF calculation.
  double r_bin_width = AppDefaults::R_BIN_WIDTH;   ///< Step size for RDF histogram.
  double q_max = AppDefaults::Q_MAX;             ///< Max momentum transfer for S(Q).
  double q_bin_width = AppDefaults::Q_BIN_WIDTH;   ///< Step size for S(Q) histogram.
  double r_int_max = AppDefaults::R_INT_MAX;     ///< Upper limit for g(r) integration.
  double angle_bin_width = AppDefaults::ANGLE_BIN_WIDTH;     ///< Step size for ADF.
  double dihedral_bin_width = AppDefaults::ANGLE_BIN_WIDTH;  ///< Step size for dihedral analysis.
  size_t max_ring_size = 8;     ///< Maximum ring size for topological analysis.

  /** @brief Map of calculator ID to its enabled state. */
  std::map<std::string, bool> active_calculators;

  double smoothing_sigma = AppDefaults::SMOOTHING_SIGMA; ///< Sigma for Gaussian kernel.
  correlation::math::KernelType smoothing_kernel = AppDefaults::SMOOTHING_KERNEL; ///< Smoothing kernel type.
  int min_frame = 0;            ///< Starting frame index.
  int max_frame = -1;           ///< Ending frame index (-1 for all).
  double time_step = AppDefaults::TIME_STEP; ///< Simulation time step in fs.

  /** @brief Bond cutoffs for S(Q) calculations. */
  std::vector<std::vector<double>> bond_cutoffs_sq;
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
   * @brief Updates the enabled state of a single calculator.
   * @param id Calculator ID (e.g., "RDF", "SQ").
   * @param enabled Whether this calculator should run.
   */
  void setCalculatorActive(const std::string &id, bool enabled) {
    options_.active_calculators[id] = enabled;
  }

  /**
   * @brief Gets a pointer to the current cell (first frame of trajectory).
   * @return Pointer to the correlation::core::Cell, or nullptr if no trajectory
   * loaded.
   */
  [[nodiscard]] const correlation::core::Cell *cell() const {
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
   *
   * @return A string containing an error message if the analysis failed, or an
   * empty string if successful.
   */
  [[nodiscard]] std::string run_analysis();

  /**
   * @brief Writes the analysis results to files (CSV, HDF5) as specified in
   * options.
   *
   * @return A string containing an error message if the writing failed, or an
   * empty string if successful.
   */
  [[nodiscard]] std::string write_files();

  /**
   * @brief Gets the atom counts for the current structure.
   * @return A map of correlation::core::Element Symbol -> Count.
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
   * @brief Calculates a recommended time step in fs based on the smallest
   * atomic mass. Uses formula: sqrt(9 * Minimal_mass / 5)
   * @return Recommended time step in femtoseconds.
   */
  [[nodiscard]] double getRecommendedTimeStep() const;

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
  [[nodiscard]] double getBondCutoff(int type1, int type2) const;

  /**
   * @brief Sets the bond cutoffs to be used in analysis.
   * @param cutoffs Matrix of cutoffs.
   */
  void setBondCutoffs(const std::vector<std::vector<double>> &cutoffs);

  /**
   * @brief Returns the names of all histograms available from the last
   * analysis.
   * @return Sorted vector of histogram names (e.g., "g(r)", "S(Q)", "PAD").
   *         Returns an empty vector if no analysis has been run yet.
   */
  [[nodiscard]] std::vector<std::string> getAvailableHistogramNames() const;

  /**
   * @brief Returns a pointer to a specific histogram from the last analysis.
   * @param name The histogram name (e.g., "g(r)").
   * @return Pointer to the Histogram, or nullptr if not found or no analysis
   * run.
   */
  [[nodiscard]] const correlation::analysis::Histogram *
  getHistogram(const std::string &name) const;

  // Callbacks
  /**
   * @brief Sets a callback function for analysis progress updates.
   * @param cb Callback: void(float progress, const std::string &message).
   */
  void setProgressCallback(std::function<void(float, const std::string &)> cb) {
    progress_callback_ = cb;
  }

private:
  /**
   * @brief Main function executed by the analysis thread.
   *
   * Responsible for orchestrating the computation of distribution functions
   * and notifying the provided progress callback upon completion or error.
   * This is typically run via std::async or std::thread.
   */
  void analysis_thread_func();

  // --- Private Data Members ---
  std::unique_ptr<correlation::core::Trajectory> trajectory_;           ///< Loaded trajectory data.
  std::unique_ptr<correlation::analysis::TrajectoryAnalyzer>
      trajectory_analyzer_;                                              ///< Analysis engine instance.
  std::unique_ptr<correlation::analysis::DistributionFunctions> df_;     ///< Combined calculation results.

  ProgramOptions options_;                                               ///< Active configuration.
  std::function<void(float, const std::string &)> progress_callback_;    ///< Progress notification hook.
};
} // namespace correlation::app
