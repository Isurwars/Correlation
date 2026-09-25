/**
 * @file AppOptions.hpp
 * @brief Configuration structures and defaults for application analysis.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/AnalysisTypes.hpp"
#include "analysis/DistributionFunctions.hpp"
#include "math/Precision.hpp"
#include "math/Smoothing.hpp"

#include <map>
#include <string>

namespace correlation::app {

/**
 * @brief Default values and messages for the application.
 */
struct AppDefaults {
  static constexpr real_t R_MAX = 20.0;       ///< Default max radius for RDF (Angstrom).
  static constexpr real_t R_BIN_WIDTH = 0.02; ///< Default bin width for RDF (Angstrom).
  static constexpr real_t Q_MAX = 20.0;       ///< Default max q for S(Q) (Angstrom^-1).
  static constexpr real_t Q_BIN_WIDTH = 0.02; ///< Default bin width for S(Q) (Angstrom^-1).
  static constexpr real_t R_INT_MAX = 10.0;   ///< Default max radius for integration (Angstrom).
  static constexpr real_t ANGLE_BIN_WIDTH = 0.25; ///< Default bin width for ADF (Degrees).
  static constexpr real_t SMOOTHING_SIGMA = 0.1;  ///< Default Gaussian smoothing sigma.
  static constexpr real_t LEF_CUTOFF = 5.0;       ///< Default cutoff for local entropy.
  static constexpr real_t LEF_SIGMA = 0.2;        ///< Default Gaussian sigma for local entropy.

  // Crystalline Defaults (1 order of magnitude smaller)
  static constexpr real_t R_BIN_WIDTH_CRYSTAL = 0.002;
  static constexpr real_t Q_BIN_WIDTH_CRYSTAL = 0.002;
  static constexpr real_t ANGLE_BIN_WIDTH_CRYSTAL = 0.1;
  static constexpr real_t SMOOTHING_SIGMA_CRYSTAL = 0.01;

  // Liquid Defaults (broad/diffuse features)
  static constexpr real_t R_BIN_WIDTH_LIQUID = 0.05;
  static constexpr real_t Q_BIN_WIDTH_LIQUID = 0.05;
  static constexpr real_t ANGLE_BIN_WIDTH_LIQUID = 0.5;
  static constexpr real_t SMOOTHING_SIGMA_LIQUID = 0.15;

  /** @brief Default smoothing kernel. */
  static constexpr decltype(correlation::math::KernelType::Gaussian) SMOOTHING_KERNEL =
      correlation::math::KernelType::Gaussian;

  static constexpr real_t TIME_STEP = 1.0; ///< Default time step (fs).

  // XRD Defaults
  static constexpr real_t XRD_LAMBDA = 1.5406;   ///< Default X-ray wavelength in Å (Cu K-alpha).
  static constexpr real_t XRD_THETA_MIN = 10.0;  ///< Default min 2-theta in degrees.
  static constexpr real_t XRD_THETA_MAX = 140.0; ///< Default max 2-theta in degrees.
  static constexpr real_t XRD_BIN_WIDTH = 0.05;  ///< Default 2-theta bin width in degrees.

  // Bond Cutoff Defaults
  static constexpr real_t BOND_MIN_FACTOR =
      0.6; ///< Default factor for minimum bond cutoff distance.
  static constexpr real_t BOND_MAX_FACTOR =
      1.2; ///< Default factor for maximum bond cutoff distance.
  static constexpr real_t BOND_GLOBAL_CUTOFF = 3.5; ///< Default global uniform bond cutoff in Å.

  // --- Status Messages ---
  static constexpr const char *MSG_RUNNING_ANALYSIS =
      "Running Analysis..."; ///< Status: Computation in progress.
  static constexpr const char *MSG_ANALYSIS_ENDED =
      "Analysis ended."; ///< Status: Successfully completed.
  static constexpr const char *MSG_SELECTING_OUTPUT =
      "Selecting output file..."; ///< UI: File picker open.
  static constexpr const char *MSG_FILE_SELECTION_CANCELLED =
      "File selection cancelled."; ///< UI: User closed picker.
  static constexpr const char *MSG_ERROR_LOADING =
      "Error loading file: "; ///< Error: IO or parsing failure.
  static constexpr const char *MSG_FILES_WRITTEN = "Files Written.";   ///< Success: Data exported.
  static constexpr const char *MSG_SAVE_CANCELLED = "Save cancelled."; ///< UI: User aborted save.
  static constexpr const char *MSG_ANALYSIS_ABORTED =
      "Analysis aborted: No trajectory loaded."; ///< Error: Missing data.
  static constexpr const char *MSG_ERROR_ANALYSIS =
      "Error during analysis: "; ///< Error: Computation failure.
  static constexpr const char *MSG_ERROR_WRITING =
      "Error during file writing: "; ///< Error: Export failure.
  static constexpr const char *MSG_NO_DATA_TO_WRITE =
      "No analysis data to write."; ///< Error: Empty results.
};

/**
 * @brief Encapsulates all configurable options for the application.
 */
struct ProgramOptions {
  std::string input_file;            ///< Path to the input trajectory file.
  std::string output_file_base;      ///< Base path/name for output files.
  bool smoothing = true;             ///< Whether to apply Gaussian smoothing to results.
  bool use_hdf5 = false;             ///< Enable HDF5 output format.
  bool use_csv = true;               ///< Enable CSV output format.
  bool use_parquet = false;          ///< Enable Parquet output format.
  real_t r_max = AppDefaults::R_MAX; ///< Max distance for RDF calculation.
  real_t r_bin_width = AppDefaults::R_BIN_WIDTH;            ///< Step size for RDF histogram.
  real_t q_max = AppDefaults::Q_MAX;                        ///< Max momentum transfer for S(Q).
  real_t q_bin_width = AppDefaults::Q_BIN_WIDTH;            ///< Step size for S(Q) histogram.
  real_t r_int_max = AppDefaults::R_INT_MAX;                ///< Upper limit for g(r) integration.
  real_t angle_bin_width = AppDefaults::ANGLE_BIN_WIDTH;    ///< Step size for ADF.
  real_t dihedral_bin_width = AppDefaults::ANGLE_BIN_WIDTH; ///< Step size for dihedral analysis.
  size_t max_ring_size = 8; ///< Maximum ring size for topological analysis.

  /** @brief Map of calculator ID to its enabled state. */
  std::map<std::string, bool> active_calculators;

  real_t smoothing_sigma = AppDefaults::SMOOTHING_SIGMA; ///< Sigma for Gaussian kernel.
  real_t lef_cutoff = AppDefaults::LEF_CUTOFF;           ///< Cutoff radius for local entropy.
  real_t lef_sigma = AppDefaults::LEF_SIGMA; ///< Gaussian standard deviation for local entropy.
  size_t hyper_samples = 10000;              ///< Number of random samples for hyperuniformity.
  correlation::math::KernelType smoothing_kernel =
      AppDefaults::SMOOTHING_KERNEL;         ///< Smoothing kernel type.
  int min_frame = 0;                         ///< Starting frame index.
  int max_frame = -1;                        ///< Ending frame index (-1 for all).
  int frame_stride = 1;                      ///< Stride between analyzed frames (>= 1).
  real_t time_step = AppDefaults::TIME_STEP; ///< Simulation time step in fs.

  int material_type = 0; ///< Material type (0: Amorphous, 1: Liquid, 2: Crystalline).

  /** @brief Parameters for X-Ray Diffraction calculation. */
  correlation::analysis::XRDParams xrd_params{
      .lambda = AppDefaults::XRD_LAMBDA,
      .theta_min = AppDefaults::XRD_THETA_MIN,
      .theta_max = AppDefaults::XRD_THETA_MAX,
      .bin_width = AppDefaults::XRD_BIN_WIDTH,
  };

  /** @brief Bond cutoff ranges for neighbor & topological calculations. */
  correlation::analysis::BondCutoffMatrix bond_cutoffs;
};

} // namespace correlation::app
