/**
 * @file CliParser.hpp
 * @brief Command-line argument parsing for correlation-cli.
 *
 * Provides the CliOptions struct and the parseArgs / printUsage helpers used by
 * the headless CLI entry point.  Extracted into its own translation unit so that
 * the parsing logic can be exercised by unit tests.
 *
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/Smoothing.hpp"

#include <span>
#include <string>

namespace correlation::cli {

/**
 * @brief All options that can be set via the command line.
 */
struct CliOptions {
  std::string input_file;
  std::string output_base;
  real_t r_max = 20.0;
  real_t r_bin_width = 0.02;
  real_t q_max = 20.0;
  real_t q_bin_width = 0.02;
  real_t angle_bin_width = 1.0;
  real_t dihedral_bin_width = 1.0;
  bool has_dihedral_bin = false;
  int min_frame = 0;
  int max_frame = -1;
  bool csv = true;
  bool hdf5 = false;
  bool parquet = false;
  bool smoothing = true;
  bool quiet = false;
  std::string disable_groups; // comma-separated list of groups to disable
  real_t time_step = 1.0;
  real_t r_int_max = 10.0;
  int max_ring_size = 8;
  real_t smoothing_sigma = 0.1;
  real_t lef_cutoff = 5.0;
  real_t lef_sigma = 0.2;
  int hyper_samples = 10000;
  correlation::math::KernelType smoothing_kernel = correlation::math::KernelType::Gaussian;
  int material_type = 0;
  bool show_version = false;
  bool show_help = false;
};

/**
 * @brief Prints the usage / help message to stderr.
 * @param program Name of the executable (argv[0]).
 */
void printUsage(const char *program);

/**
 * @brief Parses the command-line arguments into a CliOptions struct.
 * @param argv Argument values.
 * @param opts [out] Filled on success.
 * @return true if parsing succeeded, false on error (message already printed).
 */
bool parseArgs(std::span<const char *const> argv, CliOptions &opts);

/**
 * @brief C-style entry point overload forwarding to the span interface.
 */
inline bool parseArgs(int argc, const char *const *argv, CliOptions &opts) {
  return parseArgs(std::span(argv, static_cast<std::size_t>(argc)), opts);
}

} // namespace correlation::cli
