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
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "math/Smoothing.hpp"

#include <string>

namespace correlation::cli {

/**
 * @brief All options that can be set via the command line.
 */
struct CliOptions {
  std::string input_file;
  std::string output_base;
  double r_max = 20.0;
  double r_bin_width = 0.02;
  double q_max = 20.0;
  double q_bin_width = 0.02;
  double angle_bin_width = 1.0;
  double dihedral_bin_width = 1.0;
  bool has_dihedral_bin = false;
  int min_frame = 0;
  int max_frame = -1;
  bool csv = true;
  bool hdf5 = false;
  bool parquet = false;
  bool smoothing = true;
  bool quiet = false;
  std::string calculators; // comma-separated list
  double time_step = 1.0;
  double r_int_max = 10.0;
  int max_ring_size = 8;
  double smoothing_sigma = 0.1;
  correlation::math::KernelType smoothing_kernel =
      correlation::math::KernelType::Gaussian;
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
 * @param argc Argument count.
 * @param argv Argument values.
 * @param opts [out] Filled on success.
 * @return true if parsing succeeded, false on error (message already printed).
 */
bool parseArgs(int argc, char *argv[], CliOptions &opts);

} // namespace correlation::cli
