/**
 * @file main_cli.cpp
 * @brief Headless CLI entry point for Correlation.
 *
 * Provides a command-line interface to the analysis engine without requiring
 * a GUI or display server. Uses AppBackend directly, bypassing Slint.
 *
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/AppBackend.hpp"
#include "cli/CliParser.hpp"

#include "calculators/CalculatorFactory.hpp"

#include <algorithm>
#include <iostream>
#include <set>
#include <span>
#include <sstream>
#include <string>

int main(int argc, char *argv[]) {
  correlation::cli::CliOptions cli;
  if (!correlation::cli::parseArgs(std::span<char *>(argv, argc), cli)) {
    return 1;
  }

  if (cli.show_help) {
    return 0;
  }

  if (cli.show_version) {
    std::cout << "Correlation version 3.5.0\n";
    return 0;
  }

  // Map CLI options to ProgramOptions
  correlation::app::ProgramOptions opts;
  opts.input_file = cli.input_file;
  opts.output_file_base = cli.output_base;
  opts.r_max = cli.r_max;
  opts.r_bin_width = cli.r_bin_width;
  opts.q_max = cli.q_max;
  opts.q_bin_width = cli.q_bin_width;
  opts.angle_bin_width = cli.angle_bin_width;
  opts.dihedral_bin_width = cli.has_dihedral_bin ? cli.dihedral_bin_width : cli.angle_bin_width;
  opts.min_frame = cli.min_frame;
  opts.max_frame = cli.max_frame;
  opts.smoothing = cli.smoothing;
  opts.use_csv = cli.csv;
  opts.use_hdf5 = cli.hdf5;
  opts.use_parquet = cli.parquet;
  opts.time_step = cli.time_step;
  opts.r_int_max = cli.r_int_max;
  opts.max_ring_size = cli.max_ring_size;
  opts.smoothing_sigma = cli.smoothing_sigma;
  opts.smoothing_kernel = cli.smoothing_kernel;
  opts.material_type = cli.material_type;

  // Helper lambdas for string cleanup
  auto trim = [](std::string str) {
    str.erase(0, str.find_first_not_of(" \t"));
    str.erase(str.find_last_not_of(" \t") + 1);
    return str;
  };

  auto to_lower = [](std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
  };

  const auto &factory_calcs = correlation::calculators::CalculatorFactory::instance().getCalculators();

  // Parse disabled groups
  std::set<std::string> disabled_groups;
  if (!cli.disable_groups.empty()) {
    std::istringstream str_stream(cli.disable_groups);
    std::string group_name;
    while (std::getline(str_stream, group_name, ',')) {
      group_name = to_lower(trim(group_name));
      if (!group_name.empty()) {
        disabled_groups.insert(group_name);
      }
    }
  }

  for (const auto &calc : factory_calcs) {
    std::string const calc_group = to_lower(calc->getGroup());
    bool const active = (!disabled_groups.contains(calc_group));
    opts.active_calculators[calc->getName()] = active;
    opts.active_calculators[calc->getShortName()] = active;
  }

  // Create backend
  correlation::app::AppBackend backend;
  backend.setOptions(opts);

  // Set up progress callback
  if (!cli.quiet) {
    backend.setProgressCallback([](float progress, const std::string &msg) {
      std::cerr << "\r[" << static_cast<int>(progress * 100) << "%] " << msg << std::flush;
    });
  }

  // Load file
  try {
    if (!cli.quiet) {
      std::cerr << "Loading: " << cli.input_file << "\n";
    }
    std::string const msg = backend.load_file(cli.input_file);

    // Re-apply options since load_file overwrites output_file_base
    opts.output_file_base = cli.output_base;
    backend.setOptions(opts);
    if (!cli.quiet) {
      std::cerr << msg << "\n";
      std::cerr << "Frames: " << backend.getFrameCount() << "  Atoms: " << backend.getTotalAtomCount() << "\n";
    }
  } catch (const std::exception &e) {
    std::cerr << "Error loading file: " << e.what() << "\n";
    return 1;
  }

  // Run analysis
  try {
    if (!cli.quiet) {
      std::cerr << "Running analysis...\n";
    }
    std::string const err = backend.run_analysis();
    if (!err.empty()) {
      std::cerr << "\nAnalysis error: " << err << "\n";
      return 1;
    }
    if (!cli.quiet) {
      std::cerr << "\nAnalysis complete.\n";
    }
  } catch (const std::exception &e) {
    std::cerr << "\nAnalysis exception: " << e.what() << "\n";
    return 1;
  }

  // Write output files
  try {
    std::string const err = backend.write_files();
    if (!err.empty()) {
      std::cerr << "Write error: " << err << "\n";
      return 1;
    }
    if (!cli.quiet) {
      std::cerr << "Results written to: " << cli.output_base << "\n";
    }
  } catch (const std::exception &e) {
    std::cerr << "Write exception: " << e.what() << "\n";
    return 1;
  }

  return 0;
}
