/**
 * @file main_cli.cpp
 * @brief Headless CLI entry point for Correlation.
 *
 * Provides a command-line interface to the analysis engine without requiring
 * a GUI or display server. Uses AppBackend directly, bypassing Slint.
 *
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "app/AppBackend.hpp"
#include "cli/CliParser.hpp"

#include "calculators/CalculatorFactory.hpp"

#include <algorithm>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

int main(int argc, char *argv[]) {
  correlation::cli::CliOptions cli;
  if (!correlation::cli::parseArgs(argc, argv, cli)) {
    return 1;
  }

  if (cli.show_help) {
    return 0;
  }

  if (cli.show_version) {
    std::cout << "Correlation version 3.1.0\n";
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

  // Helper lambdas for string cleanup
  auto trim = [](std::string s) {
    s.erase(0, s.find_first_not_of(" \t"));
    s.erase(s.find_last_not_of(" \t") + 1);
    return s;
  };
  auto lowercase = [](std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
  };

  const auto &factory_calcs = correlation::calculators::CalculatorFactory::instance().getCalculators();

  // Populate active_calculators
  if (!cli.calculators.empty()) {
    // Disable all calculators first
    for (const auto &calc : factory_calcs) {
      opts.active_calculators[calc->getName()] = false;
      opts.active_calculators[calc->getShortName()] = false;
    }
    // Enable explicitly requested ones
    std::istringstream ss(cli.calculators);
    std::string id;
    while (std::getline(ss, id, ',')) {
      id = trim(id);
      if (!id.empty()) {
        opts.active_calculators[id] = true;
      }
    }
  } else {
    // Parse groups
    std::set<std::string> enabled_groups;
    if (!cli.groups.empty()) {
      std::istringstream ss(cli.groups);
      std::string group_name;
      while (std::getline(ss, group_name, ',')) {
        group_name = lowercase(trim(group_name));
        if (!group_name.empty()) {
          enabled_groups.insert(group_name);
        }
      }
    } else {
      // Default: radial group only
      enabled_groups.insert("radial");
    }

    bool all_groups = enabled_groups.count("all") > 0;

    for (const auto &calc : factory_calcs) {
      std::string calc_group = lowercase(calc->getGroup());
      bool active = all_groups || (enabled_groups.count(calc_group) > 0);
      opts.active_calculators[calc->getName()] = active;
      opts.active_calculators[calc->getShortName()] = active;
    }
  }

  // Create backend
  correlation::app::AppBackend backend;
  backend.setOptions(opts);

  // Set up progress callback
  if (!cli.quiet) {
    backend.setProgressCallback([](float p, const std::string &msg) {
      std::cerr << "\r[" << static_cast<int>(p * 100) << "%] " << msg << std::flush;
    });
  }

  // Load file
  try {
    if (!cli.quiet) {
      std::cerr << "Loading: " << cli.input_file << "\n";
    }
    std::string msg = backend.load_file(cli.input_file);

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
    std::string err = backend.run_analysis();
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
    std::string err = backend.write_files();
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
