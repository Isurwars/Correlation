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

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>

namespace {

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
  correlation::math::KernelType smoothing_kernel = correlation::math::KernelType::Gaussian;
  bool show_version = false;
  bool show_help = false;
};

void printUsage(const char *program) {
  std::cerr
      << "Correlation — Structural Analysis Tool (CLI Mode)\n"
      << "Usage: " << program << " <input_file> [options]\n\n"
      << "Options:\n"
      << "  -o, --output <path>       Output base path (default: input stem)\n"
      << "  --r-max <float>           Max radius for RDF (default: 20.0)\n"
      << "  --r-bin <float>           RDF bin width (default: 0.02)\n"
      << "  --q-max <float>           Max q for S(Q) (default: 20.0)\n"
      << "  --q-bin <float>           S(Q) bin width (default: 0.02)\n"
      << "  --angle-bin <float>       Angular bin width (default: 1.0)\n"
      << "  --dihedral-bin <float>    Dihedral bin width (default: copy angle-bin)\n"
      << "  --time-step <float>       Simulation time step in fs (default: 1.0)\n"
      << "  --r-int-max <float>       Max radius for g(r) integration (default: 10.0)\n"
      << "  --max-ring-size <int>     Max ring size for topology (default: 8)\n"
      << "  --smoothing-sigma <float> Bandwidth for kernel smoothing (default: 0.1)\n"
      << "  --smoothing-kernel <str>  Kernel type (gaussian, bump, triweight) (default: gaussian)\n"
      << "  --min-frame <int>         Start frame, 1-based (default: 1)\n"
      << "  --max-frame <int>         End frame, -1=all (default: -1)\n"
      << "  --csv                     Enable CSV output (default: on)\n"
      << "  --no-csv                  Disable CSV output\n"
      << "  --hdf5                    Enable HDF5 output\n"
      << "  --no-hdf5                 Disable HDF5 output (default: off)\n"
      << "  --parquet                 Enable Parquet output\n"
      << "  --no-parquet              Disable Parquet output\n"
      << "  --no-smoothing            Disable post-processing smoothing\n"
      << "  --calculators <list>      Comma-separated calculator IDs\n"
      << "  -q, --quiet               Suppress progress output\n"
      << "  -v, --version             Show version info\n"
      << "  -h, --help                Show this help message\n";
}

bool parseArgs(int argc, char *argv[], CliOptions &opts) {
  if (argc < 2) {
    printUsage(argv[0]);
    return false;
  }

  // First positional argument is the input file
  int i = 1;
  if (argv[i][0] != '-') {
    opts.input_file = argv[i];
    i++;
  }

  for (; i < argc; ++i) {
    std::string arg = argv[i];

    if (arg == "-h" || arg == "--help") {
      printUsage(argv[0]);
      opts.show_help = true;
      return true;
    } else if (arg == "-v" || arg == "--version") {
      opts.show_version = true;
      return true;
    } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
      opts.output_base = argv[++i];
    } else if (arg == "--r-max" && i + 1 < argc) {
      opts.r_max = std::stod(argv[++i]);
    } else if (arg == "--r-bin" && i + 1 < argc) {
      opts.r_bin_width = std::stod(argv[++i]);
    } else if (arg == "--q-max" && i + 1 < argc) {
      opts.q_max = std::stod(argv[++i]);
    } else if (arg == "--q-bin" && i + 1 < argc) {
      opts.q_bin_width = std::stod(argv[++i]);
    } else if (arg == "--angle-bin" && i + 1 < argc) {
      opts.angle_bin_width = std::stod(argv[++i]);
    } else if (arg == "--dihedral-bin" && i + 1 < argc) {
      opts.dihedral_bin_width = std::stod(argv[++i]);
      opts.has_dihedral_bin = true;
    } else if (arg == "--time-step" && i + 1 < argc) {
      opts.time_step = std::stod(argv[++i]);
    } else if (arg == "--r-int-max" && i + 1 < argc) {
      opts.r_int_max = std::stod(argv[++i]);
    } else if (arg == "--max-ring-size" && i + 1 < argc) {
      opts.max_ring_size = std::stoi(argv[++i]);
    } else if (arg == "--smoothing-sigma" && i + 1 < argc) {
      opts.smoothing_sigma = std::stod(argv[++i]);
    } else if (arg == "--smoothing-kernel" && i + 1 < argc) {
      std::string k_str = argv[++i];
      std::transform(k_str.begin(), k_str.end(), k_str.begin(), ::tolower);
      if (k_str == "gaussian" || k_str == "gauss") {
        opts.smoothing_kernel = correlation::math::KernelType::Gaussian;
      } else if (k_str == "bump") {
        opts.smoothing_kernel = correlation::math::KernelType::Bump;
      } else if (k_str == "triweight") {
        opts.smoothing_kernel = correlation::math::KernelType::Triweight;
      } else {
        std::cerr << "Warning: Unknown kernel type '" << k_str << "', defaulting to gaussian\n";
        opts.smoothing_kernel = correlation::math::KernelType::Gaussian;
      }
    } else if (arg == "--min-frame" && i + 1 < argc) {
      opts.min_frame = std::stoi(argv[++i]) - 1; // Convert 1-based to 0-based
      if (opts.min_frame < 0)
        opts.min_frame = 0;
    } else if (arg == "--max-frame" && i + 1 < argc) {
      opts.max_frame = std::stoi(argv[++i]);
    } else if (arg == "--csv") {
      opts.csv = true;
    } else if (arg == "--no-csv") {
      opts.csv = false;
    } else if (arg == "--hdf5") {
      opts.hdf5 = true;
    } else if (arg == "--no-hdf5") {
      opts.hdf5 = false;
    } else if (arg == "--parquet") {
      opts.parquet = true;
    } else if (arg == "--no-parquet") {
      opts.parquet = false;
    } else if (arg == "--no-smoothing") {
      opts.smoothing = false;
    } else if (arg == "--calculators" && i + 1 < argc) {
      opts.calculators = argv[++i];
    } else if (arg == "-q" || arg == "--quiet") {
      opts.quiet = true;
    } else if (opts.input_file.empty() && arg[0] != '-') {
      opts.input_file = arg;
    } else {
      std::cerr << "Unknown option: " << arg << "\n";
      printUsage(argv[0]);
      return false;
    }
  }

  if (!opts.show_help && !opts.show_version && opts.input_file.empty()) {
    std::cerr << "Error: no input file specified.\n";
    printUsage(argv[0]);
    return false;
  }

  // Default output base: same directory and stem as input
  if (opts.output_base.empty() && !opts.input_file.empty()) {
    std::filesystem::path p(opts.input_file);
    opts.output_base = (p.parent_path() / p.stem()).string();
  }

  return true;
}

} // namespace

int main(int argc, char *argv[]) {
  CliOptions cli;
  if (!parseArgs(argc, argv, cli)) {
    return 1;
  }

  if (cli.show_help) {
    return 0;
  }

  if (cli.show_version) {
    std::cout << "Correlation version 3.0.0\n";
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

  // Parse calculator list if provided
  if (!cli.calculators.empty()) {
    // First, disable all calculators
    // Then enable only the specified ones
    std::istringstream ss(cli.calculators);
    std::string id;
    while (std::getline(ss, id, ',')) {
      // Trim whitespace
      id.erase(0, id.find_first_not_of(" \t"));
      id.erase(id.find_last_not_of(" \t") + 1);
      if (!id.empty()) {
        opts.active_calculators[id] = true;
      }
    }
  }

  // Create backend
  correlation::app::AppBackend backend;
  backend.setOptions(opts);

  // Set up progress callback
  if (!cli.quiet) {
    backend.setProgressCallback([](float p, const std::string &msg) {
      std::cerr << "\r[" << static_cast<int>(p * 100) << "%] " << msg
                << std::flush;
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
      std::cerr << "Frames: " << backend.getFrameCount()
                << "  Atoms: " << backend.getTotalAtomCount() << "\n";
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
