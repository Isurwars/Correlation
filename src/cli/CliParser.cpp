/**
 * @file CliParser.cpp
 * @brief Implementation of the command-line argument parser.
 *
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "cli/CliParser.hpp"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>

namespace correlation::cli {

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

  if (!opts.show_help && !opts.show_version) {
    if (opts.r_max <= 0.0) {
      std::cerr << "Error: --r-max must be strictly positive.\n";
      return false;
    }
    if (opts.r_bin_width <= 0.0) {
      std::cerr << "Error: --r-bin must be strictly positive.\n";
      return false;
    }
    if (opts.q_max <= 0.0) {
      std::cerr << "Error: --q-max must be strictly positive.\n";
      return false;
    }
    if (opts.q_bin_width <= 0.0) {
      std::cerr << "Error: --q-bin must be strictly positive.\n";
      return false;
    }
    if (opts.angle_bin_width <= 0.0) {
      std::cerr << "Error: --angle-bin must be strictly positive.\n";
      return false;
    }
    if (opts.dihedral_bin_width <= 0.0) {
      std::cerr << "Error: --dihedral-bin must be strictly positive.\n";
      return false;
    }
    if (opts.time_step <= 0.0) {
      std::cerr << "Error: --time-step must be strictly positive.\n";
      return false;
    }
    if (opts.r_int_max <= 0.0) {
      std::cerr << "Error: --r-int-max must be strictly positive.\n";
      return false;
    }
    if (opts.max_ring_size <= 0) {
      std::cerr << "Error: --max-ring-size must be strictly positive.\n";
      return false;
    }
    if (opts.smoothing_sigma < 0.0) {
      std::cerr << "Error: --smoothing-sigma cannot be negative.\n";
      return false;
    }
  }

  // Default output base: same directory and stem as input
  if (opts.output_base.empty() && !opts.input_file.empty()) {
    std::filesystem::path p(opts.input_file);
    opts.output_base = (p.parent_path() / p.stem()).string();
  }

  return true;
}

} // namespace correlation::cli
