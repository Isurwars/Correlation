/**
 * @file CliParser.cpp
 * @brief Implementation of the command-line argument parser.
 *
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "cli/CliParser.hpp"

#include <CLI/CLI.hpp>
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>

namespace correlation::cli {

void printUsage(const char *program) {
  std::cerr << "Correlation — Structural Analysis Tool (CLI Mode)\n"
            << "Usage: " << program << " <input_file> [options]\n\n"
            << "General Options:\n"
            << "  -o, --output <path>       Output base path (default: input stem)\n"
            << "  -q, --quiet               Suppress progress output\n"
            << "  --material <str>          Material type: amorphous, liquid, crystalline (default: amorphous)\n"
            << "  -v, --version             Show version info\n"
            << "  -h, --help                Show this help message\n\n"
            << "Calculator Selection:\n"
            << "  -g, --disable-groups <list> Comma-separated groups to disable: radial, scattering,\n"
            << "                            structural, spatial, angular, dynamic, rings,\n"
            << "                            topology (default: none)\n\n"
            << "Simulation / Frame Range:\n"
            << "  --min-frame <int>         Start frame, 1-based (default: 1)\n"
            << "  --max-frame <int>         End frame, -1=all (default: -1)\n"
            << "  --time-step <float>       Simulation time step in fs (default: 1.0)\n\n"
            << "Radial Parameters:\n"
            << "  --r-max <float>           Max radius for RDF (default: 20.0)\n"
            << "  --r-bin <float>           RDF bin width (default: 0.02)\n"
            << "  --r-int-max <float>       Max radius for g(r) integration (default: 10.0)\n\n"
            << "Scattering Parameters:\n"
            << "  --q-max <float>           Max q for S(Q) (default: 20.0)\n"
            << "  --q-bin <float>           S(Q) bin width (default: 0.02)\n\n"
            << "Angular Parameters:\n"
            << "  --angle-bin <float>       Angular bin width (default: 1.0)\n"
            << "  --dihedral-bin <float>    Dihedral bin width (default: copy angle-bin)\n\n"
            << "Ring Parameters:\n"
            << "  --max-ring-size <int>     Max ring size for topology (default: 8)\n\n"
            << "Post-Processing & Output Formats:\n"
            << "  --smoothing-sigma <float> Bandwidth for kernel smoothing (default: 0.1)\n"
            << "  --smoothing-kernel <str>  Kernel type (gaussian, bump, triweight, \n"
            << "                            epanechnikov, cosine, biweight) (default: gaussian)\n"
            << "  --no-smoothing            Disable post-processing smoothing\n"
            << "  --csv                     Enable CSV output (default: on)\n"
            << "  --no-csv                  Disable CSV output\n"
            << "  --hdf5                    Enable HDF5 output\n"
            << "  --no-hdf5                 Disable HDF5 output (default: off)\n"
            << "  --parquet                 Enable Parquet output\n"
            << "  --no-parquet              Disable Parquet output\n\n"
            << "Local Entropy Parameters:\n"
            << "  --lef-cutoff <float>      Cutoff radius for local entropy (default: 5.0)\n"
            << "  --lef-sigma <float>       Gaussian standard deviation for local entropy (default: 0.2)\n";
}

namespace {

void applyMaterialDefaults(const CLI::App &app, CliOptions &opts) {
  if (opts.material_type == 2) { // Crystalline
    if (app.count("--r-bin") == 0U) {
      opts.r_bin_width = 0.002;
    }
    if (app.count("--q-bin") == 0U) {
      opts.q_bin_width = 0.002;
    }
    if (app.count("--angle-bin") == 0U) {
      opts.angle_bin_width = 0.1;
    }
    if ((app.count("--dihedral-bin") == 0U) && (app.count("--angle-bin") == 0U)) {
      opts.dihedral_bin_width = 0.1;
    }
    if (app.count("--smoothing-sigma") == 0U) {
      opts.smoothing_sigma = 0.01;
    }
  } else if (opts.material_type == 1) { // Liquid
    if (app.count("--r-bin") == 0U) {
      opts.r_bin_width = 0.05;
    }
    if (app.count("--q-bin") == 0U) {
      opts.q_bin_width = 0.05;
    }
    if (app.count("--angle-bin") == 0U) {
      opts.angle_bin_width = 2.0;
    }
    if ((app.count("--dihedral-bin") == 0U) && (app.count("--angle-bin") == 0U)) {
      opts.dihedral_bin_width = 2.0;
    }
    if (app.count("--smoothing-sigma") == 0U) {
      opts.smoothing_sigma = 0.15;
    }
  }
}

bool validateOptions(const CliOptions &opts) {
  if (opts.r_max <= 0.0) {
    std::cerr << "Error: --r-max must be strictly positive.\n";
    return false;
  }
  if (opts.r_bin_width <= 0.0) {
    std::cerr << "Error: --r-bin must be strictly positive.\n";
    return false;
  }
  if (opts.r_bin_width >= opts.r_max) {
    std::cerr << "Error: --r-bin must be strictly less than --r-max.\n";
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
  if (opts.q_bin_width >= opts.q_max) {
    std::cerr << "Error: --q-bin must be strictly less than --q-max.\n";
    return false;
  }
  if (opts.angle_bin_width <= 0.0) {
    std::cerr << "Error: --angle-bin must be strictly positive.\n";
    return false;
  }
  if (opts.angle_bin_width > 180.0) {
    std::cerr << "Error: --angle-bin must be at most 180.0 degrees.\n";
    return false;
  }
  if (opts.dihedral_bin_width <= 0.0) {
    std::cerr << "Error: --dihedral-bin must be strictly positive.\n";
    return false;
  }
  if (opts.dihedral_bin_width > 360.0) {
    std::cerr << "Error: --dihedral-bin must be at most 360.0 degrees.\n";
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
  if (opts.max_frame < -1) {
    std::cerr << "Error: --max-frame cannot be less than -1.\n";
    return false;
  }
  if (opts.max_frame >= 0 && opts.min_frame > opts.max_frame) {
    std::cerr << "Error: --min-frame cannot be greater than --max-frame.\n";
    return false;
  }
  if (opts.smoothing && opts.smoothing_sigma <= 0.0) {
    std::cerr << "Error: --smoothing-sigma must be strictly positive when smoothing is enabled.\n";
    return false;
  }
  if (opts.smoothing_sigma < 0.0) {
    std::cerr << "Error: --smoothing-sigma cannot be negative.\n";
    return false;
  }
  if (opts.lef_cutoff <= 0.0) {
    std::cerr << "Error: --lef-cutoff must be strictly positive.\n";
    return false;
  }
  if (opts.lef_sigma <= 0.0) {
    std::cerr << "Error: --lef-sigma must be strictly positive.\n";
    return false;
  }
  if (opts.hyper_samples <= 0) {
    std::cerr << "Error: --hyper-samples must be strictly positive.\n";
    return false;
  }
  if (!opts.csv && !opts.hdf5 && !opts.parquet) {
    std::cerr << "Error: At least one output format (--csv, --hdf5, or --parquet) must be enabled.\n";
    return false;
  }
  return true;
}

bool parseArgs(std::span<char *> argv, CliOptions &opts) {
  CLI::App app{"Correlation — Structural Analysis Tool (CLI Mode)"};

  // Disable Windows-style options (like /opt) so that absolute paths starting
  // with / (e.g. /data/experiment/sample.poscar) are not parsed as options on Windows.
  app.allow_windows_style_options(false);

  // Configure help and version flags in CLI11
  app.set_help_flag("-h,--help");
  app.set_version_flag("-v,--version");

  // Parse positional input file
  app.add_option("input_file", opts.input_file, "Input structure/trajectory file");

  // Options mapping
  app.add_option("-o,--output", opts.output_base);
  app.add_flag("-q,--quiet", opts.quiet);

  app.add_option("-g,--disable-groups", opts.disable_groups);

  int min_frame_temp = 1;
  app.add_option("--min-frame", min_frame_temp);
  app.add_option("--max-frame", opts.max_frame);
  app.add_option("--time-step", opts.time_step);

  app.add_option("--r-max", opts.r_max);
  app.add_option("--r-bin", opts.r_bin_width);
  app.add_option("--r-int-max", opts.r_int_max);

  app.add_option("--q-max", opts.q_max);
  app.add_option("--q-bin", opts.q_bin_width);

  app.add_option("--angle-bin", opts.angle_bin_width);
  app.add_option("--dihedral-bin", opts.dihedral_bin_width);

  app.add_option("--max-ring-size", opts.max_ring_size);

  app.add_option("--smoothing-sigma", opts.smoothing_sigma);

  app.add_option("--lef-cutoff", opts.lef_cutoff);
  app.add_option("--lef-sigma", opts.lef_sigma);

  app.add_option("--hyper-samples", opts.hyper_samples);

  std::string k_str = "gaussian";
  app.add_option("--smoothing-kernel", k_str);

  std::map<std::string, int> const mat_map = {{"amorphous", 0}, {"liquid", 1}, {"crystalline", 2}, {"crystal", 2}};
  app.add_option("--material", opts.material_type, "Material type (amorphous, liquid, crystalline)")
      ->transform(CLI::CheckedTransformer(mat_map, CLI::ignore_case));

  app.add_flag("--csv,!--no-csv", opts.csv);
  app.add_flag("--hdf5,!--no-hdf5", opts.hdf5);
  app.add_flag("--parquet,!--no-parquet", opts.parquet);
  app.add_flag("--smoothing,!--no-smoothing", opts.smoothing);

  try {
    app.parse(static_cast<int>(argv.size()), argv.data());
  } catch (const CLI::CallForHelp &e) {
    printUsage(argv[0]);
    opts.show_help = true;
    return true;
  } catch (const CLI::CallForVersion &e) {
    opts.show_version = true;
    return true;
  } catch (const CLI::ParseError &e) {
    printUsage(argv[0]);
    return false;
  }

  // Check if dihedral-bin was explicitly set
  if (app.count("--dihedral-bin") != 0U) {
    opts.has_dihedral_bin = true;
  }

  applyMaterialDefaults(app, opts);

  // Handle min-frame conversion
  if (app.count("--min-frame") != 0U) {
    opts.min_frame = min_frame_temp - 1;
    opts.min_frame = std::max(opts.min_frame, 0);
  }

  // Process smoothing kernel string (case insensitive, fallback to gaussian)
  std::ranges::transform(k_str, k_str.begin(), [](unsigned char character) { return std::tolower(character); });
  if (k_str == "gaussian" || k_str == "gauss") {
    opts.smoothing_kernel = correlation::math::KernelType::Gaussian;
  } else if (k_str == "bump") {
    opts.smoothing_kernel = correlation::math::KernelType::Bump;
  } else if (k_str == "triweight") {
    opts.smoothing_kernel = correlation::math::KernelType::Triweight;
  } else if (k_str == "epanechnikov" || k_str == "epan") {
    opts.smoothing_kernel = correlation::math::KernelType::Epanechnikov;
  } else if (k_str == "cosine" || k_str == "cos") {
    opts.smoothing_kernel = correlation::math::KernelType::Cosine;
  } else if (k_str == "biweight") {
    opts.smoothing_kernel = correlation::math::KernelType::Biweight;
  } else {
    std::cerr << "Warning: Unknown kernel type '" << k_str << "', defaulting to gaussian\n";
    opts.smoothing_kernel = correlation::math::KernelType::Gaussian;
  }

  // Missing input file check
  if (opts.input_file.empty()) {
    std::cerr << "Error: no input file specified.\n";
    printUsage(argv[0]);
    return false;
  }

  if (!validateOptions(opts)) {
    return false;
  }

  // Default output base: same directory and stem as input
  if (opts.output_base.empty()) {
    std::filesystem::path const input_path(opts.input_file);
    opts.output_base = (input_path.parent_path() / input_path.stem()).lexically_normal().string();
  }

  return true;
}
} // namespace

} // namespace correlation::cli
