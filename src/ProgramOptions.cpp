// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "../include/ProgramOptions.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>

void ProgramOptions::printHelp() {
  // Help message content remains the same...
  const std::string help_txt =
      "CORRELATION\n"
      "DESCRIPTION:\n"
      "  This program calculates the main correlation functions of a "
      "material:\n"
      "    - Radial Distribution Function (J(r)).\n"
      "    - Pair Distribution Function (g(r)).\n"
      "    - Reduced Pair Distribution Function (G(r)).\n"
      "    - Coordination Number (CN).\n"
      "    - Plane-Angle Distribution (PAD).\n\n"
      "    - Structure Factor (S(Q)).\n"
      "USAGE: correlation [OPTIONS] [input_file]\n"
      "  The minimal argument is a structure file, this program requires a "
      "file\n"
      "  that contains atom positions, crystal structure and composition.\n"
      "  Supported structure files are:\n"
      "    -*.CAR   Materials Studio structure file.\n"
      "    -*.CELL  CASTEP structure file.\n"
      "    -*.dat   ONETEP structure file.\n\n"
      "OPTIONS:\n"
      "    HELP OPTIONS:"
      "      -h, --help\n"
      "        Display this help text.\n\n"
      "    RADIAL OPTIONS:\n"
      "      -n, --normalize\n"
      "        Used to switch between weighted partials (default), or "
      "normalize\n"
      "        all the partials to 1 when r tends to infinity\n"
      "      -R, --r_cut\n"
      "        Cutoff radius in the calculation of g(r), G(r) and J(r). The "
      "default\n"
      "        radius it's set to 2 nm. The maximum recommended radius is the "
      "same as\n"
      "        shortest length of the periodic boundary conditions (PBC), "
      "anything\n"
      "        above this PBC value can be affected by periodic "
      "interactions.\n\n"
      "      -r, --r_bin_width\n"
      "        Width of the histograms for g(r) and J(r), the default is 0.05 "
      "nm.\n\n"
      "    BOND-ANGLE OPTIONS:\n"
      "      -a, --angle_bin_width\n"
      "        Width of the histograms for the PAD, default set to 1.0°.\n\n"
      "      -b, --bond_parameter\n"
      "        The ideal covalent bond length is the sum of covalent radii\n"
      "        of the two atoms. The criterion used to consider atoms as "
      "bonded\n"
      "        is the following:\n"
      "            0.6 * Sum_radii < distance < bond_parameter * Sum_radii.\n"
      "        By default the bond_parameter is set to 1.30, as a rule of "
      "thumb.\n"
      "        The default should work for most crystalline materials,\n"
      "        as well as most covalent non-crystalline materials.\n"
      "        For amorphous and liquid materials the bond_parameter should "
      "be\n"
      "        increased to match the desired distance to cut_off the bonds. \n"
      "        A bond_parameter of 1.42 is recomended for amorphous "
      "materials.\n\n"
      "        If any of the pairs is missing in the input file, the "
      "corresponding\n"
      "        bond distance will be set using the bond_parameter(1.30 by "
      "default).\n\n"
      "    OUTPUT OPTIONS:\n"
      "      -o, --out_file\n"
      "        The output file name, by default the input seed name will be "
      "used.\n\n"
      "    SMOOTHING OPTIONS:\n"
      "      -k, --kernel\n"
      "        Smoothing kernel selector (Default: 1):\n"
      "          1: Gaussian kernel.\n"
      "          2: Bump Function kernel.\n"
      "          3: Triweight kernel.\n"
      "      -K, --kernel_sigma\n"
      "        Width of the smoothing kernel, by default 0.081.\n"
      "        The recommended values for sigma are:\n"
      "        Gaussian Kernel: 0.081\n"
      "        Bump Kernel: 0.19\n"
      "      -S, --smoothing\n"
      "        Smoothing is disabled by default, this option enable "
      "smoothing\n\n"
      "CREATOR:\n"
      "  This program was created by PhD. Isaias Rodriguez Aguirre, November "
      "2020.\n"
      "  e-mail: isurwars@gmail.com\n"
      "ACKNOWLEDGMENTS:\n"
      "  This software was created during a Posdoctoral fellowship in "
      "IIM-UNAM.\n"
      "  Thanks to DGAPA-UNAM for the financial support during my fellowship.\n"
      "  And their continious support as part of the Workgroup projects with "
      "IDs:\n"
      "  IN104617 and IN116520.\n";
  std::cout << help_txt;
  exit(0);
}

ProgramOptions ProgramOptions::parse(int argc, char **argv) {
  ProgramOptions opts;
  std::vector<std::string> args(argv + 1, argv + argc);

  for (size_t i = 0; i < args.size(); ++i) {
    const std::string &arg = args[i];
    if (arg == "-h" || arg == "--help") {
      printHelp();
    } else if ((arg == "-o" || arg == "--output") && i + 1 < args.size()) {
      opts.output_file_base = args[++i];
    } else if ((arg == "-R" || arg == "--r_cut") && i + 1 < args.size()) {
      opts.r_cut = std::stod(args[++i]);
    } else if ((arg == "-r" || arg == "--r_bin") && i + 1 < args.size()) {
      opts.r_bin_width = std::stod(args[++i]);
    } else if ((arg == "-a" || arg == "--angle_bin") && i + 1 < args.size()) {
      opts.angle_bin_width = std::stod(args[++i]);
    } else if ((arg == "-b" || arg == "--bond_parameter") &&
               i + 1 < args.size()) {
      opts.bond_factor = std::stod(args[++i]);
    } else if (arg == "-S" || arg == "--smoothing") {
      opts.smoothing = true;
    } else if ((arg == "-K" || arg == "--sigma") && i + 1 < args.size()) {
      opts.smoothing_sigma = std::stod(args[++i]);
    } else if ((arg == "-k" || arg == "--kernel") && i + 1 < args.size()) {
      int k_type = std::stoi(args[++i]);
      if (k_type == 2)
        opts.smoothing_kernel = KernelType::Bump;
      else if (k_type == 3)
        opts.smoothing_kernel = KernelType::Triweight;
      else
        opts.smoothing_kernel = KernelType::Gaussian;
    } else if (arg == "-n" || arg == "--normalize") {
      opts.normalize = true;
    } else if (arg.rfind("-", 0) != 0) { // Doesn't start with '-'
      opts.input_file = arg;
    } else {
      throw std::invalid_argument("Unknown or incomplete option: " + arg);
    }
  }

  if (opts.input_file.empty()) {
    throw std::runtime_error("No input file specified.");
  }

  if (opts.output_file_base.empty()) {
    opts.output_file_base =
        std::filesystem::path(opts.input_file).stem().string();
  }

  return opts;
}
