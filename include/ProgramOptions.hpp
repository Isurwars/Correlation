#ifndef INCLUDE_PROGRAMOPTIONS_HPP_
#define INCLUDE_PROGRAMOPTIONS_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <string>

#include "Smoothing.hpp"

// Encapsulates all command-line configurable options for the application.
// This avoids global variables and makes the configuration explicit.
struct ProgramOptions {
  std::string input_file;
  std::string output_file_base;
  bool normalize = true;
  bool smoothing = false;
  double r_cut = 20.0;
  double r_bin_width = 0.02;
  double q_cut = 180.0;
  double q_bin_width = 0.157;
  double angle_bin_width = 1.0;
  double smoothing_sigma = 0.081;
  double bond_factor = 1.2;
  KernelType smoothing_kernel = KernelType::Gaussian;

  // Parses command-line arguments and populates the struct members.
  // Throws std::invalid_argument on parsing errors.
  static ProgramOptions parse(int argc, char **argv);

  // Prints the help message and exits the program.
  static void printHelp();
};

#endif // INCLUDE_PROGRAMOPTIONS_HPP_
