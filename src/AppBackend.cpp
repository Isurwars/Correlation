// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "../include/AppBackend.hpp"

#include <iostream>

#include "../include/FileIO.hpp"
#include "../include/FileWriter.hpp"

AppBackend::AppBackend() {}

void AppBackend::load_file(const std::string &path) {
  try {
    FileIO::FileType type = FileIO::determineFileType(path);
    cell_ = std::make_unique<Cell>(FileIO::readStructure(path, type));
    options_.input_file = path;
    options_.output_file_base = path;

    std::string status =
        "Loaded " + std::to_string(cell_->atomCount()) + " atoms from " + path;

  } catch (const std::exception &e) {
    std::cerr << "Error loading file: " << e.what() << std::endl;
  }
}

// Modified to accept ProgramOptions
void AppBackend::run_analysis(const ProgramOptions &options) {
  if (!cell_) {
    return;
  }

  try {
    // Create the DistributionFunctions object
    df_ = std::make_unique<DistributionFunctions>(*cell_, options.r_cut,
                                                  options.bond_factor);

    // --- Run calculations sequentially and report progress ---
    df_->calculateCoordinationNumber();
    df_->calculateRDF(options.r_cut, options.r_bin_width, options.normalize);
    df_->calculatePAD(180.0, options.angle_bin_width);
    df_->calculateSQ(20.0, 0.02);
    if (options.smoothing) {
      df_->smoothAll(options.smoothing_sigma, options.smoothing_kernel);
    }
    // --- Write results ---
    FileWriter writer(*df_);
    writer.writeAllCSVs(options.output_file_base, options.smoothing);

  } catch (const std::exception &e) {
    std::cerr << "Error during analysis: " << e.what() << std::endl;
  }
}
