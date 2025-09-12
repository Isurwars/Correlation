// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "../include/AppBackend.hpp"

#include <iostream>

#include "../include/FileIO.hpp"
#include "../include/FileWriter.hpp"

AppBackend::AppBackend() { ProgramOptions options_; }

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

void AppBackend::run_analysis() {
  if (!cell_) {
    return;
  }

  try {
    // Create the DistributionFunctions object
    df_ = std::make_unique<DistributionFunctions>(*cell_, options_.r_cut,
                                                  options_.bond_factor);

    // --- Run calculations sequentially and report progress ---
    df_->calculateRDF(options_.r_cut, options_.r_bin_width, options_.normalize);
    df_->calculatePAD(180.0, options_.angle_bin_width);
    df_->calculateSQ(20.0, 0.02);
    df_->calculateCoordinationNumber();

    // --- Write results ---
    FileWriter writer(*df_);
    writer.writeAllCSVs(options_.output_file_base, options_.smoothing);

  } catch (const std::exception &e) {
    std::cerr << "Error during analysis: " << e.what() << std::endl;
  }
}
