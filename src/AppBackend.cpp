// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "../include/AppBackend.hpp"

#include <iostream>

#include "../include/FileIO.hpp"
#include "../include/FileWriter.hpp"

//---------------------------------------------------------------------------//
//------------------------------- Constructors ------------------------------//
//---------------------------------------------------------------------------//

AppBackend::AppBackend() {}

//---------------------------------------------------------------------------//
//--------------------------------- Methods ---------------------------------//
//---------------------------------------------------------------------------//

std::string AppBackend::load_file(const std::string &path) {

  FileIO::FileType type = FileIO::determineFileType(path);
  cell_ = std::make_unique<Cell>(FileIO::readStructure(path, type));
  options_.input_file = path;
  options_.output_file_base = path;

  return "Loaded " + std::to_string(cell_->atomCount()) + " atoms from:\n" +
         path;
}

void AppBackend::run_analysis(const ProgramOptions &options) {
  if (!cell_) {
    return;
  }

  try {
    // Create the DistributionFunctions object
    df_ = std::make_unique<DistributionFunctions>(*cell_, options.r_max,
                                                  options.bond_factor);

    // --- Run calculations sequentially and report progress ---
    df_->calculateCoordinationNumber();
    df_->calculateRDF(options.r_max, options.r_bin_width);
    df_->calculatePAD(options.angle_max, options.angle_bin_width);
    df_->calculateSQ(options.q_max, options.q_bin_width, options.r_int_max);
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
