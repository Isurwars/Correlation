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

AppBackend::AppBackend() { ProgramOptions options_; }

//---------------------------------------------------------------------------//
//--------------------------------- Methods ---------------------------------//
//---------------------------------------------------------------------------//

std::string AppBackend::load_file(const std::string &path) {
  std::string display_path = path;
  std::replace(display_path.begin(), display_path.end(), '\\', '/');
  FileIO::FileType type = FileIO::determineFileType(path);
  cell_ = std::make_unique<Cell>(FileIO::readStructure(path, type));
  options_.input_file = path;
  options_.output_file_base = path;
  return "Loaded " + std::to_string(cell_->atomCount()) + " atoms from:\n" +
         display_path;
}

void AppBackend::run_analysis() {
  if (!cell_) {
    return;
  }

  try {
    // Create the DistributionFunctions object
    df_ = std::make_unique<DistributionFunctions>(*cell_, options_.r_max,
                                                  options_.bond_factor);

    // --- Run calculations sequentially and report progress ---
    df_->calculateCoordinationNumber();
    df_->calculateRDF(options_.r_max, options_.r_bin_width);
    df_->calculatePAD(options_.angle_max, options_.angle_bin_width);
    df_->calculateSQ(options_.q_max, options_.q_bin_width, options_.r_int_max);
    if (options_.smoothing) {
      df_->smoothAll(options_.smoothing_sigma, options_.smoothing_kernel);
    }
    // --- Write results ---
    FileWriter writer(*df_);
    writer.writeAllCSVs(options_.output_file_base, options_.smoothing);

  } catch (const std::exception &e) {
    std::cerr << "Error during analysis: " << e.what() << std::endl;
  }
}
