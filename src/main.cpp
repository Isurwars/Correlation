// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <iostream>

#include "../include/Cell.hpp"
#include "../include/DistributionFunctions.hpp"
#include "../include/FileIO.hpp"
#include "../include/FileWriter.hpp"
#include "../include/ProgramOptions.hpp"

int main(int argc, char **argv) {
  try {
    // 1. Parse command-line arguments into a dedicated options struct.
    ProgramOptions opts = ProgramOptions::parse(argc, argv);

    // 2. Read the input structure file to create a Cell object.
    std::cout << "Reading structure file: " << opts.input_file << std::endl;
    FileIO::FileType type = FileIO::determineFileType(opts.input_file);
    Cell main_cell = FileIO::readStructure(opts.input_file, type);
    std::cout << "File read successfully. Found " << main_cell.atomCount()
              << " atoms." << std::endl;

    // 3. Perform the correlation function calculations.
    // The DistributionFunctions class now manages its own StructureAnalyzer.
    DistributionFunctions df(main_cell, opts.r_cut, opts.bond_factor);

    std::cout << "Calculating Radial Distribution Functions..." << std::endl;
    df.calculateRDF(opts.r_cut, opts.r_bin_width, opts.normalize);

    std::cout << "Calculating Plane Angle Distribution..." << std::endl;
    df.calculatePAD(180.0, opts.angle_bin_width);

    // 4. Apply smoothing if requested by the user.
    if (opts.smoothing) {
      std::cout << "Applying smoothing with sigma = " << opts.smoothing_sigma
                << "..." << std::endl;
      df.smoothAll(opts.smoothing_sigma, opts.smoothing_kernel);
    }

    // 5. Write all calculated results to CSV files.
    std::cout << "Writing output files with base name: "
              << opts.output_file_base << std::endl;
    FileWriter writer(df);
    writer.writeAllCSVs(opts.output_file_base, opts.smoothing);

    std::cout << "\nAnalysis of " << opts.input_file
              << " completed successfully." << std::endl;

  } catch (const std::exception &e) {
    // A single, clean catch-block for all potential errors.
    std::cerr << "Error: " << e.what() << std::endl;
    return 1; // Return a non-zero exit code to indicate failure.
  }

  return 0; // Success.
}
