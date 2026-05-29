/**
 * @file FileWriter.cpp
 * @brief Implementation of the unified file writer.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "writers/FileWriter.hpp"
#include "writers/WriterFactory.hpp"

#include <fstream>

namespace correlation::writers {

FileWriter::FileWriter(const correlation::analysis::DistributionFunctions &df) : df_(df) {}

void FileWriter::write(const std::string &base_path, bool use_csv, bool use_hdf5, bool use_parquet,
                       bool smoothing) const {
  auto &factory = WriterFactory::instance();

  if (use_csv) {
    if (auto writer = factory.getWriter("CSV")) {
      writer->write(base_path, df_, smoothing);
    }
  }

  if (use_hdf5) {
    if (auto writer = factory.getWriter("HDF5")) {
      writer->write(base_path, df_, smoothing);
    }
  }

  if (use_parquet) {
    if (auto writer = factory.getWriter("Parquet")) {
      writer->write(base_path, df_, smoothing);
    }
  }

  // Write the human-readable summary file
  std::string summary_filename = base_path + "_summary.txt";
  std::ofstream summary_file(summary_filename);
  if (summary_file.is_open()) {
    summary_file << "Correlation Analysis Summary\n";
    summary_file << "============================\n\n";
    summary_file << "Calculated Dynamic Properties:\n";
    summary_file << "------------------------------\n";

    double d_msd = df_.getDiffusionCoefficientMSD();
    double d_vacf = df_.getDiffusionCoefficientVACF();
    double tau = df_.getRelaxationTime();
    double de = df_.getDeborahNumber();

    summary_file << "Self-diffusion coefficient (from MSD): ";
    if (d_msd > 0.0) {
      summary_file << d_msd << " Å²/fs\n";
    } else {
      summary_file << "N/A\n";
    }

    summary_file << "Self-diffusion coefficient (from VACF): ";
    if (d_vacf > 0.0) {
      summary_file << d_vacf << " Å²/fs\n";
    } else {
      summary_file << "N/A\n";
    }

    summary_file << "Relaxation time (from VACF): ";
    if (tau > 0.0) {
      summary_file << tau << " fs\n";
    } else {
      summary_file << "N/A\n";
    }

    summary_file << "Deborah number: ";
    if (de > 0.0) {
      summary_file << de << "\n";
    } else {
      summary_file << "N/A\n";
    }

    summary_file.close();
  }
}

} // namespace correlation::writers
