/**
 * @file FileWriter.cpp
 * @brief Implementation of the unified file writer.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "writers/FileWriter.hpp"
#include "writers/ExportBundleWriter.hpp"
#include "writers/WriterFactory.hpp"

#include <fstream>

namespace correlation::writers {

FileWriter::FileWriter(const correlation::analysis::DistributionFunctions &dists) : df_(&dists) {}

void FileWriter::write(const std::string &base_path, bool use_csv, bool use_hdf5, bool use_parquet,
                       bool smoothing) const {
  auto &factory = WriterFactory::instance();

  if (use_csv && df_ != nullptr) {
    if (auto *writer = factory.getWriter("CSV")) {
      writer->write(base_path, *df_, smoothing);
    }
  }

  if (use_hdf5 && df_ != nullptr) {
    if (auto *writer = factory.getWriter("HDF5")) {
      writer->write(base_path, *df_, smoothing);
    }
  }

  if (use_parquet && df_ != nullptr) {
    if (auto *writer = factory.getWriter("Parquet")) {
      writer->write(base_path, *df_, smoothing);
    }
  }

  writeSummaryFile(base_path);
}

std::expected<void, std::string>
FileWriter::writeBundle(const std::string &zip_path, bool use_csv, bool use_hdf5,
                        bool use_parquet, bool include_svg, bool smoothing) const {
  if (df_ == nullptr) {
    return std::unexpected("No distribution functions data available for export");
  }

  ExportBundleOptions options;
  options.include_csv = use_csv;
  options.include_hdf5 = use_hdf5;
  options.include_parquet = use_parquet;
  options.include_svg = include_svg;
  options.smoothing = smoothing;

  return ExportBundleWriter::writeBundle(zip_path, *df_, options);
}

void FileWriter::writeSummaryFile(const std::string &base_path) const {
  if (df_ == nullptr) {
    return;
  }
  const std::string summary_filename = base_path + "_summary.txt";
  std::ofstream summary_file(summary_filename);
  if (summary_file.is_open()) {
    summary_file << "Correlation Analysis Summary\n";
    summary_file << "============================\n\n";
    summary_file << "Calculated Dynamic Properties:\n";
    summary_file << "------------------------------\n";

    const real_t d_msd = df_->getDiffusionCoefficientMSD();
    const real_t d_vacf = df_->getDiffusionCoefficientVACF();
    const real_t tau = df_->getRelaxationTime();
    const real_t deb = df_->getDeborahNumber();

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
    if (deb > 0.0) {
      summary_file << deb << "\n";
    } else {
      summary_file << "N/A\n";
    }

    summary_file.close();
  }
}

} // namespace correlation::writers
