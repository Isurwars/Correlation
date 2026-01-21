// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
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

std::map<std::string, int> AppBackend::getAtomCounts() const {
  std::map<std::string, int> counts;
  if (!cell_)
    return counts;
  for (const auto &atom : cell_->atoms()) {
    counts[atom.element().symbol]++;
  }
  return counts;
}

std::vector<std::vector<double>> AppBackend::getRecommendedBondCutoffs() const {
  if (!cell_)
    return {};
  cell_->precomputeBondCutoffs();
  // We need the raw squared cutoffs or the actual distances?
  // The UI wants distances.
  const size_t num_elements = cell_->elements().size();
  std::vector<std::vector<double>> cutoffs(num_elements,
                                           std::vector<double>(num_elements));
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = 0; j < num_elements; ++j) {
      cutoffs[i][j] = cell_->getBondCutoff(i, j);
    }
  }
  return cutoffs;
}

double AppBackend::getBondCutoff(int type1, int type2) {
  if (!cell_)
    return 0.0;
  return cell_->getBondCutoff(type1, type2);
}

void AppBackend::setBondCutoffs(const std::vector<std::vector<double>> &cutoffs) {
  if (!cell_)
    return;
  // Convert distances back to squared cutoffs for the Cell
  const size_t num_elements = cell_->elements().size();
  std::vector<std::vector<double>> cutoffs_sq(
      num_elements, std::vector<double>(num_elements));
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = 0; j < num_elements; ++j) {
      cutoffs_sq[i][j] = cutoffs[i][j] * cutoffs[i][j];
    }
  }
  cell_->setBondCutoffs(cutoffs_sq);
}

void AppBackend::run_analysis() {
  if (!cell_) {
    return;
  }

  try {
    // Apply custom bond cutoffs if they were set in options
    if (!options_.bond_cutoffs_sq_.empty()) {
      cell_->setBondCutoffs(options_.bond_cutoffs_sq_);
    }

    // Create the DistributionFunctions object
    df_ = std::make_unique<DistributionFunctions>(*cell_, options_.r_max);

    // --- Run calculations sequentially and report progress ---
    df_->calculateCoordinationNumber();
    df_->calculateRDF(options_.r_max, options_.r_bin_width);
    df_->calculatePAD(180.0, options_.angle_bin_width);
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
