// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "AppBackend.hpp"

#include <iostream>

#include "FileIO.hpp"
#include "FileWriter.hpp"

//---------------------------------------------------------------------------//
//------------------------------- Constructors ------------------------------//
//---------------------------------------------------------------------------//

AppBackend::AppBackend() {}

//---------------------------------------------------------------------------//
//--------------------------------- Methods ---------------------------------//
//---------------------------------------------------------------------------//
std::string AppBackend::load_file(const std::string &path) {
  std::string display_path = path;
  std::replace(display_path.begin(), display_path.end(), '\\', '/');
  FileIO::FileType type = FileIO::determineFileType(path);

  // For now, loading a single structure file starts a new trajectory with 1 frame.
  // Ideally, FileIO::readTrajectory could handle this, but readStructure returns a Cell.
  // We can wrap it.
  
  if (type == FileIO::FileType::Arc) {
      trajectory_ = std::make_unique<Trajectory>(FileIO::readTrajectory(path, type));
  } else {
      trajectory_ = std::make_unique<Trajectory>();
      trajectory_->addFrame(FileIO::readStructure(path, type));
  }
  
  options_.input_file = path;
  options_.output_file_base = path;
  
  // Return info from the first frame
  size_t atom_count = 0;
  if (!trajectory_->getFrames().empty()) {
      atom_count = trajectory_->getFrames()[0].atomCount();
  }
  
  return "Loaded " + std::to_string(atom_count) + " atoms from:\n" +
         display_path;
}

std::map<std::string, int> AppBackend::getAtomCounts() const {
  std::map<std::string, int> counts;
  const Cell* c = cell();
  if (!c)
    return counts;
  for (const auto &atom : c->atoms()) {
    counts[atom.element().symbol]++;
  }
  return counts;
}

std::vector<std::vector<double>> AppBackend::getRecommendedBondCutoffs() const {
  if (!trajectory_ || trajectory_->getFrames().empty())
    return {};
    
  // Precompute on the trajectory (which updates its internal cache)
  trajectory_->precomputeBondCutoffs();

  const Cell& c = trajectory_->getFrames()[0];
  const size_t num_elements = c.elements().size();
  std::vector<std::vector<double>> cutoffs(num_elements,
                                           std::vector<double>(num_elements));
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = 0; j < num_elements; ++j) {
       // getBondCutoff uses indices.
       cutoffs[i][j] = trajectory_->getBondCutoff(i, j);
    }
  }
  return cutoffs;
}

double AppBackend::getBondCutoff(int type1, int type2) {
  if (!trajectory_)
      return 0.0;
  return trajectory_->getBondCutoff(type1, type2);
}

void AppBackend::setBondCutoffs(const std::vector<std::vector<double>> &cutoffs) {
  if (!trajectory_) return;
  
  // Calculate squared cutoffs
  if (trajectory_->getFrames().empty()) return;
  const size_t num_elements = trajectory_->getFrames()[0].elements().size();
  std::vector<std::vector<double>> cutoffs_sq(
      num_elements, std::vector<double>(num_elements));
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = 0; j < num_elements; ++j) {
      cutoffs_sq[i][j] = cutoffs[i][j] * cutoffs[i][j];
    }
  }
  
  trajectory_->setBondCutoffs(cutoffs_sq);
}

void AppBackend::run_analysis() {
  if (!trajectory_ || trajectory_->getFrames().empty()) {
    return;
  }

  try {
    // Apply custom bond cutoffs if they were set in options
    if (!options_.bond_cutoffs_sq_.empty()) {
       trajectory_->setBondCutoffs(options_.bond_cutoffs_sq_);
    } else {
      if (trajectory_->getBondCutoffs().empty()) {
           trajectory_->precomputeBondCutoffs();
       }
    }
    const auto& active_cutoffs = !options_.bond_cutoffs_sq_.empty() ? options_.bond_cutoffs_sq_ : trajectory_->getBondCutoffs();
    trajectory_analyzer_ = std::make_unique<TrajectoryAnalyzer>(*trajectory_, options_.r_max, active_cutoffs);
    Cell& first_frame = trajectory_->getFrames()[0];
    df_ = std::make_unique<DistributionFunctions>(first_frame, 0.0, active_cutoffs);
    if (!trajectory_analyzer_->getAnalyzers().empty()) {
        df_->setStructureAnalyzer(trajectory_analyzer_->getAnalyzers()[0].get());
    }
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
    if (options_.use_csv) {
      writer.writeAllCSVs(options_.output_file_base, options_.smoothing);
    }
    if (options_.use_hdf5) {
      writer.writeHDF(options_.output_file_base + ".h5");
    }

  } catch (const std::exception &e) {
    std::cerr << "Error during analysis: " << e.what() << std::endl;
  }
}
