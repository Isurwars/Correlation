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
  const Cell* c = cell();
  if (!c)
    return {};
    
  // We need to const_cast to call precomputeBondCutoffs because it modifies the cell.
  // Ideally this method should be const and check a flag, or mutable.
  // Or we just call it on the mutable trajectory frame if we had access.
  // Since cell() returns const Cell*, we cast.
  const_cast<Cell*>(c)->precomputeBondCutoffs();

  const size_t num_elements = c->elements().size();
  std::vector<std::vector<double>> cutoffs(num_elements,
                                           std::vector<double>(num_elements));
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = 0; j < num_elements; ++j) {
      cutoffs[i][j] = c->findElement(c->elements()[i].symbol)->id.value == i ? const_cast<Cell*>(c)->getBondCutoff(i, j) : 0.0;
    }
  }
  return cutoffs;
}

double AppBackend::getBondCutoff(int type1, int type2) {
  const Cell* c = cell();
  if (!c)
    return 0.0;
  // Again, getBondCutoff is not const in Cell.hpp (it calls precompute).
  return const_cast<Cell*>(c)->getBondCutoff(type1, type2);
}

void AppBackend::setBondCutoffs(const std::vector<std::vector<double>> &cutoffs) {
  // Apply to all frames in the trajectory? Or just the first?
  // User intention is likely global settings.
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
  
  for (auto& frame : trajectory_->getFrames()) {
      frame.setBondCutoffs(cutoffs_sq);
  }
}

void AppBackend::run_analysis() {
  if (!trajectory_ || trajectory_->getFrames().empty()) {
    return;
  }

  try {
    // Apply custom bond cutoffs if they were set in options
    if (!options_.bond_cutoffs_sq_.empty()) {
      for (auto& frame : trajectory_->getFrames()) {
         frame.setBondCutoffs(options_.bond_cutoffs_sq_);
      }
    }
    
    // 1. Instantiate TrajectoryAnalyzer
    // It will compute neighbors and angles for all frames.
    // Note: TrajectoryAnalyzer constructor takes bond_cutoffs now? 
    // It accepts them in constructor.
    trajectory_analyzer_ = std::make_unique<TrajectoryAnalyzer>(*trajectory_, options_.r_max, options_.bond_cutoffs_sq_);

    // 2. Setup DistributionFunctions for the first frame (Single Frame Analysis compatibility)
    // We pass 0.0 as cutoff to avoid re-calculation inside DF.
    Cell& first_frame = trajectory_->getFrames()[0];
    df_ = std::make_unique<DistributionFunctions>(first_frame, 0.0);
    
    // 3. Inject the pre-computed analyzer
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
