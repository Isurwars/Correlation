// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "AppBackend.hpp"

#include <iostream>

#include "FileReader.hpp"
#include "FileWriter.hpp"

//---------------------------------------------------------------------------//
//------------------------------- Constructors ------------------------------//
//---------------------------------------------------------------------------//

AppBackend::AppBackend() {}

//---------------------------------------------------------------------------//
//--------------------------------- Accessors -------------------------------//
//---------------------------------------------------------------------------//

std::map<std::string, int> AppBackend::getAtomCounts() const {
  std::map<std::string, int> counts;
  const Cell *c = cell();
  if (!c)
    return counts;
  for (const auto &atom : c->atoms()) {
    counts[atom.element().symbol]++;
  }
  return counts;
}

int AppBackend::getFrameCount() const {
  if (!trajectory_)
    return 0;
  return trajectory_->getFrames().size();
}

int AppBackend::getTotalAtomCount() const {
  if (!trajectory_ || trajectory_->getFrames().empty())
    return 0;
  return trajectory_->getFrames()[0].atomCount();
}

size_t AppBackend::getRemovedFrameCount() const {
  if (!trajectory_)
    return 0;
  return trajectory_->getRemovedFrameCount();
}

double AppBackend::getTimeStep() const {
  if (!trajectory_)
    return 1.0;
  return trajectory_->getTimeStep();
}

std::vector<std::vector<double>> AppBackend::getRecommendedBondCutoffs() const {
  if (!trajectory_ || trajectory_->getFrames().empty())
    return {};

  // Precompute on the trajectory (which updates its internal cache)
  trajectory_->precomputeBondCutoffs();

  const Cell &c = trajectory_->getFrames()[0];
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

void AppBackend::setBondCutoffs(
    const std::vector<std::vector<double>> &cutoffs) {
  if (!trajectory_)
    return;

  // Calculate squared cutoffs
  if (trajectory_->getFrames().empty())
    return;
  const size_t num_elements = trajectory_->getFrames()[0].elements().size();
  std::vector<std::vector<double>> cutoffs_sq(
      num_elements, std::vector<double>(num_elements));
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = 0; j < num_elements; ++j) {
      cutoffs_sq[i][j] = cutoffs[i][j] * cutoffs[i][j];
    }
  }

  trajectory_->setBondCutoffsSQ(cutoffs_sq);
}

//---------------------------------------------------------------------------//
//---------------------------------- Methods --------------------------------//
//---------------------------------------------------------------------------//

std::string AppBackend::load_file(const std::string &path) {
  std::string display_path = path;
  std::replace(display_path.begin(), display_path.end(), '\\', '/');
  FileReader::FileType type = FileReader::determineFileType(path);

  // For now, loading a single structure file starts a new trajectory with 1
  // frame. The determineFileType helper is used to dispatch to the correct
  // reader.

  if (type == FileReader::FileType::Arc) {
    trajectory_ =
        std::make_unique<Trajectory>(FileReader::readTrajectory(path, type));
  } else {
    trajectory_ = std::make_unique<Trajectory>();
    trajectory_->addFrame(FileReader::readStructure(path, type));
  }

  options_.input_file = path;
  options_.output_file_base = path;

  // Return info from the first frame
  size_t atom_count = 0;
  if (!trajectory_->getFrames().empty()) {
    atom_count = trajectory_->getFrames()[0].atomCount();
  }

  std::string msg = "File loaded: " + display_path;
  return msg;
}

void AppBackend::run_analysis() {
  if (!trajectory_ || trajectory_->getFrames().empty()) {
    std::cerr << "Analysis aborted: No trajectory loaded." << std::endl;
    return;
  }

  try {
    // Apply custom bond cutoffs if they were set in options
    if (!options_.bond_cutoffs_sq_.empty()) {
      trajectory_->setBondCutoffsSQ(options_.bond_cutoffs_sq_);
    } else {
      if (trajectory_->getBondCutoffsSQ().empty()) {
        trajectory_->precomputeBondCutoffs();
      }
    }
    // Determine which cutoffs to use: explicit overrides or precomputed
    // defaults.
    const auto &active_cutoffs = !options_.bond_cutoffs_sq_.empty()
                                     ? options_.bond_cutoffs_sq_
                                     : trajectory_->getBondCutoffsSQ();

    // Ensure min_frame is within bounds
    size_t start_f = options_.min_frame;
    if (start_f >= trajectory_->getFrames().size())
      start_f = 0; // Default to 0 if out of bounds

    trajectory_->setTimeStep(options_.time_step);

    if (progress_callback_)
      progress_callback_(0.0f);

    // Define progress callbacks
    auto cb_structure = [this](float p) {
      if (progress_callback_)
        progress_callback_(p * 0.7f);
    };

    auto cb_dist = [this](float p) {
      if (progress_callback_)
        progress_callback_(0.7f + p * 0.3f);
    };

    // Initialize the TrajectoryAnalyzer, which handles frame-by-frame
    // structural analysis
    trajectory_analyzer_ = std::make_unique<TrajectoryAnalyzer>(
        *trajectory_, options_.r_max, active_cutoffs, start_f,
        options_.max_frame, true, cb_structure);

    // Prepare settings
    AnalysisSettings settings;
    settings.r_max = options_.r_max;
    settings.r_bin_width = options_.r_bin_width;
    settings.q_max = options_.q_max;
    settings.q_bin_width = options_.q_bin_width;
    settings.r_int_max = options_.r_int_max;
    settings.angle_bin_width = options_.angle_bin_width;
    settings.smoothing = options_.smoothing;
    settings.smoothing_sigma = options_.smoothing_sigma;
    settings.smoothing_kernel = options_.smoothing_kernel;

    // Run parallel analysis to compute distribution functions
    // This accumulates results from all processed frames.
    df_ = DistributionFunctions::computeMean(
        *trajectory_, *trajectory_analyzer_, start_f, settings, cb_dist);

    if (df_) {
      // Ensure velocities are calculated for VACF
      if (trajectory_->getVelocities().empty()) {
        trajectory_->calculateVelocities();
      }
      df_->calculateVACF(*trajectory_, -1, start_f, options_.max_frame);
      try {
        df_->calculateVDOS();
      } catch (const std::exception &e) {
        std::cerr << "VDOS calculation failed: " << e.what() << std::endl;
      }
    }

  } catch (const std::exception &e) {
    std::cerr << "Error during analysis: " << e.what() << std::endl;
  }
}

void AppBackend::write_files() {
  if (!df_) {
    std::cerr << "No analysis data to write." << std::endl;
    return;
  }

  try {
    // --- Write results ---
    FileWriter writer(*df_);
    if (options_.use_csv) {
      writer.writeAllCSVs(options_.output_file_base, options_.smoothing);
    }
    if (options_.use_hdf5) {
      writer.writeHDF(options_.output_file_base + ".h5");
    }
    std::cout << "Files writen to: " << options_.output_file_base << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "Error during file writing: " << e.what() << std::endl;
  }
}

void AppBackend::analysis_thread_func() {
  // Implementation of analysis thread if needed, currently inline in
  // run_analysis or managed by AppController
}
