/**
 * @file AppBackend.cpp
 * @brief Implementation of the application backend.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "app/AppBackend.hpp"
#include "physics/PhysicalData.hpp"
#include "readers/FileReader.hpp"
#include "writers/FileWriter.hpp"

#include <cmath>
#include <iostream>
#include <limits>

namespace correlation::app {
//---------------------------------------------------------------------------//
//------------------------------- Constructors ------------------------------//
//---------------------------------------------------------------------------//

AppBackend::AppBackend() {}

//---------------------------------------------------------------------------//
//--------------------------------- Accessors -------------------------------//
//---------------------------------------------------------------------------//

std::map<std::string, int> AppBackend::getAtomCounts() const {
  std::map<std::string, int> counts;
  const correlation::core::Cell *c = cell();
  if (!c)
    return counts;

  const auto &elements = c->elements();
  std::vector<int> id_counts(elements.size(), 0);

  for (const auto &atom : c->atoms()) {
    id_counts[atom.element_id()]++;
  }

  for (size_t i = 0; i < elements.size(); ++i) {
    if (id_counts[i] > 0) {
      counts[elements[i].symbol] = id_counts[i];
    }
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

double AppBackend::getRecommendedTimeStep() const {
  const correlation::core::Cell *c = cell();
  if (!c || c->elements().empty()) {
    return AppDefaults::TIME_STEP;
  }

  double min_mass = std::numeric_limits<double>::max();
  bool found = false;

  for (const auto &element : c->elements()) {
    try {
      double mass = correlation::physics::getAtomicMass(element.symbol);
      if (mass < min_mass) {
        min_mass = mass;
        found = true;
      }
    } catch (const std::out_of_range &) {
      // Ignore unknown elements for this calculation
    }
  }

  if (found && min_mass > 0.0) {
    return std::sqrt(9.0 * min_mass / 5.0);
  }

  return AppDefaults::TIME_STEP;
}

std::vector<std::vector<double>> AppBackend::getRecommendedBondCutoffs() const {
  if (!trajectory_ || trajectory_->getFrames().empty())
    return {};

  // Precompute on the trajectory (which updates its internal cache)
  trajectory_->precomputeBondCutoffs();

  const correlation::core::Cell &c = trajectory_->getFrames()[0];
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

double AppBackend::getBondCutoff(int type1, int type2) const {
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

std::vector<std::string> AppBackend::getAvailableHistogramNames() const {
  if (!df_)
    return {};
  return df_->getAvailableHistograms();
}

const correlation::analysis::Histogram *
AppBackend::getHistogram(const std::string &name) const {
  if (!df_)
    return nullptr;
  try {
    return &df_->getHistogram(name);
  } catch (const std::out_of_range &) {
    return nullptr;
  }
}

//---------------------------------------------------------------------------//
//---------------------------------- Methods --------------------------------//
//---------------------------------------------------------------------------//

std::string AppBackend::load_file(const std::string &path) {
  std::string display_path = path;
  std::replace(display_path.begin(), display_path.end(), '\\', '/');
  correlation::readers::FileType type =
      correlation::readers::determineFileType(path);

  // For now, loading a single structure file starts a new trajectory with 1
  // frame. The determineFileType helper is used to dispatch to the correct
  // reader.

  if (type == correlation::readers::FileType::Arc ||
      type == correlation::readers::FileType::CastepMd ||
      type == correlation::readers::FileType::Outmol ||
      type == correlation::readers::FileType::Xdatcar) {
    trajectory_ = std::make_unique<correlation::core::Trajectory>(
        correlation::readers::readTrajectory(path, type, progress_callback_));
  } else {
    trajectory_ = std::make_unique<correlation::core::Trajectory>();
    trajectory_->addFrame(
        correlation::readers::readStructure(path, type, progress_callback_));
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

std::string AppBackend::run_analysis() {
  if (!trajectory_ || trajectory_->getFrames().empty()) {
    std::string err = AppDefaults::MSG_ANALYSIS_ABORTED;
    std::cerr << err << std::endl;
    return err;
  }

  try {
    // Apply custom bond cutoffs if they were set in options
    if (!options_.bond_cutoffs_sq.empty()) {
      trajectory_->setBondCutoffsSQ(options_.bond_cutoffs_sq);
    } else {
      if (trajectory_->getBondCutoffsSQ().empty()) {
        trajectory_->precomputeBondCutoffs();
      }
    }
    // Determine which cutoffs to use: explicit overrides or precomputed
    // defaults.
    const auto &active_cutoffs = !options_.bond_cutoffs_sq.empty()
                                     ? options_.bond_cutoffs_sq
                                     : trajectory_->getBondCutoffsSQ();

    // Ensure min_frame is within bounds
    size_t start_f = options_.min_frame;
    if (start_f >= trajectory_->getFrames().size())
      start_f = 0; // Default to 0 if out of bounds

    trajectory_->setTimeStep(options_.time_step);

    if (progress_callback_)
      progress_callback_(0.0f, "Starting analysis...");

    // Define progress callbacks
    // Since TrajectoryAnalyzer no longer does heavy lifting during
    // initialization, we give the entire progress duration to the
    // DistributionFunctions::computeMean phase.
    auto cb_structure = [this](float p, const std::string &msg) {
      if (progress_callback_)
        progress_callback_(0.0f,
                           msg); // just show message, progress is negligible
    };

    auto cb_dist = [this](float p, const std::string &msg) {
      if (progress_callback_)
        progress_callback_(p, msg);
    };

    // Initialize the TrajectoryAnalyzer, which handles frame-by-frame
    // structural analysis
    trajectory_analyzer_ =
        std::make_unique<correlation::analysis::TrajectoryAnalyzer>(
            *trajectory_, options_.r_max, active_cutoffs, start_f,
            options_.max_frame, true, cb_structure);

    // Prepare settings
    correlation::analysis::AnalysisSettings settings;
    settings.r_max = options_.r_max;
    settings.r_bin_width = options_.r_bin_width;
    settings.q_max = options_.q_max;
    settings.q_bin_width = options_.q_bin_width;
    settings.r_int_max = options_.r_int_max;
    settings.angle_bin_width = options_.angle_bin_width;
    settings.dihedral_bin_width = options_.dihedral_bin_width;
    settings.max_ring_size = options_.max_ring_size;
    settings.active_calculators = options_.active_calculators;
    settings.smoothing = options_.smoothing;
    settings.smoothing_sigma = options_.smoothing_sigma;
    settings.smoothing_kernel = options_.smoothing_kernel;

    // Run parallel analysis to compute distribution functions
    // This accumulates results from all processed frames.
    df_ = correlation::analysis::DistributionFunctions::computeMean(
        *trajectory_, *trajectory_analyzer_, start_f, settings, cb_dist);

    if (df_) {
      // Trajectory-wide calculators (VACF, VDOS)
      bool run_vacf = settings.isActive("VACF");
      bool run_vdos = settings.isActive("VDOS");
      if (run_vdos) {
        run_vacf = true; // VDOS requires VACF
      }
      if (run_vacf) {
        if (trajectory_->getVelocities().empty()) {
          trajectory_->calculateVelocities();
        }
        df_->calculateVACF(*trajectory_, -1, start_f, options_.max_frame);
        if (run_vdos) {
          try {
            df_->calculateVDOS();
          } catch (const std::exception &e) {
            std::cerr << "VDOS calculation failed: " << e.what() << std::endl;
          }
        }
      }
    }

  } catch (const std::exception &e) {
    std::string err = std::string(AppDefaults::MSG_ERROR_ANALYSIS) + e.what();
    std::cerr << "Analysis Exception: " << e.what() << std::endl;
    return err;
  } catch (...) {
    std::string err =
        std::string(AppDefaults::MSG_ERROR_ANALYSIS) + "Unknown error.";
    std::cerr << "Analysis Exception: Unknown error." << std::endl;
    return err;
  }
  return "";
}

std::string AppBackend::write_files() {
  if (!df_) {
    std::string err = AppDefaults::MSG_NO_DATA_TO_WRITE;
    std::cerr << err << std::endl;
    return err;
  }

  try {
    // --- Write results ---
    correlation::writers::FileWriter writer(*df_);
    writer.write(options_.output_file_base, options_.use_csv, options_.use_hdf5,
                 options_.use_parquet, options_.smoothing);
    std::cout << "Files writen to: " << options_.output_file_base << std::endl;
  } catch (const std::exception &e) {
    std::string err = std::string(AppDefaults::MSG_ERROR_WRITING) + e.what();
    std::cerr << err << std::endl;
    return err;
  }
  return "";
}

void AppBackend::analysis_thread_func() {
  // Implementation of analysis thread if needed, currently inline in
  // run_analysis or managed by AppController
}
} // namespace correlation::app
