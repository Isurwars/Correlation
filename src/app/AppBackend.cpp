/**
 * @file AppBackend.cpp
 * @brief Implementation of the application backend.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/AppBackend.hpp"
#include "analysis/CorrelationEngine.hpp"
#include "math/Precision.hpp"
#include "physics/PhysicalData.hpp"
#include "readers/FileReader.hpp"
#include "readers/ReaderFactory.hpp"
#include "writers/FileWriter.hpp"

#include <algorithm>
#include <filesystem>

#include <cmath>
#include <iostream>
#include <limits>

namespace correlation::app {

AppBackend::AppBackend() = default;

std::map<std::string, int> AppBackend::getAtomCounts() const {
  std::map<std::string, int> counts;
  const correlation::core::Cell *current_cell = cell();
  if (current_cell == nullptr) {
    return counts;
  }

  const auto &elements = current_cell->elements();
  std::vector<int> id_counts(elements.size(), 0);

  for (const auto &atom : current_cell->atoms()) {
    id_counts[atom.element_id()]++;
  }

  for (size_t i = 0; i < elements.size(); ++i) {
    if (id_counts[i] > 0) {
      counts[elements[i].symbol] = id_counts[i];
    }
  }

  return counts;
}

size_t AppBackend::getFrameCount() const {
  if (!trajectory_) {
    return 0;
  }
  return trajectory_->getFrameCount();
}

size_t AppBackend::getTotalAtomCount() const {
  if (!trajectory_ || trajectory_->getFrameCount() == 0) {
    return 0;
  }
  return trajectory_->firstFrame().atomCount();
}

size_t AppBackend::getRemovedFrameCount() const {
  if (!trajectory_) {
    return 0;
  }
  return trajectory_->getRemovedFrameCount();
}

real_t AppBackend::getTimeStep() const {
  if (!trajectory_) {
    return 1.0;
  }
  return trajectory_->getTimeStep();
}

real_t AppBackend::getRecommendedTimeStep() const {
  const correlation::core::Cell *current_cell = cell();
  if ((current_cell == nullptr) || current_cell->elements().empty()) {
    return AppDefaults::TIME_STEP;
  }

  real_t min_mass = std::numeric_limits<real_t>::max();
  bool found = false;

  for (const auto &element : current_cell->elements()) {
    try {
      real_t const mass = correlation::physics::getAtomicMass(element.symbol);
      if (mass < min_mass) {
        min_mass = mass;
        found = true;
      }
    } catch (const std::out_of_range &err) {
      std::cerr << "Warning: Unknown element symbol '" << element.symbol
                << "' ignored in mass calculation: " << err.what() << '\n';
    }
  }

  if (found && min_mass > 0.0) {
    return static_cast<real_t>(std::sqrt(9.0 * min_mass / 5.0));
  }

  return AppDefaults::TIME_STEP;
}

correlation::analysis::BondCutoffMatrix AppBackend::getRecommendedBondCutoffs() const {
  if (!trajectory_ || trajectory_->getFrameCount() == 0) {
    return {};
  }

  // Precompute on the trajectory (which updates its internal cache)
  trajectory_->precomputeBondCutoffs();
  return trajectory_->getBondCutoffs();
}

real_t AppBackend::getBondCutoff(size_t type1, size_t type2) const {
  if (!trajectory_) {
    return 0.0;
  }
  return trajectory_->getBondCutoff(type1, type2);
}

real_t AppBackend::getMinBondCutoff(size_t type1, size_t type2) const {
  if (!trajectory_) {
    return 0.0;
  }
  return trajectory_->getMinBondCutoff(type1, type2);
}

void AppBackend::setBondCutoffs(const correlation::analysis::BondCutoffMatrix &cutoffs) {
  if (!trajectory_) {
    return;
  }
  trajectory_->setBondCutoffs(cutoffs);
}

std::vector<std::string> AppBackend::getAvailableHistogramNames() const {
  if (!df_) {
    return {};
  }
  return df_->getAvailableHistograms();
}

const correlation::analysis::Histogram *AppBackend::getHistogram(const std::string &name) const {
  if (!df_) {
    return nullptr;
  }
  try {
    return &df_->getHistogram(name);
  } catch (const std::out_of_range &) {
    return nullptr;
  }
}

std::map<std::string, real_t> AppBackend::getAshcroftWeights() const {
  return df_ ? df_->getAshcroftWeights() : std::map<std::string, real_t>{};
}

std::string AppBackend::load_file(const std::string &path) {
  std::string display_path = path;
  std::ranges::replace(display_path, '\\', '/');
  correlation::readers::FileType const type = correlation::readers::determineFileType(path);

  // Determine whether to load as a trajectory by checking the reader's
  // isTrajectory() flag via the ReaderFactory, rather than maintaining a
  // separate hard-coded list of trajectory FileTypes.
  bool is_trajectory = false;
  std::string const ext = std::filesystem::path(path).extension().string();
  if (!ext.empty()) {
    auto *reader =
        correlation::readers::ReaderFactory::instance().getReaderForExtension({.extension = ext, .filename = path});
    if (reader != nullptr) {
      is_trajectory = reader->isTrajectory();
    }
  } else {
    // Extensionless files (e.g., bare POSCAR, XDATCAR): use the FileType
    // already determined by basename in determineFileType().
    is_trajectory = (type == correlation::readers::FileType::Xdatcar);
  }

  if (is_trajectory) {
    trajectory_ = std::make_unique<correlation::core::Trajectory>(
        correlation::readers::readTrajectory(path, type, progress_callback_));
  } else {
    trajectory_ = std::make_unique<correlation::core::Trajectory>();
    trajectory_->addFrame(correlation::readers::readStructure(path, type, progress_callback_));
  }

  options_.input_file = path;
  options_.output_file_base = path;

  std::string msg = "File loaded: " + display_path;
  return msg;
}

namespace {

[[nodiscard]] correlation::analysis::CorrelationEngineConfig toEngineConfig(const ProgramOptions &options,
                                                                            std::atomic<bool> *cancel_flag) {
  correlation::analysis::CorrelationEngineConfig config;
  config.settings.r_max = options.r_max;
  config.settings.r_bin_width = options.r_bin_width;
  config.settings.q_max = options.q_max;
  config.settings.q_bin_width = options.q_bin_width;
  config.settings.r_int_max = options.r_int_max;
  config.settings.angle_bin_width = options.angle_bin_width;
  config.settings.dihedral_bin_width = options.dihedral_bin_width;
  config.settings.max_ring_size = options.max_ring_size;
  config.settings.active_calculators = options.active_calculators;
  config.settings.smoothing = options.smoothing;
  config.settings.smoothing_sigma = options.smoothing_sigma;
  config.settings.smoothing_kernel = options.smoothing_kernel;
  config.settings.lef_cutoff = options.lef_cutoff;
  config.settings.lef_sigma = options.lef_sigma;
  config.settings.hyperuniformity_samples = options.hyper_samples;
  config.settings.cancel_flag = cancel_flag;

  config.bond_cutoffs = options.bond_cutoffs;
  config.min_frame = options.min_frame;
  config.max_frame = options.max_frame;
  config.time_step = options.time_step;
  return config;
}

} // namespace

std::string AppBackend::validateOptions() const {
  const auto config = toEngineConfig(options_, nullptr);
  return correlation::analysis::CorrelationEngine::validateConfig(config);
}

std::expected<void, std::string> AppBackend::run_analysis() {
  if (!trajectory_ || trajectory_->getFrameCount() == 0) {
    std::string const err = AppDefaults::MSG_ANALYSIS_ABORTED;
    std::cerr << err << '\n';
    return std::unexpected(err);
  }

  cancel_flag_ = false;

  std::string const validation_error = validateOptions();
  if (!validation_error.empty()) {
    return std::unexpected(validation_error);
  }

  const auto config = toEngineConfig(options_, &cancel_flag_);
  auto result = correlation::analysis::CorrelationEngine::runAnalysis(*trajectory_, config, progress_callback_);
  if (!result) {
    std::cerr << "Analysis Exception: " << result.error() << '\n';
    return std::unexpected(result.error());
  }

  df_ = std::move(result.value());
  return {};
}

std::expected<void, std::string> AppBackend::write_files() {
  if (!df_) {
    std::string const err = AppDefaults::MSG_NO_DATA_TO_WRITE;
    std::cerr << err << '\n';
    return std::unexpected(err);
  }

  try {
    // --- Write results ---
    correlation::writers::FileWriter const writer(*df_);
    writer.write(options_.output_file_base, options_.use_csv, options_.use_hdf5, options_.use_parquet,
                 options_.smoothing);
    std::cout << "Files written to: " << options_.output_file_base << '\n';
  } catch (const std::exception &e) {
    std::string const err = std::string(AppDefaults::MSG_ERROR_WRITING) + e.what();
    std::cerr << err << '\n';
    return std::unexpected(err);
  }
  return {};
}

} // namespace correlation::app
