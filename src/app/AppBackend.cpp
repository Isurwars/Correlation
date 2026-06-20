/**
 * @file AppBackend.cpp
 * @brief Implementation of the application backend.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/AppBackend.hpp"
#include "analysis/DynamicsAnalyzer.hpp"
#include "calculators/CalculatorFactory.hpp"
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
//---------------------------------------------------------------------------//
//------------------------------- Constructors ------------------------------//
//---------------------------------------------------------------------------//

AppBackend::AppBackend() = default;

//---------------------------------------------------------------------------//
//--------------------------------- Accessors -------------------------------//
//---------------------------------------------------------------------------//

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

double AppBackend::getTimeStep() const {
  if (!trajectory_) {
    return 1.0;
  }
  return trajectory_->getTimeStep();
}

double AppBackend::getRecommendedTimeStep() const {
  const correlation::core::Cell *current_cell = cell();
  if ((current_cell == nullptr) || current_cell->elements().empty()) {
    return AppDefaults::TIME_STEP;
  }

  double min_mass = std::numeric_limits<double>::max();
  bool found = false;

  for (const auto &element : current_cell->elements()) {
    try {
      double const mass = correlation::physics::getAtomicMass(element.symbol);
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
    return std::sqrt(9.0 * min_mass / 5.0);
  }

  return AppDefaults::TIME_STEP;
}

std::vector<std::vector<double>> AppBackend::getRecommendedBondCutoffs() const {
  if (!trajectory_ || trajectory_->getFrameCount() == 0) {
    return {};
  }

  // Precompute on the trajectory (which updates its internal cache)
  trajectory_->precomputeBondCutoffs();

  const correlation::core::Cell &current_cell = trajectory_->firstFrame();
  const size_t num_elements = current_cell.elements().size();
  std::vector<std::vector<double>> cutoffs(num_elements, std::vector<double>(num_elements));
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = 0; j < num_elements; ++j) {
      // getBondCutoff uses indices.
      cutoffs[i][j] = trajectory_->getBondCutoff(i, j);
    }
  }
  return cutoffs;
}

double AppBackend::getBondCutoff(size_t type1, size_t type2) const {
  if (!trajectory_) {
    return 0.0;
  }
  return trajectory_->getBondCutoff(type1, type2);
}

void AppBackend::setBondCutoffs(const std::vector<std::vector<double>> &cutoffs) {
  if (!trajectory_) {
    return;
  }

  // Calculate squared cutoffs
  if (trajectory_->getFrameCount() == 0) {
    return;
  }
  const size_t num_elements = trajectory_->firstFrame().elements().size();
  std::vector<std::vector<double>> cutoffs_sq(num_elements, std::vector<double>(num_elements));
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = 0; j < num_elements; ++j) {
      cutoffs_sq[i][j] = cutoffs[i][j] * cutoffs[i][j];
    }
  }

  trajectory_->setBondCutoffsSQ(cutoffs_sq);
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

std::map<std::string, double> AppBackend::getAshcroftWeights() const {
  return df_ ? df_->getAshcroftWeights() : std::map<std::string, double>{};
}

//---------------------------------------------------------------------------//
//---------------------------------- Methods --------------------------------//
//---------------------------------------------------------------------------//

std::string AppBackend::load_file(const std::string &path) {
  std::string display_path = path;
  std::replace(display_path.begin(), display_path.end(), '\\', '/');
  correlation::readers::FileType const type = correlation::readers::determineFileType(path);

  // Determine whether to load as a trajectory by checking the reader's
  // isTrajectory() flag via the ReaderFactory, rather than maintaining a
  // separate hard-coded list of trajectory FileTypes.
  bool is_trajectory = false;
  std::string const ext = std::filesystem::path(path).extension().string();
  if (!ext.empty()) {
    auto *reader = correlation::readers::ReaderFactory::instance().getReaderForExtension({ext, path});
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

  // Return info from the first frame
  size_t atom_count = 0;
  if (trajectory_->getFrameCount() > 0) {
    atom_count = trajectory_->firstFrame().atomCount();
  }

  std::string msg = "File loaded: " + display_path;
  return msg;
}

std::string AppBackend::validateOptions() const {
  if (options_.r_max <= 0.0) {
    return "Error: r_max must be strictly positive.";
  }
  if (options_.r_bin_width <= 0.0) {
    return "Error: r_bin_width must be strictly positive.";
  }
  if (options_.r_bin_width >= options_.r_max) {
    return "Error: r_bin_width must be strictly less than r_max.";
  }
  if (options_.q_max <= 0.0) {
    return "Error: q_max must be strictly positive.";
  }
  if (options_.q_bin_width <= 0.0) {
    return "Error: q_bin_width must be strictly positive.";
  }
  if (options_.q_bin_width >= options_.q_max) {
    return "Error: q_bin_width must be strictly less than q_max.";
  }
  if (options_.angle_bin_width <= 0.0) {
    return "Error: angle_bin_width must be strictly positive.";
  }
  if (options_.angle_bin_width > 180.0) {
    return "Error: angle_bin_width must be at most 180.0 degrees.";
  }
  if (options_.dihedral_bin_width <= 0.0) {
    return "Error: dihedral_bin_width must be strictly positive.";
  }
  if (options_.dihedral_bin_width > 360.0) {
    return "Error: dihedral_bin_width must be at most 360.0 degrees.";
  }
  if (options_.time_step <= 0.0) {
    return "Error: time_step must be strictly positive.";
  }
  if (options_.r_int_max <= 0.0) {
    return "Error: r_int_max must be strictly positive.";
  }
  if (options_.max_ring_size <= 0) {
    return "Error: max_ring_size must be strictly positive.";
  }
  if (options_.max_frame < -1) {
    return "Error: max_frame cannot be less than -1.";
  }
  if (options_.max_frame >= 0 && options_.min_frame > options_.max_frame) {
    return "Error: min_frame cannot be greater than max_frame.";
  }
  if (options_.smoothing && options_.smoothing_sigma <= 0.0) {
    return "Error: smoothing_sigma must be strictly positive when smoothing is enabled.";
  }
  if (options_.smoothing_sigma < 0.0) {
    return "Error: smoothing_sigma cannot be negative.";
  }
  return "";
}

void AppBackend::setupTrajectorySettings(size_t &start_f) {
  // Apply custom bond cutoffs if they were set in options
  if (!options_.bond_cutoffs_sq.empty()) {
    trajectory_->setBondCutoffsSQ(options_.bond_cutoffs_sq);
  } else {
    if (trajectory_->getBondCutoffsSQ().empty()) {
      trajectory_->precomputeBondCutoffs();
    }
  }

  // Ensure min_frame is within bounds
  start_f = options_.min_frame;
  if (start_f >= trajectory_->getFrameCount()) {
    start_f = 0; // Default to 0 if out of bounds
  }

  trajectory_->setTimeStep(options_.time_step);
}

void AppBackend::runTrajectoryCalculators(const correlation::analysis::AnalysisSettings &settings) {
  const auto &factory_calcs = ::correlation::calculators::CalculatorFactory::instance().getCalculators();

  // Check if we need velocities
  bool need_velocities = false;
  for (const auto &calc : factory_calcs) {
    bool const is_active = settings.isActive(calc->getName()) || settings.isActive(calc->getShortName());
    if (calc->isTrajectoryCalculator() && is_active) {
      if (calc->getName() == "VACF" || calc->getShortName() == "VACF" || calc->getName() == "vDoS" ||
          calc->getShortName() == "vDoS") {
        need_velocities = true;
        break;
      }
    }
  }

  if (need_velocities) {
    if (trajectory_->getFrameCount() == 0) {
      throw std::runtime_error("VACF requires a trajectory with frames");
    }
    trajectory_->calculateVelocities();
  }

  for (const auto &calc : factory_calcs) {
    if (!calc->isTrajectoryCalculator()) {
      continue;
    }
    if (!settings.isActive(calc->getName()) && !settings.isActive(calc->getShortName())) {
      continue;
    }

    try {
      calc->calculateTrajectory(*df_, *trajectory_, settings);
    } catch (const std::exception &e) {
      std::cerr << calc->getName() << " calculation failed: " << e.what() << '\n';
    }
  }
}

void AppBackend::calculateDynamicProperties() {
  const auto &hists = df_->getAllHistograms();

  auto it_msd = hists.find("MSD");
  if (it_msd != hists.end()) {
    const auto &hist = it_msd->second;
    auto it_total = hist.partials.find("Total");
    if (it_total != hist.partials.end()) {
      double const d_msd =
          correlation::analysis::DynamicsAnalyzer::computeDiffusionCoefficientMSD(hist.bins, it_total->second);
      df_->setDiffusionCoefficientMSD(d_msd);
    }
  }

  auto it_vacf = hists.find("VACF");
  if (it_vacf != hists.end()) {
    const auto &hist = it_vacf->second;
    auto it_total = hist.partials.find("Total");
    if (it_total != hist.partials.end()) {
      double const d_vacf =
          correlation::analysis::DynamicsAnalyzer::computeDiffusionCoefficientVACF(hist.bins, it_total->second);
      df_->setDiffusionCoefficientVACF(d_vacf);
    }
  }

  auto it_norm = hists.find("Normalized VACF");
  if (it_norm != hists.end()) {
    const auto &hist = it_norm->second;
    auto it_total = hist.partials.find("Total");
    if (it_total != hist.partials.end()) {
      double const tau = correlation::analysis::DynamicsAnalyzer::computeRelaxationTime(hist.bins, it_total->second);
      df_->setRelaxationTime(tau);
      double deborah = 0.0;
      if (!hist.bins.empty() && hist.bins.back() > 0.0) {
        deborah = tau / hist.bins.back();
      }
      df_->setDeborahNumber(deborah);
    }
  }

  // Log the calculated values to the console
  if (df_->getDiffusionCoefficientMSD() > 0.0) {
    std::cout << "Self-diffusion coefficient (from MSD): " << df_->getDiffusionCoefficientMSD() << " Å²/fs" << '\n';
  }
  if (df_->getDiffusionCoefficientVACF() > 0.0) {
    std::cout << "Self-diffusion coefficient (from VACF): " << df_->getDiffusionCoefficientVACF() << " Å²/fs" << '\n';
  }
  if (df_->getRelaxationTime() > 0.0) {
    std::cout << "Relaxation time (from VACF): " << df_->getRelaxationTime() << " fs" << '\n';
    std::cout << "Deborah number: " << df_->getDeborahNumber() << '\n';
  }
}

std::string AppBackend::run_analysis() {
  if (!trajectory_ || trajectory_->getFrameCount() == 0) {
    std::string err = AppDefaults::MSG_ANALYSIS_ABORTED;
    std::cerr << err << '\n';
    return err;
  }

  cancel_flag_ = false;

  std::string validation_error = validateOptions();
  if (!validation_error.empty()) {
    return validation_error;
  }

  try {
    size_t start_f = 0;
    setupTrajectorySettings(start_f);

    // Determine which cutoffs to use: explicit overrides or precomputed defaults.
    const auto &active_cutoffs =
        !options_.bond_cutoffs_sq.empty() ? options_.bond_cutoffs_sq : trajectory_->getBondCutoffsSQ();

    if (progress_callback_) {
      progress_callback_(0.0F, "Starting analysis...");
    }

    // Define progress callbacks
    // Since TrajectoryAnalyzer no longer does heavy lifting during
    // initialization, we give the entire progress duration to the
    // DistributionFunctions::computeMean phase.
    auto cb_structure = [this](float /*p*/, const std::string &msg) {
      if (progress_callback_) {
        progress_callback_(0.0F, msg); // just show message, progress is negligible
      }
    };

    auto cb_dist = [this](float progress, const std::string &msg) {
      if (progress_callback_) {
        progress_callback_(progress, msg);
      }
    };

    // Initialize the TrajectoryAnalyzer, which handles frame-by-frame
    // structural analysis
    trajectory_analyzer_ = std::make_unique<correlation::analysis::TrajectoryAnalyzer>(
        *trajectory_, options_.r_max, active_cutoffs, correlation::analysis::StartFrame{start_f},
        correlation::analysis::EndFrame{static_cast<size_t>(options_.max_frame)}, true, cb_structure);

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
    settings.cancel_flag = &cancel_flag_;

    // Run parallel analysis to compute distribution functions
    // This accumulates results from all processed frames.
    df_ = correlation::analysis::DistributionFunctions::computeMean(*trajectory_, *trajectory_analyzer_, start_f,
                                                                    settings, cb_dist);

    if (df_) {
      // Trajectory-wide calculators (VACF, VDOS, MSD, etc.)
      runTrajectoryCalculators(settings);

      // Check for dynamic properties and calculate/set them
      calculateDynamicProperties();
    }

  } catch (const std::exception &e) {
    std::string err = std::string(AppDefaults::MSG_ERROR_ANALYSIS) + e.what();
    std::cerr << "Analysis Exception: " << e.what() << '\n';
    return err;
  } catch (...) {
    std::string err = std::string(AppDefaults::MSG_ERROR_ANALYSIS) + "Unknown error.";
    std::cerr << "Analysis Exception: Unknown error." << '\n';
    return err;
  }
  return "";
}

std::string AppBackend::write_files() {
  if (!df_) {
    std::string err = AppDefaults::MSG_NO_DATA_TO_WRITE;
    std::cerr << err << '\n';
    return err;
  }

  try {
    // --- Write results ---
    correlation::writers::FileWriter const writer(*df_);
    writer.write(options_.output_file_base, options_.use_csv, options_.use_hdf5, options_.use_parquet,
                 options_.smoothing);
    std::cout << "Files written to: " << options_.output_file_base << '\n';
  } catch (const std::exception &e) {
    std::string err = std::string(AppDefaults::MSG_ERROR_WRITING) + e.what();
    std::cerr << err << '\n';
    return err;
  }
  return "";
}

void AppBackend::analysis_thread_func() {
  // Implementation of analysis thread if needed, currently inline in
  // run_analysis or managed by AppController
}
} // namespace correlation::app
