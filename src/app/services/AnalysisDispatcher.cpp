/**
 * @file AnalysisDispatcher.cpp
 * @brief Implementation of the AnalysisDispatcher service.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/services/AnalysisDispatcher.hpp"
#include "analysis/CorrelationEngine.hpp"
#include "writers/FileWriter.hpp"

#include <algorithm>
#include <array>
#include <iostream>

namespace correlation::app {

namespace {

[[nodiscard]] std::map<std::string, bool>
extractBondDependentCalculators(const std::map<std::string, bool> &active_calculators) {
  std::map<std::string, bool> filtered;
  if (active_calculators.empty()) {
    for (const auto &name : AnalysisDispatcher::getBondDependentCalculatorNames()) {
      filtered[name] = true;
    }
  } else {
    for (const auto &[calc_name, is_active] : active_calculators) {
      if (is_active && AnalysisDispatcher::isBondDependentCalculator(calc_name)) {
        filtered[calc_name] = true;
      }
    }
  }
  return filtered;
}

void mergeCalculatedHistograms(correlation::analysis::DistributionFunctions &target,
                               correlation::analysis::DistributionFunctions &source) {
  for (const auto &name : source.getAvailableHistograms()) {
    target.addHistogram(name, std::move(source.getMutableHistogram(name)));
  }
}

[[nodiscard]] correlation::analysis::CorrelationEngineConfig
toEngineConfig(const ProgramOptions &options, std::atomic<bool> *cancel_flag) {
  correlation::analysis::CorrelationEngineConfig config;
  config.settings.r_max = options.r_max;
  config.settings.r_bin_width = options.r_bin_width;
  config.settings.q_max = options.q_max;
  config.settings.q_bin_width = options.q_bin_width;
  config.settings.r_int_max = options.r_int_max;
  config.settings.sq_exact_direct = options.sq_exact_direct;
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
  config.settings.vacf_max_lag = options.vacf_max_lag;
  config.settings.vdos_window = options.vdos_window;
  config.settings.vdos_zero_pad = options.vdos_zero_pad;
  config.settings.cancel_flag = cancel_flag;
  config.settings.xrd_params = options.xrd_params;
  config.settings.mlip_model_path = options.mlip_model_path;
  config.settings.ldos_broadening_kernel = options.ldos_broadening_kernel;
  config.settings.ldos_broadening_sigma = options.ldos_broadening_sigma;
  config.settings.ldos_electronic_temp_k = options.ldos_electronic_temp_k;
  config.settings.ldos_gap_threshold = options.ldos_gap_threshold;
  config.settings.ldos_compute_idos = options.ldos_compute_idos;
  config.settings.eliashberg_mu_star = options.eliashberg_mu_star;
  config.settings.eliashberg_coupling_mode = options.eliashberg_coupling_mode;
  config.settings.eliashberg_disp_step = options.eliashberg_disp_step;

  config.bond_cutoffs = options.bond_cutoffs;
  config.min_frame = options.min_frame;
  config.max_frame = options.max_frame;
  config.time_step = options.time_step;
  config.settings.frame_stride = static_cast<size_t>(std::max(1, options.frame_stride));
  return config;
}

} // namespace

AnalysisDispatcher::~AnalysisDispatcher() { stopPerFramePrebake(); }

void AnalysisDispatcher::clear() noexcept {
  stopPerFramePrebake();
  const std::unique_lock lock(df_mutex_);
  df_.reset();
}

void AnalysisDispatcher::stopPerFramePrebake() noexcept {
  prebake_cancel_.store(true, std::memory_order_release);
  if (prebake_thread_.joinable()) {
    prebake_thread_.join();
  }
}

void AnalysisDispatcher::bindPerFrameTrajectory(
    std::shared_ptr<const correlation::core::Trajectory> trajectory) {
  if (!trajectory) {
    return;
  }

  std::shared_ptr<correlation::analysis::PerFrameDistributionStore> store;
  correlation::analysis::AnalysisSettings settings_copy;
  correlation::analysis::BondCutoffMatrix cutoffs_copy;
  {
    const std::shared_lock lock(df_mutex_);
    if (!df_) {
      return;
    }
    store = df_->mutablePerFrameStore();
    settings_copy = last_settings_;
    cutoffs_copy = last_bond_cutoffs_;
  }
  if (!store) {
    return;
  }
  if (cutoffs_copy.empty()) {
    cutoffs_copy = trajectory->getBondCutoffs();
  }

  using correlation::analysis::DistributionMetric;
  const std::array<bool, 4> wanted_metrics{
      !store->getMeanCurve(DistributionMetric::RDF).empty(),
      !store->getMeanCurve(DistributionMetric::PAD).empty(),
      !store->getMeanCurve(DistributionMetric::SQ).empty(),
      !store->getMeanCurve(DistributionMetric::XRD).empty(),
  };

  store->setFrameGenerator([traj = std::move(trajectory), settings = std::move(settings_copy),
                            cutoffs = std::move(cutoffs_copy),
                            wanted_metrics](std::size_t frame_idx) {
    if (!traj || frame_idx >= traj->getFrameCount()) {
      return correlation::analysis::PerFrameDistributionStore::FrameCurveArray{};
    }
    return correlation::analysis::DistributionFunctions::computeLightweightFrameCurves(
        traj->getFrame(frame_idx), settings, cutoffs, wanted_metrics);
  });
}

void AnalysisDispatcher::startPerFramePrebake(
    std::shared_ptr<const correlation::core::Trajectory> trajectory,
    std::function<void()> on_complete) {
  stopPerFramePrebake();
  if (!trajectory) {
    return;
  }
  bindPerFrameTrajectory(std::move(trajectory));

  std::shared_ptr<correlation::analysis::PerFrameDistributionStore> store;
  {
    const std::shared_lock lock(df_mutex_);
    if (df_) {
      store = df_->mutablePerFrameStore();
    }
  }
  if (!store) {
    return;
  }

  prebake_cancel_.store(false, std::memory_order_release);
  const std::size_t total_frames = store->getTotalFrames();
  prebake_thread_ = std::jthread(
      [this, store = std::move(store), total_frames, on_complete = std::move(on_complete)] {
        for (std::size_t frame_idx = 0; frame_idx < total_frames; ++frame_idx) {
          if (prebake_cancel_.load(std::memory_order_acquire)) {
            return;
          }
          if (store->needsFrameGeneration(frame_idx)) {
            store->ensureFrame(frame_idx);
          }
        }
        if (!prebake_cancel_.load(std::memory_order_acquire) && on_complete) {
          on_complete();
        }
      });
}

std::expected<void, std::string>
AnalysisDispatcher::validateOptions(const ProgramOptions &options) {
  const auto config = toEngineConfig(options, nullptr);
  return correlation::analysis::CorrelationEngine::validateConfig(config);
}

std::expected<void, std::string>
AnalysisDispatcher::runAnalysis(correlation::core::Trajectory *trajectory,
                                const ProgramOptions &options) {
  if (trajectory == nullptr || trajectory->getFrameCount() == 0) {
    std::string const err = AppDefaults::MSG_ANALYSIS_ABORTED;
    std::cerr << err << '\n';
    return std::unexpected(err);
  }

  stopPerFramePrebake();
  cancel_flag_ = false;

  if (const auto validation_res = validateOptions(options); !validation_res.has_value()) {
    return std::unexpected(validation_res.error());
  }

  const auto config = toEngineConfig(options, &cancel_flag_);
  auto result = correlation::analysis::CorrelationEngine::runAnalysis(*trajectory, config,
                                                                      progress_callback_);
  if (!result) {
    std::cerr << "Analysis Exception: " << result.error() << '\n';
    return std::unexpected(result.error());
  }

  {
    const std::unique_lock lock(df_mutex_);
    last_settings_ = config.settings;
    last_settings_.cancel_flag = nullptr;
    last_bond_cutoffs_ =
        !config.bond_cutoffs.empty() ? config.bond_cutoffs : trajectory->getBondCutoffs();
    df_ = std::move(result.value());
  }
  return {};
}

std::expected<void, std::string>
AnalysisDispatcher::runAnalysis(correlation::core::Trajectory &trajectory,
                                const ProgramOptions &options) {
  return runAnalysis(&trajectory, options);
}

std::expected<void, std::string>
AnalysisDispatcher::runFastPreviewRDF(correlation::core::Trajectory *trajectory,
                                      const ProgramOptions &options) {
  if (trajectory == nullptr || trajectory->getFrameCount() == 0) {
    return std::unexpected(AppDefaults::MSG_ANALYSIS_ABORTED);
  }

  ProgramOptions preview_opts = options;
  preview_opts.min_frame = 0;
  preview_opts.max_frame = 1;
  preview_opts.frame_stride = 1;
  preview_opts.active_calculators.clear();
  preview_opts.active_calculators["RDF"] = true;
  preview_opts.active_calculators["Radial Distribution Function"] = true;

  cancel_flag_ = false;
  if (const auto validation_res = validateOptions(preview_opts); !validation_res.has_value()) {
    return std::unexpected(validation_res.error());
  }

  auto config = toEngineConfig(preview_opts, &cancel_flag_);
  config.settings.populate_per_frame_store = false;
  auto result = correlation::analysis::CorrelationEngine::runAnalysis(*trajectory, config,
                                                                      progress_callback_);
  if (!result) {
    std::cerr << "Analysis Exception: " << result.error() << '\n';
    return std::unexpected(result.error());
  }

  {
    const std::unique_lock lock(df_mutex_);
    df_ = std::move(result.value());
  }
  return {};
}

bool AnalysisDispatcher::isBondDependentCalculator(std::string_view calc_id) noexcept {
  static constexpr std::array<std::string_view, 19> BOND_DEPENDENT = {
      "Plane Angle Distribution",
      "PAD",
      "BAD",
      "Dihedral Angle Distribution",
      "DAD",
      "Coordination Number",
      "CN",
      "Coordination",
      "Ring Size Distribution",
      "RD",
      "Rings",
      "Ring Statistics",
      "Common Neighbor Analysis",
      "CNA",
      "Cluster Size Distribution",
      "Cluster Analysis",
      "Hydrogen Bond Distribution",
      "Hydrogen Bond",
      "Motif",
  };
  return std::ranges::any_of(BOND_DEPENDENT,
                             [calc_id](std::string_view name) { return name == calc_id; });
}

std::vector<std::string> AnalysisDispatcher::getBondDependentCalculatorNames() {
  return {
      "Plane Angle Distribution",   "Dihedral Angle Distribution", "Coordination Number",
      "Ring Size Distribution",     "Common Neighbor Analysis",    "Cluster Size Distribution",
      "Hydrogen Bond Distribution",
  };
}

std::expected<void, std::string>
AnalysisDispatcher::runSelectiveRecalculation(correlation::core::Trajectory *trajectory,
                                              const ProgramOptions &options) {
  if (trajectory == nullptr || trajectory->getFrameCount() == 0) {
    std::string const err = AppDefaults::MSG_ANALYSIS_ABORTED;
    std::cerr << err << '\n';
    return std::unexpected(err);
  }

  // Fall back to complete analysis if no previous results are cached
  if (!df_) {
    return runAnalysis(trajectory, options);
  }

  cancel_flag_ = false;

  ProgramOptions selective_opts = options;
  selective_opts.active_calculators = extractBondDependentCalculators(options.active_calculators);

  // If no bond-dependent calculators are requested, return without recomputing
  if (selective_opts.active_calculators.empty()) {
    return {};
  }

  if (const auto validation_res = validateOptions(selective_opts); !validation_res.has_value()) {
    return std::unexpected(validation_res.error());
  }

  const auto config = toEngineConfig(selective_opts, &cancel_flag_);
  auto result = correlation::analysis::CorrelationEngine::runAnalysis(*trajectory, config,
                                                                      progress_callback_);
  if (!result) {
    std::cerr << "Selective Recalculation Exception: " << result.error() << '\n';
    return std::unexpected(result.error());
  }

  auto &new_df = result.value();
  if (new_df) {
    const std::unique_lock lock(df_mutex_);
    if (df_) {
      mergeCalculatedHistograms(*df_, *new_df);
    }
  }

  return {};
}

std::expected<void, std::string>
AnalysisDispatcher::runSelectiveRecalculation(correlation::core::Trajectory &trajectory,
                                              const ProgramOptions &options) {
  return runSelectiveRecalculation(&trajectory, options);
}

std::expected<void, std::string>
AnalysisDispatcher::writeFiles(const ProgramOptions &options) const {
  const std::shared_lock lock(df_mutex_);
  if (!df_) {
    std::string const err = AppDefaults::MSG_NO_DATA_TO_WRITE;
    std::cerr << err << '\n';
    return std::unexpected(err);
  }

  try {
    correlation::writers::FileWriter const writer(*df_);
    if (options.use_zip) {
      std::string zip_path = options.output_file_base;
      if (!zip_path.ends_with(".zip")) {
        zip_path += ".zip";
      }
      auto zip_res =
          writer.writeBundle(zip_path, options.use_csv, options.use_hdf5, options.use_parquet,
                             options.export_images, options.smoothing, options.export_algorithms);
      if (!zip_res) {
        return std::unexpected(zip_res.error());
      }
      std::cout << "Consolidated archive written to: " << zip_path << '\n';
    } else {
      auto folder_res = writer.writeFolder(
          options.output_file_base, options.use_csv, options.use_hdf5, options.use_parquet,
          options.export_images, options.smoothing, options.export_algorithms);
      if (!folder_res) {
        return std::unexpected(folder_res.error());
      }
      std::cout << "Categorized files written to directory: " << options.output_file_base << '\n';
    }
  } catch (const std::exception &e) {
    std::string const err = std::string(AppDefaults::MSG_ERROR_WRITING) + e.what();
    std::cerr << err << '\n';
    return std::unexpected(err);
  }
  return {};
}

std::vector<std::string> AnalysisDispatcher::getAvailableHistogramNames() const {
  const std::shared_lock lock(df_mutex_);
  if (!df_) {
    return {};
  }
  return df_->getAvailableHistograms();
}

const correlation::analysis::Histogram *
AnalysisDispatcher::getHistogram(const std::string &name) const {
  const std::shared_lock lock(df_mutex_);
  if (!df_) {
    return nullptr;
  }
  try {
    return &df_->getHistogram(name);
  } catch (const std::out_of_range &) {
    return nullptr;
  }
}

const std::map<std::string, correlation::analysis::Histogram> &
AnalysisDispatcher::getHistograms() const {
  static const std::map<std::string, correlation::analysis::Histogram> EMPTY_MAP;
  return df_ ? df_->getAllHistograms() : EMPTY_MAP;
}

std::map<std::string, real_t> AnalysisDispatcher::getAshcroftWeights() const {
  return df_ ? df_->getAshcroftWeights() : std::map<std::string, real_t>{};
}

} // namespace correlation::app
