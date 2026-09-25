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
#include <iostream>

namespace correlation::app {

namespace {

[[nodiscard]] correlation::analysis::CorrelationEngineConfig
toEngineConfig(const ProgramOptions &options, std::atomic<bool> *cancel_flag) {
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
  config.settings.xrd_params = options.xrd_params;

  config.bond_cutoffs = options.bond_cutoffs;
  config.min_frame = options.min_frame;
  config.max_frame = options.max_frame;
  config.time_step = options.time_step;
  config.settings.frame_stride = static_cast<size_t>(std::max(1, options.frame_stride));
  return config;
}

} // namespace

std::string AnalysisDispatcher::validateOptions(const ProgramOptions &options) {
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

  cancel_flag_ = false;

  std::string const validation_error = validateOptions(options);
  if (!validation_error.empty()) {
    return std::unexpected(validation_error);
  }

  const auto config = toEngineConfig(options, &cancel_flag_);
  auto result = correlation::analysis::CorrelationEngine::runAnalysis(
      *trajectory, config, progress_callback_);
  if (!result) {
    std::cerr << "Analysis Exception: " << result.error() << '\n';
    return std::unexpected(result.error());
  }

  df_ = std::move(result.value());
  return {};
}

std::expected<void, std::string>
AnalysisDispatcher::runAnalysis(correlation::core::Trajectory &trajectory,
                                const ProgramOptions &options) {
  return runAnalysis(&trajectory, options);
}

std::expected<void, std::string>
AnalysisDispatcher::writeFiles(const ProgramOptions &options) const {
  if (!df_) {
    std::string const err = AppDefaults::MSG_NO_DATA_TO_WRITE;
    std::cerr << err << '\n';
    return std::unexpected(err);
  }

  try {
    correlation::writers::FileWriter const writer(*df_);
    writer.write(options.output_file_base, options.use_csv, options.use_hdf5, options.use_parquet,
                 options.smoothing);
    std::cout << "Files written to: " << options.output_file_base << '\n';
  } catch (const std::exception &e) {
    std::string const err = std::string(AppDefaults::MSG_ERROR_WRITING) + e.what();
    std::cerr << err << '\n';
    return std::unexpected(err);
  }
  return {};
}

std::vector<std::string> AnalysisDispatcher::getAvailableHistogramNames() const {
  if (!df_) {
    return {};
  }
  return df_->getAvailableHistograms();
}

const correlation::analysis::Histogram *
AnalysisDispatcher::getHistogram(const std::string &name) const {
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
