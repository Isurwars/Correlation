/**
 * @file AppBackend.hpp
 * @brief Application backend facade coordinating TrajectoryLoader, AnalysisDispatcher, and
 * BondCutoffService.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "app/core/AppOptions.hpp"
#include "app/services/AnalysisDispatcher.hpp"
#include "app/services/BondCutoffService.hpp"
#include "app/services/TrajectoryLoader.hpp"
#include "core/Trajectory.hpp"
#include "math/Smoothing.hpp"

#include <expected>
#include <functional>
#include <map>
#include <memory>

namespace correlation::app {

/**
 * @class AppBackend
 * @brief Facade composing TrajectoryLoader, AnalysisDispatcher, and BondCutoffService.
 */
class AppBackend {
public:
  AppBackend();
  ~AppBackend() = default;

  AppBackend(const AppBackend &) = delete;
  AppBackend &operator=(const AppBackend &) = delete;
  AppBackend(AppBackend &&) = delete;
  AppBackend &operator=(AppBackend &&) = delete;

  /**
   * @brief Sets the program options.
   * @param opt The options to set.
   */
  void setOptions(const ProgramOptions &opt) { options_ = opt; }

  /**
   * @brief Gets the current program options.
   * @return Mutable reference to current options.
   */
  [[nodiscard]] ProgramOptions &options() noexcept { return options_; }

  /**
   * @brief Gets the current program options (const).
   * @return Const reference to current options.
   */
  [[nodiscard]] const ProgramOptions &options() const noexcept { return options_; }

  /**
   * @brief Updates the enabled state of a single calculator.
   */
  void setCalculatorActive(const std::string &calc_id, bool enabled) {
    options_.active_calculators[calc_id] = enabled;
  }

  /**
   * @brief Gets a pointer to the current cell (first frame of trajectory).
   */
  [[nodiscard]] const correlation::core::Cell *cell() const { return loader_.cell(); }

  /**
   * @brief Loads a file from a given path.
   */
  std::string loadFile(const std::string &path);

  /**
   * @brief Runs the analysis based on current options and loaded trajectory.
   */
  [[nodiscard]] std::expected<void, std::string> runAnalysis();

  /**
   * @brief Writes the analysis results to files (CSV, HDF5, Parquet).
   */
  [[nodiscard]] std::expected<void, std::string> writeFiles();

  /**
   * @brief Gets the atom counts for the current structure.
   */
  [[nodiscard]] std::map<std::string, int> getAtomCounts() const { return loader_.getAtomCounts(); }

  /**
   * @brief Gets the total number of frames in the trajectory.
   */
  [[nodiscard]] size_t getFrameCount() const { return loader_.getFrameCount(); }

  /**
   * @brief Gets the total number of atoms in the first frame.
   */
  [[nodiscard]] size_t getTotalAtomCount() const { return loader_.getTotalAtomCount(); }

  /**
   * @brief Gets the count of frames removed/skipped during loading.
   */
  [[nodiscard]] size_t getRemovedFrameCount() const { return loader_.getRemovedFrameCount(); }

  /**
   * @brief Gets the time step of the trajectory.
   */
  [[nodiscard]] real_t getTimeStep() const { return loader_.getTimeStep(); }

  /**
   * @brief Calculates a recommended time step in fs based on atomic masses.
   */
  [[nodiscard]] real_t getRecommendedTimeStep() const { return loader_.getRecommendedTimeStep(); }

  /**
   * @brief Calculates recommended bond cutoffs (min and max).
   */
  [[nodiscard]] correlation::analysis::BondCutoffMatrix getRecommendedBondCutoffs();

  /**
   * @brief Gets the bond cutoff for a specific pair of element types.
   */
  [[nodiscard]] real_t getBondCutoff(size_t type1, size_t type2) const {
    return cutoff_service_.getBondCutoff(loader_.trajectory(), type1, type2);
  }

  /**
   * @brief Gets the minimum bond cutoff for a specific pair of element types.
   */
  [[nodiscard]] real_t getMinBondCutoff(size_t type1, size_t type2) const {
    return cutoff_service_.getMinBondCutoff(loader_.trajectory(), type1, type2);
  }

  /**
   * @brief Sets the bond cutoffs to be used in analysis.
   */
  void setBondCutoffs(const correlation::analysis::BondCutoffMatrix &cutoffs);

  /**
   * @brief Scales recommended covalent cutoffs by a multiplicative factor.
   */
  correlation::analysis::BondCutoffMatrix applyScaledBondCutoffs(real_t scale_factor);

  /**
   * @brief Assigns uniform min and max cutoff distances across all atom pairs.
   */
  correlation::analysis::BondCutoffMatrix setUniformBondCutoff(real_t min_cutoff,
                                                               real_t max_cutoff);

  /**
   * @brief Returns names of all available histograms.
   */
  [[nodiscard]] std::vector<std::string> getAvailableHistogramNames() const {
    return dispatcher_.getAvailableHistogramNames();
  }

  /**
   * @brief Returns a pointer to a specific histogram.
   */
  [[nodiscard]] const correlation::analysis::Histogram *
  getHistogram(const std::string &name) const {
    return dispatcher_.getHistogram(name);
  }

  /**
   * @brief Returns Ashcroft-Langreth weights.
   */
  [[nodiscard]] std::map<std::string, real_t> getAshcroftWeights() const {
    return dispatcher_.getAshcroftWeights();
  }

  /**
   * @brief Returns all histograms from last analysis.
   */
  [[nodiscard]] const std::map<std::string, correlation::analysis::Histogram> &
  getHistograms() const {
    return dispatcher_.getHistograms();
  }

  /**
   * @brief Returns DistributionFunctions results.
   */
  [[nodiscard]] const correlation::analysis::DistributionFunctions *
  getDistributionFunctions() const {
    return dispatcher_.getDistributionFunctions();
  }

  /**
   * @brief Sets progress notification hook.
   */
  void setProgressCallback(std::function<void(float, const std::string &)> callback) {
    progress_callback_ = callback;
    dispatcher_.setProgressCallback(std::move(callback));
  }

  /**
   * @brief Cancels running analysis.
   */
  void cancelAnalysis() { dispatcher_.cancelAnalysis(); }

  /**
   * @brief Checks if analysis was cancelled.
   */
  [[nodiscard]] bool isCancelled() const { return dispatcher_.isCancelled(); }

  /**
   * @brief Access underlying TrajectoryLoader service.
   */
  [[nodiscard]] TrajectoryLoader &trajectoryLoader() noexcept { return loader_; }
  [[nodiscard]] const TrajectoryLoader &trajectoryLoader() const noexcept { return loader_; }
  [[nodiscard]] TrajectoryLoader &loader() noexcept { return loader_; }
  [[nodiscard]] const TrajectoryLoader &loader() const noexcept { return loader_; }

  /**
   * @brief Access underlying AnalysisDispatcher service.
   */
  [[nodiscard]] AnalysisDispatcher &analysisDispatcher() noexcept { return dispatcher_; }
  [[nodiscard]] const AnalysisDispatcher &analysisDispatcher() const noexcept {
    return dispatcher_;
  }
  [[nodiscard]] AnalysisDispatcher &dispatcher() noexcept { return dispatcher_; }
  [[nodiscard]] const AnalysisDispatcher &dispatcher() const noexcept { return dispatcher_; }

  /**
   * @brief Access underlying BondCutoffService.
   */
  [[nodiscard]] BondCutoffService &bondCutoffService() noexcept { return cutoff_service_; }
  [[nodiscard]] const BondCutoffService &bondCutoffService() const noexcept {
    return cutoff_service_;
  }
  [[nodiscard]] BondCutoffService &cutoffService() noexcept { return cutoff_service_; }
  [[nodiscard]] const BondCutoffService &cutoffService() const noexcept { return cutoff_service_; }

private:
  TrajectoryLoader loader_;
  AnalysisDispatcher dispatcher_;
  BondCutoffService cutoff_service_;
  ProgramOptions options_;
  std::function<void(float, const std::string &)> progress_callback_;
};

} // namespace correlation::app
