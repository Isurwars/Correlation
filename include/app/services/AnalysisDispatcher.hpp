/**
 * @file AnalysisDispatcher.hpp
 * @brief Service coordinating correlation analysis execution and file export.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "app/core/AppOptions.hpp"
#include "core/Trajectory.hpp"

#include <atomic>
#include <expected>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace correlation::app {

/**
 * @class AnalysisDispatcher
 * @brief Dispatches analysis workflows to CorrelationEngine, handles cancellation, and exports
 * data.
 */
class AnalysisDispatcher {
public:
  using ProgressCallback = std::function<void(float, const std::string &)>;

  AnalysisDispatcher() = default;
  ~AnalysisDispatcher() = default;

  AnalysisDispatcher(const AnalysisDispatcher &) = delete;
  AnalysisDispatcher &operator=(const AnalysisDispatcher &) = delete;
  AnalysisDispatcher(AnalysisDispatcher &&) = delete;
  AnalysisDispatcher &operator=(AnalysisDispatcher &&) = delete;

  /**
   * @brief Executes analysis for the given trajectory and options.
   * @param[in] trajectory Source trajectory with structures to analyze.
   * @param[in] options Analytical configuration parameters.
   * @return Success or error message.
   */
  std::expected<void, std::string> runAnalysis(correlation::core::Trajectory *trajectory,
                                               const ProgramOptions &options);

  /**
   * @brief Executes analysis for the given trajectory reference and options.
   * @param[in] trajectory Source trajectory reference with structures to analyze.
   * @param[in] options Analytical configuration parameters.
   * @return Success or error message.
   */
  std::expected<void, std::string> runAnalysis(correlation::core::Trajectory &trajectory,
                                               const ProgramOptions &options);

  /**
   * @brief Validates analysis configuration before execution.
   * @param[in] options Options to validate.
   * @return Empty string if valid, error message otherwise.
   */
  [[nodiscard]] static std::string validateOptions(const ProgramOptions &options);

  /**
   * @brief Writes analytical results to disk (CSV, HDF5, Parquet).
   * @param[in] options Configuration specifying output formats and destination.
   * @return Success or error message.
   */
  std::expected<void, std::string> writeFiles(const ProgramOptions &options) const;

  /**
   * @brief Signals running analysis to abort at next cancellation checkpoint.
   */
  void cancelAnalysis() noexcept { cancel_flag_ = true; }

  /**
   * @brief Checks if cancellation was requested.
   */
  [[nodiscard]] bool isCancelled() const noexcept { return cancel_flag_; }

  /**
   * @brief Sets progress notification hook.
   */
  void setProgressCallback(ProgressCallback callback) { progress_callback_ = std::move(callback); }

  /**
   * @brief Gets all available histogram names from last completed analysis.
   */
  [[nodiscard]] std::vector<std::string> getAvailableHistogramNames() const;

  /**
   * @brief Retrieves a specific analytical histogram by name.
   */
  [[nodiscard]] const correlation::analysis::Histogram *getHistogram(const std::string &name) const;

  /**
   * @brief Retrieves all analytical histograms from last run.
   */
  [[nodiscard]] const std::map<std::string, correlation::analysis::Histogram> &
  getHistograms() const;

  /**
   * @brief Retrieves Ashcroft-Langreth weights for partial structure factors.
   */
  [[nodiscard]] std::map<std::string, real_t> getAshcroftWeights() const;

  /**
   * @brief Returns active DistributionFunctions pointer.
   */
  [[nodiscard]] const correlation::analysis::DistributionFunctions *
  getDistributionFunctions() const noexcept {
    return df_.get();
  }

  /**
   * @brief Injects an existing DistributionFunctions container (useful for testing).
   */
  void setDistributionFunctions(
      std::unique_ptr<correlation::analysis::DistributionFunctions> df) noexcept {
    df_ = std::move(df);
  }

  /**
   * @brief Clears current results.
   */
  void clear() noexcept { df_.reset(); }

private:
  std::unique_ptr<correlation::analysis::DistributionFunctions> df_;
  std::atomic<bool> cancel_flag_{false};
  ProgressCallback progress_callback_;
};

} // namespace correlation::app
