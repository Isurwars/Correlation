/**
 * @file PlotExportManager.hpp
 * @brief Coordinates native file dialog, multi-dataset packaging, and background plot export.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "app/formatters/PlotSeriesManager.hpp"
#include "plotters/PlotTypes.hpp"

#include <atomic>
#include <expected>
#include <functional>
#include <map>
#include <span>
#include <string>
#include <thread>
#include <vector>

namespace correlation::app {

/**
 * @class PlotExportManager
 * @brief Coordinates OS native save dialogs, comparison datasets formatting, and export I/O.
 */
class PlotExportManager {
public:
  /**
   * @brief Callback type for reporting export results or user cancellations.
   */
  using StatusCallback = std::function<void(const std::string &)>;

  /**
   * @brief Default constructor.
   */
  PlotExportManager();

  /**
   * @brief Destructor. Joins active dialog thread if running.
   */
  ~PlotExportManager();

  PlotExportManager(const PlotExportManager &) = delete;
  PlotExportManager &operator=(const PlotExportManager &) = delete;
  PlotExportManager(PlotExportManager &&) = delete;
  PlotExportManager &operator=(PlotExportManager &&) = delete;

  /**
   * @brief Saves a plot synchronously to a destination filepath.
   * @param[in] filepath Target file path on disk (SVG or PDF).
   * @param[in] plot_name Name of the plot to save.
   * @param[in] hist Active histogram pointer.
   * @param[in] config Base plot configuration.
   * @param[in] pinned_runs Span of pinned runs for comparison exports.
   * @param[in] ashcroft_weights Map of Ashcroft weights for S(Q) partials.
   * @return Expected void on success or error description string on failure.
   */
  [[nodiscard]] static std::expected<void, std::string>
  savePlot(const std::string &filepath, const std::string &plot_name,
           const correlation::analysis::Histogram *hist, correlation::plotters::PlotConfig config,
           std::span<const PinnedRun> pinned_runs,
           const std::map<std::string, real_t> &ashcroft_weights = {});

  /**
   * @brief Launches the native OS save dialog on a background worker and executes the save on
   * completion.
   * @param[in] default_base_path Base output file path from program options.
   * @param[in] plot_name Name of the plot to export.
   * @param[in] hist Active histogram pointer.
   * @param[in] config Base plot configuration.
   * @param[in] pinned_runs Vector of pinned runs for comparison exports.
   * @param[in] ashcroft_weights Map of Ashcroft weights.
   * @param[in] on_status Callback invoked on the Slint main event loop with status text.
   */
  void exportPlotAsync(const std::string &default_base_path, const std::string &plot_name,
                       const correlation::analysis::Histogram *hist,
                       correlation::plotters::PlotConfig config, std::vector<PinnedRun> pinned_runs,
                       std::map<std::string, real_t> ashcroft_weights, StatusCallback on_status);

  /**
   * @brief Checks if a native save file dialog is currently active.
   * @return True if dialog is running, false otherwise.
   */
  [[nodiscard]] bool isDialogActive() const noexcept;

private:
  std::jthread dialog_thread_;
  std::atomic<bool> dialog_active_{false};
};

} // namespace correlation::app
