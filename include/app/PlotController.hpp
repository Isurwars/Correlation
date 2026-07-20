/**
 * @file PlotController.hpp
 * @brief Plotting and SVG rendering management.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "app/AppBackend.hpp"
#include "plotters/SvgPlotter.hpp"
#include <slint.h>

#include <atomic>
#include <chrono>
#include <map>
#include <memory>
#include <string>
#include <thread>
#include <vector>

class AppWindow;

namespace correlation::app {

struct PlotSize {
  float width;
  float height;
};

/**
 * @class PlotController
 * @brief Manages plot rendering, caching, interactions, and exports.
 */
class PlotController {
public:
  /**
   * @brief Constructs the PlotController.
   * @param[in,out] window Reference to the UI window.
   * @param[in,out] backend Reference to the application backend.
   */
  PlotController(::AppWindow &window, AppBackend &backend);

  /**
   * @brief Destructor. Joins render worker thread before destruction.
   */
  ~PlotController();

  /**
   * @brief Populates available plot names into the UI dropdown list.
   */
  void populatePlotList();

  /**
   * @brief Handles mouse movement and hover tooltip calculation on the plot canvas.
   * @param[in] mouse_x Mouse X coordinate on canvas.
   * @param[in] mouse_y Mouse Y coordinate on canvas.
   * @param[in] hover Flag indicating active mouse hover state.
   * @param[in] width Canvas layout width.
   * @param[in] height Canvas layout height.
   */
  void handleMouseMove(float mouse_x, float mouse_y, bool hover, float width, float height);

  /**
   * @brief Schedules or immediately triggers a plot redraw for the specified plot index.
   * @param[in] index Target plot selection index.
   * @param[in] immediate If true, skips debounce timer and updates immediately.
   */
  void requestPlotUpdate(int index, bool immediate = false);

  /**
   * @brief Handles saving the current plot to file (SVG/PDF).
   */
  void handleSavePlot();

  /**
   * @brief Pins the current analysis run for comparison overlay.
   */
  void handlePinRun();

  /**
   * @brief Clears all pinned comparison overlay runs.
   */
  void handleClearPinnedRuns();

  /**
   * @brief Handles plot resized callback from UI layout changes.
   * @param[in] size New plot container dimensions.
   */
  void handlePlotResized(PlotSize size);

  /**
   * @brief Timer callback for periodic plot rendering updates.
   */
  void handleUpdateTimer();

private:
  ::AppWindow &window_;
  AppBackend &backend_;

  std::thread render_thread_;

  struct RenderTaskData {
    correlation::analysis::Histogram active_hist;
    std::vector<std::pair<std::string, correlation::analysis::Histogram>> comparison_hists;
    correlation::plotters::PlotConfig config;
    correlation::plotters::HoverInfo hover;
    std::map<std::string, real_t> ashcroft_weights;
  };

  std::atomic<bool> is_rendering_{false};
  std::atomic<bool> render_pending_{false};

  std::vector<std::string> available_plot_keys_;

  struct PinnedRun {
    std::string label;
    std::map<std::string, correlation::analysis::Histogram> histograms;
  };
  std::vector<PinnedRun> pinned_runs_;

  float last_mouse_x_ = -1.0F;
  float last_mouse_y_ = -1.0F;
  bool mouse_hover_ = false;
  float last_plot_width_ = 0.0F;
  float last_plot_height_ = 0.0F;

  std::chrono::steady_clock::time_point last_replot_time_;
  slint::Timer hover_timer_;
  slint::Timer update_timer_;
  bool update_scheduled_ = false;
  int pending_plot_index_ = -1;
  bool needs_redraw_ = false;

  int last_rendered_index_ = -1;
  correlation::plotters::PlotConfig last_config_;
  correlation::plotters::HoverInfo last_hover_;
  std::size_t last_pinned_runs_count_ = 0;

  std::shared_ptr<std::string> current_svg_;

  correlation::plotters::PlotConfig buildPlotConfigFromUI();
  bool isPlotCacheHit(int index, const correlation::plotters::PlotConfig &config,
                      const correlation::plotters::HoverInfo &hover) const;
  void executePlotRender(RenderTaskData data);
  void updateTableData(const correlation::analysis::Histogram *hist);
  void executeSavePlot(const std::string &filepath, const correlation::analysis::Histogram *hist,
                       const std::string &name);
};

} // namespace correlation::app
