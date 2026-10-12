/**
 * @file PlotRenderPipeline.hpp
 * @brief Background asynchronous rendering pipeline and caching for plot visualization.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "analysis/DistributionFunctions.hpp"
#include "plotters/SvgPlotter.hpp"
#include <slint.h>

#include <atomic>
#include <condition_variable>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <vector>

namespace correlation::app {

/**
 * @struct RenderTask
 * @brief Self-contained payload consumed by the background SVG render worker.
 */
struct RenderTask {
  int index{-1};                                                         ///< Target plot index
  correlation::analysis::Histogram active_hist;                          ///< Active histogram data
  std::vector<correlation::plotters::LabeledHistogram> comparison_hists; ///< Comparison histograms
  correlation::plotters::PlotConfig config;                              ///< Output styling config
  correlation::plotters::HoverInfo hover;                                ///< Cursor hover metadata
  std::map<std::string, real_t> ashcroft_weights;                        ///< Ashcroft weights
  std::map<std::string, bool> curve_visibility;                          ///< Curve visibility flags
  std::map<std::string, std::string> custom_curve_colors;                ///< Custom curve colors
  std::size_t pinned_runs_count{0};                                      ///< Pinned runs count
};

/**
 * @class PlotRenderPipeline
 * @brief Manages background SVG rasterization/generation worker thread, task queue, and caching.
 */
class PlotRenderPipeline {
public:
  using ImageRenderedCallback = std::function<void(const slint::Image &)>;

  /**
   * @brief Constructs pipeline and starts background rendering worker thread.
   * @param[in] callback Callback invoked upon successful image generation.
   * @param[in] dispatch_to_slint_event_loop If true, dispatches via slint::invoke_from_event_loop.
   */
  explicit PlotRenderPipeline(ImageRenderedCallback callback,
                              bool dispatch_to_slint_event_loop = true);

  /**
   * @brief Destructor. Signals stop condition and joins worker thread safely.
   */
  ~PlotRenderPipeline();

  PlotRenderPipeline(const PlotRenderPipeline &) = delete;
  PlotRenderPipeline &operator=(const PlotRenderPipeline &) = delete;
  PlotRenderPipeline(PlotRenderPipeline &&) = delete;
  PlotRenderPipeline &operator=(PlotRenderPipeline &&) = delete;

  /**
   * @brief Evaluates whether a new render request matches previously rendered frame parameters.
   * @param[in] index Target plot index.
   * @param[in] config Target plot configuration.
   * @param[in] hover Target hover info.
   * @param[in] pinned_runs_count Current pinned comparison runs count.
   * @return True if previous frame exactly matches and re-render is unnecessary.
   */
  [[nodiscard]] bool isCacheHit(int index, const correlation::plotters::PlotConfig &config,
                                const correlation::plotters::HoverInfo &hover,
                                std::size_t pinned_runs_count) const noexcept;

  /**
   * @brief Invalidates cached parameters to force re-rendering on the next request.
   */
  void invalidateCache() noexcept;

  /**
   * @brief Enqueues a render task to the worker queue, discarding superseded tasks.
   * @param[in] task Self-contained render payload.
   */
  void submitTask(RenderTask task);

  /**
   * @brief Stops worker thread and cancels pending operations.
   */
  void stop() noexcept;

private:
  [[nodiscard]] std::optional<RenderTask> waitForNextTask();
  [[nodiscard]] static std::string renderSvg(RenderTask &task);
  void dispatchImage(slint::Image img);
  void renderLoop();

  ImageRenderedCallback callback_;
  bool dispatch_to_slint_event_loop_{true};
  std::jthread render_worker_;
  std::mutex render_mutex_;
  std::condition_variable render_cv_;
  std::optional<RenderTask> pending_task_;
  bool stop_render_worker_{false};
  std::shared_ptr<std::atomic<bool>> alive_{std::make_shared<std::atomic<bool>>(true)};

  mutable std::mutex cache_mutex_;
  int last_rendered_index_{-1};
  correlation::plotters::PlotConfig last_config_;
  correlation::plotters::HoverInfo last_hover_;
  std::size_t last_pinned_runs_count_{0};
};

} // namespace correlation::app
