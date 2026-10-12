/**
 * @file PlotRenderPipeline.cpp
 * @brief Background asynchronous rendering pipeline and caching for plot visualization.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/viewmodel/PlotRenderPipeline.hpp"
#include "app/formatters/PlotExportService.hpp"
#include "plotters/SvgPlotter.hpp"

#include <cmath>
#include <memory>
#include <mutex>
#include <span>
#include <utility>

namespace correlation::app {

namespace {

[[nodiscard]] bool isHoverEqual(const correlation::plotters::HoverInfo &lhs,
                                const correlation::plotters::HoverInfo &rhs) noexcept {
  if (lhs.active != rhs.active) {
    return false;
  }
  if (!lhs.active) {
    return true;
  }
  return std::abs(lhs.mouse_x - rhs.mouse_x) < 2.0F && std::abs(lhs.mouse_y - rhs.mouse_y) < 2.0F &&
         std::abs(lhs.widget_width - rhs.widget_width) < 1e-2F &&
         std::abs(lhs.widget_height - rhs.widget_height) < 1e-2F;
}

[[nodiscard]] bool isConfigEqual(const correlation::plotters::PlotConfig &lhs,
                                 const correlation::plotters::PlotConfig &rhs) noexcept {
  return lhs.theme == rhs.theme && lhs.preset_size == rhs.preset_size &&
         std::abs(lhs.width - rhs.width) < 1e-2 && std::abs(lhs.height - rhs.height) < 1e-2 &&
         lhs.palette == rhs.palette && std::abs(lhs.font_scale - rhs.font_scale) < 1e-4 &&
         std::abs(lhs.line_width - rhs.line_width) < 1e-4 &&
         std::abs(lhs.marker_size - rhs.marker_size) < 1e-4 && lhs.show_grid == rhs.show_grid &&
         lhs.show_legend == rhs.show_legend && lhs.show_markers == rhs.show_markers &&
         lhs.fill_area == rhs.fill_area && lhs.manual_x_min == rhs.manual_x_min &&
         lhs.manual_x_max == rhs.manual_x_max && lhs.manual_y_min == rhs.manual_y_min &&
         lhs.manual_y_max == rhs.manual_y_max &&
         lhs.reference_lines.size() == rhs.reference_lines.size();
}

} // namespace

PlotRenderPipeline::PlotRenderPipeline(ImageRenderedCallback callback,
                                       bool dispatch_to_slint_event_loop)
    : callback_(std::move(callback)), dispatch_to_slint_event_loop_(dispatch_to_slint_event_loop) {
  render_worker_ = std::jthread([this] { renderLoop(); });
}

PlotRenderPipeline::~PlotRenderPipeline() { stop(); }

void PlotRenderPipeline::stop() noexcept {
  alive_->store(false);
  {
    const std::scoped_lock lock(render_mutex_);
    stop_render_worker_ = true;
    pending_task_.reset();
  }
  render_cv_.notify_all();
  if (render_worker_.joinable()) {
    render_worker_.join();
  }
}

void PlotRenderPipeline::invalidateCache() noexcept {
  const std::scoped_lock lock(cache_mutex_);
  last_rendered_index_ = -1;
  last_pinned_runs_count_ = 0;
  last_config_ = {};
  last_hover_ = {};
}

bool PlotRenderPipeline::isCacheHit(int index, const correlation::plotters::PlotConfig &config,
                                    const correlation::plotters::HoverInfo &hover,
                                    std::size_t pinned_runs_count) const noexcept {
  const std::scoped_lock lock(cache_mutex_);
  return (index == last_rendered_index_ && pinned_runs_count == last_pinned_runs_count_ &&
          isConfigEqual(config, last_config_) && isHoverEqual(hover, last_hover_));
}

void PlotRenderPipeline::submitTask(RenderTask task) {
  {
    const std::scoped_lock lock(cache_mutex_);
    last_rendered_index_ = task.index;
    last_pinned_runs_count_ = task.pinned_runs_count;
    last_config_ = task.config;
    last_hover_ = task.hover;
  }
  {
    const std::scoped_lock lock(render_mutex_);
    if (stop_render_worker_) {
      return;
    }
    pending_task_ = std::move(task);
  }
  render_cv_.notify_one();
}

std::optional<RenderTask> PlotRenderPipeline::waitForNextTask() {
  std::unique_lock lock(render_mutex_);
  render_cv_.wait(lock, [this] { return stop_render_worker_ || pending_task_.has_value(); });

  if (stop_render_worker_ && !pending_task_.has_value()) {
    return std::nullopt;
  }

  if (!pending_task_.has_value()) {
    return std::nullopt;
  }

  auto task = std::move(pending_task_.value());
  pending_task_.reset();
  return task;
}

std::string PlotRenderPipeline::renderSvg(RenderTask &task) {
  if (!task.comparison_hists.empty()) {
    task.comparison_hists[0].hist = &task.active_hist;
  }

  if (task.comparison_hists.size() <= 1) {
    return correlation::plotters::renderHistogramAsSvg(task.active_hist, task.config, task.hover,
                                                       task.ashcroft_weights, task.curve_visibility,
                                                       task.custom_curve_colors);
  }

  const std::string key = PlotExportService::getComparisonKey(&task.active_hist);
  return correlation::plotters::renderComparisonSvg(task.comparison_hists, key, task.config,
                                                    task.hover);
}

void PlotRenderPipeline::dispatchImage(slint::Image img) {
  if (dispatch_to_slint_event_loop_) {
    slint::invoke_from_event_loop([this, alive = alive_, img = std::move(img)] mutable {
      if (!alive->load()) {
        return;
      }
      if (callback_) {
        callback_(img);
      }
    });
  } else {
    if (callback_ && alive_->load()) {
      callback_(img);
    }
  }
}

void PlotRenderPipeline::renderLoop() {
  while (true) {
    auto task_opt = waitForNextTask();
    if (!task_opt.has_value()) {
      break;
    }

    const std::string svg = renderSvg(task_opt.value());
    if (svg.empty()) {
      continue;
    }

    const auto *svg_bytes = reinterpret_cast<const uint8_t *>(svg.data());
    auto img = slint::Image::load_from_data(std::span<const uint8_t>(svg_bytes, svg.size()), "svg");

    dispatchImage(std::move(img));
  }
}

} // namespace correlation::app
