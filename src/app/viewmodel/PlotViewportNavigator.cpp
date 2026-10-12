/**
 * @file PlotViewportNavigator.cpp
 * @brief Implementation of PlotViewportNavigator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/viewmodel/PlotViewportNavigator.hpp"

#include <algorithm>
#include <cmath>

namespace correlation::app {

correlation::plotters::detail::NiceScale
PlotViewportNavigator::computeXScale(const correlation::analysis::Histogram &hist,
                                     const correlation::plotters::PlotConfig &config) {
  real_t raw_x_min = hist.bins.empty() ? static_cast<real_t>(0.0) : hist.bins.front();
  real_t raw_x_max = hist.bins.empty() ? static_cast<real_t>(1.0) : hist.bins.back();
  bool strict_x = false;
  if (config.manual_x_min.has_value() && config.manual_x_max.has_value()) {
    raw_x_min = *config.manual_x_min;
    raw_x_max = *config.manual_x_max;
    strict_x = true;
  }
  return correlation::plotters::detail::NiceScale(
      correlation::plotters::detail::DataRange{.min = raw_x_min, .max = raw_x_max}, 11, strict_x);
}

correlation::plotters::detail::NiceScale
PlotViewportNavigator::computeYScale(const correlation::analysis::Histogram &hist,
                                     const correlation::plotters::PlotConfig &config) {
  auto raw_y_min = static_cast<real_t>(0.0);
  auto raw_y_max = static_cast<real_t>(0.0);
  const auto &partials = hist.smoothed_partials.empty() ? hist.partials : hist.smoothed_partials;
  for (const auto &[key, value] : partials) {
    for (const real_t val : value) {
      raw_y_min = std::min(raw_y_min, val);
      raw_y_max = std::max(raw_y_max, val);
    }
  }
  const real_t y_padding = (raw_y_max - raw_y_min) * static_cast<real_t>(0.05);
  raw_y_max += y_padding;

  bool strict_y = false;
  if (config.manual_y_min.has_value() && config.manual_y_max.has_value()) {
    raw_y_min = *config.manual_y_min;
    raw_y_max = *config.manual_y_max;
    strict_y = true;
  }
  return correlation::plotters::detail::NiceScale(
      correlation::plotters::detail::DataRange{.min = raw_y_min, .max = raw_y_max}, 8, strict_y);
}

std::pair<real_t, real_t>
PlotViewportNavigator::screenToData(float screen_x, float screen_y, float canvas_w, float canvas_h,
                                    const correlation::analysis::Histogram &hist,
                                    const correlation::plotters::PlotConfig &config) {
  const auto x_scale = computeXScale(hist, config);
  const auto y_scale = computeYScale(hist, config);
  return correlation::plotters::detail::screenToData(screen_x, screen_y, canvas_w, canvas_h, config,
                                                     x_scale, y_scale);
}

bool PlotViewportNavigator::hasActiveZoom() const noexcept {
  return zoom_x_min_.has_value() || zoom_x_max_.has_value() || zoom_y_min_.has_value() ||
         zoom_y_max_.has_value();
}

void PlotViewportNavigator::applyZoomRect(float start_x, float start_y, float end_x, float end_y,
                                          float canvas_w, float canvas_h,
                                          const correlation::analysis::Histogram &hist,
                                          const correlation::plotters::PlotConfig &config) {
  if (std::abs(end_x - start_x) < 5.0F || std::abs(end_y - start_y) < 5.0F) {
    return;
  }
  if (hist.bins.empty()) {
    return;
  }

  auto [dx1, dy1] = screenToData(start_x, start_y, canvas_w, canvas_h, hist, config);
  auto [dx2, dy2] = screenToData(end_x, end_y, canvas_w, canvas_h, hist, config);

  zoom_x_min_ = std::min(dx1, dx2);
  zoom_x_max_ = std::max(dx1, dx2);
  zoom_y_min_ = std::min(dy1, dy2);
  zoom_y_max_ = std::max(dy1, dy2);
}

void PlotViewportNavigator::applyPan(float delta_x, float delta_y, float mouse_x, float mouse_y,
                                     float canvas_w, float canvas_h,
                                     const correlation::analysis::Histogram &hist,
                                     const correlation::plotters::PlotConfig &config) {
  if (std::abs(delta_x) < 0.5F && std::abs(delta_y) < 0.5F) {
    return;
  }
  if (hist.bins.empty()) {
    return;
  }

  const auto x_scale = computeXScale(hist, config);
  const auto y_scale = computeYScale(hist, config);

  const auto [d0_x, d0_y] = screenToData(mouse_x, mouse_y, canvas_w, canvas_h, hist, config);
  const auto [d1_x, d1_y] =
      screenToData(mouse_x + delta_x, mouse_y + delta_y, canvas_w, canvas_h, hist, config);

  const real_t shift_x = d1_x - d0_x;
  const real_t shift_y = d1_y - d0_y;

  const real_t cur_x_min = zoom_x_min_.value_or(x_scale.min);
  const real_t cur_x_max = zoom_x_max_.value_or(x_scale.max);
  const real_t cur_y_min = zoom_y_min_.value_or(y_scale.min);
  const real_t cur_y_max = zoom_y_max_.value_or(y_scale.max);

  zoom_x_min_ = cur_x_min - shift_x;
  zoom_x_max_ = cur_x_max - shift_x;
  zoom_y_min_ = cur_y_min - shift_y;
  zoom_y_max_ = cur_y_max - shift_y;
}

void PlotViewportNavigator::applyWheelZoom(float mouse_x, float mouse_y, float delta,
                                           float canvas_w, float canvas_h,
                                           const correlation::analysis::Histogram &hist,
                                           const correlation::plotters::PlotConfig &config) {
  if (std::abs(delta) < 0.01F) {
    return;
  }
  if (hist.bins.empty()) {
    return;
  }

  const auto x_scale = computeXScale(hist, config);
  const auto y_scale = computeYScale(hist, config);

  const auto [c_x, c_y] = screenToData(mouse_x, mouse_y, canvas_w, canvas_h, hist, config);

  const real_t zoom_factor = (delta > 0.0F) ? 0.85 : 1.18;

  const real_t cur_x_min = zoom_x_min_.value_or(x_scale.min);
  const real_t cur_x_max = zoom_x_max_.value_or(x_scale.max);
  const real_t cur_y_min = zoom_y_min_.value_or(y_scale.min);
  const real_t cur_y_max = zoom_y_max_.value_or(y_scale.max);

  zoom_x_min_ = c_x - (c_x - cur_x_min) * zoom_factor;
  zoom_x_max_ = c_x + (cur_x_max - c_x) * zoom_factor;
  zoom_y_min_ = c_y - (c_y - cur_y_min) * zoom_factor;
  zoom_y_max_ = c_y + (cur_y_max - c_y) * zoom_factor;
}

void PlotViewportNavigator::resetZoom() noexcept {
  zoom_x_min_.reset();
  zoom_x_max_.reset();
  zoom_y_min_.reset();
  zoom_y_max_.reset();
}

ViewportLimits PlotViewportNavigator::limits() const noexcept {
  return ViewportLimits{
      .x_min = zoom_x_min_,
      .x_max = zoom_x_max_,
      .y_min = zoom_y_min_,
      .y_max = zoom_y_max_,
  };
}

void PlotViewportNavigator::setLimits(ViewportLimits limits) noexcept {
  zoom_x_min_ = limits.x_min;
  zoom_x_max_ = limits.x_max;
  zoom_y_min_ = limits.y_min;
  zoom_y_max_ = limits.y_max;
}

void PlotViewportNavigator::applyToConfig(
    correlation::plotters::PlotConfig &config) const noexcept {
  if (zoom_x_min_.has_value() && zoom_x_max_.has_value()) {
    config.manual_x_min = zoom_x_min_;
    config.manual_x_max = zoom_x_max_;
  }
  if (zoom_y_min_.has_value() && zoom_y_max_.has_value()) {
    config.manual_y_min = zoom_y_min_;
    config.manual_y_max = zoom_y_max_;
  }
}

} // namespace correlation::app
