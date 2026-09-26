/**
 * @file PlotSeriesManager.cpp
 * @brief Implementation of PlotSeriesManager.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/formatters/PlotSeriesManager.hpp"
#include "plotters/SvgPlotter.hpp"

#include <algorithm>
#include <utility>

namespace correlation::app {

namespace {

[[nodiscard]] slint::Color parseHexColor(const std::string &hex) {
  if (hex.size() == 7 && hex[0] == '#') {
    auto red = static_cast<uint8_t>(std::stoul(hex.substr(1, 2), nullptr, 16));
    auto green = static_cast<uint8_t>(std::stoul(hex.substr(3, 2), nullptr, 16));
    auto blue = static_cast<uint8_t>(std::stoul(hex.substr(5, 2), nullptr, 16));
    return slint::Color::from_rgb_uint8(red, green, blue);
  }
  return slint::Color::from_rgb_uint8(128, 128, 128);
}

} // namespace

void PlotSeriesManager::setCurveVisible(int curve_id, bool visible) {
  if (curve_id >= 0 && std::cmp_less(curve_id, current_toggle_keys_.size())) {
    const std::string &key = current_toggle_keys_[curve_id];
    curve_visibility_map_[key] = visible;
  }
}

void PlotSeriesManager::setAllCurvesVisible(bool visible) {
  for (const auto &key : current_toggle_keys_) {
    curve_visibility_map_[key] = visible;
  }
}

bool PlotSeriesManager::isCurveVisible(const std::string &key, bool default_val) const {
  auto iter = curve_visibility_map_.find(key);
  if (iter != curve_visibility_map_.end()) {
    return iter->second;
  }
  return default_val;
}

void PlotSeriesManager::setCustomColor(int curve_id, const std::string &color_hex) {
  if (curve_id < 0 || std::cmp_greater_equal(curve_id, current_toggle_keys_.size())) {
    return;
  }
  const std::string &key = current_toggle_keys_[curve_id];
  custom_curve_colors_[key] = color_hex;
}

void PlotSeriesManager::pinCurrentRun(
    const std::map<std::string, correlation::analysis::Histogram> &hists) {
  if (hists.empty()) {
    return;
  }
  const std::string label = "Run " + std::to_string(pinned_runs_.size() + 1);
  pinned_runs_.push_back({.label = label, .histograms = hists});
}

void PlotSeriesManager::clearPinnedRuns() noexcept { pinned_runs_.clear(); }

void PlotSeriesManager::reset() noexcept {
  curve_visibility_map_.clear();
  custom_curve_colors_.clear();
  current_toggle_keys_.clear();
  pinned_runs_.clear();
  show_difference_curve_ = false;
}

std::shared_ptr<slint::VectorModel<CurveToggleItem>>
PlotSeriesManager::generateToggleItems(const correlation::analysis::Histogram *hist,
                                       const correlation::plotters::PlotConfig &config) {
  auto toggle_model = std::make_shared<slint::VectorModel<CurveToggleItem>>();
  current_toggle_keys_.clear();

  if (hist == nullptr) {
    return toggle_model;
  }

  const auto &partials = hist->smoothed_partials.empty() ? hist->partials : hist->smoothed_partials;

  std::vector<std::string> sorted_partial_keys;
  for (const auto &pair : partials) {
    if (pair.first != "Total") {
      sorted_partial_keys.push_back(pair.first);
    }
  }
  std::ranges::sort(sorted_partial_keys);

  const std::size_t total_curves =
      (partials.contains("Total") ? 1 : 0) + sorted_partial_keys.size() + pinned_runs_.size();

  int curve_id = 0;
  std::size_t color_idx = 0;

  auto get_hex_for_key = [&](const std::string &key) -> std::string {
    if (custom_curve_colors_.contains(key) && !custom_curve_colors_[key].empty()) {
      return custom_curve_colors_[key];
    }
    return correlation::plotters::detail::color(color_idx++, total_curves, config.palette);
  };

  if (partials.contains("Total")) {
    const bool vis =
        curve_visibility_map_.contains("Total") ? curve_visibility_map_["Total"] : true;
    curve_visibility_map_["Total"] = vis;
    current_toggle_keys_.emplace_back("Total");
    const std::string color_hex = get_hex_for_key("Total");
    toggle_model->push_back(CurveToggleItem{
        .id = curve_id++,
        .label = slint::SharedString("Total"),
        .color_hex = parseHexColor(color_hex),
        .visible = vis,
        .is_partial = false,
        .is_pinned = false,
    });
  }

  std::size_t rank = 0;
  for (const auto &p_key : sorted_partial_keys) {
    const bool default_vis = (rank < 7);
    const bool vis =
        curve_visibility_map_.contains(p_key) ? curve_visibility_map_[p_key] : default_vis;
    curve_visibility_map_[p_key] = vis;
    const std::string color_hex = get_hex_for_key(p_key);
    current_toggle_keys_.push_back(p_key);
    toggle_model->push_back(CurveToggleItem{
        .id = curve_id++,
        .label = slint::SharedString(p_key),
        .color_hex = parseHexColor(color_hex),
        .visible = vis,
        .is_partial = true,
        .is_pinned = false,
    });
    rank++;
  }

  for (const auto &pinned_run : pinned_runs_) {
    const std::string &pin_label = pinned_run.label;
    const bool vis =
        curve_visibility_map_.contains(pin_label) ? curve_visibility_map_[pin_label] : true;
    curve_visibility_map_[pin_label] = vis;
    const std::string color_hex = get_hex_for_key(pin_label);
    current_toggle_keys_.push_back(pin_label);
    toggle_model->push_back(CurveToggleItem{
        .id = curve_id++,
        .label = slint::SharedString(pin_label),
        .color_hex = parseHexColor(color_hex),
        .visible = vis,
        .is_partial = false,
        .is_pinned = true,
    });
  }

  return toggle_model;
}

} // namespace correlation::app
