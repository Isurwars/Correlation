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
#include <cmath>
#include <utility>

namespace correlation::app {

namespace {

struct RankedPartial {
  std::string key;
  real_t score{0.0};
};

[[nodiscard]] slint::Color parseHexColor(const std::string &hex) {
  if (hex.size() == 7 && hex[0] == '#') {
    const auto red = static_cast<uint8_t>(std::stoul(hex.substr(1, 2), nullptr, 16));
    const auto green = static_cast<uint8_t>(std::stoul(hex.substr(3, 2), nullptr, 16));
    const auto blue = static_cast<uint8_t>(std::stoul(hex.substr(5, 2), nullptr, 16));
    return slint::Color::from_rgb_uint8(red, green, blue);
  }
  return slint::Color::from_rgb_uint8(128, 128, 128);
}

[[nodiscard]] real_t computePartialScore(const std::string &key, const std::vector<real_t> &values,
                                         const std::map<std::string, real_t> &weights) {
  if (!weights.empty()) {
    const auto wit = weights.find(key);
    if (wit != weights.end()) {
      return wit->second;
    }
  }
  auto score = static_cast<real_t>(0.0);
  for (const real_t val : values) {
    score += std::abs(val);
  }
  return score;
}

[[nodiscard]] std::vector<RankedPartial>
rankPartials(const std::map<std::string, std::vector<real_t>> &partials,
             const std::map<std::string, real_t> &weights) {
  std::vector<RankedPartial> candidates;
  for (const auto &[key, values] : partials) {
    if (key == "Total" || key.starts_with("Frequency_")) {
      continue;
    }
    candidates.push_back({.key = key, .score = computePartialScore(key, values, weights)});
  }
  std::ranges::sort(candidates, [](const auto &lhs, const auto &rhs) {
    if (std::abs(lhs.score - rhs.score) > static_cast<real_t>(1e-12)) {
      return lhs.score > rhs.score;
    }
    return lhs.key < rhs.key;
  });
  return candidates;
}

void initDefaultVisibilities(const std::map<std::string, std::vector<real_t>> &partials,
                             const std::vector<RankedPartial> &ranked_partials,
                             std::map<std::string, bool> &curve_visibility_map) {
  if (partials.contains("Total") && !curve_visibility_map.contains("Total")) {
    curve_visibility_map["Total"] = true;
  }
  std::size_t rank = 0;
  for (const auto &partial : ranked_partials) {
    if (!curve_visibility_map.contains(partial.key)) {
      curve_visibility_map[partial.key] = (rank < 10);
    }
    rank++;
  }
}

[[nodiscard]] std::map<std::string, std::string>
assignStableCurveColors(const std::map<std::string, std::vector<real_t>> &partials,
                        const std::vector<RankedPartial> &ranked_partials,
                        const std::map<std::string, std::string> &custom_colors,
                        const correlation::plotters::PlotConfig &config) {
  std::map<std::string, std::string> assigned;
  if (partials.contains("Total")) {
    if (custom_colors.contains("Total") && !custom_colors.at("Total").empty()) {
      assigned["Total"] = custom_colors.at("Total");
    } else {
      assigned["Total"] =
          (config.theme == correlation::plotters::PlotConfig::Theme::Light) ? "#000000" : "#FFFFFF";
    }
  }
  const std::size_t count = ranked_partials.size();
  for (std::size_t i = 0; i < count; ++i) {
    const std::string &key = ranked_partials[i].key;
    if (custom_colors.contains(key) && !custom_colors.at(key).empty()) {
      assigned[key] = custom_colors.at(key);
    } else {
      assigned[key] = correlation::plotters::detail::color(i, count, config.palette);
    }
  }
  return assigned;
}

void appendTotalItem(slint::VectorModel<CurveToggleItem> &model, std::vector<std::string> &keys,
                     const std::map<std::string, std::string> &assigned_colors,
                     const std::map<std::string, bool> &visibility_map) {
  const auto iter = assigned_colors.find("Total");
  if (iter == assigned_colors.end()) {
    return;
  }
  const auto vis_it = visibility_map.find("Total");
  const bool visible = (vis_it != visibility_map.end()) ? vis_it->second : true;
  const auto curve_id = static_cast<int>(keys.size());
  keys.emplace_back("Total");
  model.push_back(CurveToggleItem{
      .id = curve_id,
      .label = slint::SharedString("Total"),
      .color_hex = parseHexColor(iter->second),
      .visible = visible,
      .is_partial = false,
      .is_pinned = false,
  });
}

void appendRankedPartials(slint::VectorModel<CurveToggleItem> &model,
                          std::vector<std::string> &keys,
                          const std::vector<RankedPartial> &ranked_partials,
                          const std::map<std::string, std::string> &assigned_colors,
                          const std::map<std::string, bool> &visibility_map, bool target_visible) {
  for (const auto &partial : ranked_partials) {
    const auto vis_it = visibility_map.find(partial.key);
    const bool visible = (vis_it != visibility_map.end()) ? vis_it->second : true;
    if (visible != target_visible) {
      continue;
    }
    const auto col_it = assigned_colors.find(partial.key);
    const std::string &color_hex = (col_it != assigned_colors.end()) ? col_it->second : "#808080";
    const auto curve_id = static_cast<int>(keys.size());
    keys.push_back(partial.key);
    model.push_back(CurveToggleItem{
        .id = curve_id,
        .label = slint::SharedString(partial.key),
        .color_hex = parseHexColor(color_hex),
        .visible = visible,
        .is_partial = true,
        .is_pinned = false,
    });
  }
}

void appendPinnedRuns(slint::VectorModel<CurveToggleItem> &model, std::vector<std::string> &keys,
                      const std::vector<PinnedRun> &pinned_runs,
                      std::map<std::string, bool> &visibility_map,
                      const std::map<std::string, std::string> &custom_colors,
                      const correlation::plotters::PlotConfig &config) {
  std::size_t pin_idx = 0;
  for (const auto &pinned_run : pinned_runs) {
    const std::string &pin_label = pinned_run.label;
    const bool vis = visibility_map.contains(pin_label) ? visibility_map[pin_label] : true;
    visibility_map[pin_label] = vis;

    std::string color_hex;
    if (custom_colors.contains(pin_label) && !custom_colors.at(pin_label).empty()) {
      color_hex = custom_colors.at(pin_label);
    } else {
      color_hex = correlation::plotters::detail::color(pin_idx + 1, config.palette);
    }

    const auto curve_id = static_cast<int>(keys.size());
    keys.push_back(pin_label);
    model.push_back(CurveToggleItem{
        .id = curve_id,
        .label = slint::SharedString(pin_label),
        .color_hex = parseHexColor(color_hex),
        .visible = vis,
        .is_partial = false,
        .is_pinned = true,
    });
    pin_idx++;
  }
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
  const auto iter = curve_visibility_map_.find(key);
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
  assigned_curve_colors_[key] = color_hex;
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
  assigned_curve_colors_.clear();
  current_toggle_keys_.clear();
  pinned_runs_.clear();
  show_difference_curve_ = false;
}

std::shared_ptr<slint::VectorModel<CurveToggleItem>>
PlotSeriesManager::generateToggleItems(const correlation::analysis::Histogram *hist,
                                       const correlation::plotters::PlotConfig &config,
                                       const std::map<std::string, real_t> &ashcroft_weights) {
  auto toggle_model = std::make_shared<slint::VectorModel<CurveToggleItem>>();
  current_toggle_keys_.clear();

  if (hist == nullptr) {
    return toggle_model;
  }

  const auto &partials = hist->smoothed_partials.empty() ? hist->partials : hist->smoothed_partials;

  const auto ranked_partials = rankPartials(partials, ashcroft_weights);
  initDefaultVisibilities(partials, ranked_partials, curve_visibility_map_);

  assigned_curve_colors_ =
      assignStableCurveColors(partials, ranked_partials, custom_curve_colors_, config);

  appendTotalItem(*toggle_model, current_toggle_keys_, assigned_curve_colors_,
                  curve_visibility_map_);
  appendRankedPartials(*toggle_model, current_toggle_keys_, ranked_partials, assigned_curve_colors_,
                       curve_visibility_map_, true);
  appendRankedPartials(*toggle_model, current_toggle_keys_, ranked_partials, assigned_curve_colors_,
                       curve_visibility_map_, false);

  appendPinnedRuns(*toggle_model, current_toggle_keys_, pinned_runs_, curve_visibility_map_,
                   custom_curve_colors_, config);

  return toggle_model;
}

} // namespace correlation::app
