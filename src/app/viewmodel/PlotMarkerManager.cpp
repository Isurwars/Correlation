/**
 * @file PlotMarkerManager.cpp
 * @brief Implementation of PlotMarkerManager.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/viewmodel/PlotMarkerManager.hpp"

#include <format>
#include <utility>

namespace correlation::app {

void PlotMarkerManager::addMarker(real_t value, bool is_vertical) {
  std::string color = is_vertical ? "#0284C7" : "#E11D48";
  std::string label = std::format("{}: {:.3f}", is_vertical ? "X" : "Y", value);

  markers_.push_back(correlation::plotters::ReferenceLine{
      .value = value,
      .is_vertical = is_vertical,
      .label = std::move(label),
      .color_hex = std::move(color),
  });
}

void PlotMarkerManager::clearMarkers() noexcept { markers_.clear(); }

std::size_t PlotMarkerManager::markersCount() const noexcept { return markers_.size(); }

std::span<const correlation::plotters::ReferenceLine> PlotMarkerManager::markers() const noexcept {
  return markers_;
}

void PlotMarkerManager::applyToConfig(correlation::plotters::PlotConfig &config) const {
  config.reference_lines = markers_;
}

} // namespace correlation::app
