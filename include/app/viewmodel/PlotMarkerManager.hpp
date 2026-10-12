/**
 * @file PlotMarkerManager.hpp
 * @brief Management and formatting of user reference marker guidelines on plots.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "plotters/SvgPlotter.hpp"

#include <span>
#include <vector>

namespace correlation::app {

/**
 * @class PlotMarkerManager
 * @brief Stores and formats vertical (+X) and horizontal (+Y) reference inspection markers.
 */
class PlotMarkerManager {
public:
  /**
   * @brief Constructs an empty PlotMarkerManager.
   */
  PlotMarkerManager() = default;

  /**
   * @brief Adds a reference guideline marker at the specified data coordinate.
   * @param[in] value Data coordinate value.
   * @param[in] is_vertical True for vertical (+X) marker, false for horizontal (+Y) marker.
   */
  void addMarker(real_t value, bool is_vertical);

  /**
   * @brief Clears all inspection reference markers.
   */
  void clearMarkers() noexcept;

  /**
   * @brief Returns the total number of placed reference markers.
   * @return Count of active markers.
   */
  [[nodiscard]] std::size_t markersCount() const noexcept;

  /**
   * @brief Returns a read-only view of all active reference markers.
   * @return Contiguous span of ReferenceLine structures.
   */
  [[nodiscard]] std::span<const correlation::plotters::ReferenceLine> markers() const noexcept;

  /**
   * @brief Populates reference markers into target PlotConfig.
   * @param[in,out] config Plot configuration to enrich.
   */
  void applyToConfig(correlation::plotters::PlotConfig &config) const;

private:
  std::vector<correlation::plotters::ReferenceLine> markers_;
};

} // namespace correlation::app
