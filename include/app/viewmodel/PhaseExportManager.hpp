/**
 * @file PhaseExportManager.hpp
 * @brief Coordinates rasterization and file export for 2D phase density heatmaps.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "app/viewmodel/PhaseSearchTypes.hpp"
#include "plotters/ColormapService.hpp"
#include <slint.h>

#include <cstdint>
#include <string>

namespace correlation::app {

/**
 * @class PhaseExportManager
 * @brief Utility manager for generating Slint heatmap textures and exporting high-resolution PNGs.
 */
class PhaseExportManager {
public:
  PhaseExportManager() = delete;

  /**
   * @brief Rasterizes a 2D phase distribution heatmap to a Slint display image.
   * @param[in] heatmap 64x64 scalar joint distribution density map.
   * @param[in] colormap Scientific colormap gradient for intensity mapping.
   * @return slint::Image containing 128x128 RGBA8 texture.
   */
  [[nodiscard]] static slint::Image renderHeatmapImage(const PhaseHeatmap2D &heatmap,
                                                       plotters::ColormapType colormap);

  /**
   * @brief Encodes and exports a 2D phase distribution heatmap as a PNG file.
   * @param[in] filepath Target destination filepath on disk.
   * @param[in] heatmap 64x64 scalar joint distribution density map.
   * @param[in] colormap Scientific colormap gradient for intensity mapping.
   * @return True if PNG file was written successfully, false on I/O or encoding error.
   */
  [[nodiscard]] static bool exportHeatmapPng(const std::string &filepath,
                                             const PhaseHeatmap2D &heatmap,
                                             plotters::ColormapType colormap);
};

} // namespace correlation::app
