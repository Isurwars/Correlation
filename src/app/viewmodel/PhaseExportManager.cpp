/**
 * @file PhaseExportManager.cpp
 * @brief Implementation of PhaseExportManager for rasterizing and exporting 2D phase heatmaps.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/viewmodel/PhaseExportManager.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <miniz.h>
#include <vector>

namespace correlation::app {

namespace {

constexpr std::size_t K_HEATMAP_GRID_BINS = 64;
constexpr std::size_t K_SLINT_PREVIEW_DIMS = 128;
constexpr std::size_t K_PNG_EXPORT_DIMS = 256;

} // namespace

slint::Image PhaseExportManager::renderHeatmapImage(const PhaseHeatmap2D &heatmap,
                                                    plotters::ColormapType colormap) {
  std::vector<slint::Rgba8Pixel> pixels(K_SLINT_PREVIEW_DIMS * K_SLINT_PREVIEW_DIMS);

  const double log_max =
      (heatmap.max_count > 0) ? std::log1p(static_cast<double>(heatmap.max_count)) : 1.0;

  for (std::size_t py = 0; py < K_SLINT_PREVIEW_DIMS; ++py) {
    // Invert Y so 0 is at bottom
    const std::size_t bin_y =
        (K_SLINT_PREVIEW_DIMS - 1 - py) * K_HEATMAP_GRID_BINS / K_SLINT_PREVIEW_DIMS;
    for (std::size_t px = 0; px < K_SLINT_PREVIEW_DIMS; ++px) {
      const std::size_t bin_x = px * K_HEATMAP_GRID_BINS / K_SLINT_PREVIEW_DIMS;
      const std::size_t cell_idx = bin_y * K_HEATMAP_GRID_BINS + bin_x;
      const std::uint32_t count =
          (cell_idx < heatmap.counts.size()) ? heatmap.counts[cell_idx] : 0U;

      const std::size_t pixel_idx = py * K_SLINT_PREVIEW_DIMS + px;
      if (count == 0) {
        pixels[pixel_idx] = slint::Rgba8Pixel{.r = 20, .g = 22, .b = 28, .a = 255};
      } else {
        const double log_c = std::log1p(static_cast<double>(count));
        const double norm = std::clamp(log_c / log_max, 0.0, 1.0);
        const auto color = plotters::ColormapService::sample(norm, colormap);
        pixels[pixel_idx] = slint::Rgba8Pixel{.r = color.r, .g = color.g, .b = color.b, .a = 255};
      }
    }
  }

  const slint::SharedPixelBuffer<slint::Rgba8Pixel> buffer(
      static_cast<std::uint32_t>(K_SLINT_PREVIEW_DIMS),
      static_cast<std::uint32_t>(K_SLINT_PREVIEW_DIMS), pixels.data());
  return slint::Image{buffer};
}

bool PhaseExportManager::exportHeatmapPng(const std::string &filepath,
                                          const PhaseHeatmap2D &heatmap,
                                          plotters::ColormapType colormap) {
  const std::string out_path = filepath.empty() ? "phase_heatmap_2d.png" : filepath;
  std::vector<std::uint8_t> raw_rgba(K_PNG_EXPORT_DIMS * K_PNG_EXPORT_DIMS * 4);

  const double log_max =
      (heatmap.max_count > 0) ? std::log1p(static_cast<double>(heatmap.max_count)) : 1.0;

  for (std::size_t py = 0; py < K_PNG_EXPORT_DIMS; ++py) {
    const std::size_t bin_y =
        (K_PNG_EXPORT_DIMS - 1 - py) * K_HEATMAP_GRID_BINS / K_PNG_EXPORT_DIMS;
    for (std::size_t px = 0; px < K_PNG_EXPORT_DIMS; ++px) {
      const std::size_t bin_x = px * K_HEATMAP_GRID_BINS / K_PNG_EXPORT_DIMS;
      const std::size_t cell_idx = bin_y * K_HEATMAP_GRID_BINS + bin_x;
      const std::uint32_t count =
          (cell_idx < heatmap.counts.size()) ? heatmap.counts[cell_idx] : 0U;

      const std::size_t offset = (py * K_PNG_EXPORT_DIMS + px) * 4;
      if (count == 0) {
        raw_rgba[offset + 0] = 20;
        raw_rgba[offset + 1] = 22;
        raw_rgba[offset + 2] = 28;
        raw_rgba[offset + 3] = 255;
      } else {
        const double log_c = std::log1p(static_cast<double>(count));
        const double norm = std::clamp(log_c / log_max, 0.0, 1.0);
        const auto color = plotters::ColormapService::sample(norm, colormap);
        raw_rgba[offset + 0] = color.r;
        raw_rgba[offset + 1] = color.g;
        raw_rgba[offset + 2] = color.b;
        raw_rgba[offset + 3] = 255;
      }
    }
  }

  std::size_t png_len = 0;
  void *png_mem =
      tdefl_write_image_to_png_file_in_memory(raw_rgba.data(), static_cast<int>(K_PNG_EXPORT_DIMS),
                                              static_cast<int>(K_PNG_EXPORT_DIMS), 4, &png_len);
  if (png_mem == nullptr || png_len == 0) {
    return false;
  }

  std::ofstream out(out_path, std::ios::binary);
  if (!out.is_open()) {
    mz_free(png_mem);
    return false;
  }

  out.write(reinterpret_cast<const char *>(png_mem), static_cast<std::streamsize>(png_len));
  mz_free(png_mem);
  return out.good();
}

} // namespace correlation::app
