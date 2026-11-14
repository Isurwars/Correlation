// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "render3d/MiniPlotRenderer.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace correlation::render3d {

namespace {

struct Bounds2D {
  real_t x_min{0.0};
  real_t x_max{1.0};
  real_t y_min{0.0};
  real_t y_max{1.0};
};

struct ViewportBox {
  float left{12.0F};
  float right{248.0F};
  float top{16.0F};
  float bottom{114.0F};
};

void blendPixel(PlotImageBuffer &buf, int pixel_x, int pixel_y, RgbaPixel color,
                float alpha_mult) noexcept {
  if (pixel_x < 0 || std::cmp_greater_equal(pixel_x, buf.width) || pixel_y < 0 ||
      std::cmp_greater_equal(pixel_y, buf.height)) {
    return;
  }
  const auto idx = static_cast<std::size_t>(pixel_y) * static_cast<std::size_t>(buf.width) +
                   static_cast<std::size_t>(pixel_x);
  auto &dst = buf.pixels[idx];

  const float src_a = (static_cast<float>(color.a) / 255.0F) * std::clamp(alpha_mult, 0.0F, 1.0F);
  const float inv_a = 1.0F - src_a;

  dst.r =
      static_cast<uint8_t>(static_cast<float>(color.r) * src_a + static_cast<float>(dst.r) * inv_a);
  dst.g =
      static_cast<uint8_t>(static_cast<float>(color.g) * src_a + static_cast<float>(dst.g) * inv_a);
  dst.b =
      static_cast<uint8_t>(static_cast<float>(color.b) * src_a + static_cast<float>(dst.b) * inv_a);
  dst.a = static_cast<uint8_t>(
      std::clamp((src_a + (static_cast<float>(dst.a) / 255.0F) * inv_a) * 255.0F, 0.0F, 255.0F));
}

void clearCanvas(PlotImageBuffer &buf, RgbaPixel bg_color) {
  std::ranges::fill(buf.pixels, bg_color);
}

void drawLineAA(PlotImageBuffer &buf, float start_x, float start_y, float end_x, float end_y,
                RgbaPixel color, float thickness) {
  const float delta_x = end_x - start_x;
  const float delta_y = end_y - start_y;
  const float len = std::sqrt(delta_x * delta_x + delta_y * delta_y);
  if (len < 0.001F) {
    blendPixel(buf, static_cast<int>(std::round(start_x)), static_cast<int>(std::round(start_y)),
               color, 1.0F);
    return;
  }

  const int steps = static_cast<int>(std::ceil(len * 2.0F));
  const float step_x = delta_x / static_cast<float>(steps);
  const float step_y = delta_y / static_cast<float>(steps);

  for (int step_idx = 0; step_idx <= steps; ++step_idx) {
    const float center_x = start_x + static_cast<float>(step_idx) * step_x;
    const float center_y = start_y + static_cast<float>(step_idx) * step_y;

    const int min_x = static_cast<int>(std::floor(center_x - thickness));
    const int max_x = static_cast<int>(std::ceil(center_x + thickness));
    const int min_y = static_cast<int>(std::floor(center_y - thickness));
    const int max_y = static_cast<int>(std::ceil(center_y + thickness));

    for (int pixel_y = min_y; pixel_y <= max_y; ++pixel_y) {
      for (int pixel_x = min_x; pixel_x <= max_x; ++pixel_x) {
        const float dist = std::sqrt(
            (static_cast<float>(pixel_x) - center_x) * (static_cast<float>(pixel_x) - center_x) +
            (static_cast<float>(pixel_y) - center_y) * (static_cast<float>(pixel_y) - center_y));
        if (dist <= thickness) {
          const float alpha = 1.0F - (dist / thickness);
          blendPixel(buf, pixel_x, pixel_y, color, alpha);
        }
      }
    }
  }
}

void drawGridLines(PlotImageBuffer &buf, const ViewportBox &box) {
  const RgbaPixel grid_color{.r = 50, .g = 56, .b = 66, .a = 160};
  const RgbaPixel base_color{.r = 75, .g = 85, .b = 100, .a = 220};

  // Baseline at bottom
  drawLineAA(buf, box.left, box.bottom, box.right, box.bottom, base_color, 1.0F);

  // Midline
  const float mid_y = (box.top + box.bottom) * 0.5F;
  drawLineAA(buf, box.left, mid_y, box.right, mid_y, grid_color, 0.8F);
}

Bounds2D computeBounds(std::span<const real_t> bins, std::span<const real_t> mean,
                       std::span<const real_t> frame, const MiniPlotConfig &config) {
  Bounds2D bounds{};
  if (!bins.empty()) {
    bounds.x_min = bins.front();
    bounds.x_max = bins.back();
  }
  if (config.manual_x_min.has_value()) {
    bounds.x_min = *config.manual_x_min;
  }
  if (config.manual_x_max.has_value()) {
    bounds.x_max = *config.manual_x_max;
  }

  real_t max_val = 0.0;
  for (const auto &val : mean) {
    max_val = std::max(max_val, val);
  }
  if (mean.empty()) {
    for (const auto &val : frame) {
      max_val = std::max(max_val, val);
    }
  }

  bounds.y_min = static_cast<real_t>(0.0);
  bounds.y_max = (max_val > static_cast<real_t>(0.0)) ? (max_val * static_cast<real_t>(1.15))
                                                      : static_cast<real_t>(1.0);
  if (config.manual_y_min.has_value()) {
    bounds.y_min = *config.manual_y_min;
  }
  if (config.manual_y_max.has_value()) {
    bounds.y_max = *config.manual_y_max;
  }
  return bounds;
}

float mapX(real_t val, const Bounds2D &bounds, const ViewportBox &box) {
  const auto denom =
      (bounds.x_max > bounds.x_min) ? (bounds.x_max - bounds.x_min) : static_cast<real_t>(1.0);
  const auto norm = static_cast<float>((val - bounds.x_min) / denom);
  return box.left + norm * (box.right - box.left);
}

float mapY(real_t val, const Bounds2D &bounds, const ViewportBox &box, bool log_scale) {
  if (log_scale) {
    constexpr auto K_FLOOR = static_cast<real_t>(1e-4);
    const real_t log_val = std::log10(std::max(val, K_FLOOR));
    const real_t log_min = std::log10(std::max(bounds.y_min, K_FLOOR));
    const real_t log_max = std::log10(std::max(bounds.y_max, K_FLOOR));
    const real_t denom = (log_max > log_min) ? (log_max - log_min) : static_cast<real_t>(1.0);
    const auto norm = static_cast<float>((log_val - log_min) / denom);
    return box.bottom - norm * (box.bottom - box.top);
  }

  const auto denom =
      (bounds.y_max > bounds.y_min) ? (bounds.y_max - bounds.y_min) : static_cast<real_t>(1.0);
  const auto norm = static_cast<float>((val - bounds.y_min) / denom);
  return box.bottom - norm * (box.bottom - box.top);
}

void renderPolyline(PlotImageBuffer &buf, std::span<const real_t> bins,
                    std::span<const real_t> values, const Bounds2D &bounds, const ViewportBox &box,
                    bool log_scale, RgbaPixel color, float thickness) {
  const std::size_t point_count = std::min(bins.size(), values.size());
  if (point_count < 2) {
    return;
  }

  for (std::size_t i = 0; i + 1 < point_count; ++i) {
    const float pt0_x = mapX(bins[i], bounds, box);
    const float pt0_y = mapY(values[i], bounds, box, log_scale);
    const float pt1_x = mapX(bins[i + 1], bounds, box);
    const float pt1_y = mapY(values[i + 1], bounds, box, log_scale);
    drawLineAA(buf, pt0_x, pt0_y, pt1_x, pt1_y, color, thickness);
  }
}

void drawPeakMarker(PlotImageBuffer &buf, const PeakInfo &peak, const Bounds2D &bounds,
                    const ViewportBox &box, bool log_scale) {
  if (!peak.valid) {
    return;
  }
  const float peak_x = mapX(peak.x, bounds, box);
  const float peak_y = mapY(peak.y, bounds, box, log_scale);

  const RgbaPixel marker_color{.r = 255, .g = 215, .b = 0, .a = 255};
  const RgbaPixel halo_color{.r = 255, .g = 215, .b = 0, .a = 90};

  // Outer halo & center dot
  for (int offset_y = -4; offset_y <= 4; ++offset_y) {
    for (int offset_x = -4; offset_x <= 4; ++offset_x) {
      const float distance =
          std::sqrt(static_cast<float>(offset_x * offset_x + offset_y * offset_y));
      if (distance <= 2.0F) {
        blendPixel(buf, static_cast<int>(peak_x) + offset_x, static_cast<int>(peak_y) + offset_y,
                   marker_color, 1.0F);
      } else if (distance <= 4.0F) {
        blendPixel(buf, static_cast<int>(peak_x) + offset_x, static_cast<int>(peak_y) + offset_y,
                   halo_color, 0.6F);
      }
    }
  }
}

} // namespace

PeakInfo MiniPlotRenderer::findPeak(std::span<const real_t> bins,
                                    std::span<const real_t> values) noexcept {
  const std::size_t point_count = std::min(bins.size(), values.size());
  if (point_count < 3) {
    return {.x = 0.0, .y = 0.0, .valid = false};
  }

  real_t max_val = values[0];
  std::size_t peak_idx = 0;
  for (std::size_t i = 1; i < point_count; ++i) {
    if (values[i] > max_val) {
      max_val = values[i];
      peak_idx = i;
    }
  }

  if (max_val <= static_cast<real_t>(0.001)) {
    return {.x = 0.0, .y = 0.0, .valid = false};
  }

  return {.x = bins[peak_idx], .y = max_val, .valid = true};
}

PlotImageBuffer MiniPlotRenderer::render(const analysis::MetricCurveData &metric_data,
                                         std::optional<std::size_t> frame_idx,
                                         const MiniPlotConfig &config, uint32_t width,
                                         uint32_t height) {
  PlotImageBuffer buffer{
      .width = width,
      .height = height,
      .pixels = std::vector<RgbaPixel>(static_cast<std::size_t>(width) *
                                       static_cast<std::size_t>(height)),
  };

  const RgbaPixel bg_color{.r = 18, .g = 21, .b = 26, .a = 230};
  clearCanvas(buffer, bg_color);

  const ViewportBox box{
      .left = 14.0F,
      .right = static_cast<float>(width) - 14.0F,
      .top = 16.0F,
      .bottom = static_cast<float>(height) - 16.0F,
  };

  drawGridLines(buffer, box);

  const auto &bins = metric_data.bins;
  const auto &mean = metric_data.mean_curve;

  std::span<const real_t> frame_span;
  if (frame_idx.has_value() && *frame_idx < metric_data.frame_curves.size()) {
    frame_span = metric_data.frame_curves[*frame_idx];
  }

  if (bins.empty()) {
    return buffer;
  }

  const Bounds2D bounds = computeBounds(bins, mean, frame_span, config);

  // 1. Draw Ensemble Mean Curve (Muted Cyan/Gray, dashed effect through thinner width)
  if (!mean.empty()) {
    const RgbaPixel mean_color{.r = 135, .g = 155, .b = 175, .a = 150};
    renderPolyline(buffer, bins, mean, bounds, box, config.y_log_scale, mean_color, 1.0F);
  }

  // 2. Draw Active Frame Curve (Vibrant Teal / Primary Accent)
  if (!frame_span.empty()) {
    const RgbaPixel frame_color{.r = 64, .g = 215, .b = 215, .a = 255};
    renderPolyline(buffer, bins, frame_span, bounds, box, config.y_log_scale, frame_color, 1.8F);

    const PeakInfo peak = findPeak(bins, frame_span);
    drawPeakMarker(buffer, peak, bounds, box, config.y_log_scale);
  } else if (!mean.empty()) {
    const PeakInfo peak = findPeak(bins, mean);
    drawPeakMarker(buffer, peak, bounds, box, config.y_log_scale);
  }

  return buffer;
}

} // namespace correlation::render3d
