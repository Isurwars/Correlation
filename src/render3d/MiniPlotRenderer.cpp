// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "render3d/MiniPlotRenderer.hpp"

#include <algorithm>
#include <cmath>

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

void blendPixel(PlotImageBuffer &buf, int x, int y, RgbaPixel color, float alpha_mult) noexcept {
  if (x < 0 || x >= static_cast<int>(buf.width) || y < 0 || y >= static_cast<int>(buf.height)) {
    return;
  }
  const auto idx = static_cast<std::size_t>(y * static_cast<int>(buf.width) + x);
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
  std::fill(buf.pixels.begin(), buf.pixels.end(), bg_color);
}

void drawLineAA(PlotImageBuffer &buf, float x0, float y0, float x1, float y1, RgbaPixel color,
                float thickness) {
  const float dx = x1 - x0;
  const float dy = y1 - y0;
  const float len = std::sqrt(dx * dx + dy * dy);
  if (len < 0.001F) {
    blendPixel(buf, static_cast<int>(std::round(x0)), static_cast<int>(std::round(y0)), color,
               1.0F);
    return;
  }

  const int steps = static_cast<int>(std::ceil(len * 2.0F));
  const float step_x = dx / static_cast<float>(steps);
  const float step_y = dy / static_cast<float>(steps);

  for (int s = 0; s <= steps; ++s) {
    const float cx = x0 + static_cast<float>(s) * step_x;
    const float cy = y0 + static_cast<float>(s) * step_y;

    const int min_x = static_cast<int>(std::floor(cx - thickness));
    const int max_x = static_cast<int>(std::ceil(cx + thickness));
    const int min_y = static_cast<int>(std::floor(cy - thickness));
    const int max_y = static_cast<int>(std::ceil(cy + thickness));

    for (int py = min_y; py <= max_y; ++py) {
      for (int px = min_x; px <= max_x; ++px) {
        const float dist = std::sqrt((static_cast<float>(px) - cx) * (static_cast<float>(px) - cx) +
                                     (static_cast<float>(py) - cy) * (static_cast<float>(py) - cy));
        if (dist <= thickness) {
          const float alpha = 1.0F - (dist / thickness);
          blendPixel(buf, px, py, color, alpha);
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
  for (const auto &v : mean) {
    max_val = std::max(max_val, v);
  }
  for (const auto &v : frame) {
    max_val = std::max(max_val, v);
  }

  bounds.y_min = 0.0;
  bounds.y_max = (max_val > 0.0) ? (max_val * 1.15) : 1.0;
  if (config.manual_y_min.has_value()) {
    bounds.y_min = *config.manual_y_min;
  }
  if (config.manual_y_max.has_value()) {
    bounds.y_max = *config.manual_y_max;
  }
  return bounds;
}

float mapX(real_t val, const Bounds2D &b, const ViewportBox &box) {
  const auto denom = (b.x_max > b.x_min) ? (b.x_max - b.x_min) : static_cast<real_t>(1.0);
  const auto norm = static_cast<float>((val - b.x_min) / denom);
  return box.left + norm * (box.right - box.left);
}

float mapY(real_t val, const Bounds2D &b, const ViewportBox &box, bool log_scale) {
  if (log_scale) {
    constexpr real_t kFloor = 1e-4;
    const real_t log_val = std::log10(std::max(val, kFloor));
    const real_t log_min = std::log10(std::max(b.y_min, kFloor));
    const real_t log_max = std::log10(std::max(b.y_max, kFloor));
    const real_t denom = (log_max > log_min) ? (log_max - log_min) : static_cast<real_t>(1.0);
    const auto norm = static_cast<float>((log_val - log_min) / denom);
    return box.bottom - norm * (box.bottom - box.top);
  }

  const auto denom = (b.y_max > b.y_min) ? (b.y_max - b.y_min) : static_cast<real_t>(1.0);
  const auto norm = static_cast<float>((val - b.y_min) / denom);
  return box.bottom - norm * (box.bottom - box.top);
}

void renderPolyline(PlotImageBuffer &buf, std::span<const real_t> bins,
                    std::span<const real_t> values, const Bounds2D &bounds, const ViewportBox &box,
                    bool log_scale, RgbaPixel color, float thickness) {
  const std::size_t n = std::min(bins.size(), values.size());
  if (n < 2) {
    return;
  }

  for (std::size_t i = 0; i + 1 < n; ++i) {
    const float x0 = mapX(bins[i], bounds, box);
    const float y0 = mapY(values[i], bounds, box, log_scale);
    const float x1 = mapX(bins[i + 1], bounds, box);
    const float y1 = mapY(values[i + 1], bounds, box, log_scale);
    drawLineAA(buf, x0, y0, x1, y1, color, thickness);
  }
}

void drawPeakMarker(PlotImageBuffer &buf, const PeakInfo &peak, const Bounds2D &bounds,
                    const ViewportBox &box, bool log_scale) {
  if (!peak.valid) {
    return;
  }
  const float px = mapX(peak.x, bounds, box);
  const float py = mapY(peak.y, bounds, box, log_scale);

  const RgbaPixel marker_color{.r = 255, .g = 215, .b = 0, .a = 255};
  const RgbaPixel halo_color{.r = 255, .g = 215, .b = 0, .a = 90};

  // Outer halo & center dot
  for (int dy = -4; dy <= 4; ++dy) {
    for (int dx = -4; dx <= 4; ++dx) {
      const float d = std::sqrt(static_cast<float>(dx * dx + dy * dy));
      if (d <= 2.0F) {
        blendPixel(buf, static_cast<int>(px) + dx, static_cast<int>(py) + dy, marker_color, 1.0F);
      } else if (d <= 4.0F) {
        blendPixel(buf, static_cast<int>(px) + dx, static_cast<int>(py) + dy, halo_color, 0.6F);
      }
    }
  }
}

} // namespace

PeakInfo MiniPlotRenderer::findPeak(std::span<const real_t> bins,
                                    std::span<const real_t> values) noexcept {
  const std::size_t n = std::min(bins.size(), values.size());
  if (n < 3) {
    return {.x = 0.0, .y = 0.0, .valid = false};
  }

  real_t max_val = values[0];
  std::size_t peak_idx = 0;
  for (std::size_t i = 1; i < n; ++i) {
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
  PlotImageBuffer buffer{.width = width,
                         .height = height,
                         .pixels =
                             std::vector<RgbaPixel>(static_cast<std::size_t>(width * height))};

  const RgbaPixel bg_color{.r = 18, .g = 21, .b = 26, .a = 230};
  clearCanvas(buffer, bg_color);

  const ViewportBox box{.left = 14.0F,
                        .right = static_cast<float>(width) - 14.0F,
                        .top = 16.0F,
                        .bottom = static_cast<float>(height) - 16.0F};

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
