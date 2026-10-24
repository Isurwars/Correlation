/**
 * @file SvgHistogramRenderer.cpp
 * @brief Implementation of single-histogram vector rendering for the SVG plotter.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "plotters/SvgHistogramRenderer.hpp"

#include "analysis/DistributionFunctions.hpp"
#include "plotters/PathFont.hpp"
#include "plotters/PlotTypes.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <format>
#include <iterator>
#include <map>
#include <optional>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace correlation::plotters {

namespace {

using namespace detail;

// ============================================================================
// Layout Metrics & Constants
// ============================================================================
namespace layout {
inline constexpr real_t MARGIN_LEFT = static_cast<real_t>(100.0);
inline constexpr real_t MARGIN_RIGHT = static_cast<real_t>(40.0);
inline constexpr real_t MARGIN_TOP = static_cast<real_t>(50.0);
inline constexpr real_t MARGIN_BOTTOM = static_cast<real_t>(90.0);

inline constexpr real_t MINOR_TICK_LENGTH = static_cast<real_t>(4.0);
inline constexpr real_t MAJOR_TICK_LENGTH = static_cast<real_t>(8.0);
inline constexpr real_t Y_TICK_LABEL_OFFSET = static_cast<real_t>(15.0);
inline constexpr real_t Y_TICK_LABEL_Y_SHIFT = static_cast<real_t>(8.0);
inline constexpr real_t X_TICK_LABEL_OFFSET = static_cast<real_t>(25.0);
inline constexpr real_t X_AXIS_TITLE_OFFSET = static_cast<real_t>(70.0);
inline constexpr real_t Y_AXIS_TITLE_OFFSET = static_cast<real_t>(35.0);

inline constexpr real_t LEGEND_RIGHT_MARGIN = static_cast<real_t>(15.0);
inline constexpr real_t LEGEND_TOP_MARGIN = static_cast<real_t>(25.0);
inline constexpr real_t LEGEND_LINE_RIGHT_PAD = static_cast<real_t>(10.0);
inline constexpr real_t LEGEND_LINE_LEFT_PAD = static_cast<real_t>(40.0);
inline constexpr real_t LEGEND_TEXT_PAD = static_cast<real_t>(45.0);
inline constexpr real_t LEGEND_TOTAL_STEP = static_cast<real_t>(26.0);
inline constexpr real_t LEGEND_HEADER_STEP = static_cast<real_t>(22.0);
inline constexpr real_t LEGEND_ITEM_STEP = static_cast<real_t>(24.0);

inline constexpr real_t SCORE_EPSILON = static_cast<real_t>(1e-12);
inline constexpr std::size_t DEFAULT_PARTIALS_LIMIT = 10;
inline constexpr std::size_t DEFAULT_BUFFER_CAPACITY = 32768;
} // namespace layout

// ============================================================================
// Value Types & Render Context
// ============================================================================

/**
 * @struct RenderContext
 * @brief Aggregate of non-owning pointers passed across pipeline stages.
 */
struct RenderContext {
  const correlation::analysis::Histogram *hist{nullptr};
  const PlotConfig *config{nullptr};
  const HoverInfo *hover{nullptr};
  const std::map<std::string, real_t> *weights{nullptr};
  const std::map<std::string, bool> *curve_visibility{nullptr};
  const std::map<std::string, std::string> *custom_colors{nullptr};
};

/**
 * @struct RenderSeries
 * @brief Processed single curve series with assigned styling color.
 */
struct RenderSeries {
  std::string key;
  std::vector<real_t> values;
  std::string color;
};

/**
 * @struct PeakCallout
 * @brief Formatted peak annotation ready for visual projection.
 */
struct PeakCallout {
  std::string label;
  real_t x{0.0};
  real_t y{0.0};
};

/**
 * @struct CanvasViewport
 * @brief Encapsulates SVG dimension bounds and coordinate projection math.
 */
struct CanvasViewport {
  real_t width{0.0};
  real_t height{0.0};
  real_t px0{0.0};
  real_t px1{0.0};
  real_t py0{0.0};
  real_t py1{0.0};

  static CanvasViewport create(const PlotConfig &config) noexcept {
    const real_t canvas_w = config.effectiveWidth();
    const real_t canvas_h = config.effectiveHeight();
    return CanvasViewport{
        .width = canvas_w,
        .height = canvas_h,
        .px0 = layout::MARGIN_LEFT,
        .px1 = canvas_w - layout::MARGIN_RIGHT,
        .py0 = layout::MARGIN_TOP,
        .py1 = canvas_h - layout::MARGIN_BOTTOM,
    };
  }

  [[nodiscard]] constexpr real_t plotWidth() const noexcept { return px1 - px0; }
  [[nodiscard]] constexpr real_t plotHeight() const noexcept { return py1 - py0; }

  [[nodiscard]] real_t mapX(real_t val, const detail::NiceScale &x_scale) const noexcept {
    return mapValue(val, x_scale.min, x_scale.max, px0, px1);
  }

  [[nodiscard]] real_t mapY(real_t val, const detail::NiceScale &y_scale) const noexcept {
    return mapValue(val, y_scale.min, y_scale.max, py1, py0);
  }

  [[nodiscard]] bool containsPixel(real_t pixel_x, real_t pixel_y) const noexcept {
    return (pixel_x >= px0 && pixel_x <= px1 && pixel_y >= py0 && pixel_y <= py1);
  }
};

// ============================================================================
// Stage 1: Data Preparation & Filtering (SeriesFilter)
// ============================================================================
class SeriesFilter {
public:
  [[nodiscard]] static bool isCurveVisible(const std::string &key,
                                           const std::map<std::string, bool> *visibility) noexcept {
    if (visibility == nullptr) {
      return true;
    }
    const auto found_it = visibility->find(key);
    return (found_it != visibility->end()) ? found_it->second : true;
  }

  [[nodiscard]] static real_t
  calculateScore(const std::string &key, const std::vector<real_t> &values,
                 const std::map<std::string, real_t> *weights) noexcept {
    if (weights != nullptr && !weights->empty()) {
      const auto wit = weights->find(key);
      if (wit != weights->end()) {
        return wit->second;
      }
    }
    real_t score = 0.0;
    for (const real_t val : values) {
      score += std::abs(val);
    }
    return score;
  }

  [[nodiscard]] static std::string
  resolveColor(const std::string &key, std::size_t partial_idx, std::size_t num_partials,
               const std::map<std::string, std::string> *custom_colors, const PlotConfig &config) {
    if (custom_colors != nullptr) {
      const auto itx = custom_colors->find(key);
      if (itx != custom_colors->end() && !itx->second.empty()) {
        return itx->second;
      }
    }
    if (key == "Total") {
      return (config.theme == PlotConfig::Theme::Light) ? "#000000" : "#FFFFFF";
    }
    return color(partial_idx, num_partials, config.palette);
  }

  [[nodiscard]] static std::vector<RenderSeries> filterSeries(const RenderContext &ctx) {
    if (ctx.hist == nullptr || ctx.config == nullptr) {
      return {};
    }
    const auto &raw_partials =
        ctx.hist->smoothed_partials.empty() ? ctx.hist->partials : ctx.hist->smoothed_partials;
    std::vector<std::pair<std::string, std::vector<real_t>>> selected;

    const auto total_it = raw_partials.find("Total");
    if (total_it != raw_partials.end() && isCurveVisible("Total", ctx.curve_visibility)) {
      selected.emplace_back("Total", total_it->second);
    }

    std::size_t non_total_count = 0;
    for (const auto &[key, ignored_val] : raw_partials) {
      if (key != "Total" && !key.starts_with("Frequency_")) {
        ++non_total_count;
      }
    }

    // Single partial dataset with Total renders only Total
    if (total_it != raw_partials.end() && non_total_count == 1) {
      return buildStyledSeries(selected, ctx);
    }

    std::vector<std::pair<std::string, real_t>> candidates;
    for (const auto &[key, value] : raw_partials) {
      if (key == "Total" || key.starts_with("Frequency_")) {
        continue;
      }
      if (!isCurveVisible(key, ctx.curve_visibility)) {
        continue;
      }
      candidates.emplace_back(key, calculateScore(key, value, ctx.weights));
    }

    std::ranges::sort(candidates, [](const auto &lhs, const auto &rhs) {
      if (std::abs(lhs.second - rhs.second) > layout::SCORE_EPSILON) {
        return lhs.second > rhs.second;
      }
      return lhs.first < rhs.first;
    });

    const bool has_custom_vis = ctx.curve_visibility != nullptr && !ctx.curve_visibility->empty();
    const std::size_t limit = has_custom_vis
                                  ? candidates.size()
                                  : std::min(candidates.size(), layout::DEFAULT_PARTIALS_LIMIT);

    for (std::size_t i = 0; i < limit; ++i) {
      const std::string &key = candidates[i].first;
      const auto iter = raw_partials.find(key);
      if (iter != raw_partials.end()) {
        selected.emplace_back(key, iter->second);
      }
    }

    return buildStyledSeries(selected, ctx);
  }

private:
  [[nodiscard]] static std::vector<RenderSeries>
  buildStyledSeries(const std::vector<std::pair<std::string, std::vector<real_t>>> &raw,
                    const RenderContext &ctx) {
    std::size_t non_total_count = 0;
    for (const auto &[key, ignored_val] : raw) {
      if (key != "Total") {
        ++non_total_count;
      }
    }
    const std::size_t num_partials = (non_total_count > 0) ? non_total_count : 1;

    std::vector<RenderSeries> styled;
    styled.reserve(raw.size());
    std::size_t partial_idx = 0;
    for (const auto &[key, vals] : raw) {
      const std::size_t idx = (key == "Total") ? 0 : partial_idx++;
      styled.push_back(RenderSeries{
          .key = key,
          .values = vals,
          .color = resolveColor(key, idx, num_partials, ctx.custom_colors, *ctx.config),
      });
    }
    return styled;
  }
};

// ============================================================================
// Stage 2: Analytics & Peak Detection (PlotAnalytics)
// ============================================================================
class PlotAnalytics {
public:
  [[nodiscard]] static std::pair<detail::NiceScale, detail::NiceScale>
  computeScales(std::span<const real_t> x_coords, std::span<const RenderSeries> series,
                const PlotConfig &config) {
    real_t raw_x_min = x_coords.empty() ? 0.0 : x_coords.front();
    real_t raw_x_max = x_coords.empty() ? 1.0 : x_coords.back();
    real_t raw_y_min = 0.0;
    real_t raw_y_max = 0.0;

    for (const auto &series_item : series) {
      for (const real_t val : series_item.values) {
        raw_y_max = std::max(raw_y_max, val);
        raw_y_min = std::min(raw_y_min, val);
      }
    }

    const real_t y_padding = (raw_y_max - raw_y_min) * static_cast<real_t>(0.05);
    raw_y_max += y_padding;

    bool strict_x = false;
    bool strict_y = false;
    if (config.manual_x_min.has_value() && config.manual_x_max.has_value()) {
      raw_x_min = *config.manual_x_min;
      raw_x_max = *config.manual_x_max;
      strict_x = true;
    }
    if (config.manual_y_min.has_value() && config.manual_y_max.has_value()) {
      raw_y_min = *config.manual_y_min;
      raw_y_max = *config.manual_y_max;
      strict_y = true;
    }

    const real_t x_step = (config.manual_ticks && config.x_major_step > static_cast<real_t>(0.0))
                              ? config.x_major_step
                              : static_cast<real_t>(0.0);
    const real_t y_step = (config.manual_ticks && config.y_major_step > static_cast<real_t>(0.0))
                              ? config.y_major_step
                              : static_cast<real_t>(0.0);

    return {
        detail::NiceScale(detail::DataRange{.min = raw_x_min, .max = raw_x_max}, 11, strict_x,
                          x_step),
        detail::NiceScale(detail::DataRange{.min = raw_y_min, .max = raw_y_max}, 8, strict_y,
                          y_step),
    };
  }

  [[nodiscard]] static std::vector<PeakCallout> detectPeaks(std::span<const real_t> x_coords,
                                                            std::span<const RenderSeries> series,
                                                            std::string_view x_label = "") {
    if (series.empty()) {
      return {};
    }
    const auto total_it = std::ranges::find_if(
        series, [](const auto &item) { return item.key == "Total" || item.key == "total"; });
    const auto &target_series = (total_it != series.end()) ? *total_it : series.front();
    const auto &y_values = target_series.values;
    const std::size_t num_pts = std::min(x_coords.size(), y_values.size());
    if (num_pts < 3) {
      return {};
    }

    const std::span<const real_t> y_span(y_values.data(), num_pts);
    const auto max_it = std::ranges::max_element(y_span);
    if (max_it == y_span.end() || *max_it <= static_cast<real_t>(0.0)) {
      return {};
    }
    const real_t max_val = *max_it;
    const real_t threshold = max_val * static_cast<real_t>(0.05);

    std::vector<std::pair<real_t, real_t>> candidates;
    for (std::size_t i = 1; i < num_pts - 1; ++i) {
      if (y_values[i] > y_values[i - 1] && y_values[i] >= y_values[i + 1] &&
          y_values[i] > threshold) {
        candidates.emplace_back(x_coords[i], y_values[i]);
      }
    }

    if (candidates.empty()) {
      const auto idx = static_cast<std::size_t>(std::distance(y_span.begin(), max_it));
      candidates.emplace_back(x_coords[idx], max_val);
    }

    std::ranges::sort(candidates,
                      [](const auto &lhs, const auto &rhs) { return lhs.second > rhs.second; });

    return {
        PeakCallout{
            .label = (x_label == "r") ? "r₁" : "Peak",
            .x = candidates[0].first,
            .y = candidates[0].second,
        },
    };
  }

  struct NearestPointHit {
    std::size_t index{0};
    std::string key;
  };

  [[nodiscard]] static std::optional<NearestPointHit>
  findNearestPoint(real_t mouse_sx, real_t mouse_sy, const CanvasViewport &viewport,
                   std::span<const real_t> x_coords, std::span<const RenderSeries> series,
                   const detail::NiceScale &x_scale, const detail::NiceScale &y_scale) {
    if (x_coords.empty() || series.empty()) {
      return std::nullopt;
    }
    auto min_dist_sq = static_cast<real_t>(1e30);
    NearestPointHit best;

    for (const auto &series_item : series) {
      const std::size_t num_points = std::min(x_coords.size(), series_item.values.size());
      for (std::size_t i = 0; i < num_points; ++i) {
        const real_t sx_i = viewport.mapX(x_coords[i], x_scale);
        const real_t sy_i = viewport.mapY(series_item.values[i], y_scale);
        const real_t dx_i = sx_i - mouse_sx;
        const real_t dy_i = sy_i - mouse_sy;
        const real_t dist_sq = dx_i * dx_i + dy_i * dy_i;
        if (dist_sq < min_dist_sq) {
          min_dist_sq = dist_sq;
          best.index = i;
          best.key = series_item.key;
        }
      }
    }
    return best;
  }
};

// ============================================================================
// Stage 3 & 4: SVG Emission (SvgEmitter)
// ============================================================================
class SvgEmitter {
public:
  explicit SvgEmitter(std::size_t capacity = layout::DEFAULT_BUFFER_CAPACITY) {
    buffer_.reserve(capacity);
  }

  void emitDefsAndHeader(const CanvasViewport &viewport, const PlotConfig &config,
                         std::span<const RenderSeries> series) {
    std::format_to(
        std::back_inserter(buffer_),
        "<svg width='{:.0f}' height='{:.0f}' xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {} {}\" "
        "shape-rendering=\"geometricPrecision\" text-rendering=\"geometricPrecision\">\n"
        "  <defs>\n"
        "    <filter id=\"tooltip-shadow\" x=\"-10%\" y=\"-10%\" width=\"120%\" height=\"120%\">\n"
        "      <feDropShadow dx=\"2\" dy=\"4\" stdDeviation=\"4\" flood-color=\"#000000\" flood-opacity=\"0.15\"/>\n"
        "    </filter>\n",
        viewport.width, viewport.height, viewport.width, viewport.height);

    for (std::size_t i = 0; i < series.size(); ++i) {
      std::format_to(
          std::back_inserter(buffer_),
          "    <linearGradient id=\"area-grad-{}\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">\n"
          "      <stop offset=\"0%\" stop-color=\"{}\" stop-opacity=\"0.35\"/>\n"
          "      <stop offset=\"100%\" stop-color=\"{}\" stop-opacity=\"0.0\"/>\n"
          "    </linearGradient>\n",
          i, series[i].color, series[i].color);
    }

    std::format_to(
        std::back_inserter(buffer_),
        "    <clipPath id=\"plot-area-clip\"><rect x=\"{:.1f}\" y=\"{:.1f}\" width=\"{:.1f}\" "
        "height=\"{:.1f}\"/></clipPath>\n"
        "  </defs>\n"
        "  <rect width=\"100%\" height=\"100%\" fill=\"{}\" rx=\"6\"/>\n",
        viewport.px0, viewport.py0, viewport.plotWidth(), viewport.plotHeight(), config.bgColor());
  }

  void emitGridAndAxes(const CanvasViewport &viewport, const detail::NiceScale &x_scale,
                       const detail::NiceScale &y_scale, const PlotConfig &config) {
    buffer_.append("  <!-- Grid & Ticks -->\n");
    emitYTicks(viewport, y_scale, config);
    emitXTicks(viewport, x_scale, config);
    emitMinorTicks(viewport, x_scale, y_scale, config);

    std::format_to(std::back_inserter(buffer_),
                   "  <rect x=\"{:.1f}\" y=\"{:.1f}\" width=\"{:.1f}\" height=\"{:.1f}\" "
                   "fill=\"none\" stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                   viewport.px0, viewport.py0, viewport.plotWidth(), viewport.plotHeight(),
                   config.axisColor());
  }

  void emitEmphasisLines(const CanvasViewport &viewport, const detail::NiceScale &y_scale,
                         std::string_view title, const PlotConfig &config) {
    std::vector<real_t> emphasis_values;
    if (title.contains("g(r)")) {
      emphasis_values.push_back(static_cast<real_t>(1.0));
    } else if (title.contains("G(r)")) {
      emphasis_values.push_back(static_cast<real_t>(0.0));
    } else if (title.contains("S(Q)") || title.contains("S(q)")) {
      emphasis_values.push_back(static_cast<real_t>(0.0));
      emphasis_values.push_back(static_cast<real_t>(1.0));
    }

    for (const real_t focus_y : emphasis_values) {
      if (focus_y >= y_scale.min && focus_y <= y_scale.max) {
        const real_t spy = viewport.mapY(focus_y, y_scale);
        const std::string_view extra =
            (std::abs(focus_y - static_cast<real_t>(1.0)) < static_cast<real_t>(1e-6))
                ? " stroke-dasharray=\"5,5\""
                : "";
        std::format_to(std::back_inserter(buffer_),
                       "  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                       "stroke=\"{}\" stroke-width=\"1.8\"{}/>\n",
                       viewport.px0, spy, viewport.px1, spy, config.axisColor(), extra);
      }
    }
  }

  void emitAreaFills(const CanvasViewport &viewport, std::span<const real_t> x_coords,
                     std::span<const RenderSeries> series, const detail::NiceScale &x_scale,
                     const detail::NiceScale &y_scale, const PlotConfig &config) {
    if (!config.fill_area || x_coords.empty()) {
      return;
    }
    for (std::size_t i = 0; i < series.size(); ++i) {
      const auto &series_item = series[i];
      const std::size_t num_points = std::min(x_coords.size(), series_item.values.size());
      if (num_points <= 1) {
        continue;
      }

      std::format_to(std::back_inserter(buffer_),
                     "  <polygon fill=\"url(#area-grad-{})\" stroke=\"none\" points=\"", i);
      const real_t sx_start = viewport.mapX(x_coords.front(), x_scale);
      std::format_to(std::back_inserter(buffer_), "{:.2f},{:.2f} ", sx_start, viewport.py1);

      for (std::size_t pt = 0; pt < num_points; ++pt) {
        const real_t screen_x = viewport.mapX(x_coords[pt], x_scale);
        const real_t screen_y = viewport.mapY(series_item.values[pt], y_scale);
        const real_t sy_clamped = std::min(screen_y, viewport.py1);
        std::format_to(std::back_inserter(buffer_), "{:.2f},{:.2f} ", screen_x, sy_clamped);
      }

      const real_t sx_end = viewport.mapX(x_coords[num_points - 1], x_scale);
      std::format_to(std::back_inserter(buffer_), "{:.2f},{:.2f}\" />\n", sx_end, viewport.py1);
    }
  }

  void emitPolylines(const CanvasViewport &viewport, std::span<const real_t> x_coords,
                     std::span<const RenderSeries> series, const detail::NiceScale &x_scale,
                     const detail::NiceScale &y_scale, const PlotConfig &config) {
    if (x_coords.empty()) {
      return;
    }
    std::string total_svg;

    for (const auto &series_item : series) {
      std::string line_svg = std::format(
          R"(  <polyline fill="none" stroke="{}" stroke-width="{:.1f}" stroke-linejoin="round" points=")",
          series_item.color, config.line_width);
      const std::size_t num_points = std::min(x_coords.size(), series_item.values.size());
      for (std::size_t pt = 0; pt < num_points; ++pt) {
        const real_t screen_x = viewport.mapX(x_coords[pt], x_scale);
        const real_t screen_y = viewport.mapY(series_item.values[pt], y_scale);
        line_svg += std::format("{:.2f},{:.2f} ", screen_x, screen_y);
      }
      line_svg += "\" />\n";

      if (series_item.key == "Total") {
        total_svg = std::move(line_svg);
      } else {
        buffer_.append(line_svg);
      }
    }

    if (!total_svg.empty()) {
      buffer_.append(total_svg);
    }
  }

  void emitMarkers(const CanvasViewport &viewport, std::span<const real_t> x_coords,
                   std::span<const RenderSeries> series, const detail::NiceScale &x_scale,
                   const detail::NiceScale &y_scale, const PlotConfig &config) {
    if (!config.show_markers || x_coords.empty()) {
      return;
    }
    for (const auto &series_item : series) {
      const std::size_t num_points = std::min(x_coords.size(), series_item.values.size());
      for (std::size_t pt = 0; pt < num_points; ++pt) {
        const real_t screen_x = viewport.mapX(x_coords[pt], x_scale);
        const real_t screen_y = viewport.mapY(series_item.values[pt], y_scale);
        std::format_to(std::back_inserter(buffer_),
                       "  <circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"{:.1f}\" fill=\"{}\" "
                       "stroke=\"none\"/>\n",
                       screen_x, screen_y, config.marker_size, series_item.color);
      }
    }
  }

  void emitReferenceLines(const CanvasViewport &viewport, const detail::NiceScale &x_scale,
                          const detail::NiceScale &y_scale, const PlotConfig &config) {
    if (config.reference_lines.empty()) {
      return;
    }
    for (const auto &ref : config.reference_lines) {
      const std::string_view color =
          ref.color_hex.empty() ? config.axisColor() : std::string_view(ref.color_hex);
      if (ref.is_vertical) {
        if (ref.value >= x_scale.min && ref.value <= x_scale.max) {
          const real_t spx = viewport.mapX(ref.value, x_scale);
          std::format_to(std::back_inserter(buffer_),
                         "  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"1.8\" stroke-dasharray=\"4,4\"/>\n",
                         spx, viewport.py0, spx, viewport.py1, color);
          if (!ref.label.empty()) {
            buffer_.append(renderTextAsPath(ref.label, spx + static_cast<real_t>(5.0),
                                            viewport.py0 + static_cast<real_t>(18.0),
                                            static_cast<real_t>(16.0) * config.font_scale,
                                            TextAnchor::Start, color, config.use_native_text));
          }
        }
      } else {
        if (ref.value >= y_scale.min && ref.value <= y_scale.max) {
          const real_t spy = viewport.mapY(ref.value, y_scale);
          std::format_to(std::back_inserter(buffer_),
                         "  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"1.8\" stroke-dasharray=\"4,4\"/>\n",
                         viewport.px0, spy, viewport.px1, spy, color);
          if (!ref.label.empty()) {
            buffer_.append(renderTextAsPath(ref.label, viewport.px0 + static_cast<real_t>(8.0),
                                            spy - static_cast<real_t>(5.0),
                                            static_cast<real_t>(16.0) * config.font_scale,
                                            TextAnchor::Start, color, config.use_native_text));
          }
        }
      }
    }
  }

  void emitPeakCallouts(const CanvasViewport &viewport, const detail::NiceScale &x_scale,
                        const detail::NiceScale &y_scale, std::span<const PeakCallout> peaks,
                        const PlotConfig &config) {
    if (!config.show_peak_callouts || peaks.empty()) {
      return;
    }
    for (const auto &peak : peaks) {
      const real_t screen_x = viewport.mapX(peak.x, x_scale);
      const real_t screen_y = viewport.mapY(peak.y, y_scale);
      if (!viewport.containsPixel(screen_x, screen_y)) {
        continue;
      }

      std::format_to(
          std::back_inserter(buffer_),
          "  <polygon points=\"{:.1f},{:.1f} {:.1f},{:.1f} {:.1f},{:.1f} {:.1f},{:.1f}\" "
          "fill=\"#ffffff\" stroke=\"#70b0ff\" stroke-width=\"2.0\"/>\n",
          screen_x, screen_y - static_cast<real_t>(6.0), screen_x + static_cast<real_t>(6.0),
          screen_y, screen_x, screen_y + static_cast<real_t>(6.0),
          screen_x - static_cast<real_t>(6.0), screen_y);

      const std::string text = std::format("{} {:.2f} ({:.2f})", peak.label, peak.x, peak.y);
      const real_t badge_y =
          std::max(viewport.py0 + static_cast<real_t>(20.0), screen_y - static_cast<real_t>(18.0));

      std::format_to(std::back_inserter(buffer_),
                     "  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                     "stroke=\"#70b0ff\" stroke-width=\"1.0\" stroke-dasharray=\"2,2\"/>\n",
                     screen_x, screen_y - static_cast<real_t>(6.0), screen_x,
                     badge_y + static_cast<real_t>(6.0));
      buffer_.append(renderTextAsPath(
          text, screen_x, badge_y, static_cast<real_t>(16.0) * config.font_scale,
          TextAnchor::Middle, (config.theme == PlotConfig::Theme::Light) ? "#0072b2" : "#70b0ff",
          config.use_native_text));
    }
  }

  void emitLegend(const CanvasViewport &viewport, std::span<const RenderSeries> series,
                  const PlotConfig &config) {
    if (!config.show_legend || series.empty()) {
      return;
    }
    const real_t legend_x = viewport.px1 - layout::LEGEND_RIGHT_MARGIN;
    real_t legend_y = viewport.py0 + layout::LEGEND_TOP_MARGIN;

    // 1. Total entry
    const auto total_it =
        std::ranges::find_if(series, [](const auto &item) { return item.key == "Total"; });
    if (total_it != series.end()) {
      std::format_to(std::back_inserter(buffer_),
                     "  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                     "stroke=\"{}\" stroke-width=\"4.0\"/>\n",
                     legend_x - layout::LEGEND_LINE_LEFT_PAD, legend_y,
                     legend_x - layout::LEGEND_LINE_RIGHT_PAD, legend_y, total_it->color);
      buffer_.append(renderTextAsPath("Total", legend_x - layout::LEGEND_TEXT_PAD,
                                      legend_y + static_cast<real_t>(6.0),
                                      static_cast<real_t>(20.0) * config.font_scale,
                                      TextAnchor::End, config.textColor(), config.use_native_text));
      legend_y += layout::LEGEND_TOTAL_STEP;
    }

    // 2. Partials
    const bool has_partials =
        std::ranges::any_of(series, [](const auto &item) { return item.key != "Total"; });
    if (!has_partials) {
      return;
    }

    buffer_.append(renderTextAsPath(
        "Partials:", legend_x - layout::LEGEND_LINE_RIGHT_PAD, legend_y + static_cast<real_t>(4.0),
        static_cast<real_t>(16.0) * config.font_scale, TextAnchor::End,
        (config.theme == PlotConfig::Theme::Light) ? "#555555" : "#cbd5e1",
        config.use_native_text));
    legend_y += layout::LEGEND_HEADER_STEP;

    for (const auto &series_item : series) {
      if (series_item.key == "Total") {
        continue;
      }
      std::format_to(std::back_inserter(buffer_),
                     "  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                     "stroke=\"{}\" stroke-width=\"3.0\"/>\n",
                     legend_x - layout::LEGEND_LINE_LEFT_PAD, legend_y,
                     legend_x - layout::LEGEND_LINE_RIGHT_PAD, legend_y, series_item.color);
      buffer_.append(renderTextAsPath(series_item.key, legend_x - layout::LEGEND_TEXT_PAD,
                                      legend_y + static_cast<real_t>(5.0),
                                      static_cast<real_t>(18.0) * config.font_scale,
                                      TextAnchor::End, config.textColor(), config.use_native_text));
      legend_y += layout::LEGEND_ITEM_STEP;
    }
  }

  void emitTitlesAndLabels(const CanvasViewport &viewport, const std::string &x_label,
                           const std::string &y_label, const PlotConfig &config) {
    buffer_.append(renderTextAsPath(
        x_label, (viewport.px0 + viewport.px1) / static_cast<real_t>(2.0),
        viewport.py1 + layout::X_AXIS_TITLE_OFFSET, static_cast<real_t>(32.0) * config.font_scale,
        TextAnchor::Middle, config.textColor(), config.use_native_text));

    std::format_to(
        std::back_inserter(buffer_), "  <g transform=\"translate({:.1f}, {:.1f}) rotate(-90)\">\n",
        layout::Y_AXIS_TITLE_OFFSET, (viewport.py0 + viewport.py1) / static_cast<real_t>(2.0));
    buffer_.append(renderTextAsPath(y_label, static_cast<real_t>(0.0), static_cast<real_t>(0.0),
                                    static_cast<real_t>(32.0) * config.font_scale,
                                    TextAnchor::Middle, config.textColor(),
                                    config.use_native_text));
    buffer_.append("  </g>\n");
  }

  void emitHover(const CanvasViewport &viewport, const HoverInfo &hover,
                 std::span<const real_t> x_coords, std::span<const RenderSeries> series,
                 const detail::NiceScale &x_scale, const detail::NiceScale &y_scale,
                 const PlotConfig &config) {
    if (!hover.active || hover.widget_width <= static_cast<real_t>(0.0) ||
        hover.widget_height <= static_cast<real_t>(0.0) || x_coords.empty()) {
      return;
    }

    const real_t svg_aspect = viewport.width / viewport.height;
    const real_t widget_aspect = hover.widget_width / hover.widget_height;
    real_t scale = 1.0;
    real_t offset_x = 0.0;
    real_t offset_y = 0.0;

    if (widget_aspect > svg_aspect) {
      scale = hover.widget_height / viewport.height;
      offset_x = (hover.widget_width - viewport.width * scale) / static_cast<real_t>(2.0);
    } else {
      scale = hover.widget_width / viewport.width;
      offset_y = (hover.widget_height - viewport.height * scale) / static_cast<real_t>(2.0);
    }

    const real_t mouse_sx = (hover.mouse_x - offset_x) / scale;
    const real_t mouse_sy = (hover.mouse_y - offset_y) / scale;

    if (!viewport.containsPixel(mouse_sx, mouse_sy)) {
      return;
    }

    const auto hit = PlotAnalytics::findNearestPoint(mouse_sx, mouse_sy, viewport, x_coords, series,
                                                     x_scale, y_scale);
    if (!hit) {
      return;
    }

    const std::size_t idx = hit->index;
    const real_t sx_data = viewport.mapX(x_coords[idx], x_scale);

    std::format_to(std::back_inserter(buffer_),
                   "  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                   "stroke=\"{}\" stroke-width=\"1.5\" stroke-dasharray=\"4,4\"/>\n",
                   sx_data, viewport.py0, sx_data, viewport.py1, config.axisColor());

    auto snapped_sy_data = static_cast<real_t>(-1.0);
    for (const auto &series_item : series) {
      if (idx < series_item.values.size()) {
        const real_t sy_data = viewport.mapY(series_item.values[idx], y_scale);
        if (series_item.key == hit->key) {
          snapped_sy_data = sy_data;
        }
        std::format_to(
            std::back_inserter(buffer_),
            "  <circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"6\" fill=\"{}\" stroke=\"{}\" stroke-width=\"2\"/>\n",
            sx_data, sy_data, series_item.color, config.bgColor());
      }
    }

    if (snapped_sy_data >= viewport.py0 && snapped_sy_data <= viewport.py1) {
      std::format_to(std::back_inserter(buffer_),
                     "  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                     "stroke=\"{}\" stroke-width=\"1.5\" stroke-dasharray=\"4,4\"/>\n",
                     viewport.px0, snapped_sy_data, sx_data, snapped_sy_data, config.axisColor());
    }
  }

  void appendRaw(std::string_view str) { buffer_.append(str); }

  [[nodiscard]] std::string finalize() {
    buffer_.append("</svg>\n");
    return std::move(buffer_);
  }

private:
  std::string buffer_;

  void emitYTicks(const CanvasViewport &viewport, const detail::NiceScale &y_scale,
                  const PlotConfig &config) {
    for (const real_t y_val : y_scale.ticks) {
      const real_t spy = viewport.mapY(y_val, y_scale);
      if (config.show_grid) {
        std::format_to(std::back_inserter(buffer_),
                       "  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                       "stroke=\"{}\" stroke-width=\"0.8\" stroke-dasharray=\"3,3\"/>\n",
                       viewport.px0, spy, viewport.px1, spy, config.gridColor());
      }
      std::format_to(std::back_inserter(buffer_),
                     "  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                     "stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                     viewport.px0 - layout::MAJOR_TICK_LENGTH, spy, viewport.px0, spy,
                     config.axisColor());
      buffer_.append(renderTextAsPath(
          fmtScientific(y_val), viewport.px0 - layout::Y_TICK_LABEL_OFFSET,
          spy + layout::Y_TICK_LABEL_Y_SHIFT, static_cast<real_t>(23.0) * config.font_scale,
          TextAnchor::End, config.textColor(), config.use_native_text));
    }
  }

  void emitXTicks(const CanvasViewport &viewport, const detail::NiceScale &x_scale,
                  const PlotConfig &config) {
    for (const real_t x_val : x_scale.ticks) {
      const real_t spx = viewport.mapX(x_val, x_scale);
      if (config.show_grid) {
        std::format_to(std::back_inserter(buffer_),
                       "  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                       "stroke=\"{}\" stroke-width=\"0.8\" stroke-dasharray=\"3,3\"/>\n",
                       spx, viewport.py0, spx, viewport.py1, config.gridColor());
      }
      std::format_to(std::back_inserter(buffer_),
                     "  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                     "stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                     spx, viewport.py1, spx, viewport.py1 + layout::MAJOR_TICK_LENGTH,
                     config.axisColor());
      buffer_.append(
          renderTextAsPath(fmtScientific(x_val), spx, viewport.py1 + layout::X_TICK_LABEL_OFFSET,
                           static_cast<real_t>(23.0) * config.font_scale, TextAnchor::Middle,
                           config.textColor(), config.use_native_text));
    }
  }

  void emitMinorTicks(const CanvasViewport &viewport, const detail::NiceScale &x_scale,
                      const detail::NiceScale &y_scale, const PlotConfig &config) {
    if (config.minor_subdivisions <= 1) {
      return;
    }
    const auto subs = static_cast<real_t>(config.minor_subdivisions);

    // Minor Y
    if (y_scale.ticks.size() >= 2) {
      for (std::size_t i = 0; i < y_scale.ticks.size() - 1; ++i) {
        const real_t tick_start = y_scale.ticks[i];
        const real_t delta = (y_scale.ticks[i + 1] - tick_start) / subs;
        for (int sub_idx = 1; sub_idx < config.minor_subdivisions; ++sub_idx) {
          const real_t spy =
              viewport.mapY(tick_start + static_cast<real_t>(sub_idx) * delta, y_scale);
          if (spy >= viewport.py0 && spy <= viewport.py1) {
            std::format_to(std::back_inserter(buffer_),
                           "  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                           "stroke=\"{}\" stroke-width=\"1.0\"/>\n",
                           viewport.px0 - layout::MINOR_TICK_LENGTH, spy, viewport.px0, spy,
                           config.axisColor());
          }
        }
      }
    }

    // Minor X
    if (x_scale.ticks.size() >= 2) {
      for (std::size_t i = 0; i < x_scale.ticks.size() - 1; ++i) {
        const real_t tick_start = x_scale.ticks[i];
        const real_t delta = (x_scale.ticks[i + 1] - tick_start) / subs;
        for (int sub_idx = 1; sub_idx < config.minor_subdivisions; ++sub_idx) {
          const real_t spx =
              viewport.mapX(tick_start + static_cast<real_t>(sub_idx) * delta, x_scale);
          if (spx >= viewport.px0 && spx <= viewport.px1) {
            std::format_to(std::back_inserter(buffer_),
                           "  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                           "stroke=\"{}\" stroke-width=\"1.0\"/>\n",
                           spx, viewport.py1, spx, viewport.py1 + layout::MINOR_TICK_LENGTH,
                           config.axisColor());
          }
        }
      }
    }
  }
};

[[nodiscard]] std::string renderNoDataPlaceholder(const CanvasViewport &viewport,
                                                  const PlotConfig &config) {
  if (config.use_native_text) {
    return std::format(
        "<svg width='{0:.0f}' height='{1:.0f}' xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {0:.0f} {1:.0f}\">"
        "<rect width=\"100%\" height=\"100%\" fill=\"{2}\"/>"
        "<text x=\"{3:.1f}\" y=\"{4:.1f}\" font-family=\"'Outfit', 'Plus Jakarta Sans', 'Inter', 'Roboto', 'Helvetica Neue', "
        "sans-serif\" font-size=\"{5:.1f}\" text-anchor=\"middle\" fill=\"{6}\">No data available</text></svg>",
        viewport.width, viewport.height, config.bgColor(),
        viewport.width / static_cast<real_t>(2.0),
        viewport.height / static_cast<real_t>(2.0) + static_cast<real_t>(8.0),
        static_cast<real_t>(24.0) * config.font_scale, config.textColor());
  }

  const std::string no_data_path = Roboto::instance().render(TextRenderParameters{
      .text = "No data available",
      .start_x = viewport.width / static_cast<real_t>(2.0),
      .start_y = viewport.height / static_cast<real_t>(2.0) + static_cast<real_t>(8.0),
      .font_size = static_cast<real_t>(24.0) * config.font_scale,
      .anchor = "middle",
  });

  return std::format(
      "<svg width='{0:.0f}' height='{1:.0f}' xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {0:.0f} {1:.0f}\">"
      "<rect width=\"100%\" height=\"100%\" fill=\"{2}\"/>"
      "<path d=\"{3}\" fill=\"{4}\" fill-rule=\"evenodd\" stroke=\"none\"/></svg>",
      viewport.width, viewport.height, config.bgColor(), no_data_path, config.textColor());
}

[[nodiscard]] std::pair<std::string, std::string>
formatAxisLabels(const correlation::analysis::Histogram &hist) {
  std::string x_label = hist.x_label.empty() ? "x" : hist.x_label;
  std::string y_label = hist.y_label.empty() ? "y" : hist.y_label;

  if (!hist.x_unit.empty()) {
    x_label += std::format(" ({})", hist.x_unit);
  }
  if (!hist.y_unit.empty()) {
    std::string y_unit = hist.y_unit;
    if (y_unit == "Å^-1") {
      y_unit = "Å⁻¹";
    }
    y_label += std::format(" ({})", y_unit);
  }
  return {x_label, y_label};
}

} // namespace

// ============================================================================
// Public Entry Point
// ============================================================================
std::string renderHistogramAsSvg(const correlation::analysis::Histogram &hist,
                                 const PlotConfig &config, const HoverInfo &hover,
                                 const std::map<std::string, real_t> &weights,
                                 const std::map<std::string, bool> &curve_visibility,
                                 const std::map<std::string, std::string> &custom_colors) {
  const RenderContext ctx{
      .hist = &hist,
      .config = &config,
      .hover = &hover,
      .weights = &weights,
      .curve_visibility = &curve_visibility,
      .custom_colors = &custom_colors,
  };

  const CanvasViewport viewport = CanvasViewport::create(config);

  // Stage 1: Filter & prepare series
  const std::vector<RenderSeries> series = SeriesFilter::filterSeries(ctx);
  if (hist.bins.empty() || series.empty()) {
    return renderNoDataPlaceholder(viewport, config);
  }

  // Stage 2: Analytics, scales & peak detection
  const auto [x_scale, y_scale] = PlotAnalytics::computeScales(hist.bins, series, config);
  const std::vector<PeakCallout> peaks =
      PlotAnalytics::detectPeaks(hist.bins, series, hist.x_label);
  const auto [x_label, y_label] = formatAxisLabels(hist);
  const std::string_view title = hist.title.empty() ? "Histogram" : std::string_view(hist.title);

  // Stage 3 & 4: SVG Emission
  SvgEmitter emitter;
  emitter.emitDefsAndHeader(viewport, config, series);
  emitter.emitGridAndAxes(viewport, x_scale, y_scale, config);
  emitter.emitEmphasisLines(viewport, y_scale, title, config);

  emitter.appendRaw("  <g clip-path=\"url(#plot-area-clip)\">\n");
  emitter.emitAreaFills(viewport, hist.bins, series, x_scale, y_scale, config);
  emitter.emitPolylines(viewport, hist.bins, series, x_scale, y_scale, config);
  emitter.emitMarkers(viewport, hist.bins, series, x_scale, y_scale, config);
  emitter.emitPeakCallouts(viewport, x_scale, y_scale, peaks, config);
  emitter.emitReferenceLines(viewport, x_scale, y_scale, config);
  emitter.appendRaw("  </g>\n");

  emitter.emitLegend(viewport, series, config);
  emitter.emitTitlesAndLabels(viewport, x_label, y_label, config);
  emitter.emitHover(viewport, hover, hist.bins, series, x_scale, y_scale, config);

  return emitter.finalize();
}

} // namespace correlation::plotters
