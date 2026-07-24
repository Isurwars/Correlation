/**
 * @file SvgHistogramRenderer.hpp
 * @brief Single histogram rendering logic for the SVG plotter.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "plotters/PlotTypes.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <format>
#include <map>
#include <ranges>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace correlation::plotters {
namespace detail {

/**
 * @struct SvgHistogramRenderer
 * @brief Internal state and rendering logic for single-histogram SVG generation.
 */
struct SvgHistogramRenderer {
  const correlation::analysis::Histogram *hist = nullptr;
  const PlotConfig *config = nullptr;
  const HoverInfo *hover = nullptr;
  const std::map<std::string, real_t> *weights = nullptr;

  std::string title;
  std::string x_label;
  std::string y_label;

  real_t kW = static_cast<real_t>(0.0);
  real_t kH = static_cast<real_t>(0.0);
  real_t px0 = static_cast<real_t>(0.0);
  real_t px1 = static_cast<real_t>(0.0);
  real_t py0 = static_cast<real_t>(0.0);
  real_t py1 = static_cast<real_t>(0.0);

  std::map<std::string, std::vector<real_t>> partials;
  const std::vector<real_t> *xs = nullptr;

  detail::NiceScale xScale;
  detail::NiceScale yScale;

  std::ostringstream svg;

  SvgHistogramRenderer(const correlation::analysis::Histogram &histogram, const PlotConfig &cfg, const HoverInfo &hov,
                       const std::map<std::string, real_t> &component_weights)
      : hist(&histogram), config(&cfg), hover(&hov), weights(&component_weights), xs(&histogram.bins) {
    initializeLayoutAndLabels();
    filterPartials();
  }

  void initializeLayoutAndLabels() {
    title = hist->title.empty() ? "Histogram" : hist->title;
    x_label = hist->x_label.empty() ? "x" : hist->x_label;
    y_label = hist->y_label.empty() ? "y" : hist->y_label;

    if (!hist->x_unit.empty()) {
      x_label += std::format(" ({})", hist->x_unit);
    }
    if (!hist->y_unit.empty()) {
      std::string y_unit = hist->y_unit;
      if (y_unit == "Å^-1") {
        y_unit = "Å⁻¹";
      }
      y_label += std::format(" ({})", y_unit);
    }

    kW = config->effective_width();
    kH = config->effective_height();
    const real_t kLeft = static_cast<real_t>(100.0);
    const real_t kRight = static_cast<real_t>(40.0);
    const real_t kTop = static_cast<real_t>(50.0);
    const real_t kBot = static_cast<real_t>(90.0);

    px0 = kLeft;
    px1 = kW - kRight;
    py0 = kTop;
    py1 = kH - kBot;
  }

  void filterPartials() {
    const auto &raw_partials = hist->smoothed_partials.empty() ? hist->partials : hist->smoothed_partials;
    auto total_it = raw_partials.find("Total");
    if (total_it != raw_partials.end()) {
      partials["Total"] = total_it->second;
    }

    std::vector<std::pair<std::string, real_t>> candidates;
    for (const auto &[key, value] : raw_partials) {
      if (key == "Total" || key.starts_with("Frequency_")) {
        continue;
      }
      real_t score = static_cast<real_t>(0.0);
      if (weights != nullptr && !weights->empty()) {
        auto wit = weights->find(key);
        if (wit != weights->end()) {
          score = wit->second;
        }
      } else {
        for (real_t val : value) {
          score += std::abs(val);
        }
      }
      candidates.emplace_back(key, score);
    }

    std::ranges::sort(candidates, [](const auto &lhs, const auto &rhs) { return lhs.second > rhs.second; });

    std::size_t limit = std::min(candidates.size(), std::size_t(10));
    for (std::size_t i = 0; i < limit; ++i) {
      const std::string &key = candidates.at(i).first;
      auto iter = raw_partials.find(key);
      if (iter != raw_partials.end()) {
        partials[key] = iter->second;
      }
    }
  }

  bool hasNoData() const { return xs == nullptr || xs->empty() || partials.empty(); }

  std::string renderNoData() const {
    if (config->use_native_text) {
      return std::format(
          "<svg width='{0:.0f}' height='{1:.0f}' xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {0:.0f} {1:.0f}\">"
          "<rect width=\"100%\" height=\"100%\" fill=\"{2}\"/>"
          "<text x=\"{3:.1f}\" y=\"{4:.1f}\" font-family=\"'Outfit', 'Plus Jakarta Sans', 'Inter', 'Roboto', 'Helvetica Neue', "
          "sans-serif\" font-size=\"{5:.1f}\" text-anchor=\"middle\" fill=\"{6}\">No data available</text></svg>",
          kW, kH, config->bg_color(), kW / static_cast<real_t>(2.0),
          kH / static_cast<real_t>(2.0) + static_cast<real_t>(8.0), static_cast<real_t>(24.0) * config->font_scale,
          config->text_color());
    }
    std::string no_data_path = Roboto::instance().render(TextRenderParameters{
        .text = "No data available",
        .start_x = kW / static_cast<real_t>(2.0),
        .start_y = kH / static_cast<real_t>(2.0) + static_cast<real_t>(8.0),
        .font_size = static_cast<real_t>(24.0) * config->font_scale,
        .anchor = "middle",
    });
    return std::format(
        "<svg width='{0:.0f}' height='{1:.0f}' xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {0:.0f} {1:.0f}\">"
        "<rect width=\"100%\" height=\"100%\" fill=\"{2}\"/>"
        "<path d=\"{3}\" fill=\"{4}\" fill-rule=\"evenodd\" stroke=\"none\"/></svg>",
        kW, kH, config->bg_color(), no_data_path, config->text_color());
  }

  void computeScales() {
    if (xs == nullptr || xs->empty()) {
      return;
    }
    real_t raw_x_min = xs->front();
    real_t raw_x_max = xs->back();
    real_t raw_y_min = static_cast<real_t>(0.0);
    real_t raw_y_max = static_cast<real_t>(0.0);
    for (const auto &[key, value] : partials) {
      for (real_t val : value) {
        raw_y_max = std::max(raw_y_max, val);
        raw_y_min = std::min(raw_y_min, val);
      }
    }
    real_t y_padding = (raw_y_max - raw_y_min) * static_cast<real_t>(0.05);
    raw_y_max += y_padding;

    xScale = detail::NiceScale(
        detail::DataRange{
            .min = raw_x_min,
            .max = raw_x_max,
        },
        11);
    yScale = detail::NiceScale(
        detail::DataRange{
            .min = raw_y_min,
            .max = raw_y_max,
        },
        8);
  }

  void writeHeader() {
    svg << std::format("<svg width='{:.0f}' height='{:.0f}' xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {} {}\" "
                       "shape-rendering=\"geometricPrecision\" text-rendering=\"geometricPrecision\">\n",
                       kW, kH, kW, kH);
    svg << "  <defs>\n"
        << "    <filter id=\"tooltip-shadow\" x=\"-10%\" y=\"-10%\" width=\"120%\" height=\"120%\">\n"
        << "      <feDropShadow dx=\"2\" dy=\"4\" stdDeviation=\"4\" flood-color=\"#000000\" flood-opacity=\"0.15\"/>\n"
        << "    </filter>\n";

    std::size_t color_idx = 0;
    for (const auto &[key, value] : partials) {
      std::string col = color(color_idx, config->palette);
      std::string grad_id = std::format("area-grad-{}", color_idx);
      svg << std::format("    <linearGradient id=\"{}\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">\n"
                         "      <stop offset=\"0%\" stop-color=\"{}\" stop-opacity=\"0.35\"/>\n"
                         "      <stop offset=\"100%\" stop-color=\"{}\" stop-opacity=\"0.0\"/>\n"
                         "    </linearGradient>\n",
                         grad_id, col, col);
      color_idx++;
    }
    svg << "  </defs>\n";
    svg << std::format("  <rect width=\"100%\" height=\"100%\" fill=\"{}\" rx=\"6\"/>\n", config->bg_color());
  }

  void drawGridAndAxes() {
    svg << "  <!-- Grid & Ticks -->\n";
    for (real_t y_val : yScale.ticks) {
      real_t spy = mapValue(y_val, yScale.min, yScale.max, py1, py0);
      if (config->show_grid) {
        svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                           "stroke=\"{}\" stroke-width=\"0.8\" stroke-dasharray=\"3,3\"/>\n",
                           px0, spy, px1, spy, config->grid_color());
      }
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                         px0 - static_cast<real_t>(8.0), spy, px0, spy, config->axis_color());
      svg << renderTextAsPath(fmtScientific(y_val), px0 - static_cast<real_t>(15.0), spy + static_cast<real_t>(7.0),
                              static_cast<real_t>(20.0) * config->font_scale, TextAnchor::End, config->text_color(),
                              config->use_native_text);
    }

    for (real_t x_val : xScale.ticks) {
      real_t spx = mapValue(x_val, xScale.min, xScale.max, px0, px1);
      if (config->show_grid) {
        svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                           "stroke=\"{}\" stroke-width=\"0.8\" stroke-dasharray=\"3,3\"/>\n",
                           spx, py0, spx, py1, config->grid_color());
      }
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                         spx, py1, spx, py1 + static_cast<real_t>(8.0), config->axis_color());
      svg << renderTextAsPath(fmtScientific(x_val), spx, py1 + static_cast<real_t>(25.0),
                              static_cast<real_t>(20.0) * config->font_scale, TextAnchor::Middle, config->text_color(),
                              config->use_native_text);
    }

    svg << std::format("  <rect x=\"{:.1f}\" y=\"{:.1f}\" width=\"{:.1f}\" height=\"{:.1f}\" "
                       "fill=\"none\" stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                       px0, py0, px1 - px0, py1 - py0, config->axis_color());
  }

  void drawEmphasisLines() {
    std::vector<real_t> emphasis_values;
    if (title.contains("g(r)")) {
      emphasis_values.push_back(static_cast<real_t>(1.0));
    } else if (title.contains("G(r)")) {
      emphasis_values.push_back(static_cast<real_t>(0.0));
    } else if (title.contains("S(Q)") || title.contains("S(q)")) {
      emphasis_values.push_back(static_cast<real_t>(0.0));
      emphasis_values.push_back(static_cast<real_t>(1.0));
    }

    for (real_t focus_y : emphasis_values) {
      if (focus_y >= yScale.min && focus_y <= yScale.max) {
        real_t spy = mapValue(focus_y, yScale.min, yScale.max, py1, py0);
        std::string extra = (std::abs(focus_y - static_cast<real_t>(1.0)) < static_cast<real_t>(1e-6))
                                ? " stroke-dasharray=\"5,5\""
                                : "";
        svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                           "stroke=\"{}\" stroke-width=\"1.8\"{}/>\n",
                           px0, spy, px1, spy, config->axis_color(), extra);
      }
    }
  }

  void drawAreaFills() {
    if (config == nullptr || !config->fill_area || xs == nullptr) {
      return;
    }
    std::size_t fill_ci = 0;
    for (const auto &[key, value] : partials) {
      std::string grad_id = std::format("area-grad-{}", fill_ci);
      std::size_t num_points = std::min(xs->size(), value.size());
      if (num_points > 1) {
        svg << std::format("  <polygon fill=\"url(#{})\" stroke=\"none\" points=\"", grad_id);
        real_t sx_start = mapValue(xs->at(0), xScale.min, xScale.max, px0, px1);
        svg << std::format("{:.2f},{:.2f} ", sx_start, py1);
        for (std::size_t point_idx = 0; point_idx < num_points; ++point_idx) {
          real_t screen_x = mapValue(xs->at(point_idx), xScale.min, xScale.max, px0, px1);
          real_t screen_y = mapValue(value.at(point_idx), yScale.min, yScale.max, py1, py0);
          real_t sy_clamped = std::min(screen_y, py1);
          svg << std::format("{:.2f},{:.2f} ", screen_x, sy_clamped);
        }
        real_t sx_end = mapValue(xs->at(num_points - 1), xScale.min, xScale.max, px0, px1);
        svg << std::format("{:.2f},{:.2f}", sx_end, py1);
        svg << "\" />\n";
      }
      fill_ci++;
    }
  }

  std::vector<std::pair<std::string, std::string>> drawPolylines() {
    std::vector<std::pair<std::string, std::string>> legend;
    if (xs == nullptr) {
      return legend;
    }
    std::size_t color_idx = 0;
    for (const auto &[key, value] : partials) {
      const std::string col = color(color_idx++, config->palette);
      svg << std::format(
          R"(  <polyline fill="none" stroke="{}" stroke-width="{:.1f}" stroke-linejoin="round" points=")", col,
          config->line_width);
      std::size_t num_points = std::min(xs->size(), value.size());
      for (std::size_t point_idx = 0; point_idx < num_points; ++point_idx) {
        real_t screen_x = mapValue(xs->at(point_idx), xScale.min, xScale.max, px0, px1);
        real_t screen_y = mapValue(value.at(point_idx), yScale.min, yScale.max, py1, py0);
        svg << std::format("{:.2f},{:.2f} ", screen_x, screen_y);
      }
      svg << "\" />\n";
      legend.emplace_back(key, col);
    }
    return legend;
  }

  void drawMarkers() {
    if (!config->show_markers || xs == nullptr) {
      return;
    }
    std::size_t marker_ci = 0;
    for (const auto &[key, value] : partials) {
      const std::string col = color(marker_ci++, config->palette);
      std::size_t num_points = std::min(xs->size(), value.size());
      for (std::size_t point_idx = 0; point_idx < num_points; ++point_idx) {
        real_t screen_x = mapValue(xs->at(point_idx), xScale.min, xScale.max, px0, px1);
        real_t screen_y = mapValue(value.at(point_idx), yScale.min, yScale.max, py1, py0);
        svg << std::format("  <circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"3.5\" fill=\"{}\" stroke=\"none\"/>\n", screen_x,
                           screen_y, col);
      }
    }
  }

  void drawLegend(const std::vector<std::pair<std::string, std::string>> &legend) {
    if (!config->show_legend) {
      return;
    }
    real_t legend_x = px1 - static_cast<real_t>(15.0);
    real_t legend_y = py0 + static_cast<real_t>(25.0);
    for (const auto &iter : std::views::reverse(legend)) {
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"4.0\"/>\n",
                         legend_x - static_cast<real_t>(40.0), legend_y, legend_x - static_cast<real_t>(10.0), legend_y,
                         iter.second);
      svg << renderTextAsPath(iter.first, legend_x - static_cast<real_t>(45.0), legend_y + static_cast<real_t>(6.0),
                              static_cast<real_t>(18.0) * config->font_scale, TextAnchor::End, config->text_color(),
                              config->use_native_text);
      legend_y += static_cast<real_t>(28.0);
    }
  }

  void drawTitlesAndLabels() {
    svg << renderTextAsPath(x_label, (px0 + px1) / static_cast<real_t>(2.0), py1 + static_cast<real_t>(75.0),
                            static_cast<real_t>(28.0) * config->font_scale, TextAnchor::Middle, config->text_color(),
                            config->use_native_text);

    svg << std::format("  <g transform=\"translate({:.1f}, {:.1f}) rotate(-90)\">\n", static_cast<real_t>(40.0),
                       (py0 + py1) / static_cast<real_t>(2.0));
    svg << renderTextAsPath(y_label, static_cast<real_t>(0.0), static_cast<real_t>(0.0),
                            static_cast<real_t>(28.0) * config->font_scale, TextAnchor::Middle, config->text_color(),
                            config->use_native_text);
    svg << "  </g>\n";
  }

  /**
   * @struct NearestPoint
   * @brief Holds index and partial key for the nearest hover point in single-histogram rendering.
   */
  struct NearestPoint {
    std::size_t index = 0; ///< Bin index of nearest data point.
    std::string key;       ///< Partial key of nearest dataset line.
  };

  NearestPoint findNearestPoint(real_t mouse_sx, real_t mouse_sy) const {
    real_t min_dist_sq = static_cast<real_t>(1e30);
    NearestPoint best;
    for (const auto &[key, value] : partials) {
      std::size_t num_points = std::min(xs->size(), value.size());
      for (std::size_t point_idx = 0; point_idx < num_points; ++point_idx) {
        real_t sx_i = mapValue(xs->at(point_idx), xScale.min, xScale.max, px0, px1);
        real_t sy_i = mapValue(value.at(point_idx), yScale.min, yScale.max, py1, py0);
        real_t dx_i = sx_i - mouse_sx;
        real_t dy_i = sy_i - mouse_sy;
        real_t dist_sq = dx_i * dx_i + dy_i * dy_i;
        if (dist_sq < min_dist_sq) {
          min_dist_sq = dist_sq;
          best.index = point_idx;
          best.key = key;
        }
      }
    }
    return best;
  }

  /**
   * @struct TooltipPosition
   * @brief Position coordinates and snapped data values for single-histogram tooltips.
   */
  struct TooltipPosition {
    real_t sx_data;         ///< Canvas X coordinate.
    real_t snapped_sy_data; ///< Canvas Y coordinate of target point.
    real_t target_x;        ///< Physical x-axis data coordinate.
  };

  void drawTooltipBox(const TooltipPosition &pos, const std::string &best_key,
                      const std::vector<std::tuple<std::string, real_t, std::string>> &hover_values) {
    if (hover_values.empty()) {
      return;
    }
    real_t tooltip_w = static_cast<real_t>(200.0);
    real_t tooltip_h = static_cast<real_t>(35.0) + static_cast<real_t>(22.0) * static_cast<real_t>(hover_values.size());

    real_t tooltip_x = (pos.sx_data < kW / static_cast<real_t>(2.0))
                           ? pos.sx_data + static_cast<real_t>(15.0)
                           : pos.sx_data - tooltip_w - static_cast<real_t>(15.0);
    real_t tooltip_y = (pos.snapped_sy_data >= static_cast<real_t>(0.0))
                           ? pos.snapped_sy_data - tooltip_h / static_cast<real_t>(2.0)
                           : py0 + static_cast<real_t>(15.0);
    tooltip_y =
        std::max(py0 + static_cast<real_t>(8.0), std::min(py1 - tooltip_h - static_cast<real_t>(8.0), tooltip_y));

    std::string card_bg = (config->theme == PlotConfig::Theme::Light) ? "#FFFFFF" : "#181825";
    std::string card_border = (config->theme == PlotConfig::Theme::Light) ? "#dddddd" : "#45475a";
    std::string fg_color = (config->theme == PlotConfig::Theme::Light) ? "#333333" : "#cdd6f4";

    svg << std::format(
        "  <rect x=\"{:.1f}\" y=\"{:.1f}\" width=\"{:.1f}\" height=\"{:.1f}\" rx=\"6\" "
        "fill=\"{}\" fill-opacity=\"0.92\" stroke=\"{}\" stroke-width=\"1.5\" filter=\"url(#tooltip-shadow)\"/>\n",
        tooltip_x, tooltip_y, tooltip_w, tooltip_h, card_bg, card_border);

    std::string x_unit_str = hist->x_unit;
    std::string header_txt =
        std::format("{} = {:.4f}{}", x_label, pos.target_x, x_unit_str.empty() ? "" : " " + x_unit_str);
    auto paren = x_label.find(" (");
    if (paren != std::string::npos) {
      header_txt = std::format("{} = {:.4f}", x_label.substr(0, paren), pos.target_x);
    }

    svg << renderTextAsPath(header_txt, tooltip_x + static_cast<real_t>(12.0), tooltip_y + static_cast<real_t>(20.0),
                            static_cast<real_t>(14.0) * config->font_scale, TextAnchor::Start, fg_color,
                            config->use_native_text);

    real_t cur_y = tooltip_y + static_cast<real_t>(42.0);
    for (const auto &[name, val, col] : hover_values) {
      svg << std::format("  <circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"5\" fill=\"{}\"/>\n",
                         tooltip_x + static_cast<real_t>(18.0), cur_y - static_cast<real_t>(4.0), col);
      std::string line_txt = std::format("{}: {:.4f}", name, val);
      if (name == best_key) {
        line_txt += " (nearest)";
      }
      svg << renderTextAsPath(line_txt, tooltip_x + static_cast<real_t>(30.0), cur_y,
                              static_cast<real_t>(13.0) * config->font_scale, TextAnchor::Start, fg_color,
                              config->use_native_text);
      cur_y += static_cast<real_t>(22.0);
    }
  }

  void drawHoverTooltip() {
    if (hover == nullptr || !hover->active || hover->widget_width <= static_cast<real_t>(0.0) ||
        hover->widget_height <= static_cast<real_t>(0.0) || xs == nullptr) {
      return;
    }
    real_t svg_aspect = kW / kH;
    real_t widget_aspect = hover->widget_width / hover->widget_height;
    real_t scale = static_cast<real_t>(1.0);
    real_t offset_x = static_cast<real_t>(0.0);
    real_t offset_y = static_cast<real_t>(0.0);

    if (widget_aspect > svg_aspect) {
      scale = hover->widget_height / kH;
      offset_x = (hover->widget_width - kW * scale) / static_cast<real_t>(2.0);
    } else {
      scale = hover->widget_width / kW;
      offset_y = (hover->widget_height - kH * scale) / static_cast<real_t>(2.0);
    }

    real_t mouse_sx = (hover->mouse_x - offset_x) / scale;
    real_t mouse_sy = (hover->mouse_y - offset_y) / scale;

    if (mouse_sx < px0 || mouse_sx > px1 || mouse_sy < py0 || mouse_sy > py1) {
      return;
    }

    NearestPoint best = findNearestPoint(mouse_sx, mouse_sy);

    std::size_t idx = best.index;
    real_t target_x = xs->at(idx);
    real_t sx_data = mapValue(target_x, xScale.min, xScale.max, px0, px1);

    svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                       "stroke=\"{}\" stroke-width=\"1.5\" stroke-dasharray=\"4,4\"/>\n",
                       sx_data, py0, sx_data, py1, config->axis_color());

    std::vector<std::tuple<std::string, real_t, std::string>> hover_values;
    std::size_t color_idx = 0;
    real_t snapped_sy_data = static_cast<real_t>(-1.0);
    for (const auto &[key, value] : partials) {
      if (idx < value.size()) {
        const std::string col = color(color_idx++, config->palette);
        real_t y_val = value.at(idx);
        real_t sy_data = mapValue(y_val, yScale.min, yScale.max, py1, py0);

        if (key == best.key) {
          snapped_sy_data = sy_data;
        }

        svg << std::format(
            "  <circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"6\" fill=\"{}\" stroke=\"{}\" stroke-width=\"2\"/>\n", sx_data,
            sy_data, col, config->bg_color());

        hover_values.emplace_back(key, y_val, col);
      }
    }

    if (snapped_sy_data >= py0 && snapped_sy_data <= py1) {
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"1.5\" stroke-dasharray=\"4,4\"/>\n",
                         px0, snapped_sy_data, sx_data, snapped_sy_data, config->axis_color());
    }

    drawTooltipBox({sx_data, snapped_sy_data, target_x}, best.key, hover_values);
  }

  std::string getResult() {
    svg << "</svg>\n";
    return svg.str();
  }
};

} // namespace detail

/**
 * @brief Renders a `Histogram` as a self-contained SVG string.
 */
inline std::string renderHistogramAsSvg(const correlation::analysis::Histogram &hist, const PlotConfig &config = {},
                                        const HoverInfo &hover = {},
                                        const std::map<std::string, real_t> &weights = {}) {
  detail::SvgHistogramRenderer renderer(hist, config, hover, weights);
  if (renderer.hasNoData()) {
    return renderer.renderNoData();
  }
  renderer.computeScales();
  renderer.writeHeader();
  renderer.drawGridAndAxes();
  renderer.drawEmphasisLines();
  renderer.drawAreaFills();
  auto legend = renderer.drawPolylines();
  renderer.drawMarkers();
  renderer.drawLegend(legend);
  renderer.drawTitlesAndLabels();
  renderer.drawHoverTooltip();
  return renderer.getResult();
}

} // namespace correlation::plotters
