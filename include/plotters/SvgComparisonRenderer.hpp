/**
 * @file SvgComparisonRenderer.hpp
 * @brief Multi-histogram comparison rendering logic for the SVG plotter.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "plotters/PlotTypes.hpp"

#include <algorithm>
#include <cstddef>
#include <format>
#include <ranges>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace correlation::plotters {
namespace detail {

struct SvgComparisonRenderer {
  const std::vector<LabeledHistogram> *datasets = nullptr;
  const std::string *partial_key = nullptr;
  const PlotConfig *config = nullptr;
  const HoverInfo *hover = nullptr;

  real_t kWidth = static_cast<real_t>(0.0);
  real_t kHeight = static_cast<real_t>(0.0);
  real_t px0 = static_cast<real_t>(0.0);
  real_t px1 = static_cast<real_t>(0.0);
  real_t py0 = static_cast<real_t>(0.0);
  real_t py1 = static_cast<real_t>(0.0);

  real_t raw_x_min = static_cast<real_t>(1e30);
  real_t raw_x_max = static_cast<real_t>(-1e30);
  real_t raw_y_min = static_cast<real_t>(0.0);
  real_t raw_y_max = static_cast<real_t>(-1e30);

  detail::NiceScale xScale;
  detail::NiceScale yScale;

  std::string x_label;
  std::string y_label;

  std::ostringstream svg;
  std::vector<std::pair<std::string, std::string>> legend;

  SvgComparisonRenderer(const std::vector<LabeledHistogram> &datasets_ref, const std::string &partial_key_ref,
                        const PlotConfig &config_ref, const HoverInfo &hover_ref)
      : datasets(&datasets_ref), partial_key(&partial_key_ref), config(&config_ref), hover(&hover_ref) {}

  void initializeLabelsAndLayout() {
    const auto &ref = *datasets->front().hist;
    x_label = ref.x_label.empty() ? "x" : ref.x_label;
    y_label = ref.y_label.empty() ? "y" : ref.y_label;

    if (!ref.x_unit.empty()) {
      x_label += std::format(" ({})", ref.x_unit);
    }
    if (!ref.y_unit.empty()) {
      std::string y_unit = ref.y_unit;
      if (y_unit == "Å^-1") {
        y_unit = "Å⁻¹";
      }
      y_label += std::format(" ({})", y_unit);
    }

    kWidth = config->effective_width();
    kHeight = config->effective_height();
    const real_t kLeft = static_cast<real_t>(100.0);
    const real_t kRight = static_cast<real_t>(40.0);
    const real_t kTop = static_cast<real_t>(50.0);
    const real_t kBot = static_cast<real_t>(90.0);

    px0 = kLeft;
    px1 = kWidth - kRight;
    py0 = kTop;
    py1 = kHeight - kBot;
  }

  void computeGlobalRanges() {
    for (const auto &dataset : *datasets) {
      const auto &partials =
          dataset.hist->smoothed_partials.empty() ? dataset.hist->partials : dataset.hist->smoothed_partials;
      auto iterator = partials.find(*partial_key);
      if (iterator == partials.end()) {
        continue;
      }
      const auto &x_values = dataset.hist->bins;
      const auto &y_values = iterator->second;
      if (!x_values.empty()) {
        raw_x_min = std::min(raw_x_min, x_values.front());
        raw_x_max = std::max(raw_x_max, x_values.back());
      }
      for (real_t y_value : y_values) {
        raw_y_max = std::max(raw_y_max, y_value);
        raw_y_min = std::min(raw_y_min, y_value);
      }
    }
  }

  bool hasValidRanges() const { return raw_x_max > raw_x_min && raw_y_max > raw_y_min; }

  std::string renderNoComparisonData() const {
    if (config->use_native_text) {
      return std::format(
          "<svg width='{0:.0f}' height='{1:.0f}' xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {0:.0f} {1:.0f}\">"
          "<rect width=\"100%\" height=\"100%\" fill=\"{2}\"/>"
          "<text x=\"{3:.1f}\" y=\"{4:.1f}\" font-family=\"'Outfit', 'Plus Jakarta Sans', 'Inter', 'Roboto', 'Helvetica Neue', "
          "sans-serif\" font-size=\"{5:.1f}\" text-anchor=\"middle\" fill=\"{6}\">No comparison data</text></svg>",
          kWidth, kHeight, config->bg_color(), kWidth / static_cast<real_t>(2.0),
          kHeight / static_cast<real_t>(2.0) + static_cast<real_t>(8.0), static_cast<real_t>(24.0) * config->font_scale,
          config->text_color());
    }
    std::string no_data_path = Roboto::instance().render("No comparison data", kWidth / static_cast<real_t>(2.0),
                                                         kHeight / static_cast<real_t>(2.0) + static_cast<real_t>(8.0),
                                                         static_cast<real_t>(24.0) * config->font_scale, "middle");
    return std::format(
        "<svg width='{0:.0f}' height='{1:.0f}' xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {0:.0f} {1:.0f}\">"
        "<rect width=\"100%\" height=\"100%\" fill=\"{2}\"/>"
        "<path d=\"{3}\" fill=\"{4}\" fill-rule=\"evenodd\" stroke=\"none\"/></svg>",
        kWidth, kHeight, config->bg_color(), no_data_path, config->text_color());
  }

  void computeScales() {
    real_t y_padding = (raw_y_max - raw_y_min) * static_cast<real_t>(0.05);
    raw_y_max += y_padding;

    xScale = detail::NiceScale(detail::DataRange{.min = raw_x_min, .max = raw_x_max}, 11);
    yScale = detail::NiceScale(detail::DataRange{.min = raw_y_min, .max = raw_y_max}, 8);
  }

  void writeHeaderAndDefs() {
    svg << std::format("<svg width='{:.0f}' height='{:.0f}' xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {} {}\" "
                       "shape-rendering=\"geometricPrecision\" text-rendering=\"geometricPrecision\">\n",
                       kWidth, kHeight, kWidth, kHeight);
    svg << "  <defs>\n"
        << "    <filter id=\"tooltip-shadow\" x=\"-10%\" y=\"-10%\" width=\"120%\" height=\"120%\">\n"
        << "      <feDropShadow dx=\"2\" dy=\"4\" stdDeviation=\"4\" flood-color=\"#000000\" flood-opacity=\"0.15\"/>\n"
        << "    </filter>\n";

    std::size_t color_idx = 0;
    for (const auto &dataset : *datasets) {
      std::string col = detail::color(color_idx, config->palette);
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
    for (real_t y_value : yScale.ticks) {
      real_t spy = detail::mapValue(y_value, yScale.min, yScale.max, py1, py0);
      if (config->show_grid) {
        svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                           "stroke=\"{}\" stroke-width=\"0.8\" stroke-dasharray=\"3,3\"/>\n",
                           px0, spy, px1, spy, config->grid_color());
      }
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                         px0 - static_cast<real_t>(8.0), spy, px0, spy, config->axis_color());
      svg << renderTextAsPath(detail::fmtScientific(y_value), px0 - static_cast<real_t>(15.0),
                              spy + static_cast<real_t>(7.0), static_cast<real_t>(20.0) * config->font_scale,
                              TextAnchor::End, config->text_color(), config->use_native_text);
    }

    for (real_t x_value : xScale.ticks) {
      real_t spx = detail::mapValue(x_value, xScale.min, xScale.max, px0, px1);
      if (config->show_grid) {
        svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                           "stroke=\"{}\" stroke-width=\"0.8\" stroke-dasharray=\"3,3\"/>\n",
                           spx, py0, spx, py1, config->grid_color());
      }
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                         spx, py1, spx, py1 + static_cast<real_t>(8.0), config->axis_color());
      svg << renderTextAsPath(detail::fmtScientific(x_value), spx, py1 + static_cast<real_t>(25.0),
                              static_cast<real_t>(20.0) * config->font_scale, TextAnchor::Middle, config->text_color(),
                              config->use_native_text);
    }

    svg << std::format("  <rect x=\"{:.1f}\" y=\"{:.1f}\" width=\"{:.1f}\" height=\"{:.1f}\" "
                       "fill=\"none\" stroke=\"{}\" stroke-width=\"1.5\"/>\n",
                       px0, py0, px1 - px0, py1 - py0, config->axis_color());
  }

  void drawAreaFills() {
    if (config->fill_area) {
      std::size_t fill_ci = 0;
      for (const auto &dataset : *datasets) {
        const auto &partials =
            dataset.hist->smoothed_partials.empty() ? dataset.hist->partials : dataset.hist->smoothed_partials;
        auto iterator = partials.find(*partial_key);
        if (iterator != partials.end()) {
          const auto &x_values = dataset.hist->bins;
          const auto &y_values = iterator->second;
          std::string grad_id = std::format("area-grad-{}", fill_ci);
          std::size_t n_points = std::min(x_values.size(), y_values.size());
          if (n_points > 1) {
            svg << std::format("  <polygon fill=\"url(#{})\" stroke=\"none\" points=\"", grad_id);
            real_t sp_x_start = detail::mapValue(x_values[0], xScale.min, xScale.max, px0, px1);
            svg << std::format("{:.2f},{:.2f} ", sp_x_start, py1);
            for (std::size_t i = 0; i < n_points; ++i) {
              real_t sp_x = detail::mapValue(x_values[i], xScale.min, xScale.max, px0, px1);
              real_t sp_y = detail::mapValue(y_values[i], yScale.min, yScale.max, py1, py0);
              real_t sp_y_clamped = std::min(sp_y, py1);
              svg << std::format("{:.2f},{:.2f} ", sp_x, sp_y_clamped);
            }
            real_t sp_x_end = detail::mapValue(x_values[n_points - 1], xScale.min, xScale.max, px0, px1);
            svg << std::format("{:.2f},{:.2f}", sp_x_end, py1);
            svg << "\" />\n";
          }
        }
        fill_ci++;
      }
    }
  }

  void drawPolylines() {
    std::size_t color_idx = 0;
    for (const auto &dataset : *datasets) {
      const auto &partials =
          dataset.hist->smoothed_partials.empty() ? dataset.hist->partials : dataset.hist->smoothed_partials;
      auto partial_iter = partials.find(*partial_key);
      if (partial_iter == partials.end()) {
        continue;
      }

      const auto &x_values = dataset.hist->bins;
      const auto &y_values = partial_iter->second;
      const std::string col = detail::color(color_idx++, config->palette);

      svg << std::format(
          R"(  <polyline fill="none" stroke="{}" stroke-width="{:.1f}" stroke-linejoin="round" points=")", col,
          config->line_width);
      std::size_t n_points = std::min(x_values.size(), y_values.size());
      for (std::size_t i = 0; i < n_points; ++i) {
        real_t sp_x = detail::mapValue(x_values[i], xScale.min, xScale.max, px0, px1);
        real_t sp_y = detail::mapValue(y_values[i], yScale.min, yScale.max, py1, py0);
        svg << std::format("{:.2f},{:.2f} ", sp_x, sp_y);
      }
      svg << "\" />\n";
      legend.emplace_back(dataset.label, col);
    }
  }

  void drawMarkers() {
    if (config->show_markers) {
      std::size_t marker_color_idx = 0;
      for (const auto &dataset : *datasets) {
        const auto &partials =
            dataset.hist->smoothed_partials.empty() ? dataset.hist->partials : dataset.hist->smoothed_partials;
        auto partial_iter = partials.find(*partial_key);
        if (partial_iter != partials.end()) {
          const auto &x_values = dataset.hist->bins;
          const auto &y_values = partial_iter->second;
          const std::string col = detail::color(marker_color_idx++, config->palette);
          std::size_t n_points = std::min(x_values.size(), y_values.size());
          for (std::size_t i = 0; i < n_points; ++i) {
            real_t sp_x = detail::mapValue(x_values[i], xScale.min, xScale.max, px0, px1);
            real_t sp_y = detail::mapValue(y_values[i], yScale.min, yScale.max, py1, py0);
            svg << std::format("  <circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"3.5\" fill=\"{}\" stroke=\"none\"/>\n", sp_x,
                               sp_y, col);
          }
        } else {
          marker_color_idx++;
        }
      }
    }
  }

  void drawLegend() {
    if (config->show_legend) {
      real_t legend_x = px1 - static_cast<real_t>(15.0);
      real_t legend_y = py0 + static_cast<real_t>(25.0);
      for (auto [label, color] : std::views::reverse(legend)) {
        svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                           "stroke=\"{}\" stroke-width=\"4.0\"/>\n",
                           legend_x - static_cast<real_t>(40.0), legend_y, legend_x - static_cast<real_t>(10.0),
                           legend_y, color);
        svg << renderTextAsPath(label, legend_x - static_cast<real_t>(45.0), legend_y + static_cast<real_t>(6.0),
                                static_cast<real_t>(18.0) * config->font_scale, TextAnchor::End, config->text_color(),
                                config->use_native_text);
        legend_y += static_cast<real_t>(28.0);
      }
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

  struct NearestPoint {
    std::size_t index = 0;
    std::string label;
  };

  NearestPoint findNearestPoint(real_t sp_x, real_t sp_y) const {
    real_t min_dist_sq = static_cast<real_t>(1e30);
    NearestPoint best;

    for (const auto &dataset : *datasets) {
      const auto &partials =
          dataset.hist->smoothed_partials.empty() ? dataset.hist->partials : dataset.hist->smoothed_partials;
      auto pit = partials.find(*partial_key);
      if (pit == partials.end()) {
        continue;
      }
      const auto &x_values = dataset.hist->bins;
      const auto &y_values = pit->second;
      std::size_t n_points = std::min(x_values.size(), y_values.size());
      for (std::size_t i = 0; i < n_points; ++i) {
        real_t sp_x_i = detail::mapValue(x_values[i], xScale.min, xScale.max, px0, px1);
        real_t sp_y_i = detail::mapValue(y_values[i], yScale.min, yScale.max, py1, py0);
        real_t dx_i = sp_x_i - sp_x;
        real_t dy_i = sp_y_i - sp_y;
        real_t dist_sq = dx_i * dx_i + dy_i * dy_i;
        if (dist_sq < min_dist_sq) {
          min_dist_sq = dist_sq;
          best.index = i;
          best.label = dataset.label;
        }
      }
    }
    return best;
  }

  std::vector<std::tuple<std::string, real_t, std::string>>
  collectHoverValuesAndDrawMarkers(std::size_t idx, const std::string &best_label, real_t sx_data,
                                   real_t &snapped_sy_data) {
    std::vector<std::tuple<std::string, real_t, std::string>> hover_values;
    std::size_t color_idx = 0;
    for (const auto &dataset : *datasets) {
      const auto &partials =
          dataset.hist->smoothed_partials.empty() ? dataset.hist->partials : dataset.hist->smoothed_partials;
      auto pit = partials.find(*partial_key);
      if (pit != partials.end() && idx < pit->second.size()) {
        const std::string col = detail::color(color_idx++, config->palette);
        real_t y_val = pit->second[idx];
        real_t sy_data = detail::mapValue(y_val, yScale.min, yScale.max, py1, py0);

        if (dataset.label == best_label) {
          snapped_sy_data = sy_data;
        }

        svg << std::format(
            "  <circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"6\" fill=\"{}\" stroke=\"{}\" stroke-width=\"2\"/>\n", sx_data,
            sy_data, col, config->bg_color());

        hover_values.emplace_back(dataset.label, y_val, col);
      } else {
        color_idx++;
      }
    }
    return hover_values;
  }

  struct TooltipPosition {
    real_t sx_data = static_cast<real_t>(0.0);
    real_t snapped_sy_data = static_cast<real_t>(0.0);
    real_t target_x = static_cast<real_t>(0.0);
  };

  void drawTooltipBox(const TooltipPosition &pos, const std::string &best_label,
                      const std::vector<std::tuple<std::string, real_t, std::string>> &hover_values) {
    if (hover_values.empty()) {
      return;
    }
    real_t tooltip_w = static_cast<real_t>(200.0);
    real_t tooltip_h = static_cast<real_t>(35.0) + static_cast<real_t>(22.0) * static_cast<real_t>(hover_values.size());

    real_t tooltip_x = (pos.sx_data < kWidth / static_cast<real_t>(2.0))
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

    const auto &ref = *datasets->front().hist;
    std::string x_unit_str = ref.x_unit;
    std::string header_txt =
        std::format("{} = {:.4f}{}", x_label, pos.target_x, x_unit_str.empty() ? "" : " " + x_unit_str);
    auto paren = x_label.find(" (");
    if (paren != std::string::npos) {
      header_txt = std::format("{} = {:.4f}", x_label.substr(0, paren), pos.target_x);
    }

    svg << renderTextAsPath(header_txt, tooltip_x + static_cast<real_t>(12.0), tooltip_y + static_cast<real_t>(24.0),
                            static_cast<real_t>(14.0) * config->font_scale, TextAnchor::Start, fg_color,
                            config->use_native_text);

    real_t cur_y = tooltip_y + static_cast<real_t>(46.0);
    for (const auto &[name, val, col] : hover_values) {
      svg << std::format("  <circle cx=\"{:.1f}\" cy=\"{:.1f}\" r=\"5\" fill=\"{}\"/>\n",
                         tooltip_x + static_cast<real_t>(18.0), cur_y - static_cast<real_t>(5.0), col);
      std::string line_txt = std::format("{}: {:.4f}", name, val);
      if (name == best_label) {
        line_txt += " (nearest)";
      }
      svg << renderTextAsPath(line_txt, tooltip_x + static_cast<real_t>(30.0), cur_y,
                              static_cast<real_t>(13.0) * config->font_scale, TextAnchor::Start, fg_color,
                              config->use_native_text);
      cur_y += static_cast<real_t>(22.0);
    }
  }

  void drawHoverTooltip() {
    if (!hover->active || hover->widget_width <= static_cast<real_t>(0.0) ||
        hover->widget_height <= static_cast<real_t>(0.0)) {
      return;
    }
    real_t svg_aspect = kWidth / kHeight;
    real_t widget_aspect = hover->widget_width / hover->widget_height;
    real_t scale = static_cast<real_t>(1.0);
    real_t offset_x = static_cast<real_t>(0.0);
    real_t offset_y = static_cast<real_t>(0.0);

    if (widget_aspect > svg_aspect) {
      scale = hover->widget_height / kHeight;
      offset_x = (hover->widget_width - kWidth * scale) / static_cast<real_t>(2.0);
    } else {
      scale = hover->widget_width / kWidth;
      offset_y = (hover->widget_height - kHeight * scale) / static_cast<real_t>(2.0);
    }

    real_t sp_x = (hover->mouse_x - offset_x) / scale;
    real_t sp_y = (hover->mouse_y - offset_y) / scale;

    if (sp_x < px0 || sp_x > px1 || sp_y < py0 || sp_y > py1) {
      return;
    }

    auto [best_idx, best_label] = findNearestPoint(sp_x, sp_y);

    const auto &x_bins = datasets->front().hist->bins;
    if (x_bins.empty() || best_idx >= x_bins.size()) {
      return;
    }

    real_t target_x = x_bins[best_idx];
    real_t sx_data = detail::mapValue(target_x, xScale.min, xScale.max, px0, px1);

    svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                       "stroke=\"{}\" stroke-width=\"1.5\" stroke-dasharray=\"4,4\"/>\n",
                       sx_data, py0, sx_data, py1, config->axis_color());

    real_t snapped_sy_data = static_cast<real_t>(-1.0);
    auto hover_values = collectHoverValuesAndDrawMarkers(best_idx, best_label, sx_data, snapped_sy_data);

    if (snapped_sy_data >= py0 && snapped_sy_data <= py1) {
      svg << std::format("  <line x1=\"{:.1f}\" y1=\"{:.1f}\" x2=\"{:.1f}\" y2=\"{:.1f}\" "
                         "stroke=\"{}\" stroke-width=\"1.5\" stroke-dasharray=\"4,4\"/>\n",
                         px0, snapped_sy_data, sx_data, snapped_sy_data, config->axis_color());
    }

    drawTooltipBox({sx_data, snapped_sy_data, target_x}, best_label, hover_values);
  }

  std::string getResult() {
    svg << "</svg>\n";
    return svg.str();
  }
};

} // namespace detail

/**
 * @brief Renders multiple histograms overlaid on a single plot.
 */
inline std::string renderComparisonSvg(const std::vector<LabeledHistogram> &datasets,
                                       const std::string &partial_key = "Total", const PlotConfig &config = {},
                                       const HoverInfo &hover = {}) {
  if (datasets.empty()) {
    return "";
  }
  detail::SvgComparisonRenderer renderer(datasets, partial_key, config, hover);
  renderer.initializeLabelsAndLayout();
  renderer.computeGlobalRanges();
  if (!renderer.hasValidRanges()) {
    return renderer.renderNoComparisonData();
  }
  renderer.computeScales();
  renderer.writeHeaderAndDefs();
  renderer.drawGridAndAxes();
  renderer.drawAreaFills();
  renderer.drawPolylines();
  renderer.drawMarkers();
  renderer.drawLegend();
  renderer.drawTitlesAndLabels();
  renderer.drawHoverTooltip();
  return renderer.getResult();
}

} // namespace correlation::plotters
