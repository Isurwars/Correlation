/**
 * @file PlotController.cpp
 * @brief Implementation of PlotController.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/viewmodel/PlotController.hpp"
#include "AppWindow.h"
#include "app/formatters/PlotExportService.hpp"
#include "app/formatters/PlotTableFormatter.hpp"
#include <nfd.h>

#include <algorithm>
#include <filesystem>
#include <format>
#include <map>
#include <span>
#include <utility>

namespace correlation::app {

namespace {
template <typename T> T safeParse(const slint::SharedString &str, T default_value) {
  try {
    if constexpr (std::is_same_v<T, float>) {
      return std::stof(str.data());
    } else if constexpr (std::is_same_v<T, real_t>) {
      return std::stod(str.data());
    } else {
      return default_value;
    }
  } catch (const std::exception &e) {
    return default_value;
  }
}

correlation::plotters::detail::NiceScale
computeXScale(const correlation::analysis::Histogram &hist,
              const correlation::plotters::PlotConfig &config) {
  real_t raw_x_min = hist.bins.empty() ? static_cast<real_t>(0.0) : hist.bins.front();
  real_t raw_x_max = hist.bins.empty() ? static_cast<real_t>(1.0) : hist.bins.back();
  bool strict_x = false;
  if (config.manual_x_min.has_value() && config.manual_x_max.has_value()) {
    raw_x_min = *config.manual_x_min;
    raw_x_max = *config.manual_x_max;
    strict_x = true;
  }
  return correlation::plotters::detail::NiceScale(
      correlation::plotters::detail::DataRange{.min = raw_x_min, .max = raw_x_max}, 11, strict_x);
}

correlation::plotters::detail::NiceScale
computeYScale(const correlation::analysis::Histogram &hist,
              const correlation::plotters::PlotConfig &config) {
  auto raw_y_min = static_cast<real_t>(0.0);
  auto raw_y_max = static_cast<real_t>(0.0);
  const auto &partials = hist.smoothed_partials.empty() ? hist.partials : hist.smoothed_partials;
  for (const auto &[key, value] : partials) {
    for (const real_t val : value) {
      raw_y_max = std::max(raw_y_max, val);
      raw_y_min = std::min(raw_y_min, val);
    }
  }
  const real_t y_padding = (raw_y_max - raw_y_min) * static_cast<real_t>(0.05);
  raw_y_max += y_padding;

  bool strict_y = false;
  if (config.manual_y_min.has_value() && config.manual_y_max.has_value()) {
    raw_y_min = *config.manual_y_min;
    raw_y_max = *config.manual_y_max;
    strict_y = true;
  }
  return correlation::plotters::detail::NiceScale(
      correlation::plotters::detail::DataRange{.min = raw_y_min, .max = raw_y_max}, 8, strict_y);
}

} // namespace

PlotController::PlotController(::AppWindow &window, AnalysisDispatcher &dispatcher,
                               const ProgramOptions &options)
    : window_(window), dispatcher_(dispatcher), options_(options) {
  update_timer_.start(slint::TimerMode::Repeated, std::chrono::milliseconds(1000),
                      [this]() { handleUpdateTimer(); });
}

PlotController::~PlotController() {
  update_timer_.stop();
  hover_timer_.stop();
  if (dialog_thread_.joinable()) {
    dialog_thread_.join();
  }
  if (render_thread_.joinable()) {
    render_thread_.join();
  }
}

void PlotController::handlePlotResized(PlotSize size) {
  if (std::abs(size.width - last_plot_width_) < 1.0F &&
      std::abs(size.height - last_plot_height_) < 1.0F) {
    return;
  }
  last_plot_width_ = size.width;
  last_plot_height_ = size.height;
  const int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    slint::invoke_from_event_loop([this, current_idx] { requestPlotUpdate(current_idx, true); });
  }
}

void PlotController::handleUpdateTimer() {
  const int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    requestPlotUpdate(current_idx, false);
  }
}

correlation::plotters::PlotConfig PlotController::buildPlotConfigFromUI() {
  correlation::plotters::PlotConfig config;
  config.theme = window_.get_is_dark() ? correlation::plotters::PlotConfig::Theme::Dark
                                       : correlation::plotters::PlotConfig::Theme::Light;

  const int size_preset_val = window_.get_export_config().size_preset;
  if (size_preset_val == 1) {
    config.preset_size = correlation::plotters::PlotConfig::PresetSize::SingleColumn;
  } else if (size_preset_val == 2) {
    config.preset_size = correlation::plotters::PlotConfig::PresetSize::DoubleColumn;
  } else if (size_preset_val == 3) {
    config.preset_size = correlation::plotters::PlotConfig::PresetSize::Presentation;
  } else {
    config.preset_size = correlation::plotters::PlotConfig::PresetSize::Default;
    if (last_plot_width_ > 1.0F && last_plot_height_ > 1.0F) {
      config.width = static_cast<real_t>(last_plot_width_);
      config.height = static_cast<real_t>(last_plot_height_);
    }
  }

  const int palette_val = window_.get_export_config().palette;
  switch (palette_val) {
  case 1:
    config.palette = correlation::plotters::PlotConfig::Palette::Grayscale;
    break;
  case 2:
    config.palette = correlation::plotters::PlotConfig::Palette::Viridis;
    break;
  case 3:
    config.palette = correlation::plotters::PlotConfig::Palette::Magma;
    break;
  case 4:
    config.palette = correlation::plotters::PlotConfig::Palette::Heatmap;
    break;
  case 5:
    config.palette = correlation::plotters::PlotConfig::Palette::Rainbow;
    break;
  case 6:
    config.palette = correlation::plotters::PlotConfig::Palette::Turbo;
    break;
  case 7:
    config.palette = correlation::plotters::PlotConfig::Palette::Plasma;
    break;
  case 8:
    config.palette = correlation::plotters::PlotConfig::Palette::Inferno;
    break;
  case 9:
    config.palette = correlation::plotters::PlotConfig::Palette::Cividis;
    break;
  default:
    config.palette = correlation::plotters::PlotConfig::Palette::OkabeIto;
    break;
  }

  config.font_scale = safeParse(window_.get_export_config().font_scale, 1.0F);
  if (config.font_scale <= 0.0F) {
    config.font_scale = 1.0F;
  }

  config.line_width = safeParse(window_.get_export_config().line_width, 3.0F);
  if (config.line_width <= 0.0F) {
    config.line_width = 3.0F;
  }

  config.marker_size = safeParse(window_.get_export_config().marker_size, 3.5F);
  if (config.marker_size <= 0.0F) {
    config.marker_size = 3.5F;
  }

  config.show_grid = window_.get_export_config().show_grid;
  config.show_legend = window_.get_export_config().show_legend;
  config.show_markers = window_.get_export_config().show_markers;
  config.fill_area = window_.get_export_config().fill_area;

  config.manual_x_min = zoom_x_min_;
  config.manual_x_max = zoom_x_max_;
  config.manual_y_min = zoom_y_min_;
  config.manual_y_max = zoom_y_max_;
  config.reference_lines = reference_markers_;

  return config;
}

void PlotController::populatePlotList() {
  last_rendered_index_ = -1;
  auto names = dispatcher_.getAvailableHistogramNames();

  std::map<std::string, int> priority = {
      {"g_r", 0},       {"g_r_unweighted", 1}, {"H_r", 2},      {"G_r", 3},
      {"J_r", 4},       {"S_q", 10},           {"XRD", 11},     {"PAD", 20},
      {"PAD_raw", 21},  {"DAD", 22},           {"DAD_raw", 23}, {"CN", 24},
      {"RD", 25},       {"MSD", 30},           {"VACF", 31},    {"VDOS", 32},
      {"sigma2_N", 40}, {"chi_H", 41}};

  std::ranges::sort(names, [&](const std::string &lhs, const std::string &rhs) {
    const int prio_a = priority.contains(lhs) ? priority.at(lhs) : 100;
    const int prio_b = priority.contains(rhs) ? priority.at(rhs) : 100;
    if (prio_a != prio_b) {
      return prio_a < prio_b;
    }
    return lhs < rhs;
  });

  available_plot_keys_ = names;

  auto menu_model = std::make_shared<slint::VectorModel<MenuItem>>();
  for (const auto &name : names) {
    MenuItem item;
    const correlation::analysis::Histogram *hist = dispatcher_.getHistogram(name);
    const std::string display_text = (hist != nullptr && !hist->title.empty()) ? hist->title : name;
    item.text = slint::SharedString(display_text);
    item.enabled = true;
    menu_model->push_back(item);
  }
  window_.set_plot_items(menu_model);

  // Update dynamic properties
  const auto *distribution_functions = dispatcher_.getDistributionFunctions();
  if (distribution_functions != nullptr) {
    real_t msd = distribution_functions->getDiffusionCoefficientMSD();
    if (msd > 0.0) {
      window_.set_diff_msd(slint::SharedString(std::format("{:.6f} Å²/fs", msd)));
    } else {
      window_.set_diff_msd("");
    }

    real_t vacf = distribution_functions->getDiffusionCoefficientVACF();
    if (vacf > 0.0) {
      window_.set_diff_vacf(slint::SharedString(std::format("{:.6f} Å²/fs", vacf)));
    } else {
      window_.set_diff_vacf("");
    }

    real_t tau = distribution_functions->getRelaxationTime();
    if (tau > 0.0) {
      window_.set_relaxation_time(slint::SharedString(std::format("{:.4f} fs", tau)));
    } else {
      window_.set_relaxation_time("");
    }

    real_t deb = distribution_functions->getDeborahNumber();
    if (deb > 0.0) {
      window_.set_deborah_number(slint::SharedString(std::format("{:.4f}", deb)));
    } else {
      window_.set_deborah_number("");
    }
  } else {
    window_.set_diff_msd("");
    window_.set_diff_vacf("");
    window_.set_relaxation_time("");
    window_.set_deborah_number("");
  }

  window_.set_selected_plot_index(names.empty() ? -1 : 0);
}

void PlotController::handleMouseMove(float mouse_x, float mouse_y, bool hover, float width,
                                     float height) {
  const bool actual_hover = hover;

  if (std::abs(mouse_x - last_mouse_x_) < 0.5F && std::abs(mouse_y - last_mouse_y_) < 0.5F &&
      actual_hover == mouse_hover_ && std::abs(width - last_plot_width_) < 1e-2F &&
      std::abs(height - last_plot_height_) < 1e-2F) {
    return;
  }

  const bool hover_changed = (actual_hover != mouse_hover_);

  last_mouse_x_ = mouse_x;
  last_mouse_y_ = mouse_y;
  mouse_hover_ = actual_hover;
  last_plot_width_ = width;
  last_plot_height_ = height;

  if (!actual_hover) {
    window_.set_hover_coord_text("");
  } else {
    const int current_idx = window_.get_selected_plot_index();
    if (current_idx >= 0 && std::cmp_less(current_idx, available_plot_keys_.size())) {
      const std::string &name = available_plot_keys_[current_idx];
      const auto *hist = dispatcher_.getHistogram(name);
      if (hist != nullptr && !hist->bins.empty()) {
        auto config = buildPlotConfigFromUI();
        auto geom = correlation::plotters::detail::getViewportGeometry(config);
        auto [svg_x, svg_y] =
            correlation::plotters::detail::screenToSvg(mouse_x, mouse_y, width, height, config);
        if (svg_x >= geom.px0 && svg_x <= geom.px1 && svg_y >= geom.py0 && svg_y <= geom.py1) {
          auto x_scale = computeXScale(*hist, config);
          auto y_scale = computeYScale(*hist, config);
          auto [d_x, d_y] = correlation::plotters::detail::screenToData(
              mouse_x, mouse_y, width, height, config, x_scale, y_scale);
          const std::string text = std::format("X: {:.4g}  Y: {:.4g}", d_x, d_y);
          window_.set_hover_coord_text(slint::SharedString(text));
        } else {
          window_.set_hover_coord_text("");
        }
      } else {
        window_.set_hover_coord_text("");
      }
    } else {
      window_.set_hover_coord_text("");
    }
  }

  const int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    requestPlotUpdate(current_idx, hover_changed || !actual_hover);
  }
}

void PlotController::requestPlotUpdate(int index, bool immediate) {
  if (index < 0 || std::cmp_greater_equal(index, available_plot_keys_.size())) {
    return;
  }

  const std::string &name = available_plot_keys_[index];
  const correlation::analysis::Histogram *hist = dispatcher_.getHistogram(name);
  if (hist == nullptr) {
    return;
  }

  updateTableData(hist);
  updateCurveToggleItems(hist);

  if (immediate) {
    needs_redraw_ = true;
  }

  const correlation::plotters::PlotConfig config = buildPlotConfigFromUI();
  correlation::plotters::HoverInfo hover;
  hover.active = mouse_hover_;
  hover.mouse_x = last_mouse_x_;
  hover.mouse_y = last_mouse_y_;
  hover.widget_width = last_plot_width_;
  hover.widget_height = last_plot_height_;

  if (isPlotCacheHit(index, config, hover) && !needs_redraw_) {
    return;
  }

  if (is_rendering_) {
    render_pending_ = true;
    pending_plot_index_ = index;
    return;
  }

  is_rendering_ = true;
  needs_redraw_ = false;
  last_rendered_index_ = index;
  last_pinned_runs_count_ = series_manager_.getPinnedRunsCount();
  last_config_ = config;
  last_hover_ = hover;

  RenderTaskData data;
  data.active_hist = *hist;
  data.config = config;
  data.config.show_difference_curve = series_manager_.shouldShowDifference();
  data.hover = hover;
  data.ashcroft_weights = dispatcher_.getAshcroftWeights();
  data.curve_visibility = series_manager_.getCurveVisibilityMap();
  data.custom_curve_colors = series_manager_.getCustomColors();

  data.comparison_hists.push_back({.label = "Current", .hist = &data.active_hist});
  for (const auto &pinned_run : series_manager_.getPinnedRuns()) {
    auto hist_it = pinned_run.histograms.find(name);
    if (hist_it != pinned_run.histograms.end()) {
      const bool vis = series_manager_.isCurveVisible(pinned_run.label, true);
      data.comparison_hists.push_back(
          {.label = pinned_run.label, .hist = &hist_it->second, .style = {.visible = vis}});
    }
  }

  if (render_thread_.joinable()) {
    render_thread_.join();
  }

  if (window_.get_selected_plot_index() != index) {
    window_.set_selected_plot_index(index);
  }

  executePlotRender(std::move(data));
}

void PlotController::updateTableData(const correlation::analysis::Histogram *hist) {
  const auto model = PlotTableFormatter::formatTable(hist);
  window_.set_table_headers(model.headers);
  window_.set_table_rows(model.rows);
}

void PlotController::handleSetCurveColor(int curve_id, const slint::SharedString &color_hex) {
  series_manager_.setCustomColor(curve_id, std::string(color_hex.data()));
  const int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    requestPlotUpdate(current_idx, true);
  }
}

void PlotController::updateCurveToggleItems(const correlation::analysis::Histogram *hist) {
  const auto active_config = buildPlotConfigFromUI();
  window_.set_curve_toggle_items(series_manager_.generateToggleItems(hist, active_config));
}

bool PlotController::isPlotCacheHit(int index, const correlation::plotters::PlotConfig &config,
                                    const correlation::plotters::HoverInfo &hover) const {
  return (index == last_rendered_index_ &&
          series_manager_.getPinnedRunsCount() == last_pinned_runs_count_ &&
          config.theme == last_config_.theme && config.preset_size == last_config_.preset_size &&
          std::abs(config.width - last_config_.width) < 1e-2 &&
          std::abs(config.height - last_config_.height) < 1e-2 &&
          config.palette == last_config_.palette &&
          std::abs(config.font_scale - last_config_.font_scale) < 1e-4 &&
          std::abs(config.line_width - last_config_.line_width) < 1e-4 &&
          std::abs(config.marker_size - last_config_.marker_size) < 1e-4 &&
          config.show_grid == last_config_.show_grid &&
          config.show_legend == last_config_.show_legend &&
          config.show_markers == last_config_.show_markers &&
          config.fill_area == last_config_.fill_area &&
          config.manual_x_min == last_config_.manual_x_min &&
          config.manual_x_max == last_config_.manual_x_max &&
          config.manual_y_min == last_config_.manual_y_min &&
          config.manual_y_max == last_config_.manual_y_max &&
          config.reference_lines.size() == last_config_.reference_lines.size() &&
          hover.active == last_hover_.active &&
          std::abs(hover.mouse_x - last_hover_.mouse_x) < 1e-2 &&
          std::abs(hover.mouse_y - last_hover_.mouse_y) < 1e-2 &&
          std::abs(hover.widget_width - last_hover_.widget_width) < 1e-2 &&
          std::abs(hover.widget_height - last_hover_.widget_height) < 1e-2);
}

void PlotController::executePlotRender(RenderTaskData data) {
  if (render_thread_.joinable()) {
    render_thread_.join();
  }
  render_thread_ = std::jthread([this, data = std::move(data)]() mutable {
    std::string svg;
    if (data.comparison_hists.size() <= 1) {
      svg = correlation::plotters::renderHistogramAsSvg(
          data.active_hist, data.config, data.hover, data.ashcroft_weights, data.curve_visibility,
          data.custom_curve_colors);
    } else {
      const std::string key = PlotExportService::getComparisonKey(&data.active_hist);
      svg = correlation::plotters::renderComparisonSvg(data.comparison_hists, key, data.config,
                                                       data.hover);
    }

    slint::invoke_from_event_loop([this, svg = std::move(svg)]() {
      const auto *svg_bytes = reinterpret_cast<const uint8_t *>(svg.data());
      auto img = slint::private_api::load_image_from_embedded_data(
          std::span<const uint8_t>(svg_bytes, svg.size()), "svg");
      window_.set_preview_plot(img);

      is_rendering_ = false;

      if (render_pending_) {
        render_pending_ = false;
        requestPlotUpdate(pending_plot_index_, true);
      }
    });
  });
}

void PlotController::handleSavePlot() {
  if (dialog_active_.exchange(true)) {
    return;
  }

  const int index = window_.get_selected_plot_index();
  if (index < 0 || std::cmp_greater_equal(index, available_plot_keys_.size())) {
    dialog_active_.store(false);
    return;
  }
  const std::string name = available_plot_keys_[index];
  const correlation::analysis::Histogram *hist = dispatcher_.getHistogram(name);
  if (hist == nullptr) {
    dialog_active_.store(false);
    return;
  }

  std::string default_path = options_.output_file_base;
  if (!default_path.empty()) {
    default_path += "_" + name;
  } else {
    default_path = name;
  }

  const std::filesystem::path def_path(default_path);
  std::string default_dir = def_path.parent_path().string();
  std::string default_name = def_path.filename().string();

  if (dialog_thread_.joinable()) {
    dialog_thread_.join();
  }

  dialog_thread_ = std::jthread([this, default_dir = std::move(default_dir),
                                 default_name = std::move(default_name), hist, name]() {
    std::array<nfdfilteritem_t, 2> filter_list = {{{
                                                       .name = "SVG Image",
                                                       .spec = "svg",
                                                   },
                                                   {
                                                       .name = "PDF Document",
                                                       .spec = "pdf",
                                                   }}};
    const nfdfiltersize_t filter_count = filter_list.size();

    nfdchar_t *out_path = nullptr;
    const nfdresult_t result =
        NFD_SaveDialogU8(&out_path, filter_list.data(), filter_count,
                         default_dir.empty() ? nullptr : default_dir.c_str(),
                         default_name.empty() ? nullptr : default_name.c_str());

    if (result == NFD_OKAY) {
      std::string filepath(out_path);
      NFD_FreePathU8(out_path);
      slint::invoke_from_event_loop([this, filepath = std::move(filepath), hist, name]() {
        dialog_active_.store(false);
        executeSavePlot(filepath, hist, name);
      });
    } else if (result == NFD_CANCEL) {
      slint::invoke_from_event_loop([this]() {
        dialog_active_.store(false);
        window_.set_analysis_status_text(slint::SharedString(AppDefaults::MSG_SAVE_CANCELLED));
      });
    } else {
      std::string error_msg = "Error: ";
      error_msg += NFD_GetError();
      slint::invoke_from_event_loop([this, error_msg = std::move(error_msg)]() {
        dialog_active_.store(false);
        window_.set_analysis_status_text(slint::SharedString(error_msg));
      });
    }
  });
}

void PlotController::executeSavePlot(const std::string &filepath,
                                     const correlation::analysis::Histogram *hist,
                                     const std::string &name) {
  if (hist == nullptr) {
    return;
  }

  correlation::plotters::PlotConfig config = buildPlotConfigFromUI();
  config.use_native_text = true;

  if (config.preset_size == correlation::plotters::PlotConfig::PresetSize::Default) {
    config.width = 1200.0;
    config.height = 900.0;
  }

  std::expected<void, std::string> result;
  if (series_manager_.getPinnedRuns().empty()) {
    result = PlotExportService::exportHistogram(filepath, *hist, config,
                                                dispatcher_.getAshcroftWeights());
  } else {
    std::vector<correlation::plotters::LabeledHistogram> datasets;
    datasets.reserve(series_manager_.getPinnedRunsCount() + 1);
    datasets.push_back({.label = "Current", .hist = hist});
    for (const auto &pinned_run : series_manager_.getPinnedRuns()) {
      auto hist_it = pinned_run.histograms.find(name);
      if (hist_it != pinned_run.histograms.end()) {
        datasets.push_back({.label = pinned_run.label, .hist = &hist_it->second});
      }
    }
    const std::string comp_key = PlotExportService::getComparisonKey(hist);
    result = PlotExportService::exportComparison(filepath, datasets, comp_key, config);
  }

  if (result.has_value()) {
    window_.set_analysis_status_text(slint::SharedString(name + " plot saved successfully."));
  } else {
    window_.set_analysis_status_text(slint::SharedString("Failed to save plot: " + result.error()));
  }
}

void PlotController::handlePinRun() {
  const auto &hists = dispatcher_.getHistograms();
  if (hists.empty()) {
    return;
  }

  series_manager_.pinCurrentRun(hists);

  slint::invoke_from_event_loop([this]() {
    window_.set_pinned_runs_count(static_cast<int>(series_manager_.getPinnedRunsCount()));

    const int current_idx = window_.get_selected_plot_index();
    if (current_idx >= 0) {
      requestPlotUpdate(current_idx, true);
    }
  });
}

void PlotController::handleClearPinnedRuns() {
  series_manager_.clearPinnedRuns();

  slint::invoke_from_event_loop([this]() {
    window_.set_pinned_runs_count(0);

    const int current_idx = window_.get_selected_plot_index();
    if (current_idx >= 0) {
      requestPlotUpdate(current_idx, true);
    }
  });
}

void PlotController::handleToggleCurveVisibility(int curve_id, bool visible) {
  series_manager_.setCurveVisible(curve_id, visible);
  slint::invoke_from_event_loop([this]() {
    const int current_idx = window_.get_selected_plot_index();
    if (current_idx >= 0) {
      requestPlotUpdate(current_idx, true);
    }
  });
}

void PlotController::handleToggleAllCurves(bool visible) {
  series_manager_.setAllCurvesVisible(visible);
  slint::invoke_from_event_loop([this]() {
    const int current_idx = window_.get_selected_plot_index();
    if (current_idx >= 0) {
      requestPlotUpdate(current_idx, true);
    }
  });
}

void PlotController::handleToggleDifferencePlot(bool show_difference) {
  series_manager_.setShowDifference(show_difference);
  slint::invoke_from_event_loop([this]() {
    const int current_idx = window_.get_selected_plot_index();
    if (current_idx >= 0) {
      requestPlotUpdate(current_idx, true);
    }
  });
}

void PlotController::handleZoomRect(float start_x, float start_y, float end_x, float end_y) {
  if (std::abs(end_x - start_x) < 5.0F || std::abs(end_y - start_y) < 5.0F) {
    return;
  }
  const int current_idx = window_.get_selected_plot_index();
  if (current_idx < 0 || std::cmp_greater_equal(current_idx, available_plot_keys_.size())) {
    return;
  }
  const std::string &name = available_plot_keys_[current_idx];
  const auto *hist = dispatcher_.getHistogram(name);
  if (hist == nullptr || hist->bins.empty()) {
    return;
  }

  auto config = buildPlotConfigFromUI();
  auto x_scale = computeXScale(*hist, config);
  auto y_scale = computeYScale(*hist, config);

  auto [dx1, dy1] = correlation::plotters::detail::screenToData(
      start_x, start_y, last_plot_width_, last_plot_height_, config, x_scale, y_scale);
  auto [dx2, dy2] = correlation::plotters::detail::screenToData(
      end_x, end_y, last_plot_width_, last_plot_height_, config, x_scale, y_scale);

  zoom_x_min_ = std::min(dx1, dx2);
  zoom_x_max_ = std::max(dx1, dx2);
  zoom_y_min_ = std::min(dy1, dy2);
  zoom_y_max_ = std::max(dy1, dy2);

  window_.set_has_active_zoom(true);
  requestPlotUpdate(current_idx, true);
}

void PlotController::handlePan(float delta_x, float delta_y) {
  if (std::abs(delta_x) < 0.5F && std::abs(delta_y) < 0.5F) {
    return;
  }
  const int current_idx = window_.get_selected_plot_index();
  if (current_idx < 0 || std::cmp_greater_equal(current_idx, available_plot_keys_.size())) {
    return;
  }
  const std::string &name = available_plot_keys_[current_idx];
  const auto *hist = dispatcher_.getHistogram(name);
  if (hist == nullptr || hist->bins.empty()) {
    return;
  }

  const auto config = buildPlotConfigFromUI();
  const auto x_scale = computeXScale(*hist, config);
  const auto y_scale = computeYScale(*hist, config);

  const auto [d0_x, d0_y] = correlation::plotters::detail::screenToData(
      last_mouse_x_, last_mouse_y_, last_plot_width_, last_plot_height_, config, x_scale, y_scale);
  const auto [d1_x, d1_y] = correlation::plotters::detail::screenToData(
      last_mouse_x_ + delta_x, last_mouse_y_ + delta_y, last_plot_width_, last_plot_height_, config,
      x_scale, y_scale);

  const real_t shift_x = d1_x - d0_x;
  const real_t shift_y = d1_y - d0_y;

  const real_t cur_x_min = zoom_x_min_.value_or(x_scale.min);
  const real_t cur_x_max = zoom_x_max_.value_or(x_scale.max);
  const real_t cur_y_min = zoom_y_min_.value_or(y_scale.min);
  const real_t cur_y_max = zoom_y_max_.value_or(y_scale.max);

  zoom_x_min_ = cur_x_min - shift_x;
  zoom_x_max_ = cur_x_max - shift_x;
  zoom_y_min_ = cur_y_min - shift_y;
  zoom_y_max_ = cur_y_max - shift_y;

  window_.set_has_active_zoom(true);
  requestPlotUpdate(current_idx, true);
}

void PlotController::handleWheelZoom(float mouse_x, float mouse_y, float delta) {
  if (std::abs(delta) < 0.01F) {
    return;
  }
  const int current_idx = window_.get_selected_plot_index();
  if (current_idx < 0 || std::cmp_greater_equal(current_idx, available_plot_keys_.size())) {
    return;
  }
  const std::string &name = available_plot_keys_[current_idx];
  const auto *hist = dispatcher_.getHistogram(name);
  if (hist == nullptr || hist->bins.empty()) {
    return;
  }

  const auto config = buildPlotConfigFromUI();
  const auto x_scale = computeXScale(*hist, config);
  const auto y_scale = computeYScale(*hist, config);

  const auto [c_x, c_y] = correlation::plotters::detail::screenToData(
      mouse_x, mouse_y, last_plot_width_, last_plot_height_, config, x_scale, y_scale);

  const real_t zoom_factor = (delta > 0.0F) ? 0.85 : 1.18;

  const real_t cur_x_min = zoom_x_min_.value_or(x_scale.min);
  const real_t cur_x_max = zoom_x_max_.value_or(x_scale.max);
  const real_t cur_y_min = zoom_y_min_.value_or(y_scale.min);
  const real_t cur_y_max = zoom_y_max_.value_or(y_scale.max);

  zoom_x_min_ = c_x - (c_x - cur_x_min) * zoom_factor;
  zoom_x_max_ = c_x + (cur_x_max - c_x) * zoom_factor;
  zoom_y_min_ = c_y - (c_y - cur_y_min) * zoom_factor;
  zoom_y_max_ = c_y + (cur_y_max - c_y) * zoom_factor;

  window_.set_has_active_zoom(true);
  requestPlotUpdate(current_idx, true);
}

void PlotController::handleResetZoom() {
  zoom_x_min_.reset();
  zoom_x_max_.reset();
  zoom_y_min_.reset();
  zoom_y_max_.reset();

  window_.set_has_active_zoom(false);
  const int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    requestPlotUpdate(current_idx, true);
  }
}

void PlotController::handleAddMarkerAt(float click_x, float click_y, bool is_vertical) {
  const int current_idx = window_.get_selected_plot_index();
  if (current_idx < 0 || std::cmp_greater_equal(current_idx, available_plot_keys_.size())) {
    return;
  }
  const std::string &name = available_plot_keys_[current_idx];
  const auto *hist = dispatcher_.getHistogram(name);
  if (hist == nullptr || hist->bins.empty()) {
    return;
  }

  auto config = buildPlotConfigFromUI();
  auto x_scale = computeXScale(*hist, config);
  auto y_scale = computeYScale(*hist, config);

  auto [d_x, d_y] = correlation::plotters::detail::screenToData(
      click_x, click_y, last_plot_width_, last_plot_height_, config, x_scale, y_scale);

  real_t val = is_vertical ? d_x : d_y;
  std::string color = is_vertical ? "#0284C7" : "#E11D48";
  std::string label = std::format("{}: {:.3f}", is_vertical ? "X" : "Y", val);

  reference_markers_.push_back(correlation::plotters::ReferenceLine{
      .value = val,
      .is_vertical = is_vertical,
      .label = std::move(label),
      .color_hex = std::move(color),
  });

  window_.set_markers_count(static_cast<int>(reference_markers_.size()));
  requestPlotUpdate(current_idx, true);
}

void PlotController::handleAddMarker(bool is_vertical) {
  handleAddMarkerAt(last_mouse_x_, last_mouse_y_, is_vertical);
}

void PlotController::handleClearMarkers() {
  reference_markers_.clear();
  window_.set_markers_count(0);
  const int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    requestPlotUpdate(current_idx, true);
  }
}

void PlotController::resetPlotState() noexcept {
  series_manager_.reset();
  zoom_x_min_.reset();
  zoom_x_max_.reset();
  zoom_y_min_.reset();
  zoom_y_max_.reset();
  reference_markers_.clear();
  slint::invoke_from_event_loop([this]() {
    window_.set_pinned_runs_count(0);
    window_.set_has_active_zoom(false);
    window_.set_markers_count(0);
    window_.set_hover_coord_text("");
    const int current_idx = window_.get_selected_plot_index();
    if (current_idx >= 0) {
      requestPlotUpdate(current_idx, true);
    }
  });
}

} // namespace correlation::app
