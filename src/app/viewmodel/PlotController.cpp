/**
 * @file PlotController.cpp
 * @brief Implementation of PlotController.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/viewmodel/PlotController.hpp"
#include "AppWindow.h"
#include "app/formatters/PlotTableFormatter.hpp"
#include "app/viewmodel/Structure3DController.hpp"

#include <algorithm>
#include <format>
#include <map>
#include <span>
#include <utility>

namespace correlation::app {

namespace {
template <typename T> T safeParse(const slint::SharedString &str, T default_value) {
  const std::string_view raw_str_view(
      str.data(), str.data() != nullptr ? std::char_traits<char>::length(str.data()) : 0);
  if (raw_str_view.empty()) {
    return default_value;
  }
  try {
    if constexpr (std::is_floating_point_v<T>) {
      return static_cast<T>(std::stod(std::string(raw_str_view)));
    } else if constexpr (std::is_integral_v<T>) {
      return static_cast<T>(std::stoi(std::string(raw_str_view)));
    } else {
      return default_value;
    }
  } catch (...) {
    return default_value;
  }
}

correlation::plotters::detail::NiceScale
computeXScale(const correlation::analysis::Histogram &hist,
              const correlation::plotters::PlotConfig &config) {
  return PlotViewportNavigator::computeXScale(hist, config);
}

correlation::plotters::detail::NiceScale
computeYScale(const correlation::analysis::Histogram &hist,
              const correlation::plotters::PlotConfig &config) {
  return PlotViewportNavigator::computeYScale(hist, config);
}

void sortPlotNames(std::vector<std::string> &names) {
  static const std::map<std::string, int> PRIORITY = {
      // Real-space pair & radial distributions
      {"g_r", 0},
      {"g_r_unweighted", 1},
      {"H_r", 2},
      {"G_r", 3},
      {"J_r", 4},
      {"SDF", 5},

      // Reciprocal-space & scattering
      {"S_q", 10},
      {"XRD", 11},

      // Angular & torsional distributions
      {"PAD", 20},
      {"PAD_raw", 21},
      {"BAD", 22},
      {"DAD", 23},
      {"DAD_raw", 24},

      // Coordination, topology & networks
      {"CN", 30},
      {"CNA", 31},
      {"voronoi_cn", 32},
      {"Voronoi Coordination Number", 32},
      {"voronoi_volume", 33},
      {"Voronoi Volume", 33},
      {"voronoi_sphericity", 34},
      {"Voronoi Sphericity", 34},
      {"voronoi_signatures", 35},
      {"Voronoi Signatures", 35},
      {"RD", 36},
      {"cluster_size", 37},
      {"Cluster Size", 37},
      {"HBond", 38},

      // Order parameters & local structure
      {"Q4", 40},
      {"Q6", 41},
      {"W4_hat", 42},
      {"W6_hat", 43},
      {"Q4_bar", 44},
      {"Q6_bar", 45},
      {"G6_r", 46},
      {"TWF_Nconn", 47},
      {"Q_tet", 48},
      {"COP", 49},
      {"LEF", 50},
      {"sigma2_N", 51},
      {"chi_H", 52},

      // Dynamics & temporal correlations
      {"MSD", 50},
      {"D_eff", 51},
      {"VACF", 52},
      {"vacf_normalized", 53},
      {"Normalized VACF", 53},
      {"VDOS", 54},

      // Electronic structure & MLIP projections
      {"TDOS", 60},
      {"tdos", 61},
      {"MotifProjectedTDOS_CNA", 62},
      {"MotifProjectedTDOS_Steinhardt", 63},
  };

  std::ranges::sort(names, [](const std::string &lhs, const std::string &rhs) {
    const auto it_a = PRIORITY.find(lhs);
    const auto it_b = PRIORITY.find(rhs);
    const int prio_a = (it_a != PRIORITY.end()) ? it_a->second : 100;
    const int prio_b = (it_b != PRIORITY.end()) ? it_b->second : 100;
    if (prio_a != prio_b) {
      return prio_a < prio_b;
    }
    return lhs < rhs;
  });
}

std::shared_ptr<slint::VectorModel<MenuItem>>
buildPlotMenuModel(std::span<const std::string> names, const AnalysisDispatcher &dispatcher) {
  auto menu_model = std::make_shared<slint::VectorModel<MenuItem>>();
  for (const auto &name : names) {
    MenuItem item;
    const correlation::analysis::Histogram *hist = dispatcher.getHistogram(name);
    const bool has_custom_title = (hist != nullptr) && (!hist->title.empty());
    const std::string display_text = has_custom_title ? hist->title : name;
    item.text = slint::SharedString(display_text);
    item.enabled = true;
    menu_model->push_back(item);
  }
  return menu_model;
}

void updateDynamicProperties(
    const ::AppWindow &window,
    const correlation::analysis::DistributionFunctions *distribution_functions) {
  if (distribution_functions == nullptr) {
    window.set_diff_msd("");
    window.set_diff_vacf("");
    window.set_relaxation_time("");
    window.set_deborah_number("");
    return;
  }

  const real_t msd = distribution_functions->getDiffusionCoefficientMSD();
  window.set_diff_msd(msd > 0.0 ? slint::SharedString(std::format("{:.6f} Å²/fs", msd)) : "");

  const real_t vacf = distribution_functions->getDiffusionCoefficientVACF();
  window.set_diff_vacf(vacf > 0.0 ? slint::SharedString(std::format("{:.6f} Å²/fs", vacf)) : "");

  const real_t tau = distribution_functions->getRelaxationTime();
  window.set_relaxation_time(tau > 0.0 ? slint::SharedString(std::format("{:.4f} fs", tau)) : "");

  const real_t deb = distribution_functions->getDeborahNumber();
  window.set_deborah_number(deb > 0.0 ? slint::SharedString(std::format("{:.4f}", deb)) : "");
}

} // namespace

PlotController::PlotController(::AppWindow &window, AnalysisDispatcher &dispatcher,
                               const ProgramOptions &options)
    : window_(window), dispatcher_(dispatcher), options_(options),
      pipeline_([this](const slint::Image &img) { window_.set_preview_plot(img); }) {}

PlotController::~PlotController() = default;

void PlotController::handlePlotResized(PlotSize size) {
  if (std::abs(size.width - last_plot_width_) < 4.0F &&
      std::abs(size.height - last_plot_height_) < 4.0F) {
    return;
  }
  last_plot_width_ = size.width;
  last_plot_height_ = size.height;
  if (on_plot_geometry_changed_) {
    on_plot_geometry_changed_(last_plot_width_, last_plot_height_);
  }
  const int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    requestPlotUpdate(current_idx, true);
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
    config.width = 1200.0;
    config.height = 900.0;
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

  config.manual_ticks = window_.get_export_config().manual_ticks;
  config.x_major_step = safeParse(window_.get_export_config().x_major_step, 0.0F);
  config.y_major_step = safeParse(window_.get_export_config().y_major_step, 0.0F);
  config.minor_subdivisions = window_.get_export_config().minor_subdivisions;

  const int target_idx =
      (active_plot_index_ >= 0) ? active_plot_index_ : window_.get_selected_plot_index();
  if (target_idx >= 0 && metric_settings_.contains(target_idx)) {
    const auto &settings = metric_settings_.at(target_idx);
    config.show_grid = settings.show_major_grid;
    config.show_minor_grid = settings.show_minor_grid;
    config.y_log_scale = settings.y_log_scale;
    if (!settings.auto_x) {
      config.manual_x_min = settings.manual_x_min;
      config.manual_x_max = settings.manual_x_max;
    } else {
      config.manual_x_min = std::nullopt;
      config.manual_x_max = std::nullopt;
    }
    if (!settings.auto_y) {
      config.manual_y_min = settings.manual_y_min;
      config.manual_y_max = settings.manual_y_max;
    } else {
      config.manual_y_min = std::nullopt;
      config.manual_y_max = std::nullopt;
    }
    if (!settings.auto_ticks) {
      config.manual_ticks = true;
      config.x_major_step = settings.x_major_step;
      config.y_major_step = settings.y_major_step;
      config.minor_subdivisions = settings.minor_subdivisions;
    } else {
      config.manual_ticks = false;
      config.x_major_step = static_cast<real_t>(0.0);
      config.y_major_step = static_cast<real_t>(0.0);
      config.minor_subdivisions = 0;
    }
    if (settings.line_width > static_cast<real_t>(0.0)) {
      config.line_width = settings.line_width;
    }
  }

  viewport_navigator_.applyToConfig(config);
  marker_manager_.applyToConfig(config);

  return config;
}

void PlotController::populatePlotList() {
  pipeline_.invalidateCache();
  auto names = dispatcher_.getAvailableHistogramNames();

  sortPlotNames(names);
  available_plot_keys_ = names;

  window_.set_plot_items(buildPlotMenuModel(names, dispatcher_));
  updateDynamicProperties(window_, dispatcher_.getDistributionFunctions());

  window_.set_selected_plot_index(names.empty() ? -1 : 0);
  if (!names.empty()) {
    requestPlotUpdate(0, true);
  }
}

void PlotController::handleMouseMove(float mouse_x, float mouse_y, bool hover, float width,
                                     float height) {
  const bool actual_hover = hover;

  if (std::abs(mouse_x - last_mouse_x_) < 0.5F && std::abs(mouse_y - last_mouse_y_) < 0.5F &&
      actual_hover == mouse_hover_ && std::abs(width - last_plot_width_) < 1e-2F &&
      std::abs(height - last_plot_height_) < 1e-2F) {
    return;
  }

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
        const auto config = buildPlotConfigFromUI();
        const auto geom = correlation::plotters::detail::getViewportGeometry(config);
        auto [svg_x, svg_y] =
            correlation::plotters::detail::screenToSvg(mouse_x, mouse_y, width, height, config);
        if (svg_x >= geom.px0 && svg_x <= geom.px1 && svg_y >= geom.py0 && svg_y <= geom.py1) {
          const auto x_scale = computeXScale(*hist, config);
          const auto y_scale = computeYScale(*hist, config);
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
}

void PlotController::requestPlotUpdate(int index, bool immediate) {
  if (index < 0 || std::cmp_greater_equal(index, available_plot_keys_.size())) {
    return;
  }

  if (active_plot_index_ != index) {
    active_plot_index_ = index;
    syncMetricSettingsToUI(index);
  }

  const std::string &name = available_plot_keys_[index];
  const correlation::analysis::Histogram *hist = dispatcher_.getHistogram(name);
  if (hist == nullptr) {
    return;
  }

  const correlation::plotters::PlotConfig config = buildPlotConfigFromUI();
  correlation::plotters::HoverInfo hover;
  hover.active =
      false; // Retired in interactive mode: Slint renders hardware-accelerated crosshairs
  hover.mouse_x = last_mouse_x_;
  hover.mouse_y = last_mouse_y_;
  hover.widget_width = last_plot_width_;
  hover.widget_height = last_plot_height_;

  const std::size_t pinned_runs_count = series_manager_.getPinnedRunsCount();
  if (!immediate && pipeline_.isCacheHit(index, config, hover, pinned_runs_count)) {
    return;
  }

  updateTableData(hist);
  updateCurveToggleItems(hist);

  RenderTask task;
  task.index = index;
  task.active_hist = *hist;
  task.config = config;
  task.config.show_difference_curve = series_manager_.shouldShowDifference();
  task.hover = hover;
  task.ashcroft_weights = dispatcher_.getAshcroftWeights();
  task.curve_visibility = series_manager_.getCurveVisibilityMap();
  task.custom_curve_colors = series_manager_.getAssignedCurveColors();
  task.pinned_runs_count = pinned_runs_count;

  correlation::plotters::LabeledHistogram current_entry;
  current_entry.label = "Current";
  task.comparison_hists.push_back(std::move(current_entry));
  for (const auto &pinned_run : series_manager_.getPinnedRuns()) {
    const auto hist_it = pinned_run.histograms.find(name);
    if (hist_it != pinned_run.histograms.end()) {
      const bool vis = series_manager_.isCurveVisible(pinned_run.label, true);
      correlation::plotters::LabeledHistogram pinned_entry;
      pinned_entry.label = pinned_run.label;
      pinned_entry.hist = &hist_it->second;
      pinned_entry.style.visible = vis;
      task.comparison_hists.push_back(std::move(pinned_entry));
    }
  }

  if (window_.get_selected_plot_index() != index) {
    window_.set_selected_plot_index(index);
  }

  pipeline_.submitTask(std::move(task));
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
  window_.set_curve_toggle_items(
      series_manager_.generateToggleItems(hist, active_config, dispatcher_.getAshcroftWeights()));
}

void PlotController::handleSavePlot() {
  const int index = window_.get_selected_plot_index();
  if (index < 0 || std::cmp_greater_equal(index, available_plot_keys_.size())) {
    return;
  }
  const std::string name = available_plot_keys_[index];
  const correlation::analysis::Histogram *hist = dispatcher_.getHistogram(name);
  if (hist == nullptr) {
    return;
  }

  const auto config = buildPlotConfigFromUI();
  const auto pinned_runs = series_manager_.getPinnedRuns();
  const std::vector<PinnedRun> pinned_vec(pinned_runs.begin(), pinned_runs.end());

  export_manager_.exportPlotAsync(options_.output_file_base, name, hist, config, pinned_vec,
                                  dispatcher_.getAshcroftWeights(),
                                  [this](const std::string &status) {
                                    window_.set_analysis_status_text(slint::SharedString(status));
                                  });
}

void PlotController::handlePinRun() {
  const auto &hists = dispatcher_.getHistograms();
  if (hists.empty()) {
    return;
  }

  series_manager_.pinCurrentRun(hists);
  window_.set_pinned_runs_count(static_cast<int>(series_manager_.getPinnedRunsCount()));

  const int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    requestPlotUpdate(current_idx, true);
  }
}

void PlotController::handleClearPinnedRuns() {
  series_manager_.clearPinnedRuns();
  window_.set_pinned_runs_count(0);

  const int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    requestPlotUpdate(current_idx, true);
  }
}

void PlotController::handleToggleCurveVisibility(int curve_id, bool visible) {
  series_manager_.setCurveVisible(curve_id, visible);
  const int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    requestPlotUpdate(current_idx, true);
  }
}

void PlotController::handleToggleAllCurves(bool visible) {
  series_manager_.setAllCurvesVisible(visible);
  const int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    requestPlotUpdate(current_idx, true);
  }
}

void PlotController::handleToggleDifferencePlot(bool show_difference) {
  series_manager_.setShowDifference(show_difference);
  const int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    requestPlotUpdate(current_idx, true);
  }
}

void PlotController::handleZoomRect(float start_x, float start_y, float end_x, float end_y) {
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
  viewport_navigator_.applyZoomRect(start_x, start_y, end_x, end_y, last_plot_width_,
                                    last_plot_height_, *hist, config);

  syncZoomLimitsToSettings(current_idx);
  window_.set_has_active_zoom(viewport_navigator_.hasActiveZoom());
  requestPlotUpdate(current_idx, true);
}

void PlotController::handlePan(float delta_x, float delta_y) {
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
  viewport_navigator_.applyPan(delta_x, delta_y, last_mouse_x_, last_mouse_y_, last_plot_width_,
                               last_plot_height_, *hist, config);

  syncZoomLimitsToSettings(current_idx);
  window_.set_has_active_zoom(viewport_navigator_.hasActiveZoom());
  requestPlotUpdate(current_idx, true);
}

void PlotController::handleWheelZoom(float mouse_x, float mouse_y, float delta) {
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
  viewport_navigator_.applyWheelZoom(mouse_x, mouse_y, delta, last_plot_width_, last_plot_height_,
                                     *hist, config);

  syncZoomLimitsToSettings(current_idx);
  window_.set_has_active_zoom(viewport_navigator_.hasActiveZoom());
  requestPlotUpdate(current_idx, true);
}

void PlotController::handleResetZoom() {
  viewport_navigator_.resetZoom();
  window_.set_has_active_zoom(false);
  const int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    auto &settings = metric_settings_[current_idx];
    settings.auto_x = true;
    settings.manual_x_min = std::nullopt;
    settings.manual_x_max = std::nullopt;
    settings.auto_y = true;
    settings.manual_y_min = std::nullopt;
    settings.manual_y_max = std::nullopt;
    syncMetricSettingsToUI(current_idx);
    syncSettingsTo3DHUD(current_idx, settings);
    pipeline_.invalidateCache();
    requestPlotUpdate(current_idx, true);
  }
}

void PlotController::syncZoomLimitsToSettings(int current_idx) {
  if (current_idx < 0) {
    return;
  }
  const auto limits = viewport_navigator_.limits();
  auto &settings = metric_settings_[current_idx];
  if (limits.x_min.has_value() || limits.x_max.has_value()) {
    settings.auto_x = false;
    settings.manual_x_min = limits.x_min;
    settings.manual_x_max = limits.x_max;
  }
  if (limits.y_min.has_value() || limits.y_max.has_value()) {
    settings.auto_y = false;
    settings.manual_y_min = limits.y_min;
    settings.manual_y_max = limits.y_max;
  }
  syncMetricSettingsToUI(current_idx);
  syncSettingsTo3DHUD(current_idx, settings);
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

  const auto config = buildPlotConfigFromUI();
  const auto [d_x, d_y] = PlotViewportNavigator::screenToData(click_x, click_y, last_plot_width_,
                                                              last_plot_height_, *hist, config);

  const real_t val = is_vertical ? d_x : d_y;
  marker_manager_.addMarker(val, is_vertical);

  window_.set_markers_count(static_cast<int>(marker_manager_.markersCount()));
  requestPlotUpdate(current_idx, true);
}

void PlotController::handleAddMarker(bool is_vertical) {
  handleAddMarkerAt(last_mouse_x_, last_mouse_y_, is_vertical);
}

void PlotController::handleClearMarkers() {
  marker_manager_.clearMarkers();
  window_.set_markers_count(0);
  const int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    requestPlotUpdate(current_idx, true);
  }
}

void PlotController::resetPlotState() noexcept {
  series_manager_.reset();
  viewport_navigator_.resetZoom();
  marker_manager_.clearMarkers();
  pipeline_.invalidateCache();
  slint::invoke_from_event_loop([this] {
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

void PlotController::syncMetricSettingsToUI(int index) {
  if (index < 0) {
    return;
  }
  auto &settings = metric_settings_[index];
  ::PlotAxisSettings ui_settings;
  ui_settings.auto_x = settings.auto_x;
  ui_settings.x_min = settings.manual_x_min.has_value()
                          ? slint::SharedString(std::format("{:.4g}", *settings.manual_x_min))
                          : slint::SharedString("");
  ui_settings.x_max = settings.manual_x_max.has_value()
                          ? slint::SharedString(std::format("{:.4g}", *settings.manual_x_max))
                          : slint::SharedString("");
  ui_settings.auto_y = settings.auto_y;
  ui_settings.y_min = settings.manual_y_min.has_value()
                          ? slint::SharedString(std::format("{:.4g}", *settings.manual_y_min))
                          : slint::SharedString("");
  ui_settings.y_max = settings.manual_y_max.has_value()
                          ? slint::SharedString(std::format("{:.4g}", *settings.manual_y_max))
                          : slint::SharedString("");
  ui_settings.auto_ticks = settings.auto_ticks;
  ui_settings.x_major_step = settings.x_major_step > static_cast<real_t>(0.0)
                                 ? slint::SharedString(std::format("{:.4g}", settings.x_major_step))
                                 : slint::SharedString("");
  ui_settings.y_major_step = settings.y_major_step > static_cast<real_t>(0.0)
                                 ? slint::SharedString(std::format("{:.4g}", settings.y_major_step))
                                 : slint::SharedString("");
  ui_settings.minor_subdivisions = settings.minor_subdivisions;
  ui_settings.y_log_scale = settings.y_log_scale;
  ui_settings.show_major_grid = settings.show_major_grid;
  ui_settings.show_minor_grid = settings.show_minor_grid;
  ui_settings.line_width = slint::SharedString(std::format("{:.1f}", settings.line_width));

  window_.set_axis_settings(ui_settings);
}

void PlotController::syncSettingsTo3DHUD(
    int index, const correlation::plotters::MetricPlotSettings &settings) {
  if (structure_3d_controller_ == nullptr || index < 0 ||
      std::cmp_greater_equal(index, available_plot_keys_.size())) {
    return;
  }
  const std::string &name = available_plot_keys_[index];
  std::optional<analysis::DistributionMetric> metric;
  if (name == "g_r" || name.contains("g_r") || name.contains("Pair Distribution") ||
      name.contains("Radial Distribution") || name.contains("g(r)") || name.contains("RDF")) {
    metric = analysis::DistributionMetric::RDF;
  } else if (name == "PAD" || name.contains("PAD") || name.contains("Planar Angle") ||
             name.contains("theta")) {
    metric = analysis::DistributionMetric::PAD;
  } else if (name == "S_q" || name.contains("S_q") || name.contains("Structure Factor") ||
             name.contains("S(Q)") || name.contains("SQ")) {
    metric = analysis::DistributionMetric::SQ;
  }

  if (metric.has_value()) {
    const render3d::MiniPlotConfig mini_cfg{
        .manual_x_min = settings.auto_x ? std::nullopt : settings.manual_x_min,
        .manual_x_max = settings.auto_x ? std::nullopt : settings.manual_x_max,
        .manual_y_min = settings.auto_y ? std::nullopt : settings.manual_y_min,
        .manual_y_max = settings.auto_y ? std::nullopt : settings.manual_y_max,
        .y_log_scale = settings.y_log_scale,
    };
    structure_3d_controller_->setMetricPlotLimits(*metric, mini_cfg);
  }
}

namespace {

void parseXAxisSettings(const ::PlotAxisSettings &ui_settings,
                        correlation::plotters::MetricPlotSettings &settings) {
  const std::string_view x_min_sv(ui_settings.x_min.data(),
                                  ui_settings.x_min.data() != nullptr
                                      ? std::char_traits<char>::length(ui_settings.x_min.data())
                                      : 0);
  const std::string_view x_max_sv(ui_settings.x_max.data(),
                                  ui_settings.x_max.data() != nullptr
                                      ? std::char_traits<char>::length(ui_settings.x_max.data())
                                      : 0);
  const bool has_manual_x = !x_min_sv.empty() || !x_max_sv.empty();
  settings.auto_x = has_manual_x ? false : ui_settings.auto_x;

  if (!settings.auto_x) {
    settings.manual_x_min =
        !x_min_sv.empty()
            ? std::optional<real_t>(safeParse<real_t>(ui_settings.x_min, static_cast<real_t>(0.0)))
            : std::nullopt;
    settings.manual_x_max =
        !x_max_sv.empty()
            ? std::optional<real_t>(safeParse<real_t>(ui_settings.x_max, static_cast<real_t>(10.0)))
            : std::nullopt;
    if (settings.manual_x_min.has_value() && settings.manual_x_max.has_value() &&
        *settings.manual_x_min >= *settings.manual_x_max) {
      std::swap(settings.manual_x_min, settings.manual_x_max);
    }
  } else {
    settings.manual_x_min = std::nullopt;
    settings.manual_x_max = std::nullopt;
  }
}

void parseYAxisSettings(const ::PlotAxisSettings &ui_settings,
                        correlation::plotters::MetricPlotSettings &settings) {
  const std::string_view y_min_sv(ui_settings.y_min.data(),
                                  ui_settings.y_min.data() != nullptr
                                      ? std::char_traits<char>::length(ui_settings.y_min.data())
                                      : 0);
  const std::string_view y_max_sv(ui_settings.y_max.data(),
                                  ui_settings.y_max.data() != nullptr
                                      ? std::char_traits<char>::length(ui_settings.y_max.data())
                                      : 0);
  const bool has_manual_y = !y_min_sv.empty() || !y_max_sv.empty();
  settings.auto_y = has_manual_y ? false : ui_settings.auto_y;

  if (!settings.auto_y) {
    settings.manual_y_min =
        !y_min_sv.empty()
            ? std::optional<real_t>(safeParse<real_t>(ui_settings.y_min, static_cast<real_t>(0.0)))
            : std::nullopt;
    settings.manual_y_max =
        !y_max_sv.empty()
            ? std::optional<real_t>(safeParse<real_t>(ui_settings.y_max, static_cast<real_t>(1.0)))
            : std::nullopt;
    if (settings.manual_y_min.has_value() && settings.manual_y_max.has_value() &&
        *settings.manual_y_min >= *settings.manual_y_max) {
      std::swap(settings.manual_y_min, settings.manual_y_max);
    }
  } else {
    settings.manual_y_min = std::nullopt;
    settings.manual_y_max = std::nullopt;
  }
}

void parseTickAndGridSettings(const ::PlotAxisSettings &ui_settings,
                              correlation::plotters::MetricPlotSettings &settings) {
  const std::string_view x_step_sv(ui_settings.x_major_step.data(),
                                   ui_settings.x_major_step.data() != nullptr
                                       ? std::char_traits<char>::length(
                                             ui_settings.x_major_step.data())
                                       : 0);
  const std::string_view y_step_sv(ui_settings.y_major_step.data(),
                                   ui_settings.y_major_step.data() != nullptr
                                       ? std::char_traits<char>::length(
                                             ui_settings.y_major_step.data())
                                       : 0);
  const bool has_manual_steps = !x_step_sv.empty() || !y_step_sv.empty();
  settings.auto_ticks = has_manual_steps ? false : ui_settings.auto_ticks;

  if (!settings.auto_ticks) {
    settings.x_major_step = safeParse<real_t>(ui_settings.x_major_step, static_cast<real_t>(0.0));
    settings.y_major_step = safeParse<real_t>(ui_settings.y_major_step, static_cast<real_t>(0.0));
  } else {
    settings.x_major_step = static_cast<real_t>(0.0);
    settings.y_major_step = static_cast<real_t>(0.0);
  }

  settings.minor_subdivisions = std::clamp(ui_settings.minor_subdivisions, 0, 5);
  settings.y_log_scale = ui_settings.y_log_scale;
  settings.show_major_grid = ui_settings.show_major_grid;
  settings.show_minor_grid = ui_settings.show_minor_grid;
}

} // namespace

void PlotController::handleApplyAxisSettings() {
  if (active_plot_index_ < 0) {
    active_plot_index_ = window_.get_selected_plot_index();
  }
  if (active_plot_index_ < 0) {
    return;
  }
  const auto ui_settings = window_.get_axis_settings();
  auto &settings = metric_settings_[active_plot_index_];

  parseXAxisSettings(ui_settings, settings);
  parseYAxisSettings(ui_settings, settings);
  parseTickAndGridSettings(ui_settings, settings);

  viewport_navigator_.setLimits({
      .x_min = settings.manual_x_min,
      .x_max = settings.manual_x_max,
      .y_min = settings.manual_y_min,
      .y_max = settings.manual_y_max,
  });

  syncMetricSettingsToUI(active_plot_index_);
  syncSettingsTo3DHUD(active_plot_index_, settings);
  pipeline_.invalidateCache();
  requestPlotUpdate(active_plot_index_, true);
}

void PlotController::handleResetAxisSettings() {
  if (active_plot_index_ < 0) {
    return;
  }
  metric_settings_[active_plot_index_] = correlation::plotters::MetricPlotSettings{};
  syncMetricSettingsToUI(active_plot_index_);
  viewport_navigator_.resetZoom();
  syncSettingsTo3DHUD(active_plot_index_, metric_settings_[active_plot_index_]);
  pipeline_.invalidateCache();
  requestPlotUpdate(active_plot_index_, true);
}

void PlotController::handleSetCurveWidth(const slint::SharedString &width_str) {
  const auto curve_width = safeParse<real_t>(width_str, static_cast<real_t>(3.0));
  if (curve_width <= static_cast<real_t>(0.0)) {
    return;
  }
  if (active_plot_index_ >= 0) {
    metric_settings_[active_plot_index_].line_width = curve_width;
  }
  pipeline_.invalidateCache();
  if (active_plot_index_ >= 0) {
    requestPlotUpdate(active_plot_index_, true);
  }
}

} // namespace correlation::app
