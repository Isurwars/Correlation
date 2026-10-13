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

#include <algorithm>
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
      {"COP", 46},
      {"LEF", 47},
      {"sigma2_N", 48},
      {"chi_H", 49},

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

  task.comparison_hists.push_back({.label = "Current", .hist = nullptr});
  for (const auto &pinned_run : series_manager_.getPinnedRuns()) {
    const auto hist_it = pinned_run.histograms.find(name);
    if (hist_it != pinned_run.histograms.end()) {
      const bool vis = series_manager_.isCurveVisible(pinned_run.label, true);
      task.comparison_hists.push_back(
          {.label = pinned_run.label, .hist = &hist_it->second, .style = {.visible = vis}});
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
                                  dispatcher_.getAshcroftWeights(), [this](std::string status) {
                                    window_.set_analysis_status_text(
                                        slint::SharedString(std::move(status)));
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

  window_.set_has_active_zoom(viewport_navigator_.hasActiveZoom());
  requestPlotUpdate(current_idx, true);
}

void PlotController::handleResetZoom() {
  viewport_navigator_.resetZoom();
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

} // namespace correlation::app
