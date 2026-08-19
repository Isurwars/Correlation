/**
 * @file PlotController.cpp
 * @brief Implementation of PlotController.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/PlotController.hpp"
#include "AppWindow.h"
#include "plotters/PdfPlotter.hpp"
#include <nfd.h>

#include <algorithm>
#include <filesystem>
#include <format>
#include <fstream>
#include <map>
#include <span>

namespace correlation::app {

namespace {
std::string getComparisonKey(const correlation::analysis::Histogram *hist) {
  std::string key = "Total";
  const auto &partials = hist->smoothed_partials.empty() ? hist->partials : hist->smoothed_partials;
  if (!partials.empty() && !partials.contains(key)) {
    key = partials.begin()->first;
  }
  return key;
}

template <typename T> T safe_parse(const slint::SharedString &str, T default_value) {
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

slint::Color parseHexColor(const std::string &hex) {
  if (hex.size() == 7 && hex[0] == '#') {
    auto red = static_cast<uint8_t>(std::stoul(hex.substr(1, 2), nullptr, 16));
    auto green = static_cast<uint8_t>(std::stoul(hex.substr(3, 2), nullptr, 16));
    auto blue = static_cast<uint8_t>(std::stoul(hex.substr(5, 2), nullptr, 16));
    return slint::Color::from_rgb_uint8(red, green, blue);
  }
  return slint::Color::from_rgb_uint8(128, 128, 128);
}
} // namespace

PlotController::PlotController(::AppWindow &window, AppBackend &backend) : window_(window), backend_(backend) {
  update_timer_.start(slint::TimerMode::Repeated, std::chrono::milliseconds(1000), [this]() { handleUpdateTimer(); });
}

PlotController::~PlotController() {
  if (render_thread_.joinable()) {
    render_thread_.join();
  }
}

void PlotController::handlePlotResized(PlotSize size) {
  if (std::abs(size.width - last_plot_width_) < 1.0F && std::abs(size.height - last_plot_height_) < 1.0F) {
    return;
  }
  last_plot_width_ = size.width;
  last_plot_height_ = size.height;
  int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    slint::invoke_from_event_loop([this, current_idx] { requestPlotUpdate(current_idx, true); });
  }
}

void PlotController::handleUpdateTimer() {
  int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    requestPlotUpdate(current_idx, false);
  }
}

correlation::plotters::PlotConfig PlotController::buildPlotConfigFromUI() {
  correlation::plotters::PlotConfig config;
  config.theme = window_.get_is_dark() ? correlation::plotters::PlotConfig::Theme::Dark
                                       : correlation::plotters::PlotConfig::Theme::Light;

  int size_preset_val = window_.get_export_config().size_preset;
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

  int palette_val = window_.get_export_config().palette;
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

  config.font_scale = safe_parse(window_.get_export_config().font_scale, 1.0F);
  if (config.font_scale <= 0.0F) {
    config.font_scale = 1.0F;
  }

  config.line_width = safe_parse(window_.get_export_config().line_width, 3.0F);
  if (config.line_width <= 0.0F) {
    config.line_width = 3.0F;
  }

  config.show_grid = window_.get_export_config().show_grid;
  config.show_legend = window_.get_export_config().show_legend;
  config.show_markers = window_.get_export_config().show_markers;
  config.fill_area = window_.get_export_config().fill_area;

  return config;
}

void PlotController::populatePlotList() {
  last_rendered_index_ = -1;
  auto names = backend_.getAvailableHistogramNames();

  std::map<std::string, int> priority = {{"g_r", 0},  {"G_r", 1},   {"J_r", 2},   {"S_q", 10},      {"XRD", 11},
                                         {"BAD", 20}, {"PAD", 21},  {"DAD", 22},  {"CN", 23},       {"RD", 24},
                                         {"MSD", 30}, {"VACF", 31}, {"VDOS", 32}, {"sigma2_N", 40}, {"chi_H", 41}};

  std::ranges::sort(names, [&](const std::string &lhs, const std::string &rhs) {
    int prio_a = priority.contains(lhs) ? priority.at(lhs) : 100;
    int prio_b = priority.contains(rhs) ? priority.at(rhs) : 100;
    if (prio_a != prio_b) {
      return prio_a < prio_b;
    }
    return lhs < rhs;
  });

  available_plot_keys_ = names;

  auto menu_model = std::make_shared<slint::VectorModel<MenuItem>>();
  for (const auto &name : names) {
    MenuItem item;
    const correlation::analysis::Histogram *hist = backend_.getHistogram(name);
    std::string display_text = (hist != nullptr && !hist->title.empty()) ? hist->title : name;
    item.text = slint::SharedString(display_text);
    item.enabled = true;
    menu_model->push_back(item);
  }
  window_.set_plot_items(menu_model);

  // Update dynamic properties
  const auto *distribution_functions = backend_.getDistributionFunctions();
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

void PlotController::handleMouseMove(float mouse_x, float mouse_y, bool hover, float width, float height) {
  bool actual_hover = hover;

  if (std::abs(mouse_x - last_mouse_x_) < 0.5F && std::abs(mouse_y - last_mouse_y_) < 0.5F &&
      actual_hover == mouse_hover_ && std::abs(width - last_plot_width_) < 1e-2F &&
      std::abs(height - last_plot_height_) < 1e-2F) {
    return;
  }

  bool hover_changed = (actual_hover != mouse_hover_);

  last_mouse_x_ = mouse_x;
  last_mouse_y_ = mouse_y;
  mouse_hover_ = actual_hover;
  last_plot_width_ = width;
  last_plot_height_ = height;

  int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    requestPlotUpdate(current_idx, hover_changed || !actual_hover);
  }
}

void PlotController::requestPlotUpdate(int index, bool immediate) {
  if (index < 0 || index >= static_cast<int>(available_plot_keys_.size())) {
    return;
  }

  const std::string &name = available_plot_keys_[index];
  const correlation::analysis::Histogram *hist = backend_.getHistogram(name);
  if (hist == nullptr) {
    return;
  }

  updateTableData(hist);
  updateCurveToggleItems(hist);

  if (immediate) {
    needs_redraw_ = true;
  }

  correlation::plotters::PlotConfig config = buildPlotConfigFromUI();
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
  last_pinned_runs_count_ = pinned_runs_.size();
  last_config_ = config;
  last_hover_ = hover;

  RenderTaskData data;
  data.active_hist = *hist;
  data.config = config;
  data.hover = hover;
  data.ashcroft_weights = backend_.getAshcroftWeights();
  data.curve_visibility = curve_visibility_map_;

  data.comparison_hists.push_back({.label = "Current", .hist = &data.active_hist});
  for (const auto &pinned_run : pinned_runs_) {
    auto hist_it = pinned_run.histograms.find(name);
    if (hist_it != pinned_run.histograms.end()) {
      bool vis = curve_visibility_map_.contains(pinned_run.label) ? curve_visibility_map_[pinned_run.label] : true;
      data.comparison_hists.push_back({.label = pinned_run.label, .hist = &hist_it->second, .style = {.visible = vis}});
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
  // Update table data
  const auto &partials = hist->smoothed_partials.empty() ? hist->partials : hist->smoothed_partials;
  std::vector<std::string> keys;
  for (const auto &pair : partials) {
    if (pair.first != "Total") {
      keys.push_back(pair.first);
    }
  }
  std::sort(keys.begin(), keys.end());
  if (partials.contains("Total")) {
    keys.insert(keys.begin(), "Total");
  }

  auto slint_headers = std::make_shared<slint::VectorModel<slint::SharedString>>();
  std::string x_header = hist->x_label;
  if (!hist->x_unit.empty()) {
    x_header += " (" + hist->x_unit + ")";
  }
  slint_headers->push_back(slint::SharedString(x_header));
  for (const auto &key : keys) {
    slint_headers->push_back(slint::SharedString(key));
  }
  window_.set_table_headers(slint_headers);

  auto slint_rows = std::make_shared<slint::VectorModel<TableRow>>();
  size_t num_bins = hist->bins.size();
  for (size_t i = 0; i < num_bins; ++i) {
    TableRow row;
    auto row_values = std::make_shared<slint::VectorModel<slint::SharedString>>();

    // First column: bin value
    row_values->push_back(slint::SharedString(std::format("{:.4f}", hist->bins[i])));

    // Other columns: partial values
    for (const auto &key : keys) {
      real_t val = 0.0;
      auto iterator = partials.find(key);
      if (iterator != partials.end() && i < iterator->second.size()) {
        val = iterator->second[i];
      }
      row_values->push_back(slint::SharedString(std::format("{:.6g}", val)));
    }

    row.values = row_values;
    slint_rows->push_back(row);
  }
  window_.set_table_rows(slint_rows);
}

void PlotController::handleSetCurveColor(int curve_id, const slint::SharedString &color_hex) {
  if (curve_id < 0 || curve_id >= static_cast<int>(current_toggle_keys_.size())) {
    return;
  }
  const std::string &key = current_toggle_keys_[curve_id];
  custom_curve_colors_[key] = std::string(color_hex.data());
  int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    requestPlotUpdate(current_idx, true);
  }
}

void PlotController::updateCurveToggleItems(const correlation::analysis::Histogram *hist) {
  auto toggle_model = std::make_shared<slint::VectorModel<CurveToggleItem>>();
  current_toggle_keys_.clear();

  const auto &partials = hist->smoothed_partials.empty() ? hist->partials : hist->smoothed_partials;
  correlation::plotters::PlotConfig active_config = buildPlotConfigFromUI();

  std::vector<std::string> sorted_partial_keys;
  for (const auto &pair : partials) {
    if (pair.first != "Total") {
      sorted_partial_keys.push_back(pair.first);
    }
  }
  std::sort(sorted_partial_keys.begin(), sorted_partial_keys.end());

  std::size_t total_curves = (partials.contains("Total") ? 1 : 0) + sorted_partial_keys.size() + pinned_runs_.size();

  int curve_id = 0;
  std::size_t color_idx = 0;

  auto get_hex_for_key = [&](const std::string &key) -> std::string {
    if (custom_curve_colors_.contains(key) && !custom_curve_colors_[key].empty()) {
      return custom_curve_colors_[key];
    }
    return correlation::plotters::detail::color(color_idx++, total_curves, active_config.palette);
  };

  if (partials.contains("Total")) {
    bool vis = curve_visibility_map_.contains("Total") ? curve_visibility_map_["Total"] : true;
    curve_visibility_map_["Total"] = vis;
    current_toggle_keys_.push_back("Total");
    std::string color_hex = get_hex_for_key("Total");
    toggle_model->push_back(CurveToggleItem{
        .id = curve_id++,
        .label = slint::SharedString("Total"),
        .color_hex = parseHexColor(color_hex),
        .visible = vis,
        .is_partial = false,
        .is_pinned = false,
    });
  }

  std::size_t rank = 0;
  for (const auto &p_key : sorted_partial_keys) {
    bool default_vis = (rank < 7);
    bool vis = curve_visibility_map_.contains(p_key) ? curve_visibility_map_[p_key] : default_vis;
    curve_visibility_map_[p_key] = vis;
    std::string color_hex = get_hex_for_key(p_key);
    current_toggle_keys_.push_back(p_key);
    toggle_model->push_back(CurveToggleItem{
        .id = curve_id++,
        .label = slint::SharedString(p_key),
        .color_hex = parseHexColor(color_hex),
        .visible = vis,
        .is_partial = true,
        .is_pinned = false,
    });
    rank++;
  }

  for (std::size_t i = 0; i < pinned_runs_.size(); ++i) {
    std::string pin_label = pinned_runs_[i].label;
    bool vis = curve_visibility_map_.contains(pin_label) ? curve_visibility_map_[pin_label] : true;
    curve_visibility_map_[pin_label] = vis;
    std::string color_hex = get_hex_for_key(pin_label);
    current_toggle_keys_.push_back(pin_label);
    toggle_model->push_back(CurveToggleItem{
        .id = curve_id++,
        .label = slint::SharedString(pin_label),
        .color_hex = parseHexColor(color_hex),
        .visible = vis,
        .is_partial = false,
        .is_pinned = true,
    });
  }

  window_.set_curve_toggle_items(toggle_model);
}

bool PlotController::isPlotCacheHit(int index, const correlation::plotters::PlotConfig &config,
                                    const correlation::plotters::HoverInfo &hover) const {
  return (index == last_rendered_index_ && pinned_runs_.size() == last_pinned_runs_count_ &&
          config.theme == last_config_.theme && config.preset_size == last_config_.preset_size &&
          std::abs(config.width - last_config_.width) < 1e-2 && std::abs(config.height - last_config_.height) < 1e-2 &&
          config.palette == last_config_.palette && std::abs(config.font_scale - last_config_.font_scale) < 1e-4 &&
          std::abs(config.line_width - last_config_.line_width) < 1e-4 && config.show_grid == last_config_.show_grid &&
          config.show_legend == last_config_.show_legend && config.show_markers == last_config_.show_markers &&
          config.fill_area == last_config_.fill_area && hover.active == last_hover_.active &&
          std::abs(hover.mouse_x - last_hover_.mouse_x) < 1e-2 &&
          std::abs(hover.mouse_y - last_hover_.mouse_y) < 1e-2 &&
          std::abs(hover.widget_width - last_hover_.widget_width) < 1e-2 &&
          std::abs(hover.widget_height - last_hover_.widget_height) < 1e-2);
}

void PlotController::executePlotRender(RenderTaskData data) {
  if (render_thread_.joinable()) {
    render_thread_.join();
  }
  render_thread_ = std::thread([this, data = std::move(data)]() mutable {
    data.config.show_difference_curve = show_difference_curve_;
    std::string svg;
    if (data.comparison_hists.size() <= 1) {
      svg =
          correlation::plotters::renderHistogramAsSvg(data.active_hist, data.config, data.hover, data.ashcroft_weights,
                                                      data.curve_visibility, custom_curve_colors_);
    } else {
      std::string key = "Total";
      const auto &partials =
          data.active_hist.smoothed_partials.empty() ? data.active_hist.partials : data.active_hist.smoothed_partials;
      if (!partials.empty() && !partials.contains(key)) {
        key = partials.begin()->first;
      }
      svg = correlation::plotters::renderComparisonSvg(data.comparison_hists, key, data.config, data.hover);
    }

    static std::atomic<uint64_t> file_counter{0};
    auto temp_dir = std::filesystem::temp_directory_path();
    auto temp_path = temp_dir / ("correlation_preview_" + std::to_string(file_counter++) + ".svg");

    std::string final_temp_path;
    std::ofstream out(temp_path);
    if (out) {
      out << svg;
      out.close();
      final_temp_path = temp_path.string();
    }

    slint::invoke_from_event_loop([this, final_temp_path, svg]() {
      if (!final_temp_path.empty()) {
        auto img = slint::Image::load_from_path(slint::SharedString(final_temp_path));
        window_.set_preview_plot(img);

        std::error_code error_code;
        std::filesystem::remove(final_temp_path, error_code);
      } else {
        const auto *svg_bytes = std::bit_cast<const uint8_t *>(svg.data());
        auto img =
            slint::private_api::load_image_from_embedded_data(std::span<const uint8_t>(svg_bytes, svg.size()), "svg");

        window_.set_preview_plot(img);
      }

      is_rendering_ = false;

      if (render_pending_) {
        render_pending_ = false;
        requestPlotUpdate(pending_plot_index_, true);
      }
    });
  });
}

void PlotController::handleSavePlot() {
  int index = window_.get_selected_plot_index();
  if (index < 0 || index >= static_cast<int>(available_plot_keys_.size())) {
    return;
  }
  const std::string &name = available_plot_keys_[index];
  const correlation::analysis::Histogram *hist = backend_.getHistogram(name);
  if (hist == nullptr) {
    return;
  }

  std::string default_path = backend_.options().output_file_base;
  if (!default_path.empty()) {
    default_path += "_" + name;
  } else {
    default_path = name;
  }

  std::filesystem::path def_path(default_path);
  std::string default_dir = def_path.parent_path().string();
  std::string default_name = def_path.filename().string();

  std::array<nfdfilteritem_t, 2> filterList = {{{
                                                    .name = "SVG Image",
                                                    .spec = "svg",
                                                },
                                                {
                                                    .name = "PDF Document",
                                                    .spec = "pdf",
                                                }}};
  nfdfiltersize_t filterCount = filterList.size();

  nfdchar_t *outPath = nullptr;
  nfdresult_t result =
      NFD_SaveDialogU8(&outPath, filterList.data(), filterCount, default_dir.empty() ? nullptr : default_dir.c_str(),
                       default_name.empty() ? nullptr : default_name.c_str());

  if (result == NFD_OKAY) {
    std::string filepath(outPath);
    NFD_FreePathU8(outPath);
    executeSavePlot(filepath, hist, name);
  } else if (result == NFD_CANCEL) {
    window_.set_analysis_status_text(slint::SharedString(AppDefaults::MSG_SAVE_CANCELLED));
  } else {
    std::string error_msg = "Error: ";
    error_msg += NFD_GetError();
    window_.set_analysis_status_text(slint::SharedString(error_msg));
  }
}

void PlotController::executeSavePlot(const std::string &filepath, const correlation::analysis::Histogram *hist,
                                     const std::string &name) {
  correlation::plotters::PlotConfig config = buildPlotConfigFromUI();
  config.use_native_text = true;

  if (config.preset_size == correlation::plotters::PlotConfig::PresetSize::Default) {
    config.width = 1200.0;
    config.height = 900.0;
  }

  std::string ext = std::filesystem::path(filepath).extension().string();
  bool save_as_pdf = (ext == ".pdf" || ext == ".PDF");

  auto build_datasets = [&]() {
    std::vector<correlation::plotters::LabeledHistogram> datasets;
    datasets.push_back({"Current", hist});
    for (const auto &pinned_run : pinned_runs_) {
      auto hist_it = pinned_run.histograms.find(name);
      if (hist_it != pinned_run.histograms.end()) {
        datasets.push_back({pinned_run.label, &hist_it->second});
      }
    }
    return datasets;
  };

  if (save_as_pdf) {
    if (pinned_runs_.empty()) {
      correlation::plotters::renderHistogramAsPdf(*hist, filepath, config);
    } else {
      correlation::plotters::renderComparisonPdf(build_datasets(), {getComparisonKey(hist), filepath}, config);
    }
    window_.set_analysis_status_text(slint::SharedString(name + " plot saved as PDF successfully."));
  } else {
    std::string svg;
    if (pinned_runs_.empty()) {
      svg = correlation::plotters::renderHistogramAsSvg(*hist, config, {}, backend_.getAshcroftWeights());
    } else {
      svg = correlation::plotters::renderComparisonSvg(build_datasets(), getComparisonKey(hist), config);
    }
    std::ofstream out(filepath);
    if (out.is_open()) {
      out << svg;
      out.close();
      window_.set_analysis_status_text(slint::SharedString(name + " plot saved successfully."));
    } else {
      window_.set_analysis_status_text(slint::SharedString("Failed to open file for saving plot."));
    }
  }
}

void PlotController::handlePinRun() {
  const auto &hists = backend_.getHistograms();
  if (hists.empty()) {
    return;
  }

  std::string label = "Run " + std::to_string(pinned_runs_.size() + 1);
  pinned_runs_.push_back({label, hists});

  slint::invoke_from_event_loop([this]() {
    window_.set_pinned_runs_count(static_cast<int>(pinned_runs_.size()));

    int current_idx = window_.get_selected_plot_index();
    if (current_idx >= 0) {
      requestPlotUpdate(current_idx, true);
    }
  });
}

void PlotController::handleClearPinnedRuns() {
  pinned_runs_.clear();

  slint::invoke_from_event_loop([this]() {
    window_.set_pinned_runs_count(0);

    int current_idx = window_.get_selected_plot_index();
    if (current_idx >= 0) {
      requestPlotUpdate(current_idx, true);
    }
  });
}

void PlotController::handleToggleCurveVisibility(int curve_id, bool visible) {
  if (curve_id >= 0 && curve_id < static_cast<int>(current_toggle_keys_.size())) {
    const std::string &key = current_toggle_keys_[curve_id];
    curve_visibility_map_[key] = visible;
  }
  slint::invoke_from_event_loop([this]() {
    int current_idx = window_.get_selected_plot_index();
    if (current_idx >= 0) {
      requestPlotUpdate(current_idx, true);
    }
  });
}

void PlotController::handleToggleAllCurves(bool visible) {
  for (const auto &key : current_toggle_keys_) {
    curve_visibility_map_[key] = visible;
  }
  slint::invoke_from_event_loop([this]() {
    int current_idx = window_.get_selected_plot_index();
    if (current_idx >= 0) {
      requestPlotUpdate(current_idx, true);
    }
  });
}

void PlotController::handleToggleDifferencePlot(bool show_difference) {
  show_difference_curve_ = show_difference;
  slint::invoke_from_event_loop([this]() {
    int current_idx = window_.get_selected_plot_index();
    if (current_idx >= 0) {
      requestPlotUpdate(current_idx, true);
    }
  });
}

} // namespace correlation::app
