/**
 * @file AppController.cpp
 * @brief Implementation of the application controller.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif

#include "AppWindow.h"
#include "app/AnalysisRunner.hpp"
#include "app/AppController.hpp"

#include "app/BondCutoffMapper.hpp"
#include "app/FileIOHandler.hpp"
#include "app/InputValidator.hpp"
#include "app/PlotController.hpp"
#include "app/PresetController.hpp"
#include "app/UpdateChecker.hpp"
#include "calculators/CalculatorFactory.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <map>
#include <string>
#include <vector>

namespace correlation::app {

AppController::AppController(::AppWindow &window, AppBackend &backend)
    : window_(window), backend_(backend) {
  // Initialize Native File Dialog
  NFD_Init();

  // Populate the calculator groups from the factory
  populateCalculatorGroups();

  // Set export formats availability based on compile-time definitions
#ifdef CORRELATION_USE_HDF5
  window_.set_hdf5_available(true);
#else
  window_.set_hdf5_available(false);
#endif

#ifdef CORRELATION_USE_ARROW
  window_.set_parquet_available(true);
#else
  window_.set_parquet_available(false);
#endif

  // Set application version string for UI display
  window_.set_app_version(slint::SharedString(std::string("v") + CORRELATION_VERSION_STRING));

  analysis_runner_ = std::make_unique<AnalysisRunner>(window_, backend_, *this);
  file_io_handler_ = std::make_unique<FileIOHandler>(window_, backend_, *this);
  input_validator_ = std::make_unique<InputValidator>(window_, backend_, *this);
  plot_controller_ = std::make_unique<PlotController>(window_, backend_);
  preset_controller_ = std::make_unique<PresetController>(window_, backend_, *this);

  // default options to UI
  handleOptionstoUI();

  // Connect the UI signals to the controller's member functions.
  // We use lambdas to capture 'this' and call the appropriate method.
  window_.on_run_analysis([this]() { analysis_runner_->handleRunAnalysis(); });
  window_.on_cancel_analysis([this]() { backend_.cancel_analysis(); });
  window_.on_browse_file([this]() { file_io_handler_->handleBrowseFile(); });
  window_.on_reload_file([this]() { file_io_handler_->handleReloadFile(); });
  window_.on_write_files([this]() { file_io_handler_->handleWriteFiles(); });
  window_.on_validate_inputs([this]() {
    slint::invoke_from_event_loop(
        [this]() { static_cast<void>(input_validator_->validateInputs()); });
  });

  // Handle calculator toggle: update backend options and refresh the UI model
  window_.on_toggle_calculator([this](const slint::SharedString &calc_id, bool enabled) {
    backend_.setCalculatorActive(std::string(calc_id.data()), enabled);
    populateCalculatorGroups();
    updateActiveGroupFlags();
  });

  // Handle plot selection: generate SVG and push to UI
  window_.on_select_plot([this](int index) {
    slint::invoke_from_event_loop(
        [this, index]() { plot_controller_->requestPlotUpdate(index, true); });
  });

  // Handle curve visibility toggle from PreviewCard checklist
  window_.on_toggle_curve_visibility([this](int curve_id, bool visible) {
    plot_controller_->handleToggleCurveVisibility(curve_id, visible);
  });

  window_.on_toggle_all_curves(
      [this](bool visible) { plot_controller_->handleToggleAllCurves(visible); });

  window_.on_set_curve_color([this](int curve_id, const slint::SharedString &color_hex) {
    plot_controller_->handleSetCurveColor(curve_id, color_hex);
  });

  // Handle mouse move on preview plot
  window_.on_mouse_move(
      [this](float mouse_x, float mouse_y, bool hover, float width, float height) {
        plot_controller_->handleMouseMove(mouse_x, mouse_y, hover, width, height);
      });

  // Handle save plot request (SVG or PDF)
  window_.on_save_plot([this]() { plot_controller_->handleSavePlot(); });

  // Handle pin run request
  window_.on_pin_run([this]() { plot_controller_->handlePinRun(); });

  // Handle clear pinned runs request
  window_.on_clear_pinned_runs([this]() { plot_controller_->handleClearPinnedRuns(); });

  // Handle difference plot toggle
  window_.on_toggle_difference_plot(
      [this](bool show) { plot_controller_->handleToggleDifferencePlot(show); });

  // Handle preset load, save, delete requests
  window_.on_load_preset([this](int index) {
    slint::invoke_from_event_loop([this, index]() { preset_controller_->handleLoadPreset(index); });
  });
  window_.on_save_preset([this](const slint::SharedString &name) {
    slint::invoke_from_event_loop([this, sname = std::string(name.data())]() {
      preset_controller_->handleSavePreset(sname);
    });
  });
  window_.on_delete_preset([this](int index) {
    slint::invoke_from_event_loop(
        [this, index]() { preset_controller_->handleDeletePreset(index); });
  });
  window_.on_material_type_changed([this](int type) {
    slint::invoke_from_event_loop(
        [this, type]() { preset_controller_->handleMaterialTypeChanged(type); });
  });

  // Handle plot resized callback from UI
  window_.on_plot_resized([this](float width, float height) {
    plot_controller_->handlePlotResized({
        .width = width,
        .height = height,
    });
  });

  // Handle layout geometry changed callback from UI
  window_.on_layout_geometry_changed([this](float left_w, float mid_w) {
    settings_.left_col_width = left_w;
    settings_.middle_col_width = mid_w;
    saveSettings();
  });

  // Handle reset bond cutoffs request from UI
  window_.on_reset_bond_cutoffs([this]() { setBondCutoffs(); });
  window_.on_reset_rdf_options([this]() { handleResetRDFOptions(); });
  window_.on_reset_angle_options([this]() { handleResetAngleOptions(); });
  window_.on_reset_sq_options([this]() { handleResetSQOptions(); });
  window_.on_reset_xrd_options([this]() { handleResetXRDOptions(); });
  window_.on_xrd_radiation_preset_changed([this](int idx) { handleXRDPresetChanged(idx); });
  window_.on_apply_scaled_cutoffs([this](float factor) { handleApplyScaledCutoffs(factor); });
  window_.on_set_uniform_cutoff([this](float max_val) { handleSetUniformCutoff(max_val); });
  window_.on_reset_rings_options([this]() { handleResetRingsOptions(); });
  window_.on_reset_smoothing_options([this]() { handleResetSmoothingOptions(); });
  window_.on_reset_advanced_options([this]() { handleResetAdvancedOptions(); });
  window_.on_reset_trajectory_options([this]() { handleResetTrajectoryOptions(); });
  window_.on_reset_export_settings([this]() { handleResetExportSettings(); });
  window_.on_reset_analyses_selection([this]() { handleResetAnalysesSelection(); });
  window_.on_reset_material_type([this]() { handleResetMaterialType(); });
  window_.on_reset_presets([this]() { window_.set_selected_preset(-1); });
  window_.on_clear_comparison_curves([this]() { handleClearComparisonCurves(); });

  // Handle open external URL (e.g. download update from browser)
  window_.on_open_url([](const slint::SharedString &url) {
    UpdateChecker::openUrlInBrowser(std::string(url.data()));
  });

  // Initial load of settings and preset list
  loadSettings();
  preset_controller_->refreshPresetList();

  // Asynchronously check for GitHub releases in the background without interrupting the user
  UpdateChecker::checkForUpdatesAsync(window_, CORRELATION_VERSION_STRING);
}

AppController::~AppController() {
  saveSettings();
  // Quit Native File Dialog
  NFD_Quit();
}

void AppController::loadSettings() {
  settings_ = SettingsManager::load();
  window_.set_left_col_width(settings_.left_col_width);
  window_.set_last_col_width(settings_.left_col_width);
  window_.set_middle_col_width(settings_.middle_col_width);
  if (settings_.window_width > 0 && settings_.window_height > 0) {
    window_.window().set_size(
        slint::PhysicalSize({.width = settings_.window_width, .height = settings_.window_height}));
  }
}

void AppController::saveSettings() const {
  AppSettings settings = settings_;
  const float last_w = window_.get_last_col_width();
  settings.left_col_width = (last_w >= 180.0F) ? last_w : window_.get_left_col_width();
  settings.middle_col_width = window_.get_middle_col_width();
  const auto phys_size = window_.window().size();
  if (phys_size.width > 0 && phys_size.height > 0) {
    settings.window_width = phys_size.width;
    settings.window_height = phys_size.height;
  }
  SettingsManager::save(settings);
}

// Safe conversion helper
namespace {
/**
 * @brief Safely converts a Slint SharedString to a numeric type with a default
 * fallback.
 *
 * @tparam T The numeric type to return (e.g., float, real_t).
 * @param str The Slint string to parse.
 * @param default_value The value to return if parsing fails.
 * @return The parsed value or default_value on error.
 */
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
    // Optionally, log the error or update a UI status message
    return default_value;
  }
}
} // namespace

void AppController::handleOptionstoUI() {
  ProgramOptions opt = backend_.options();
  window_.set_in_file_text(slint::SharedString(opt.input_file));
  {
    auto opts = window_.get_analysis_options();
    opts.smoothing_enabled = opt.smoothing;
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.r_max = slint::SharedString(std::format("{:.2f}", opt.r_max));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.r_bin_width = slint::SharedString(std::format("{:.2f}", opt.r_bin_width));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.q_max = slint::SharedString(std::format("{:.2f}", opt.q_max));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.q_bin_width = slint::SharedString(std::format("{:.2f}", opt.q_bin_width));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.r_int_max = slint::SharedString(std::format("{:.2f}", opt.r_int_max));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.angle_bin_width = slint::SharedString(std::format("{:.2f}", opt.angle_bin_width));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.dihedral_bin_width = slint::SharedString(std::format("{:.2f}", opt.dihedral_bin_width));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.max_ring_size = slint::SharedString(std::to_string(opt.max_ring_size));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.smoothing_sigma = slint::SharedString(std::format("{:.2f}", opt.smoothing_sigma));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.smoothing_kernel = static_cast<int>(opt.smoothing_kernel);
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.material_type = opt.material_type;
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.lef_cutoff = slint::SharedString(std::format("{:.2f}", opt.lef_cutoff));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.lef_sigma = slint::SharedString(std::format("{:.2f}", opt.lef_sigma));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.hyper_samples = slint::SharedString(std::to_string(opt.hyper_samples));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.xrd_radiation_preset = 0;
    opts.xrd_lambda = slint::SharedString(std::format("{:.4f}", opt.xrd_params.lambda));
    opts.xrd_theta_min = slint::SharedString(std::format("{:.1f}", opt.xrd_params.theta_min));
    opts.xrd_theta_max = slint::SharedString(std::format("{:.1f}", opt.xrd_params.theta_max));
    opts.xrd_bin_width = slint::SharedString(std::format("{:.2f}", opt.xrd_params.bin_width));
    window_.set_analysis_options(opts);
  }

  {
    auto opts = window_.get_analysis_options();
    opts.min_frame = slint::SharedString(std::to_string(opt.min_frame + 1));
    window_.set_analysis_options(opts);
  } // UI is 1-based
  if (opt.max_frame == -1) {
    {
      auto opts = window_.get_analysis_options();
      opts.max_frame = "End";
      window_.set_analysis_options(opts);
    }
  } else {
    {
      auto opts = window_.get_analysis_options();
      opts.max_frame = slint::SharedString(std::to_string(opt.max_frame));
      window_.set_analysis_options(opts);
    }
  }
  {
    auto opts = window_.get_analysis_options();
    opts.time_step = slint::SharedString(std::format("{:.2f}", opt.time_step));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.frame_stride = slint::SharedString(std::to_string(opt.frame_stride));
    window_.set_analysis_options(opts);
  }
  updateActiveGroupFlags();
};

void AppController::updateActiveGroupFlags() {
  const auto &calculators =
      ::correlation::calculators::CalculatorFactory::instance().getCalculators();
  const auto &opts = backend_.options();

  bool has_radial = false;
  bool has_scattering = false;
  bool has_angular = false;
  bool has_rings = false;

  for (const auto &calc : calculators) {
    std::string_view const grp = calc->getGroup();
    bool enabled = true; // default on
    auto calc_it = opts.active_calculators.find(std::string(calc->getName()));
    if (calc_it != opts.active_calculators.end()) {
      enabled = calc_it->second;
    }
    if (enabled) {
      if (grp == "Radial") {
        has_radial = true;
      } else if (grp == "Scattering") {
        has_scattering = true;
      } else if (grp == "Angular") {
        has_angular = true;
      } else if (grp == "Rings") {
        has_rings = true;
      }
    }
  }

  window_.set_has_radial_active(has_radial);
  window_.set_has_scattering_active(has_scattering);
  window_.set_has_angular_active(has_angular);
  window_.set_has_rings_active(has_rings);
}

ProgramOptions AppController::handleOptionsfromUI() {
  ProgramOptions opt;
  const std::string input_path_str = window_.get_in_file_text().data();
  const std::filesystem::path full_path(input_path_str);
  std::filesystem::path output_path = full_path.parent_path() / full_path.stem();
  opt.input_file = input_path_str;
  opt.output_file_base = output_path.make_preferred().string();
  opt.smoothing = true;
  opt.r_max = safeParse(window_.get_analysis_options().r_max, opt.r_max);
  opt.r_bin_width = safeParse(window_.get_analysis_options().r_bin_width, opt.r_bin_width);
  opt.q_max = safeParse(window_.get_analysis_options().q_max, opt.q_max);
  opt.q_bin_width = safeParse(window_.get_analysis_options().q_bin_width, opt.q_bin_width);
  opt.r_int_max = safeParse(window_.get_analysis_options().r_int_max, opt.r_int_max);
  opt.angle_bin_width =
      safeParse(window_.get_analysis_options().angle_bin_width, opt.angle_bin_width);
  opt.dihedral_bin_width =
      safeParse(window_.get_analysis_options().dihedral_bin_width, opt.dihedral_bin_width);
  opt.max_ring_size = static_cast<size_t>(safeParse(window_.get_analysis_options().max_ring_size,
                                                    static_cast<real_t>(opt.max_ring_size)));
  opt.hyper_samples = static_cast<size_t>(safeParse(window_.get_analysis_options().hyper_samples,
                                                    static_cast<real_t>(opt.hyper_samples)));

  // Collect active_calculators from the UI model
  const auto groups = window_.get_calculator_groups();
  for (size_t gi = 0; gi < groups->row_count(); ++gi) {
    const auto maybe_group = groups->row_data(gi);
    if (!maybe_group.has_value()) {
      continue;
    }
    const auto &group = maybe_group.value();
    for (size_t ci = 0; ci < group.calculators->row_count(); ++ci) {
      const auto maybe_calc = group.calculators->row_data(ci);
      if (!maybe_calc.has_value()) {
        continue;
      }
      const auto &calc = maybe_calc.value();
      opt.active_calculators[std::string(calc.id.data())] = calc.enabled;
    }
  }

  opt.smoothing_sigma =
      safeParse(window_.get_analysis_options().smoothing_sigma, opt.smoothing_sigma);
  opt.smoothing_kernel =
      static_cast<correlation::math::KernelType>(window_.get_analysis_options().smoothing_kernel);
  opt.material_type = window_.get_analysis_options().material_type;
  opt.lef_cutoff = safeParse(window_.get_analysis_options().lef_cutoff, opt.lef_cutoff);
  opt.lef_sigma = safeParse(window_.get_analysis_options().lef_sigma, opt.lef_sigma);
  opt.xrd_params.lambda =
      safeParse(window_.get_analysis_options().xrd_lambda, opt.xrd_params.lambda);
  opt.xrd_params.theta_min =
      safeParse(window_.get_analysis_options().xrd_theta_min, opt.xrd_params.theta_min);
  opt.xrd_params.theta_max =
      safeParse(window_.get_analysis_options().xrd_theta_max, opt.xrd_params.theta_max);
  opt.xrd_params.bin_width =
      safeParse(window_.get_analysis_options().xrd_bin_width, opt.xrd_params.bin_width);

  // Parse Frame Selection
  // - Handles string presets "start" and "end" case-insensitively.
  // - Numeric values are 1-based in UI, converted to 0-based for backend.

  // Helper lambda for case-insensitive comparison
  auto to_lower = [](const std::string &str) -> std::string {
    std::string data = str;
    std::ranges::transform(data, data.begin(),
                           [](unsigned char chr) { return static_cast<char>(std::tolower(chr)); });
    return data;
  };

  // Frame Selection
  try {
    const std::string min_s = window_.get_analysis_options().min_frame.data();
    const std::string min_s_lower = to_lower(min_s);

    if (min_s_lower == "start") {
      opt.min_frame = 0;
    } else if (min_s_lower == "end") {
      opt.min_frame = std::max(0, static_cast<int>(backend_.getFrameCount()) - 1);
    } else {
      opt.min_frame = std::max(0, std::stoi(min_s) - 1); // UI is 1-based
    }
  } catch (const std::exception &) {
    opt.min_frame = 0;
  }

  try {
    const std::string max_s = window_.get_analysis_options().max_frame.data();
    const std::string max_s_lower = to_lower(max_s);

    if (max_s_lower == "end" || max_s.empty()) {
      opt.max_frame = -1;
    } else if (max_s_lower == "start") {
      opt.max_frame = 1; // 1-based index 1 -> implies reading only the first
                         // frame. In TrajectoryAnalyzer loop: i < max_frame. So
                         // max_frame=1 means process frame 0 only.
      opt.max_frame = 1;
    } else {
      opt.max_frame = std::stoi(max_s);
    }
  } catch (const std::exception &) {
    opt.max_frame = -1;
  }

  opt.time_step = safeParse(window_.get_analysis_options().time_step, opt.time_step);

  try {
    const std::string stride_s = window_.get_analysis_options().frame_stride.data();
    opt.frame_stride = std::max(1, std::stoi(stride_s));
  } catch (const std::exception &) {
    opt.frame_stride = 1;
  }

  // Handle Bond Cutoffs
  opt.bond_cutoffs = getBondCutoffs();

  return opt;
};

void AppController::setBondCutoffs() {
  if (backend_.cell() == nullptr) {
    return;
  }

  const auto &elements = backend_.cell()->elements();
  const auto entries = BondCutoffMapper::createDefaultCutoffEntries(elements);

  auto slint_cutoffs = std::make_shared<slint::VectorModel<BondCutoff>>();
  for (const auto &entry : entries) {
    slint_cutoffs->push_back({
        .element1 = slint::SharedString(entry.element1),
        .element2 = slint::SharedString(entry.element2),
        .min_distance = slint::SharedString(entry.min_distance),
        .max_distance = slint::SharedString(entry.max_distance),
    });
  }

  window_.set_bond_cutoffs(slint_cutoffs);
  window_.set_bond_cutoffs_reset_trigger(window_.get_bond_cutoffs_reset_trigger() + 1);
  backend_.setBondCutoffs(getBondCutoffs());
}

correlation::analysis::BondCutoffMatrix AppController::getBondCutoffs() {
  auto slint_cutoffs = window_.get_bond_cutoffs();
  if (backend_.cell() == nullptr || slint_cutoffs == nullptr) {
    return {};
  }

  const auto &elements = backend_.cell()->elements();
  std::vector<CutoffEntry> entries;
  entries.reserve(slint_cutoffs->row_count());

  for (size_t k = 0; k < slint_cutoffs->row_count(); ++k) {
    auto maybe_item = slint_cutoffs->row_data(k);
    if (!maybe_item.has_value()) {
      continue;
    }
    const auto &item = maybe_item.value();
    entries.push_back(CutoffEntry{
        .element1 = std::string(item.element1.data()),
        .element2 = std::string(item.element2.data()),
        .min_distance = std::string(item.min_distance.data()),
        .max_distance = std::string(item.max_distance.data()),
    });
  }

  return BondCutoffMapper::parseCutoffMatrix(entries, elements);
}

void AppController::populateCalculatorGroups() {
  const auto &calculators =
      ::correlation::calculators::CalculatorFactory::instance().getCalculators();
  const auto &opts = backend_.options();

  // Collect group names in insertion order
  std::vector<std::string> group_order;
  std::map<std::string, std::vector<const correlation::calculators::BaseCalculator *>> groups_map;
  for (const auto &calc : calculators) {
    std::string const grp = std::string(calc->getGroup());
    if (!groups_map.contains(grp)) {
      group_order.push_back(grp);
    }
    groups_map[grp].push_back(calc.get());
  }

  auto groups_model = std::make_shared<slint::VectorModel<CalculatorGroup>>();
  for (const auto &grp_name : group_order) {
    auto calcs_model = std::make_shared<slint::VectorModel<CalculatorInfo>>();
    for (const auto *calc : groups_map.at(grp_name)) {
      bool enabled = true; // default on
      auto calc_it = opts.active_calculators.find(std::string(calc->getName()));
      if (calc_it != opts.active_calculators.end()) {
        enabled = calc_it->second;
      }
      CalculatorInfo info;
      info.id = slint::SharedString(calc->getName());
      info.name = slint::SharedString(calc->getName());
      info.description = slint::SharedString(calc->getDescription());
      info.enabled = enabled;
      calcs_model->push_back(info);
    }
    CalculatorGroup group;
    group.name = slint::SharedString(grp_name);
    group.calculators = calcs_model;
    groups_model->push_back(group);
  }
  window_.set_calculator_groups(groups_model);
  window_.set_total_calculator_count(static_cast<int>(calculators.size()));
}

void AppController::handleResetRDFOptions() {
  auto opts = window_.get_analysis_options();
  opts.r_max = slint::SharedString(std::format("{:.2f}", AppDefaults::R_MAX));
  if (opts.material_type == 2) {
    opts.r_bin_width = slint::SharedString(std::format("{:.3f}", AppDefaults::R_BIN_WIDTH_CRYSTAL));
  } else if (opts.material_type == 1) {
    opts.r_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::R_BIN_WIDTH_LIQUID));
  } else {
    opts.r_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::R_BIN_WIDTH));
  }
  window_.set_analysis_options(opts);
  static_cast<void>(input_validator_->validateInputs());
}

void AppController::handleResetAngleOptions() {
  auto opts = window_.get_analysis_options();
  if (opts.material_type == 2) {
    opts.angle_bin_width =
        slint::SharedString(std::format("{:.2f}", AppDefaults::ANGLE_BIN_WIDTH_CRYSTAL));
    opts.dihedral_bin_width =
        slint::SharedString(std::format("{:.2f}", AppDefaults::ANGLE_BIN_WIDTH_CRYSTAL));
  } else if (opts.material_type == 1) {
    opts.angle_bin_width =
        slint::SharedString(std::format("{:.2f}", AppDefaults::ANGLE_BIN_WIDTH_LIQUID));
    opts.dihedral_bin_width =
        slint::SharedString(std::format("{:.2f}", AppDefaults::ANGLE_BIN_WIDTH_LIQUID));
  } else {
    opts.angle_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::ANGLE_BIN_WIDTH));
    opts.dihedral_bin_width =
        slint::SharedString(std::format("{:.2f}", AppDefaults::ANGLE_BIN_WIDTH));
  }
  window_.set_analysis_options(opts);
  static_cast<void>(input_validator_->validateInputs());
}

void AppController::handleResetSQOptions() {
  auto opts = window_.get_analysis_options();
  opts.q_max = slint::SharedString(std::format("{:.2f}", AppDefaults::Q_MAX));
  opts.r_int_max = slint::SharedString(std::format("{:.2f}", AppDefaults::R_INT_MAX));
  if (opts.material_type == 2) {
    opts.q_bin_width = slint::SharedString(std::format("{:.3f}", AppDefaults::Q_BIN_WIDTH_CRYSTAL));
  } else if (opts.material_type == 1) {
    opts.q_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::Q_BIN_WIDTH_LIQUID));
  } else {
    opts.q_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::Q_BIN_WIDTH));
  }
  window_.set_analysis_options(opts);
  static_cast<void>(input_validator_->validateInputs());
}

void AppController::handleResetRingsOptions() {
  auto opts = window_.get_analysis_options();
  opts.max_ring_size = slint::SharedString(std::to_string(ProgramOptions{}.max_ring_size));
  window_.set_analysis_options(opts);
  static_cast<void>(input_validator_->validateInputs());
}

void AppController::handleResetSmoothingOptions() {
  auto opts = window_.get_analysis_options();
  opts.smoothing_enabled = ProgramOptions{}.smoothing;
  opts.smoothing_kernel = static_cast<int>(AppDefaults::SMOOTHING_KERNEL);
  if (opts.material_type == 2) {
    opts.smoothing_sigma =
        slint::SharedString(std::format("{:.2f}", AppDefaults::SMOOTHING_SIGMA_CRYSTAL));
  } else if (opts.material_type == 1) {
    opts.smoothing_sigma =
        slint::SharedString(std::format("{:.2f}", AppDefaults::SMOOTHING_SIGMA_LIQUID));
  } else {
    opts.smoothing_sigma = slint::SharedString(std::format("{:.2f}", AppDefaults::SMOOTHING_SIGMA));
  }
  window_.set_analysis_options(opts);
  static_cast<void>(input_validator_->validateInputs());
}

void AppController::handleResetAdvancedOptions() {
  auto opts = window_.get_analysis_options();
  opts.lef_cutoff = slint::SharedString(std::format("{:.2f}", AppDefaults::LEF_CUTOFF));
  opts.lef_sigma = slint::SharedString(std::format("{:.2f}", AppDefaults::LEF_SIGMA));
  opts.hyper_samples = slint::SharedString(std::to_string(ProgramOptions{}.hyper_samples));
  window_.set_analysis_options(opts);
  static_cast<void>(input_validator_->validateInputs());
}

void AppController::handleResetTrajectoryOptions() {
  auto opts = window_.get_analysis_options();
  if (backend_.getFrameCount() > 0) {
    opts.time_step = slint::SharedString(std::format("{:.2f}", backend_.getRecommendedTimeStep()));
    opts.min_frame = "1";
    opts.max_frame = slint::SharedString(std::to_string(backend_.getFrameCount()));
  } else {
    opts.time_step = slint::SharedString(std::format("{:.2f}", AppDefaults::TIME_STEP));
    opts.min_frame = "1";
    opts.max_frame = "End";
  }
  opts.frame_stride = "1";
  window_.set_analysis_options(opts);
  static_cast<void>(input_validator_->validateInputs());
}

void AppController::handleResetExportSettings() {
  ExportConfig cfg;
  cfg.size_preset = 0;
  cfg.palette = 0;
  cfg.font_scale = "1.0";
  cfg.line_width = "3.0";
  cfg.marker_size = "3.5";
  cfg.show_legend = true;
  cfg.show_grid = true;
  cfg.show_markers = false;
  cfg.fill_area = false;
  window_.set_export_config(cfg);
  if (plot_controller_) {
    plot_controller_->requestPlotUpdate(window_.get_selected_plot_index(), true);
  }
}

void AppController::handleResetAnalysesSelection() {
  populateCalculatorGroups();
  updateActiveGroupFlags();
}

void AppController::handleResetMaterialType() {
  auto opts = window_.get_analysis_options();
  opts.material_type = 0;
  window_.set_analysis_options(opts);
  if (preset_controller_) {
    preset_controller_->handleMaterialTypeChanged(0);
  }
}

void AppController::handleClearComparisonCurves() {
  auto empty_model = std::make_shared<slint::VectorModel<ComparisonCurveData>>();
  window_.set_comparison_curves(empty_model);
  if (plot_controller_) {
    plot_controller_->handleClearPinnedRuns();
  }
}

void AppController::handleResetXRDOptions() {
  auto opts = window_.get_analysis_options();
  opts.xrd_radiation_preset = 0;
  opts.xrd_lambda = slint::SharedString(std::format("{:.4f}", AppDefaults::XRD_LAMBDA));
  opts.xrd_theta_min = slint::SharedString(std::format("{:.1f}", AppDefaults::XRD_THETA_MIN));
  opts.xrd_theta_max = slint::SharedString(std::format("{:.1f}", AppDefaults::XRD_THETA_MAX));
  opts.xrd_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::XRD_BIN_WIDTH));
  window_.set_analysis_options(opts);
  static_cast<void>(input_validator_->validateInputs());
}

void AppController::handleXRDPresetChanged(int preset_idx) {
  auto opts = window_.get_analysis_options();
  opts.xrd_radiation_preset = preset_idx;
  switch (preset_idx) {
  case 0: // Cu-Ka
    opts.xrd_lambda = "1.5406";
    break;
  case 1: // Mo-Ka
    opts.xrd_lambda = "0.7107";
    break;
  case 2: // Co-Ka
    opts.xrd_lambda = "1.7890";
    break;
  case 3: // Cr-Ka
    opts.xrd_lambda = "2.2897";
    break;
  case 4: // Custom
  default:
    break;
  }
  window_.set_analysis_options(opts);
  static_cast<void>(input_validator_->validateInputs());
}

void AppController::handleApplyScaledCutoffs(float scale_factor) {
  if (backend_.cell() == nullptr || scale_factor <= 0.0F) {
    return;
  }
  const auto scaled_cutoffs = backend_.applyScaledBondCutoffs(static_cast<real_t>(scale_factor));
  const auto &elements = backend_.cell()->elements();
  const auto num_elements = elements.size();

  auto slint_cutoffs = std::make_shared<slint::VectorModel<BondCutoff>>();
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      const real_t min_d = std::sqrt(scaled_cutoffs[i][j].min_sq);
      const real_t max_d = std::sqrt(scaled_cutoffs[i][j].max_sq);
      slint_cutoffs->push_back(BondCutoff{
          .element1 = slint::SharedString(elements[i].symbol),
          .element2 = slint::SharedString(elements[j].symbol),
          .min_distance = slint::SharedString(std::format("{:.2f}", min_d)),
          .max_distance = slint::SharedString(std::format("{:.2f}", max_d)),
      });
    }
  }

  window_.set_bond_cutoffs(slint_cutoffs);
  window_.set_bond_cutoffs_reset_trigger(window_.get_bond_cutoffs_reset_trigger() + 1);
}

void AppController::handleSetUniformCutoff(float max_cutoff) {
  if (backend_.cell() == nullptr || max_cutoff <= 0.0F) {
    return;
  }
  const auto uniform_cutoffs = backend_.setUniformBondCutoff(0.0, static_cast<real_t>(max_cutoff));
  const auto &elements = backend_.cell()->elements();
  const auto num_elements = elements.size();

  auto slint_cutoffs = std::make_shared<slint::VectorModel<BondCutoff>>();
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = i; j < num_elements; ++j) {
      const real_t min_d = std::sqrt(uniform_cutoffs[i][j].min_sq);
      const real_t max_d = std::sqrt(uniform_cutoffs[i][j].max_sq);
      slint_cutoffs->push_back(BondCutoff{
          .element1 = slint::SharedString(elements[i].symbol),
          .element2 = slint::SharedString(elements[j].symbol),
          .min_distance = slint::SharedString(std::format("{:.2f}", min_d)),
          .max_distance = slint::SharedString(std::format("{:.2f}", max_d)),
      });
    }
  }

  window_.set_bond_cutoffs(slint_cutoffs);
  window_.set_bond_cutoffs_reset_trigger(window_.get_bond_cutoffs_reset_trigger() + 1);
}

} // namespace correlation::app
