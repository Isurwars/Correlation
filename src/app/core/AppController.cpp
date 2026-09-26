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
#include "app/core/AppController.hpp"
#include "app/services/AnalysisRunner.hpp"

#include "app/formatters/BondCutoffMapper.hpp"
#include "app/services/FileIOHandler.hpp"
#include "app/services/InputValidator.hpp"
#include "app/services/OptionsResetService.hpp"
#include "app/services/OptionsSyncService.hpp"
#include "app/services/UpdateChecker.hpp"
#include "app/viewmodel/BondCutoffController.hpp"
#include "app/viewmodel/PlotController.hpp"
#include "app/viewmodel/PresetController.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "physics/PhysicalData.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <map>
#include <string>
#include <vector>

namespace correlation::app {

AppController::AppController(::AppWindow &window, TrajectoryLoader &loader,
                             AnalysisDispatcher &dispatcher, ProgramOptions &options)
    : window_(window), loader_(loader), dispatcher_(dispatcher), options_(options),
      bond_cutoff_controller_(window, loader, options) {
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

  analysis_runner_ =
      std::make_unique<AnalysisRunner>(window_, loader_, dispatcher_, options_, *this);
  file_io_handler_ =
      std::make_unique<FileIOHandler>(window_, loader_, dispatcher_, options_, *this);
  input_validator_ = std::make_unique<InputValidator>(window_, *this);
  plot_controller_ = std::make_unique<PlotController>(window_, dispatcher_, options_);
  preset_controller_ = std::make_unique<PresetController>(window_, options_, *this);

  // default options to UI
  handleOptionstoUI();

  // Connect the UI signals to the controller's member functions.
  // We use lambdas to capture 'this' and call the appropriate method.
  window_.on_run_analysis([this]() { analysis_runner_->handleRunAnalysis(); });
  window_.on_cancel_analysis([this]() { dispatcher_.cancelAnalysis(); });
  window_.on_browse_file([this]() { file_io_handler_->handleBrowseFile(); });
  window_.on_reload_file([this]() { file_io_handler_->handleReloadFile(); });
  window_.on_write_files([this]() { file_io_handler_->handleWriteFiles(); });
  window_.on_validate_inputs([this]() {
    slint::invoke_from_event_loop(
        [this]() { static_cast<void>(input_validator_->validateInputs()); });
  });

  // Handle calculator toggle: update options and refresh the UI model
  window_.on_toggle_calculator([this](const slint::SharedString &calc_id, bool enabled) {
    options_.active_calculators[std::string(calc_id.data())] = enabled;
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
  window_.on_apply_min_factor([this](float factor) { handleApplyMinFactor(factor); });
  window_.on_apply_max_factor([this](float factor) { handleApplyMaxFactor(factor); });
  window_.on_apply_global_cutoff([this](float cutoff) { handleApplyGlobalCutoff(cutoff); });
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

AppController::AppController(::AppWindow &window, TrajectoryLoader &loader,
                             AnalysisDispatcher &dispatcher,
                             [[maybe_unused]] BondCutoffService &cutoff_service,
                             ProgramOptions &options)
    : AppController(window, loader, dispatcher, options) {}

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
void AppController::handleOptionstoUI() { OptionsSyncService::writeToUI(window_, options_); }

void AppController::updateActiveGroupFlags() {
  OptionsSyncService::updateActiveGroupFlags(window_, options_);
}

ProgramOptions AppController::handleOptionsfromUI() {
  const auto cutoffs = bond_cutoff_controller_.getBondCutoffs();
  auto opt_expected = OptionsSyncService::readFromUI(window_, loader_.getFrameCount(), cutoffs);
  if (opt_expected) {
    return *opt_expected;
  }
  ProgramOptions default_opts;
  default_opts.bond_cutoffs = cutoffs;
  return default_opts;
}

void AppController::setBondCutoffs() { bond_cutoff_controller_.setBondCutoffs(); }

correlation::analysis::BondCutoffMatrix AppController::getBondCutoffs() {
  return bond_cutoff_controller_.getBondCutoffs();
}

void AppController::populateCalculatorGroups() {
  const auto &calculators =
      ::correlation::calculators::CalculatorFactory::instance().getCalculators();
  const auto &opts = options_;

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
  OptionsResetService::resetRDF(window_);
  static_cast<void>(input_validator_->validateInputs());
}

void AppController::handleResetAngleOptions() {
  OptionsResetService::resetAngle(window_);
  static_cast<void>(input_validator_->validateInputs());
}

void AppController::handleResetSQOptions() {
  OptionsResetService::resetSQ(window_);
  static_cast<void>(input_validator_->validateInputs());
}

void AppController::handleResetRingsOptions() {
  OptionsResetService::resetRings(window_);
  static_cast<void>(input_validator_->validateInputs());
}

void AppController::handleResetSmoothingOptions() {
  OptionsResetService::resetSmoothing(window_);
  static_cast<void>(input_validator_->validateInputs());
}

void AppController::handleResetAdvancedOptions() {
  OptionsResetService::resetAdvanced(window_);
  static_cast<void>(input_validator_->validateInputs());
}

void AppController::handleResetTrajectoryOptions() {
  OptionsResetService::resetTrajectory(window_, loader_);
  static_cast<void>(input_validator_->validateInputs());
}

void AppController::handleResetExportSettings() {
  OptionsResetService::resetExportSettings(window_);
  if (plot_controller_) {
    plot_controller_->requestPlotUpdate(window_.get_selected_plot_index(), true);
  }
}

void AppController::handleResetAnalysesSelection() {
  populateCalculatorGroups();
  updateActiveGroupFlags();
}

void AppController::handleResetMaterialType() {
  OptionsResetService::resetMaterialType(window_);
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
  OptionsResetService::resetXRD(window_);
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
  bond_cutoff_controller_.applyScaledCutoffs(scale_factor);
}

void AppController::handleSetUniformCutoff(float max_cutoff) {
  bond_cutoff_controller_.setUniformCutoff(max_cutoff);
}

void AppController::handleApplyMinFactor(float min_factor) {
  bond_cutoff_controller_.applyCovalentFactor(min_factor, FactorBound::Min);
}

void AppController::handleApplyMaxFactor(float max_factor) {
  bond_cutoff_controller_.applyCovalentFactor(max_factor, FactorBound::Max);
}

void AppController::handleApplyGlobalCutoff(float global_cutoff) {
  bond_cutoff_controller_.applyGlobalCutoff(global_cutoff);
}

} // namespace correlation::app
