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

#include "app/services/FileIOHandler.hpp"
#include "app/services/InputValidator.hpp"
#include "app/services/OptionsResetService.hpp"
#include "app/services/OptionsSyncService.hpp"
#include "app/services/UpdateChecker.hpp"
#include "app/viewmodel/BondCutoffController.hpp"
#include "app/viewmodel/PlotController.hpp"
#include "app/viewmodel/PresetController.hpp"
#include "app/viewmodel/Structure3DController.hpp"
#include "calculators/CalculatorFactory.hpp"

#include <cstdlib>
#include <map>
#include <string>
#include <vector>

namespace correlation::app {

AppController::AppController(::AppWindow &window, TrajectoryLoader &loader,
                             AnalysisDispatcher &dispatcher, ProgramOptions &options)
    : window_(window), loader_(loader), dispatcher_(dispatcher), options_(options),
      bond_cutoff_controller_(window, loader, options, &dispatcher) {
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
  structure_3d_controller_ = std::make_unique<Structure3DController>(window_, options_);

  bond_cutoff_controller_.setPlotController(plot_controller_.get());
  plot_controller_->setOnPlotGeometryChanged([this](float width, float height) {
    bond_cutoff_controller_.updateCalibrationPixels(width, height);
  });

  // default options to UI
  handleOptionstoUI();

  // Connect the UI signals to the controller's member functions.
  registerActionCallbacks();
  registerPlotCallbacks();
  registerStructure3DCallbacks();
  registerPresetAndOptionCallbacks();
  registerExportCallbacks();

  // Handle open external URL (e.g. download update from browser)
  window_.on_open_url([](const slint::SharedString &url) {
    UpdateChecker::openUrlInBrowser(std::string(url.data()));
  });

  // Initial load of settings and preset list
  loadSettings();
  preset_controller_->refreshPresetList();
  preset_controller_->initDomainPresets();
  bond_cutoff_controller_.evaluateSanityState();

  // Asynchronously check for GitHub releases in the background without interrupting the user
  UpdateChecker::checkForUpdatesAsync(window_, CORRELATION_VERSION_STRING);
}

void AppController::registerActionCallbacks() {
  window_.on_run_analysis([this] { analysis_runner_->handleRunAnalysis(); });
  window_.on_cancel_analysis([this] { dispatcher_.cancelAnalysis(); });
  window_.on_browse_file([this] { file_io_handler_->handleBrowseFile(); });
  window_.on_reload_file([this] { file_io_handler_->handleReloadFile(); });
  window_.on_write_files([this] { file_io_handler_->handleWriteFiles(); });
  window_.on_validate_inputs([this] {
    slint::invoke_from_event_loop(
        [this] { static_cast<void>(input_validator_->validateInputs()); });
  });

  window_.on_toggle_calculator([this](const slint::SharedString &calc_id, bool enabled) {
    options_.active_calculators[std::string(calc_id.data())] = enabled;
    populateCalculatorGroups();
    updateActiveGroupFlags();
  });

  window_.on_layout_geometry_changed([this](float left_w, float mid_w) {
    settings_.left_col_width = left_w;
    settings_.middle_col_width = mid_w;
    saveSettings();
  });
}

void AppController::registerPlotCallbacks() {
  window_.on_select_plot([this](int index) { plot_controller_->requestPlotUpdate(index, true); });

  window_.on_toggle_curve_visibility([this](int curve_id, bool visible) {
    plot_controller_->handleToggleCurveVisibility(curve_id, visible);
  });

  window_.on_toggle_all_curves(
      [this](bool visible) { plot_controller_->handleToggleAllCurves(visible); });

  window_.on_set_curve_color([this](int curve_id, const slint::SharedString &color_hex) {
    plot_controller_->handleSetCurveColor(curve_id, color_hex);
  });

  window_.on_mouse_move(
      [this](float mouse_x, float mouse_y, bool hover, float width, float height) {
        plot_controller_->handleMouseMove(mouse_x, mouse_y, hover, width, height);
      });

  window_.on_save_plot([this] { plot_controller_->handleSavePlot(); });
  window_.on_pin_run([this] { plot_controller_->handlePinRun(); });
  window_.on_clear_pinned_runs([this] { plot_controller_->handleClearPinnedRuns(); });
  window_.on_toggle_difference_plot(
      [this](bool show) { plot_controller_->handleToggleDifferencePlot(show); });

  window_.on_plot_resized([this](float width, float height) {
    plot_controller_->handlePlotResized({
        .width = width,
        .height = height,
    });
  });

  window_.on_zoom_to_rect([this](float start_x, float start_y, float end_x, float end_y) {
    plot_controller_->handleZoomRect(start_x, start_y, end_x, end_y);
  });
  window_.on_plot_pan(
      [this](float delta_x, float delta_y) { plot_controller_->handlePan(delta_x, delta_y); });
  window_.on_plot_wheel_zoom([this](float mouse_x, float mouse_y, float delta) {
    plot_controller_->handleWheelZoom(mouse_x, mouse_y, delta);
  });
  window_.on_reset_zoom([this] { plot_controller_->handleResetZoom(); });
  window_.on_add_marker_at([this](float click_x, float click_y, bool is_vert) {
    plot_controller_->handleAddMarkerAt(click_x, click_y, is_vert);
  });
  window_.on_add_marker([this](bool is_vert) { plot_controller_->handleAddMarker(is_vert); });
  window_.on_clear_markers([this] { plot_controller_->handleClearMarkers(); });

  window_.on_select_calibration_pair([this](int idx) {
    window_.set_active_calibration_pair_index(idx);
    bond_cutoff_controller_.updateCalibrationPixels(plot_controller_->getLastPlotWidth(),
                                                    plot_controller_->getLastPlotHeight());
  });

  window_.on_auto_detect_cutoff(
      [this](int idx) { bond_cutoff_controller_.autoDetectCutoffs(idx); });

  window_.on_reset_calibration_pair(
      [this](int idx) { bond_cutoff_controller_.resetCalibrationPair(idx); });

  window_.on_drag_cutoff_handle([this](int idx, bool is_max, float mouse_x) {
    bond_cutoff_controller_.dragCutoffHandle(idx, is_max, mouse_x);
  });
}

void AppController::registerStructure3DCallbacks() {
  if (!structure_3d_controller_) {
    return;
  }

  window_.on_rotate_3d([this](float delta_x, float delta_y) {
    structure_3d_controller_->handleRotate(delta_x, delta_y);
  });

  window_.on_pan_3d([this](float delta_x, float delta_y) {
    structure_3d_controller_->handlePan(delta_x, delta_y);
  });

  window_.on_zoom_3d([this](float delta) { structure_3d_controller_->handleZoom(delta); });

  window_.on_reset_view_3d([this] { structure_3d_controller_->handleResetView(); });

  window_.on_toggle_bonds_3d(
      [this](bool enabled) { structure_3d_controller_->handleToggleBonds(enabled); });

  window_.on_toggle_wireframe_3d(
      [this](bool enabled) { structure_3d_controller_->handleToggleWireframe(enabled); });

  window_.on_toggle_in_cell_3d(
      [this](bool enabled) { structure_3d_controller_->handleToggleInCell(enabled); });

  window_.on_seek_frame_3d(
      [this](int frame_idx) { structure_3d_controller_->handleSeekFrame(frame_idx); });

  window_.on_step_frame_3d(
      [this](bool forward) { structure_3d_controller_->handleStepFrame(forward); });

  window_.on_viewport_3d_resized([this](float width, float height) {
    structure_3d_controller_->handleViewportResized(width, height);
  });

  window_.on_preview_tab_changed([this](int tab) {
    const bool is_3d = (tab == 2);
    structure_3d_controller_->set3DModeActive(is_3d);
  });

  window_.on_toggle_phase_search_3d(
      [this](bool enabled) { structure_3d_controller_->setPhaseSearchEnabled(enabled); });

  window_.on_set_phase_preset_3d([this](int preset_idx) {
    structure_3d_controller_->setCrystalPreset(static_cast<CrystalPreset>(preset_idx));
  });

  window_.on_set_phase_descriptor_3d([this](int desc_idx) {
    structure_3d_controller_->setPhaseDescriptor(static_cast<PhaseDescriptor>(desc_idx));
  });

  window_.on_set_phase_colormap_3d([this](int cmap_idx) {
    structure_3d_controller_->setPhaseColormap(static_cast<plotters::ColormapType>(cmap_idx));
  });

  window_.on_set_phase_thresholds_3d([this](float min_val, float max_val) {
    structure_3d_controller_->setPhaseThresholds(min_val, max_val);
  });

  window_.on_set_phase_intensity_3d(
      [this](float intensity) { structure_3d_controller_->setPhaseHaloIntensity(intensity); });

  window_.on_set_phase_radius_ratio_3d(
      [this](float ratio) { structure_3d_controller_->setPhaseHaloRadiusRatio(ratio); });

  window_.on_set_phase_dim_non_matching_3d(
      [this](bool dim) { structure_3d_controller_->setPhaseDimNonMatching(dim); });

  window_.on_set_phase_hide_non_matching_3d(
      [this](bool hide) { structure_3d_controller_->setPhaseHideNonMatching(hide); });

  window_.on_set_phase_2d_mode_3d(
      [this](bool is_2d) { structure_3d_controller_->setPhase2DMode(is_2d); });

  window_.on_set_phase_descriptor_y_3d([this](int desc_idx) {
    structure_3d_controller_->setPhaseDescriptorY(static_cast<PhaseDescriptor>(desc_idx));
  });

  window_.on_set_phase_thresholds_2d_3d([this](float min_x, float max_x, float min_y, float max_y) {
    structure_3d_controller_->setPhaseThresholds2D(min_x, max_x, min_y, max_y);
  });

  window_.on_export_phase_heatmap_png_3d([this]() {
    // Generate default export filename in current directory or user pictures
    structure_3d_controller_->exportHeatmapPng("phase_heatmap_2d.png");
  });
}

void AppController::registerPresetAndOptionCallbacks() {
  window_.on_load_preset([this](int index) {
    slint::invoke_from_event_loop([this, index] { preset_controller_->handleLoadPreset(index); });
  });
  window_.on_save_preset([this](const slint::SharedString &name) {
    slint::invoke_from_event_loop(
        [this, sname = std::string(name.data())] { preset_controller_->handleSavePreset(sname); });
  });
  window_.on_delete_preset([this](int index) {
    slint::invoke_from_event_loop([this, index] { preset_controller_->handleDeletePreset(index); });
  });
  window_.on_material_type_changed([this](int type) {
    slint::invoke_from_event_loop(
        [this, type] { preset_controller_->handleMaterialTypeChanged(type); });
  });

  window_.on_apply_domain_preset([this](int preset_id) {
    slint::invoke_from_event_loop([this, preset_id] {
      preset_controller_->handleApplyDomainPreset(preset_id);
      bond_cutoff_controller_.evaluateSanityState();
    });
  });

  window_.on_fix_sanity_issue(
      [this](int action) { bond_cutoff_controller_.handleFixSanityIssue(action); });

  window_.on_reset_bond_cutoffs([this] {
    setBondCutoffs();
    if (structure_3d_controller_) {
      structure_3d_controller_->updateBondCutoffScale(1.15);
    }
  });
  window_.on_reset_rdf_options([this] {
    OptionsResetService::resetRDF(window_);
    static_cast<void>(input_validator_->validateInputs());
  });
  window_.on_reset_angle_options([this] {
    OptionsResetService::resetAngle(window_);
    static_cast<void>(input_validator_->validateInputs());
  });
  window_.on_reset_sq_options([this] {
    OptionsResetService::resetSQ(window_);
    static_cast<void>(input_validator_->validateInputs());
  });
  window_.on_reset_xrd_options([this] {
    OptionsResetService::resetXRD(window_);
    static_cast<void>(input_validator_->validateInputs());
  });
  window_.on_xrd_radiation_preset_changed([this](int idx) {
    OptionsResetService::setXRDPreset(window_, idx);
    static_cast<void>(input_validator_->validateInputs());
  });
  window_.on_apply_scaled_cutoffs([this](float factor) {
    bond_cutoff_controller_.applyScaledCutoffs(factor);
    if (structure_3d_controller_ && factor > 0.0F) {
      structure_3d_controller_->updateBondCutoffScale(static_cast<double>(factor));
    }
  });
  window_.on_set_uniform_cutoff(
      [this](float max_val) { bond_cutoff_controller_.setUniformCutoff(max_val); });
  window_.on_apply_min_factor([this](float factor) {
    bond_cutoff_controller_.applyCovalentFactor(factor, FactorBound::Min);
  });
  window_.on_apply_max_factor([this](float factor) {
    bond_cutoff_controller_.applyCovalentFactor(factor, FactorBound::Max);
  });
  window_.on_apply_global_cutoff(
      [this](float cutoff) { bond_cutoff_controller_.applyGlobalCutoff(cutoff); });
  window_.on_reset_rings_options([this] {
    OptionsResetService::resetRings(window_);
    static_cast<void>(input_validator_->validateInputs());
  });
  window_.on_reset_smoothing_options([this] {
    OptionsResetService::resetSmoothing(window_);
    static_cast<void>(input_validator_->validateInputs());
  });
  window_.on_reset_advanced_options([this] {
    OptionsResetService::resetAdvanced(window_);
    static_cast<void>(input_validator_->validateInputs());
  });
  window_.on_reset_trajectory_options([this] {
    OptionsResetService::resetTrajectory(window_, loader_);
    static_cast<void>(input_validator_->validateInputs());
  });
  window_.on_reset_export_settings([this] {
    OptionsResetService::resetExportSettings(window_);
    if (plot_controller_) {
      plot_controller_->requestPlotUpdate(window_.get_selected_plot_index(), true);
    }
  });
  window_.on_reset_material_type([this] {
    OptionsResetService::resetMaterialType(window_);
    if (preset_controller_) {
      preset_controller_->handleMaterialTypeChanged(0);
    }
  });
  window_.on_reset_analyses_selection([this] {
    populateCalculatorGroups();
    updateActiveGroupFlags();
  });
  window_.on_reset_presets([this] { window_.set_selected_preset(-1); });
  window_.on_clear_comparison_curves([this] {
    const auto empty_model = std::make_shared<slint::VectorModel<ComparisonCurveData>>();
    window_.set_comparison_curves(empty_model);
    if (plot_controller_) {
      plot_controller_->handleClearPinnedRuns();
    }
  });
}

void AppController::registerExportCallbacks() {
  window_.on_toggle_export_analysis([this](const slint::SharedString &analysis_id, bool enabled) {
    const auto model = window_.get_export_analysis_items();
    if (!model) {
      return;
    }
    for (size_t i = 0; i < model->row_count(); ++i) {
      auto maybe_item = model->row_data(i);
      if (maybe_item && maybe_item->id == analysis_id) {
        auto item = *maybe_item;
        item.enabled = enabled;
        model->set_row_data(i, item);
        break;
      }
    }
  });

  window_.on_toggle_all_export_analyses([this](bool enabled) {
    const auto model = window_.get_export_analysis_items();
    if (!model) {
      return;
    }
    for (size_t i = 0; i < model->row_count(); ++i) {
      auto maybe_item = model->row_data(i);
      if (maybe_item) {
        auto item = *maybe_item;
        item.enabled = enabled;
        model->set_row_data(i, item);
      }
    }
  });

  window_.on_reset_export_options([this] {
    ExportFormatOptions opts{};
    opts.export_csv = true;
    opts.export_images = true;
    opts.export_hdf5 = false;
    opts.export_arrow = false;
    opts.use_zip = false;
    window_.set_export_format_options(opts);
    window_.set_export_all_analyses(true);
    const auto model = window_.get_export_analysis_items();
    if (model) {
      for (size_t i = 0; i < model->row_count(); ++i) {
        auto maybe_item = model->row_data(i);
        if (maybe_item) {
          auto item = *maybe_item;
          item.enabled = true;
          model->set_row_data(i, item);
        }
      }
    }
  });
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

void AppController::loadTrajectory(const std::string &filepath,
                                   std::function<void(bool)> on_complete) {
  if (file_io_handler_) {
    file_io_handler_->startLoadingTrajectory(filepath, std::move(on_complete));
  }
}

void AppController::runAnalysis() {
  if (analysis_runner_) {
    analysis_runner_->handleRunAnalysis();
  }
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

namespace {
bool isDefaultDisabledCalculator(const correlation::calculators::BaseCalculator *calc) noexcept {
  if (calc == nullptr) {
    return false;
  }
  const auto name = calc->getName();
  const auto short_name = calc->getShortName();
  if (short_name == "RD" || name == "Ring Size Distribution") {
    return true;
  }
  if (short_name == "DAD" || name == "Dihedral Angle Distribution") {
    return true;
  }
  if (name.contains("(GPU)") || short_name.ends_with("_GPU")) {
    return true;
  }
  return false;
}
} // namespace

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

  const auto groups_model = std::make_shared<slint::VectorModel<CalculatorGroup>>();
  for (const auto &grp_name : group_order) {
    const auto calcs_model = std::make_shared<slint::VectorModel<CalculatorInfo>>();
    for (const auto *calc : groups_map.at(grp_name)) {
      bool enabled = !isDefaultDisabledCalculator(calc);
      if (const auto calc_it = opts.active_calculators.find(std::string(calc->getName()));
          calc_it != opts.active_calculators.end()) {
        enabled = calc_it->second;
      } else if (const auto short_it =
                     opts.active_calculators.find(std::string(calc->getShortName()));
                 short_it != opts.active_calculators.end()) {
        enabled = short_it->second;
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

void AppController::populateExportAnalyses() {
  const auto names = dispatcher_.getAvailableHistogramNames();
  const auto model = std::make_shared<slint::VectorModel<ExportAnalysisItem>>();
  for (const auto &name : names) {
    if (name.ends_with("_raw")) {
      continue;
    }
    const auto *hist = dispatcher_.getHistogram(name);
    const std::string label = (hist != nullptr && !hist->title.empty()) ? hist->title : name;
    ExportAnalysisItem item;
    item.id = slint::SharedString(name);
    item.name = slint::SharedString(label);
    item.enabled = true;
    model->push_back(item);
  }
  window_.set_export_analysis_items(model);
  window_.set_export_all_analyses(true);
}

} // namespace correlation::app
