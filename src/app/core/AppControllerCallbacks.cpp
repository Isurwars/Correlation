/**
 * @file AppControllerCallbacks.cpp
 * @brief UI Action, plot, analysis, and preset callback registrations for AppController.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "AppWindow.h"
#include "app/core/AppController.hpp"
#include "app/services/AnalysisRunner.hpp"
#include "app/services/FileIOHandler.hpp"
#include "app/services/InputValidator.hpp"
#include "app/services/OptionsResetService.hpp"
#include "app/viewmodel/BondCutoffController.hpp"
#include "app/viewmodel/PlotController.hpp"
#include "app/viewmodel/PresetController.hpp"
#include "app/viewmodel/Structure3DController.hpp"
#include "calculators/CalculatorFactory.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace correlation::app {

namespace {

struct CommandSpec {
  std::string_view id;
  std::string_view title;
  std::string_view category;
  std::string_view shortcut;
};

constexpr std::array<CommandSpec, 18> COMMAND_CATALOG{{
    {.id = "open_file",
     .title = "Open Structure / Trajectory File",
     .category = "File",
     .shortcut = "Ctrl+O"},
    {.id = "run_analysis",
     .title = "Run Enabled Analyses",
     .category = "Analysis",
     .shortcut = "Ctrl+R"},
    {.id = "cancel_analysis",
     .title = "Cancel Running Analysis",
     .category = "Analysis",
     .shortcut = "Esc"},
    {.id = "export_data",
     .title = "Export All Data Files",
     .category = "Export",
     .shortcut = "Ctrl+S"},
    {.id = "export_plot",
     .title = "Export Current Plot Image",
     .category = "Export",
     .shortcut = "Ctrl+Shift+S"},
    {.id = "export_structure",
     .title = "Export Structure File",
     .category = "Export",
     .shortcut = "Ctrl+E"},
    {.id = "copy_svg",
     .title = "Copy Plot SVG to Clipboard",
     .category = "Plot",
     .shortcut = "Ctrl+C"},
    {.id = "view_2d", .title = "Switch to 2D Plot View", .category = "View", .shortcut = "1"},
    {.id = "view_3d", .title = "Switch to 3D Structure View", .category = "View", .shortcut = "2"},
    {.id = "view_dual", .title = "Switch to Dual Split View", .category = "View", .shortcut = "3"},
    {.id = "view_data",
     .title = "Switch to Data Explorer View",
     .category = "View",
     .shortcut = "4"},
    {.id = "reset_view",
     .title = "Reset Plot Zoom & 3D Camera",
     .category = "View",
     .shortcut = "H"},
    {.id = "toggle_zen",
     .title = "Toggle Zen Mode (Hide Sidebars)",
     .category = "Workspace",
     .shortcut = "Ctrl+\\"},
    {.id = "toggle_theme",
     .title = "Toggle Dark / Light Theme",
     .category = "Workspace",
     .shortcut = "Ctrl+D"},
    {.id = "show_shortcuts",
     .title = "Show Keyboard Shortcuts",
     .category = "Help",
     .shortcut = "?"},
    {.id = "sample_pdsi",
     .title = "Load Sample: a-PdSi Metallic Glass",
     .category = "Sample",
     .shortcut = ""},
    {.id = "sample_cellulose",
     .title = "Load Sample: Cellulose Polymer",
     .category = "Sample",
     .shortcut = ""},
    {.id = "sample_si",
     .title = "Load Sample: Silicon Crystal",
     .category = "Sample",
     .shortcut = ""},
}};

[[nodiscard]] bool matchesQuery(std::string_view haystack, std::string_view needle) {
  if (needle.empty()) {
    return true;
  }
  const auto match =
      std::ranges::search(haystack, needle, [](unsigned char lhs, unsigned char rhs) {
        return std::tolower(lhs) == std::tolower(rhs);
      });
  return !match.empty();
}

[[nodiscard]] bool isCommandEnabled(std::string_view cmd_id, const AppWindow &window) {
  if (cmd_id == "run_analysis") {
    return window.get_file_loaded() && !window.get_has_validation_errors() &&
           !window.get_analysis_running();
  }
  if (cmd_id == "cancel_analysis") {
    return window.get_analysis_running();
  }
  if (cmd_id == "export_data") {
    return window.get_analysis_done();
  }
  if (cmd_id == "export_plot" || cmd_id == "copy_svg") {
    return window.get_selected_plot_index() >= 0;
  }
  if (cmd_id == "export_structure") {
    return window.get_file_loaded();
  }
  return true;
}

[[nodiscard]] std::shared_ptr<slint::VectorModel<CommandItem>>
buildFilteredCommands(const AppWindow &window, std::string_view query) {
  std::vector<CommandItem> items;
  items.reserve(COMMAND_CATALOG.size());
  for (const auto &spec : COMMAND_CATALOG) {
    if (!matchesQuery(spec.title, query) && !matchesQuery(spec.category, query)) {
      continue;
    }
    CommandItem item{};
    item.id = slint::SharedString(spec.id);
    item.title = slint::SharedString(spec.title);
    item.category = slint::SharedString(spec.category);
    item.shortcut = slint::SharedString(spec.shortcut);
    item.enabled = isCommandEnabled(spec.id, window);
    items.push_back(std::move(item));
  }
  return std::make_shared<slint::VectorModel<CommandItem>>(std::move(items));
}

[[nodiscard]] std::string resolveSampleRelativePath(std::string_view sample_name) {
  if (sample_name == "a-PdSi") {
    return "examples/a-PdSi/a-PdSi.car";
  }
  if (sample_name == "Cellulose") {
    return "examples/Cellulose/Cellulose.md";
  }
  if (sample_name == "Si") {
    return "examples/Si/Si.car";
  }
  return {};
}

[[nodiscard]] std::string resolveSampleFilePath(std::string_view sample_name) {
  const std::string rel = resolveSampleRelativePath(sample_name);
  if (rel.empty()) {
    return {};
  }
  std::error_code ec;
  for (const std::string_view prefix : {"", "../", "../../", "../../../"}) {
    const std::filesystem::path candidate = std::filesystem::path(prefix) / rel;
    if (std::filesystem::exists(candidate, ec) && !ec) {
      return candidate.string();
    }
  }
#ifndef _WIN32
  const auto exe_path = std::filesystem::read_symlink("/proc/self/exe", ec);
  if (!ec && !exe_path.empty()) {
    const auto exe_dir = exe_path.parent_path();
    for (const std::string_view prefix : {"", "../", "../../", "../share/correlation/"}) {
      const std::filesystem::path candidate = exe_dir / prefix / rel;
      if (std::filesystem::exists(candidate, ec) && !ec) {
        return candidate.string();
      }
    }
  }
#endif
  return {};
}

} // namespace

void AppController::registerActionCallbacks() {
  window_.set_command_items(buildFilteredCommands(window_, ""));
  window_.on_filter_commands([this](const slint::SharedString &query) {
    window_.set_command_items(buildFilteredCommands(window_, std::string_view(query.data())));
  });
  window_.on_run_analysis([this] { analysis_runner_->handleRunAnalysis(); });
  window_.on_cancel_analysis([this] { dispatcher_.cancelAnalysis(); });
  window_.on_browse_file([this] { file_io_handler_->handleBrowseFile(); });
  window_.on_load_sample([this](const slint::SharedString &sample_name) {
    const std::string_view name(sample_name.data());
    const std::string resolved = resolveSampleFilePath(name);
    if (!resolved.empty()) {
      window_.invoke_show_toast(slint::SharedString("Loading " + std::string(name) + " sample..."),
                                0);
      loadTrajectory(resolved);
    } else {
      const std::string msg = "Sample dataset '" + std::string(name) + "' not found on disk.";
      window_.set_file_status_text(slint::SharedString(msg));
      window_.invoke_show_toast(slint::SharedString(msg), 3);
    }
  });
  window_.on_load_recent_file([this](const slint::SharedString &recent_path) {
    const std::string path(recent_path.data());
    std::error_code ec;
    if (!path.empty() && std::filesystem::exists(path, ec) && !ec) {
      loadTrajectory(path);
    } else {
      const std::string fname = std::filesystem::path(path).filename().string();
      const std::string msg = "Recent file not found: " + (fname.empty() ? path : fname);
      window_.set_file_status_text(slint::SharedString(msg));
      window_.invoke_show_toast(slint::SharedString(msg), 3);
    }
  });
  window_.on_reload_file([this] { file_io_handler_->handleReloadFile(); });
  window_.on_rename_element([this](const slint::SharedString &old_sym,
                                   const slint::SharedString &new_sym) {
    file_io_handler_->handleRenameElement(std::string(old_sym.data()), std::string(new_sym.data()));
  });
  window_.on_write_files([this] { file_io_handler_->handleWriteFiles(); });
  window_.on_export_structure([this] { file_io_handler_->handleExportStructure(); });
  window_.on_validate_inputs([this] {
    slint::invoke_from_event_loop(
        [this] { static_cast<void>(input_validator_->validateInputs()); });
  });

  window_.on_toggle_calculator([this](const slint::SharedString &calc_id, bool enabled) {
    const std::string id_str(calc_id.data());
    options_.active_calculators[id_str] = enabled;
    const auto *calc =
        ::correlation::calculators::CalculatorFactory::instance().getCalculator(id_str);
    if (calc != nullptr) {
      options_.active_calculators[std::string(calc->getName())] = enabled;
      if (!calc->getShortName().empty()) {
        options_.active_calculators[std::string(calc->getShortName())] = enabled;
      }
      const bool is_vdos =
          (calc->getShortName() == "VDOS" || calc->getName() == "Vibrational Density of States");
      const bool is_vacf = (calc->getShortName() == "VACF" ||
                            calc->getName() == "Velocity Autocorrelation Function");
      if (is_vdos && enabled) {
        options_.active_calculators["Velocity Autocorrelation Function"] = true;
        options_.active_calculators["VACF"] = true;
      } else if (is_vacf && !enabled) {
        options_.active_calculators["Vibrational Density of States"] = false;
        options_.active_calculators["VDOS"] = false;
      }
    }
    populateCalculatorGroups();
    updateActiveGroupFlags();
    markResultsStaleIfDone();
  });

  registerAnalysesCallbacks();

  window_.on_layout_geometry_changed([this](float left_w, float mid_w) {
    settings_.left_col_width = left_w;
    settings_.middle_col_width = mid_w;
    saveSettings();
  });
  window_.on_right_layout_changed([this](float right_w) {
    settings_.right_col_width = std::clamp(right_w, 200.0F, 520.0F);
    saveSettings();
  });
  window_.on_dual_split_changed([this](float ratio) {
    settings_.dual_split_ratio = std::clamp(ratio, 0.15F, 0.85F);
    saveSettings();
  });
  window_.on_theme_changed([this](bool dark) {
    settings_.dark_theme = dark;
    saveSettings();
  });
}

void AppController::registerAnalysesCallbacks() {
  window_.on_toggle_calculator_group([this](const slint::SharedString &group_name, bool enabled) {
    setCalculatorsEnabled(enabled, std::string_view(group_name.data()));
    markResultsStaleIfDone();
  });

  window_.on_toggle_all_calculators([this](bool enabled) {
    setCalculatorsEnabled(enabled);
    markResultsStaleIfDone();
  });

  window_.on_set_calculator_group_expanded(
      [this](const slint::SharedString &group_name, bool expanded) {
        setCalculatorGroupsExpanded(expanded, std::string_view(group_name.data()));
      });

  window_.on_set_all_calculator_groups_expanded(
      [this](bool expanded) { setCalculatorGroupsExpanded(expanded); });
}

void AppController::registerPlotCallbacks() {
  window_.on_select_plot([this](int index) { plot_controller_->requestPlotUpdate(index, true); });

  window_.on_toggle_curve_visibility([this](int curve_id, bool visible) {
    plot_controller_->handleToggleCurveVisibility(curve_id, visible);
  });

  window_.on_toggle_all_curves(
      [this](bool visible) { plot_controller_->handleToggleAllCurves(visible); });

  window_.on_filter_curves(
      [this](const slint::SharedString &filter) { plot_controller_->handleFilterCurves(filter); });

  window_.on_invert_curves([this] { plot_controller_->handleInvertCurveVisibility(); });

  window_.on_set_curve_color([this](int curve_id, const slint::SharedString &color_hex) {
    plot_controller_->handleSetCurveColor(curve_id, color_hex);
  });

  window_.on_mouse_move(
      [this](float mouse_x, float mouse_y, bool hover, float width, float height) {
        plot_controller_->handleMouseMove(mouse_x, mouse_y, hover, width, height);
      });

  window_.on_save_plot([this] { plot_controller_->handleSavePlot(); });
  window_.on_copy_plot_clipboard([this] { plot_controller_->handleCopyPlotToClipboard(); });
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

  window_.on_tag_peak_at(
      [this](float peak_x, float peak_y) { plot_controller_->handleTagPeakAt(peak_x, peak_y); });
  window_.on_tag_peak_range(
      [this](float x_start, float x_end) { plot_controller_->handleTagPeakRange(x_start, x_end); });
  window_.on_remove_tagged_peak([this](int idx) { plot_controller_->handleRemoveTaggedPeak(idx); });
  window_.on_clear_tagged_peaks([this] { plot_controller_->handleClearTaggedPeaks(); });
  window_.on_export_tagged_peaks([this] { plot_controller_->handleExportTaggedPeaks(); });

  window_.on_apply_axis_settings([this] { plot_controller_->handleApplyAxisSettings(); });
  window_.on_reset_axis_settings([this] { plot_controller_->handleResetAxisSettings(); });
  window_.on_set_curve_width([this](const slint::SharedString &shared_str) {
    plot_controller_->handleSetCurveWidth(shared_str);
  });

  window_.on_select_calibration_pair([this](int idx) {
    window_.set_active_calibration_pair_index(idx);
    bond_cutoff_controller_.updateCalibrationPixels(plot_controller_->getLastPlotWidth(),
                                                    plot_controller_->getLastPlotHeight());
  });

  window_.on_auto_detect_cutoff([this](int idx) {
    bond_cutoff_controller_.autoDetectCutoffs(idx);
    markResultsStaleIfDone();
  });

  window_.on_reset_calibration_pair([this](int idx) {
    bond_cutoff_controller_.resetCalibrationPair(idx);
    markResultsStaleIfDone();
  });

  window_.on_drag_cutoff_handle([this](int idx, bool is_max, float mouse_x) {
    bond_cutoff_controller_.dragCutoffHandle(idx, is_max, mouse_x);
    markResultsStaleIfDone();
  });
}

void AppController::registerPresetAndOptionCallbacks() {
  window_.on_load_preset([this](int index) {
    slint::invoke_from_event_loop([this, index] {
      preset_controller_->handleLoadPreset(index);
      markResultsStaleIfDone();
    });
  });
  window_.on_save_preset([this](const slint::SharedString &name) {
    slint::invoke_from_event_loop(
        [this, sname = std::string(name.data())] { preset_controller_->handleSavePreset(sname); });
  });
  window_.on_delete_preset([this](int index) {
    slint::invoke_from_event_loop([this, index] { preset_controller_->handleDeletePreset(index); });
  });
  window_.on_material_type_changed([this](int type) {
    slint::invoke_from_event_loop([this, type] {
      preset_controller_->handleMaterialTypeChanged(type);
      markResultsStaleIfDone();
    });
  });

  window_.on_apply_domain_preset([this](int preset_id) {
    slint::invoke_from_event_loop([this, preset_id] {
      preset_controller_->handleApplyDomainPreset(preset_id);
      bond_cutoff_controller_.evaluateSanityState();
      markResultsStaleIfDone();
    });
  });

  window_.on_fix_sanity_issue([this](int action) {
    bond_cutoff_controller_.handleFixSanityIssue(action);
    markResultsStaleIfDone();
  });

  window_.on_reset_bond_cutoffs([this] {
    setBondCutoffs();
    if (structure_3d_controller_) {
      structure_3d_controller_->updateBondCutoffScale(1.15);
    }
    markResultsStaleIfDone();
  });
  window_.on_reset_rdf_options([this] {
    OptionsResetService::resetRDF(window_);
    static_cast<void>(input_validator_->validateInputs());
    markResultsStaleIfDone();
  });
  window_.on_reset_angle_options([this] {
    OptionsResetService::resetAngle(window_);
    static_cast<void>(input_validator_->validateInputs());
    markResultsStaleIfDone();
  });
  window_.on_reset_sq_options([this] {
    OptionsResetService::resetSQ(window_);
    static_cast<void>(input_validator_->validateInputs());
    markResultsStaleIfDone();
  });
  window_.on_reset_xrd_options([this] {
    OptionsResetService::resetXRD(window_);
    static_cast<void>(input_validator_->validateInputs());
    markResultsStaleIfDone();
  });
  window_.on_xrd_radiation_preset_changed([this](int idx) {
    OptionsResetService::setXRDPreset(window_, idx);
    static_cast<void>(input_validator_->validateInputs());
    markResultsStaleIfDone();
  });
  window_.on_apply_scaled_cutoffs([this](float factor) {
    bond_cutoff_controller_.applyScaledCutoffs(factor);
    if (structure_3d_controller_ && factor > 0.0F) {
      structure_3d_controller_->updateBondCutoffScale(static_cast<double>(factor));
    }
    markResultsStaleIfDone();
  });
  window_.on_set_uniform_cutoff([this](float max_val) {
    bond_cutoff_controller_.setUniformCutoff(max_val);
    markResultsStaleIfDone();
  });
  window_.on_apply_min_factor([this](float factor) {
    bond_cutoff_controller_.applyCovalentFactor(factor, FactorBound::Min);
    markResultsStaleIfDone();
  });
  window_.on_apply_max_factor([this](float factor) {
    bond_cutoff_controller_.applyCovalentFactor(factor, FactorBound::Max);
    markResultsStaleIfDone();
  });
  window_.on_apply_global_cutoff([this](float cutoff) {
    bond_cutoff_controller_.applyGlobalCutoff(cutoff);
    markResultsStaleIfDone();
  });
  window_.on_reset_rings_options([this] {
    OptionsResetService::resetRings(window_);
    static_cast<void>(input_validator_->validateInputs());
    markResultsStaleIfDone();
  });
  window_.on_reset_smoothing_options([this] {
    OptionsResetService::resetSmoothing(window_);
    static_cast<void>(input_validator_->validateInputs());
    markResultsStaleIfDone();
  });
  window_.on_reset_advanced_options([this] {
    OptionsResetService::resetAdvanced(window_);
    static_cast<void>(input_validator_->validateInputs());
    markResultsStaleIfDone();
  });
  window_.on_reset_trajectory_options([this] {
    OptionsResetService::resetTrajectory(window_, loader_);
    static_cast<void>(input_validator_->validateInputs());
    markResultsStaleIfDone();
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
    markResultsStaleIfDone();
  });
  window_.on_reset_analyses_selection([this] {
    options_.active_calculators.clear();
    populateCalculatorGroups();
    updateActiveGroupFlags();
    markResultsStaleIfDone();
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

} // namespace correlation::app
