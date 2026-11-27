/**
 * @file AppController.cpp
 * @brief Core lifecycle and service orchestration for the application controller.
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
#include "app/services/OptionsSyncService.hpp"
#include "app/services/UpdateChecker.hpp"
#include "app/viewmodel/BondCutoffController.hpp"
#include "app/viewmodel/PlotController.hpp"
#include "app/viewmodel/PresetController.hpp"
#include "app/viewmodel/Structure3DController.hpp"

#include <nfd.h>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <span>
#include <string>
#include <utility>

namespace correlation::app {

namespace {

void syncRecentFilesToWindow(::AppWindow &window, std::span<const std::string> files) {
  const auto recent_model = std::make_shared<slint::VectorModel<slint::SharedString>>();
  const auto items_model = std::make_shared<slint::VectorModel<RecentFileItem>>();
  for (const auto &file : files) {
    recent_model->push_back(slint::SharedString(file));
    std::string base = std::filesystem::path(file).filename().string();
    if (base.empty()) {
      base = file;
    }
    items_model->push_back(RecentFileItem{
        .basename = slint::SharedString(base),
        .full_path = slint::SharedString(file),
    });
  }
  window.set_recent_files(recent_model);
  window.set_recent_file_items(items_model);
}

} // namespace

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
  plot_controller_->setTrajectoryLoader(&loader_);
  plot_controller_->setStructure3DController(structure_3d_controller_.get());
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

AppController::AppController(::AppWindow &window, TrajectoryLoader &loader,
                             AnalysisDispatcher &dispatcher,
                             [[maybe_unused]] BondCutoffService &cutoff_service,
                             ProgramOptions &options)
    : AppController(window, loader, dispatcher, options) {}

AppController::~AppController() {
  alive_->store(false, std::memory_order_release);
  if (export_3d_thread_.joinable()) {
    export_3d_thread_.join();
  }
  saveSettings();
  // Quit Native File Dialog
  NFD_Quit();
}

void AppController::loadSettings() {
  settings_ = SettingsManager::load();
  window_.set_left_col_width(settings_.left_col_width);
  window_.set_last_col_width(settings_.left_col_width);
  window_.set_middle_col_width(settings_.middle_col_width);
  window_.set_right_col_width(settings_.right_col_width);
  window_.set_last_right_col_width(settings_.right_col_width);
  window_.set_preview_dual_split_ratio(settings_.dual_split_ratio);
  window_.set_is_dark_theme(settings_.dark_theme);
  if (settings_.window_width > 0 && settings_.window_height > 0) {
    window_.window().set_size(
        slint::PhysicalSize({.width = settings_.window_width, .height = settings_.window_height}));
  }

  syncRecentFilesToWindow(window_, settings_.recent_files);
}

void AppController::addRecentFile(const std::string &filepath) {
  if (filepath.empty()) {
    return;
  }
  auto &recents = settings_.recent_files;
  const auto [first, last] = std::ranges::remove(recents, filepath);
  recents.erase(first, last);
  recents.insert(recents.begin(), filepath);
  if (recents.size() > 5) {
    recents.resize(5);
  }
  saveSettings();
  syncRecentFilesToWindow(window_, recents);
}

void AppController::saveSettings() const {
  AppSettings settings = settings_;
  const float last_w = window_.get_last_col_width();
  settings.left_col_width = (last_w >= 180.0F) ? last_w : window_.get_left_col_width();
  settings.middle_col_width = window_.get_middle_col_width();
  const float last_rw = window_.get_last_right_col_width();
  settings.right_col_width = (last_rw >= 200.0F) ? last_rw : window_.get_right_col_width();
  settings.dual_split_ratio = std::clamp(window_.get_preview_dual_split_ratio(), 0.15F, 0.85F);
  settings.dark_theme = window_.get_is_dark_theme();
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

void AppController::markResultsStaleIfDone() noexcept {
  if (window_.get_analysis_done() && !window_.get_analysis_running()) {
    window_.set_results_stale(true);
  }
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

} // namespace correlation::app
