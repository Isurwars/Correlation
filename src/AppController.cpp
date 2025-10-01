// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#if defined(_WIN32)
#define NOMINMAX
#include <Windows.h>
#endif

#include "../include/AppController.hpp"

#include <filesystem>
#include <string>
#include <vector>

//---------------------------------------------------------------------------//
//------------------------------- Constructors ------------------------------//
//---------------------------------------------------------------------------//

AppController::AppController(AppWindow &ui, AppBackend &backend)
    : ui_(ui), backend_(backend) {

  // default options to UI
  handleOptionstoUI(ui_);

  // Connect the UI signals to the controller's member functions
  ui_.on_run_analysis([this]() { handleRunAnalysis(); });
  ui_.on_browse_file([this]() { handleBrowseFile(); });
  ui_.on_check_file_dialog_status([this]() { handleCheckFileDialogStatus(); });
}

//---------------------------------------------------------------------------//
//--------------------------------- Helpers ---------------------------------//
//---------------------------------------------------------------------------//

void AppController::handleOptionstoUI(AppWindow &ui) {
  ProgramOptions opt;
  ui.set_in_file_text(slint::SharedString(opt.input_file));
  ui.set_normalize(opt.normalize);
  ui.set_smoothing(opt.smoothing);
  ui.set_r_max(slint::SharedString(std::format("{:.2f}", opt.r_max)));
  ui.set_r_bin_width(
      slint::SharedString(std::format("{:.2f}", opt.r_bin_width)));
  ui.set_q_max(slint::SharedString(std::format("{:.2f}", opt.q_max)));
  ui.set_q_bin_width(
      slint::SharedString(std::format("{:.2f}", opt.q_bin_width)));
  ui.set_r_int_max(slint::SharedString(std::format("{:.2f}", opt.r_int_max)));
  ui.set_angle_max(slint::SharedString(std::format("{:.2f}", opt.angle_max)));
  ui.set_angle_bin_width(
      slint::SharedString(std::format("{:.2f}", opt.angle_bin_width)));
  ui.set_bond_factor(
      slint::SharedString(std::format("{:.2f}", opt.bond_factor)));
  ui.set_smoothing_sigma(
      slint::SharedString(std::format("{:.2f}", opt.smoothing_sigma)));
  ui.set_smoothing_kernel(static_cast<int>(opt.smoothing_kernel));
};

ProgramOptions AppController::handleOptionsfromUI(AppWindow &ui) {
  ProgramOptions opt;
  const std::string input_path_str = ui_.get_in_file_text().data();
  std::filesystem::path full_path(input_path_str);
  const std::string output_path_base =
      full_path.parent_path().string() + "/" + full_path.stem().string();
  opt.input_file = input_path_str;
  opt.output_file_base = output_path_base;
  opt.normalize = ui_.get_normalize();
  opt.smoothing = ui_.get_smoothing();
  opt.r_max = std::stof(ui_.get_r_max().data());
  opt.r_bin_width = std::stof(ui_.get_r_bin_width().data());
  opt.q_max = std::stof(ui_.get_q_max().data());
  opt.q_bin_width = std::stof(ui_.get_q_bin_width().data());
  opt.r_int_max = std::stof(ui_.get_r_int_max().data());
  opt.angle_max = std::stof(ui_.get_angle_max().data());
  opt.angle_bin_width = std::stof(ui_.get_angle_bin_width().data());
  opt.bond_factor = std::stof(ui_.get_bond_factor().data());
  opt.smoothing_sigma = std::stof(ui_.get_smoothing_sigma().data());
  opt.smoothing_kernel = static_cast<KernelType>(ui_.get_smoothing_kernel());

  return opt;
};

//---------------------------------------------------------------------------//
//--------------------------------- Methods ---------------------------------//
//---------------------------------------------------------------------------//

void AppController::handleRunAnalysis() {
  // Create a ProgramOptions object from the UI properties
  ProgramOptions opt = handleOptionsfromUI(ui_);
  // run analysis
  backend_.run_analysis(opt);
}

void AppController::handleBrowseFile() {
  std::vector<std::string> filters = {"Supported Structure Files",
                                      "*.car *.cell *.cif *.dat",
                                      "Materials Studio CAR",
                                      "*.car",
                                      "CASTEP CELL",
                                      "*.cell",
                                      "CIF files",
                                      "*.cif",
                                      "ONETEP DAT",
                                      "*.dat",
                                      "All Files",
                                      "*"};

  current_file_dialog_ = std::make_unique<pfd::open_file>(
      "Select a structure file", "", filters, pfd::opt::multiselect);

  ui_.set_timer_running(true);
  ui_.set_text_timer_running(false);
  ui_.set_text_opacity(true);
}

void AppController::handleCheckFileDialogStatus() {
  if (current_file_dialog_ && current_file_dialog_->ready(0)) {
    auto selection = current_file_dialog_->result();
    if (!selection.empty()) {
      ui_.set_in_file_text(slint::SharedString(selection[0]));
      backend_.load_file(selection[0]);
      ui_.set_status_text("File: " + slint::SharedString(selection[0]) +
                          " loaded successfully.");
    } else {
      ui_.set_status_text("File selection cancelled.");
    }
    ui_.set_timer_running(false);
    current_file_dialog_.reset();
  }
}
