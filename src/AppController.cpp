// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "../include/AppController.hpp"

#include <filesystem>
#include <string>
#include <vector>

//---------------------------------------------------------------------------//
//------------------------------- Constructors ------------------------------//
//---------------------------------------------------------------------------//

AppController::AppController(AppWindow &ui, AppBackend &backend)
    : ui_(ui), backend_(backend) {

  // Connect the UI signals to the controller's member functions
  ui_.on_run_analysis([this]() { handleRunAnalysis(); });
  ui_.on_browse_file([this]() { handleBrowseFile(); });
  ui_.on_check_file_dialog_status([this]() { handleCheckFileDialogStatus(); });
}

//---------------------------------------------------------------------------//
//--------------------------------- Methods ---------------------------------//
//---------------------------------------------------------------------------//

void AppController::handleRunAnalysis() {
  // Create a ProgramOptions object from the UI properties
  ProgramOptions options;
  const std::string input_path_str = ui_.get_in_file_text().data();
  std::filesystem::path full_path(input_path_str);
  const std::string output_path_base =
      full_path.parent_path().string() + "/" + full_path.stem().string();
  options.input_file = input_path_str;
  options.output_file_base = output_path_base;
  options.r_max = ui_.get_r_max();
  options.r_bin_width = ui_.get_r_bin_width();
  options.q_max = ui_.get_q_max();
  options.q_bin_width = ui_.get_q_bin_width();
  options.angle_max = ui_.get_angle_max();
  options.angle_bin_width = ui_.get_angle_bin_width();
  options.bond_factor = ui_.get_bond_factor();
  options.smoothing = ui_.get_smoothing();
  options.smoothing_sigma = ui_.get_smoothing_sigma();
  options.smoothing_kernel =
      static_cast<KernelType>(ui_.get_smoothing_kernel());

  backend_.run_analysis(options);
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
