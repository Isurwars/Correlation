// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
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

// Safe conversion helper
float safe_stof(const slint::SharedString &s, float default_value) {
  try {
    return std::stof(s.data());
  } catch (const std::exception &e) {
    // Optionally, log the error or update a UI status message
    return default_value;
  }
}

void AppController::handleOptionstoUI(AppWindow &ui) {
  ProgramOptions opt = backend_.options();
  ui.set_in_file_text(slint::SharedString(opt.input_file));
  ui.set_smoothing(opt.smoothing);
  ui.set_use_hdf5(opt.use_hdf5);
  ui.set_use_csv(opt.use_csv);
  ui.set_r_max(slint::SharedString(std::format("{:.2f}", opt.r_max)));
  ui.set_r_bin_width(
      slint::SharedString(std::format("{:.2f}", opt.r_bin_width)));
  ui.set_q_max(slint::SharedString(std::format("{:.2f}", opt.q_max)));
  ui.set_q_bin_width(
      slint::SharedString(std::format("{:.2f}", opt.q_bin_width)));
  ui.set_r_int_max(slint::SharedString(std::format("{:.2f}", opt.r_int_max)));
  ui.set_angle_bin_width(
      slint::SharedString(std::format("{:.2f}", opt.angle_bin_width)));
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
  opt.smoothing = true; // Smoothing is now always on by default
  opt.use_hdf5 = ui_.get_use_hdf5();
  opt.use_csv = ui_.get_use_csv();
  opt.r_max = safe_stof(ui_.get_r_max(), opt.r_max);
  opt.r_bin_width = safe_stof(ui_.get_r_bin_width(), opt.r_bin_width);
  opt.q_max = safe_stof(ui_.get_q_max(), opt.q_max);
  opt.q_bin_width = safe_stof(ui_.get_q_bin_width(), opt.q_bin_width);
  opt.r_int_max = safe_stof(ui_.get_r_int_max(), opt.r_int_max);
  opt.angle_bin_width =
      safe_stof(ui_.get_angle_bin_width(), opt.angle_bin_width);
  opt.smoothing_sigma =
      safe_stof(ui_.get_smoothing_sigma(), opt.smoothing_sigma);
  opt.smoothing_kernel = static_cast<KernelType>(ui_.get_smoothing_kernel());

  // Handle Bond Cutoffs
  auto cutoffs = getBondCutoffs(ui_);
  size_t n = cutoffs.size();
  opt.bond_cutoffs_sq_.resize(n, std::vector<double>(n));
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      opt.bond_cutoffs_sq_[i][j] = cutoffs[i][j] * cutoffs[i][j];
    }
  }

  return opt;
};

//---------------------------------------------------------------------------//
//--------------------------------- Methods ---------------------------------//
//---------------------------------------------------------------------------//

void AppController::handleRunAnalysis() {
  // Create a ProgramOptions object from the UI properties
  backend_.setOptions(handleOptionsfromUI(ui_));

  // run analysis
  backend_.run_analysis();
  ui_.set_analysis_status_text("Analysis ended.");
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
  ui_.set_text_opacity(true);
}

void AppController::handleCheckFileDialogStatus() {
  std::string message = "File selection cancelled.";
  if (current_file_dialog_ && current_file_dialog_->ready(0)) {
    auto selection = current_file_dialog_->result();
    if (!selection.empty()) {
      ui_.set_in_file_text(slint::SharedString(selection[0]));
      try {
        message = backend_.load_file(selection[0]);
      } catch (const std::exception &e) {
        message = "Error loading file: " + std::string(e.what());
      }
    }
    ui_.set_file_status_text(slint::SharedString(message));
    ui_.set_timer_running(false);
    
    // Populate atom counts and bond cutoffs in UI
    if (!selection.empty()) {
      ui_.set_file_loaded(true);
      
      // Atom Counts
      auto atom_counts_map = backend_.getAtomCounts();
      auto slint_atom_counts = std::make_shared<slint::VectorModel<AtomCount>>();
      for (const auto& [symbol, count] : atom_counts_map) {
          slint_atom_counts->push_back({slint::SharedString(symbol), count});
      }
      ui_.set_atom_counts(slint_atom_counts);

      // Bond Cutoffs
      setBondCutoffs(ui_);
    }

    current_file_dialog_.reset();
  }
}

void AppController::setBondCutoffs(AppWindow &ui) {
  auto recommended = backend_.getRecommendedBondCutoffs();
  auto elements = backend_.cell()->elements();
  auto slint_cutoffs = std::make_shared<slint::VectorModel<BondCutoff>>();

  for (size_t i = 0; i < elements.size(); ++i) {
    for (size_t j = i; j < elements.size(); ++j) {
      slint_cutoffs->push_back({
          slint::SharedString(elements[i].symbol),
          slint::SharedString(elements[j].symbol),
          slint::SharedString(std::format("{:.2f}", recommended[i][j]))
      });
    }
  }
  ui.set_bond_cutoffs(slint_cutoffs);
}

std::vector<std::vector<double>> AppController::getBondCutoffs(AppWindow &ui) {
  auto slint_cutoffs = ui.get_bond_cutoffs();
  if (!backend_.cell()) return {};
  auto elements = backend_.cell()->elements();
  size_t num_elements = elements.size();
  std::vector<std::vector<double>> cutoffs(num_elements, std::vector<double>(num_elements, 0.0));

  for (size_t k = 0; k < slint_cutoffs->row_count(); ++k) {
    auto item = slint_cutoffs->row_data(k).value();
    std::string s1 = item.element1.data();
    std::string s2 = item.element2.data();
    double dist = std::stod(item.distance.data());

    int i = -1, j = -1;
    for (size_t e = 0; e < num_elements; ++e) {
      if (elements[e].symbol == s1) i = (int)e;
      if (elements[e].symbol == s2) j = (int)e;
    }

    if (i != -1 && j != -1) {
      cutoffs[i][j] = dist;
      cutoffs[j][i] = dist;
    }
  }
  return cutoffs;
}
