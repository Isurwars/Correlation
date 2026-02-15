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

  // Connect the UI signals to the controller's member functions.
  // We use lambdas to capture 'this' and call the appropriate method.
  ui_.on_run_analysis([this]() { handleRunAnalysis(); });
  ui_.on_write_files([this]() { handleWriteFiles(); });
  ui_.on_browse_file([this]() { handleBrowseFile(); });
  ui_.on_check_file_dialog_status([this]() { handleCheckFileDialogStatus(); });
}

AppController::~AppController() {
  if (analysis_thread_.joinable()) {
    analysis_thread_.join();
  }
}

//---------------------------------------------------------------------------//
//----------------------------- Private Methods -----------------------------//
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

void AppController::updateProgress(float p) {
  if (p < 0.0f)
    p = 0.0f;
  if (p > 1.0f)
    p = 1.0f;
  slint::invoke_from_event_loop([=, this]() { ui_.set_progress(p); });
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

  ui.set_min_frame(
      slint::SharedString(std::to_string(opt.min_frame + 1))); // UI is 1-based
  if (opt.max_frame == -1) {
    ui.set_max_frame("End");
  } else {
    ui.set_max_frame(slint::SharedString(std::to_string(opt.max_frame)));
  }
  ui.set_time_step(slint::SharedString(std::format("{:.2f}", opt.time_step)));
};

ProgramOptions AppController::handleOptionsfromUI(AppWindow &ui) {
  ProgramOptions opt;
  const std::string input_path_str = ui_.get_in_file_text().data();
  std::filesystem::path full_path(input_path_str);
  std::filesystem::path output_path =
      full_path.parent_path() / full_path.stem();
  opt.input_file = input_path_str;
  opt.output_file_base = output_path.make_preferred().string();
  opt.smoothing = true;
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

  // Frame Selection Logic:
  // - "Start" maps to 0
  // - "End" maps to the last frame index (or -1 for max_frame to indicate
  // 'all')
  // - Numeric values are 1-based in UI, converted to 0-based for backend.

  // Helper lambda for case-insensitive comparison
  auto to_lower = [](const std::string &s) {
    std::string data = s;
    std::transform(data.begin(), data.end(), data.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return data;
  };

  // Frame Selection
  try {
    std::string min_s = ui_.get_min_frame().data();
    std::string min_s_lower = to_lower(min_s);

    if (min_s_lower == "start") {
      opt.min_frame = 0;
    } else if (min_s_lower == "end") {
      opt.min_frame = backend_.getFrameCount() - 1;
      if (opt.min_frame < 0)
        opt.min_frame = 0;
    } else {
      opt.min_frame = std::stoi(min_s) - 1; // UI is 1-based
    }
    if (opt.min_frame < 0)
      opt.min_frame = 0;
  } catch (...) {
    opt.min_frame = 0;
  }

  try {
    std::string max_s = ui_.get_max_frame().data();
    std::string max_s_lower = to_lower(max_s);

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
  } catch (...) {
    opt.max_frame = -1;
  }

  opt.time_step = safe_stof(ui_.get_time_step(), opt.time_step);

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

void AppController::setBondCutoffs(AppWindow &ui) {
  auto recommended = backend_.getRecommendedBondCutoffs();
  auto elements = backend_.cell()->elements();
  auto slint_cutoffs = std::make_shared<slint::VectorModel<BondCutoff>>();

  for (size_t i = 0; i < elements.size(); ++i) {
    for (size_t j = i; j < elements.size(); ++j) {
      slint_cutoffs->push_back(
          {slint::SharedString(elements[i].symbol),
           slint::SharedString(elements[j].symbol),
           slint::SharedString(std::format("{:.2f}", recommended[i][j]))});
    }
  }
  ui.set_bond_cutoffs(slint_cutoffs);
}

std::vector<std::vector<double>> AppController::getBondCutoffs(AppWindow &ui) {
  auto slint_cutoffs = ui.get_bond_cutoffs();
  if (!backend_.cell())
    return {};
  auto elements = backend_.cell()->elements();
  size_t num_elements = elements.size();
  std::vector<std::vector<double>> cutoffs(
      num_elements, std::vector<double>(num_elements, 0.0));

  for (size_t k = 0; k < slint_cutoffs->row_count(); ++k) {
    auto item = slint_cutoffs->row_data(k).value();
    std::string s1 = item.element1.data();
    std::string s2 = item.element2.data();
    double dist = std::stod(item.distance.data());

    int i = -1, j = -1;
    for (size_t e = 0; e < num_elements; ++e) {
      if (elements[e].symbol == s1)
        i = (int)e;
      if (elements[e].symbol == s2)
        j = (int)e;
    }

    if (i != -1 && j != -1) {
      cutoffs[i][j] = dist;
      cutoffs[j][i] = dist;
    }
  }
  return cutoffs;
}

//---------------------------------------------------------------------------//
//---------------------------------- Methods --------------------------------//
//---------------------------------------------------------------------------//

void AppController::handleRunAnalysis() {
  ui_.set_analysis_done(false); // Reset done state
  ui_.set_analysis_running(true);
  ui_.set_progress(0.0f);
  ui_.set_analysis_status_text("Running Analysis...");

  // Create a ProgramOptions object from the UI properties
  backend_.setOptions(handleOptionsfromUI(ui_));

  // Set the progress callback
  backend_.setProgressCallback([this](float p) { updateProgress(p); });

  // run analysis in a separate thread
  if (analysis_thread_.joinable()) {
    analysis_thread_.join();
  }

  analysis_thread_ = std::thread([this]() {
    backend_.run_analysis();

    slint::invoke_from_event_loop([this]() {
      ui_.set_analysis_running(false);
      ui_.set_analysis_status_text("Analysis ended.");
      ui_.set_analysis_done(true);
      ui_.set_progress(1.0f);
    });
  });

  // Detach or move is not enough, we need to keep the thread object alive.
  // We keep it as a member variable.
}

void AppController::handleWriteFiles() {
  backend_.setOptions(handleOptionsfromUI(ui_));

  std::string default_path = backend_.options().output_file_base;

  current_save_dialog_ = std::make_unique<pfd::save_file>(
      "Select Output File Name", default_path,
      std::vector<std::string>{"All Files", "*"}, pfd::opt::none);

  ui_.set_timer_running(true);
  ui_.set_text_opacity(true);
  ui_.set_analysis_status_text("Selecting output file...");
}

void AppController::handleBrowseFile() {
  std::vector<std::string> filters = {"Supported Structure Files",
                                      "*arc *.car *.cell *.cif *.dat",
                                      "Materials Studio CAR",
                                      "*.car",
                                      "Materials Studio ARC",
                                      "*.arc",
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
  if (current_file_dialog_ && current_file_dialog_->ready(0)) {
    std::string message = "File selection cancelled.";
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
      auto slint_atom_counts =
          std::make_shared<slint::VectorModel<AtomCount>>();
      for (const auto &[symbol, count] : atom_counts_map) {
        slint_atom_counts->push_back({slint::SharedString(symbol), count});
      }
      ui_.set_atom_counts(slint_atom_counts);

      // Bond Cutoffs
      setBondCutoffs(ui_);

      // File Info
      ui_.set_num_frames(backend_.getFrameCount());
      ui_.set_total_atoms(backend_.getTotalAtomCount());
      ui_.set_removed_frames_count(
          static_cast<int>(backend_.getRemovedFrameCount()));
      ui_.set_time_step(
          slint::SharedString(std::format("{:.2f}", backend_.getTimeStep())));
    }

    current_file_dialog_.reset();
  }

  if (current_save_dialog_ && current_save_dialog_->ready(0)) {
    std::string result = current_save_dialog_->result();
    if (!result.empty()) {
      std::filesystem::path p(result);
      if (p.has_extension()) {
        p.replace_extension("");
      }

      ProgramOptions opts = handleOptionsfromUI(ui_);
      opts.output_file_base = p.string();
      backend_.setOptions(opts);
      backend_.write_files();
      ui_.set_analysis_status_text("Files Written.");
    } else {
      ui_.set_analysis_status_text("Save cancelled.");
    }
    ui_.set_timer_running(false);
    current_save_dialog_.reset();
  }
}
