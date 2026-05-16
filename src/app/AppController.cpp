/**
 * @file AppController.cpp
 * @brief Implementation of the application controller.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#if defined(_WIN32)
#define NOMINMAX
#include <Windows.h>
#endif

#include "app/AppController.hpp"
#include "calculators/CalculatorFactory.hpp"
#include "plotters/SvgPlotter.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <format>
#include <map>
#include <string>
#include <thread>
#include <vector>

namespace correlation::app {

//---------------------------------------------------------------------------//
//------------------------------- Constructors ------------------------------//
//---------------------------------------------------------------------------//

AppController::AppController(AppWindow &ui, AppBackend &backend)
    : ui_(ui), backend_(backend) {

  // Populate the calculator groups from the factory
  populateCalculatorGroups(ui_);

  // default options to UI
  handleOptionstoUI(ui_);

  // Connect the UI signals to the controller's member functions.
  // We use lambdas to capture 'this' and call the appropriate method.
  ui_.on_run_analysis([this]() { handleRunAnalysis(); });
  ui_.on_write_files([this]() { handleWriteFiles(); });
  ui_.on_browse_file([this]() { handleBrowseFile(); });
  ui_.on_check_file_dialog_status([this]() { handleCheckFileDialogStatus(); });

  // Handle calculator toggle: update backend options and refresh the UI model
  ui_.on_toggle_calculator([this](slint::SharedString id, bool enabled) {
    backend_.setCalculatorActive(std::string(id.data()), enabled);
  });

  // Handle plot selection: generate SVG and push to UI
  ui_.on_select_plot([this](int index) { handleSelectPlot(index); });

  // Handle save plot request
  ui_.on_save_plot([this]() { handleSavePlot(); });
}

AppController::~AppController() {
  if (analysis_thread_.joinable()) {
    analysis_thread_.join();
  }
  if (load_thread_.joinable()) {
    load_thread_.join();
  }
}

//---------------------------------------------------------------------------//
//----------------------------- Private Methods -----------------------------//
//---------------------------------------------------------------------------//

// Safe conversion helper
/**
 * @brief Safely converts a Slint SharedString to a float with a default
 * fallback.
 *
 * @param s The Slint string to parse.
 * @param default_value The value to return if parsing fails.
 * @return The parsed float or default_value on error.
 */
float safe_stof(const slint::SharedString &s, float default_value) {
  try {
    return std::stof(s.data());
  } catch (const std::exception &e) {
    // Optionally, log the error or update a UI status message
    return default_value;
  }
}

void AppController::updateProgress(float p, const std::string &msg) {
  if (p < 0.0f)
    p = 0.0f;
  if (p > 1.0f)
    p = 1.0f;
  slint::invoke_from_event_loop([=, this]() {
    ui_.set_progress(p);
    if (!msg.empty()) {
      ui_.set_analysis_status_text(slint::SharedString(msg));
    }
  });
}

void AppController::handleOptionstoUI(AppWindow &ui) {
  ProgramOptions opt = backend_.options();
  ui.set_in_file_text(slint::SharedString(opt.input_file));
  ui.set_smoothing(opt.smoothing);
  ui.set_use_hdf5(opt.use_hdf5);
  ui.set_use_csv(opt.use_csv);
  ui.set_use_parquet(opt.use_parquet);
  ui.set_r_max(slint::SharedString(std::format("{:.2f}", opt.r_max)));
  ui.set_r_bin_width(
      slint::SharedString(std::format("{:.2f}", opt.r_bin_width)));
  ui.set_q_max(slint::SharedString(std::format("{:.2f}", opt.q_max)));
  ui.set_q_bin_width(
      slint::SharedString(std::format("{:.2f}", opt.q_bin_width)));
  ui.set_r_int_max(slint::SharedString(std::format("{:.2f}", opt.r_int_max)));
  ui.set_angle_bin_width(
      slint::SharedString(std::format("{:.2f}", opt.angle_bin_width)));
  ui.set_dihedral_bin_width(
      slint::SharedString(std::format("{:.2f}", opt.dihedral_bin_width)));
  ui.set_max_ring_size(slint::SharedString(std::to_string(opt.max_ring_size)));
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
  opt.use_parquet = ui_.get_use_parquet();
  opt.r_max = safe_stof(ui_.get_r_max(), opt.r_max);
  opt.r_bin_width = safe_stof(ui_.get_r_bin_width(), opt.r_bin_width);
  opt.q_max = safe_stof(ui_.get_q_max(), opt.q_max);
  opt.q_bin_width = safe_stof(ui_.get_q_bin_width(), opt.q_bin_width);
  opt.r_int_max = safe_stof(ui_.get_r_int_max(), opt.r_int_max);
  opt.angle_bin_width =
      safe_stof(ui_.get_angle_bin_width(), opt.angle_bin_width);
  opt.dihedral_bin_width =
      safe_stof(ui_.get_dihedral_bin_width(), opt.dihedral_bin_width);
  opt.max_ring_size = static_cast<size_t>(safe_stof(
      ui_.get_max_ring_size(), static_cast<float>(opt.max_ring_size)));

  // Collect active_calculators from the UI model
  auto groups = ui_.get_calculator_groups();
  for (size_t gi = 0; gi < groups->row_count(); ++gi) {
    auto group = groups->row_data(gi).value();
    for (size_t ci = 0; ci < group.calculators->row_count(); ++ci) {
      auto calc = group.calculators->row_data(ci).value();
      opt.active_calculators[std::string(calc.id.data())] = calc.enabled;
    }
  }

  opt.smoothing_sigma =
      safe_stof(ui_.get_smoothing_sigma(), opt.smoothing_sigma);
  opt.smoothing_kernel =
      static_cast<correlation::math::KernelType>(ui_.get_smoothing_kernel());

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
  opt.bond_cutoffs_sq.resize(n, std::vector<double>(n));
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      opt.bond_cutoffs_sq[i][j] = cutoffs[i][j] * cutoffs[i][j];
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

void AppController::populateCalculatorGroups(AppWindow &ui) {
  const auto &calculators =
      ::correlation::calculators::CalculatorFactory::instance().getCalculators();
  const auto &opts = backend_.options();

  // Collect group names in insertion order
  std::vector<std::string> group_order;
  std::map<std::string,
           std::vector<const correlation::calculators::BaseCalculator *>>
      groups_map;
  for (const auto &calc : calculators) {
    const std::string &grp = calc->getGroup();
    if (groups_map.find(grp) == groups_map.end()) {
      group_order.push_back(grp);
    }
    groups_map[grp].push_back(calc.get());
  }

  auto groups_model = std::make_shared<slint::VectorModel<CalculatorGroup>>();
  for (const auto &grp_name : group_order) {
    auto calcs_model = std::make_shared<slint::VectorModel<CalculatorInfo>>();
    for (const auto *calc : groups_map.at(grp_name)) {
      bool enabled = true; // default on
      auto it = opts.active_calculators.find(calc->getName());
      if (it != opts.active_calculators.end()) {
        enabled = it->second;
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
  ui.set_calculator_groups(groups_model);
  ui.set_total_calculator_count(static_cast<int>(calculators.size()));
}

//---------------------------------------------------------------------------//
//---------------------------------- Methods --------------------------------//
//---------------------------------------------------------------------------//

void AppController::handleRunAnalysis() {
  ui_.set_analysis_done(false); // Reset done state
  ui_.set_analysis_running(true);
  ui_.set_progress(0.0f);
  ui_.set_analysis_status_text(
      slint::SharedString(AppDefaults::MSG_RUNNING_ANALYSIS));

  // Create a ProgramOptions object from the UI properties
  backend_.setOptions(handleOptionsfromUI(ui_));

  // Set the progress callback
  backend_.setProgressCallback(
      [this](float p, const std::string &msg) { updateProgress(p, msg); });

  // run analysis in a separate thread
  if (analysis_thread_.joinable()) {
    analysis_thread_.join();
  }

  analysis_thread_ = std::thread([this]() {
    std::string err = backend_.run_analysis();

    slint::invoke_from_event_loop([this, err]() {
      ui_.set_analysis_running(false);
      if (err.empty()) {
        ui_.set_analysis_status_text(
            slint::SharedString(AppDefaults::MSG_ANALYSIS_ENDED));
      } else {
        ui_.set_analysis_status_text(slint::SharedString(err));
      }
      ui_.set_analysis_done(true);
      ui_.set_progress(1.0f);

      // Populate the plot dropdown and auto-preview the first histogram
      populatePlotList();
      if (!available_plot_keys_.empty()) {
        handleSelectPlot(0);
      }
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
  ui_.set_analysis_status_text(
      slint::SharedString(AppDefaults::MSG_SELECTING_OUTPUT));
}

void AppController::handleBrowseFile() {
  std::vector<std::string> filters = {
      "Supported Structure Files",
      "*arc *.car *.cell *.cif *.dat *.md *.outmol *.poscar *.contcar *.vasp *.xdatcar",
      "Materials Studio CAR",
      "*.car",
      "Materials Studio ARC",
      "*.arc",
      "CASTEP CELL",
      "*.cell",
      "CASTEP MD",
      "*.md",
      "CIF files",
      "*.cif",
      "ONETEP DAT",
      "*.dat",
      "DMol3 Outmol",
      "*.outmol",
      "VASP POSCAR/CONTCAR",
      "POSCAR CONTCAR *.poscar *.contcar *.vasp",
      "VASP XDATCAR",
      "XDATCAR *.xdatcar",
      "All Files",
      "*"};

  current_file_dialog_ = std::make_unique<pfd::open_file>(
      "Select a structure file", "", filters, pfd::opt::multiselect);

  ui_.set_timer_running(true);
  ui_.set_text_opacity(true);
}

void AppController::handleCheckFileDialogStatus() {
  if (current_file_dialog_ && current_file_dialog_->ready(0)) {
    auto selection = current_file_dialog_->result();
    if (!selection.empty()) {
      std::string filepath = selection[0];
      ui_.set_in_file_text(slint::SharedString(filepath));
      ui_.set_file_status_text(slint::SharedString("Loading file..."));
      ui_.set_timer_running(true);
      ui_.set_text_opacity(true);
      ui_.set_progress(0.0f);

      if (load_thread_.joinable()) {
        load_thread_.join();
      }

      load_thread_ = std::thread([this, filepath]() {
        backend_.setProgressCallback([this](float p, const std::string &msg) {
          slint::invoke_from_event_loop([=, this]() {
            if (!msg.empty()) {
              ui_.set_file_status_text(slint::SharedString(msg));
            }
          });
        });

        std::string message;
        bool success = false;
        try {
          message = backend_.load_file(filepath);
          success = true;
        } catch (const std::exception &e) {
          message = std::string(AppDefaults::MSG_ERROR_LOADING) +
                    std::string(e.what());
        }

        slint::invoke_from_event_loop([this, message, success]() {
          ui_.set_file_status_text(slint::SharedString(message));
          ui_.set_timer_running(false);
          ui_.set_text_opacity(false);

          if (success && backend_.cell()) {
            ui_.set_file_loaded(true);

            // Atom Counts
            auto atom_counts_map = backend_.getAtomCounts();
            auto slint_atom_counts =
                std::make_shared<slint::VectorModel<AtomCount>>();
            for (const auto &[symbol, count] : atom_counts_map) {
              slint_atom_counts->push_back(
                  {slint::SharedString(symbol), count});
            }
            ui_.set_atom_counts(slint_atom_counts);

            // Bond Cutoffs
            setBondCutoffs(ui_);

            // File Info
            ui_.set_num_frames(backend_.getFrameCount());
            ui_.set_total_atoms(backend_.getTotalAtomCount());
            ui_.set_removed_frames_count(
                static_cast<int>(backend_.getRemovedFrameCount()));
            ui_.set_time_step(slint::SharedString(
                std::format("{:.2f}", backend_.getRecommendedTimeStep())));

            // Update Run Analysis Card Frame Info
            ui_.set_min_frame("1");
            ui_.set_max_frame(
                slint::SharedString(std::to_string(backend_.getFrameCount())));
          }
        });
      });
    } else {
      std::string message = AppDefaults::MSG_FILE_SELECTION_CANCELLED;
      ui_.set_file_status_text(slint::SharedString(message));
      slint::invoke_from_event_loop([this]() {
        ui_.set_timer_running(false);
        ui_.set_text_opacity(false);
      });
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
      std::string err = backend_.write_files();
      if (err.empty()) {
        ui_.set_analysis_status_text(
            slint::SharedString(AppDefaults::MSG_FILES_WRITTEN));
      } else {
        ui_.set_analysis_status_text(slint::SharedString(err));
      }
    } else {
      ui_.set_analysis_status_text(
          slint::SharedString(AppDefaults::MSG_SAVE_CANCELLED));
    }
    slint::invoke_from_event_loop([this]() {
      ui_.set_timer_running(false);
      ui_.set_text_opacity(false);
    });
    current_save_dialog_.reset();
  }

  if (current_plot_save_dialog_ && current_plot_save_dialog_->ready(0)) {
    std::string result = current_plot_save_dialog_->result();
    if (!result.empty()) {
      int index = ui_.get_selected_plot_index();
      if (index >= 0 && index < static_cast<int>(available_plot_keys_.size())) {
        const std::string &name = available_plot_keys_[index];
        const correlation::analysis::Histogram *hist = backend_.getHistogram(name);
        if (hist) {
          correlation::plotters::PlotConfig config;
          config.theme = ui_.get_is_dark()
                             ? correlation::plotters::PlotConfig::Theme::Dark
                             : correlation::plotters::PlotConfig::Theme::Light;
          std::string svg = correlation::plotters::renderHistogramAsSvg(*hist, config);
          
          std::ofstream out(result);
          if (out.is_open()) {
            out << svg;
            out.close();
            ui_.set_analysis_status_text(slint::SharedString(name + " plot saved successfully."));
          } else {
            ui_.set_analysis_status_text(slint::SharedString("Failed to open file for saving plot."));
          }
        }
      }
    } else {
      ui_.set_analysis_status_text(
          slint::SharedString(AppDefaults::MSG_SAVE_CANCELLED));
    }
    slint::invoke_from_event_loop([this]() {
      ui_.set_timer_running(false);
      ui_.set_text_opacity(false);
    });
    current_plot_save_dialog_.reset();
  }
}

//---------------------------------------------------------------------------//
//---------------------------- Plot Preview Methods -------------------------//
//---------------------------------------------------------------------------//

void AppController::populatePlotList() {
  auto names = backend_.getAvailableHistogramNames();

  // Custom ordering based on user request:
  // 1. Radial: g(r), G(r), J(r)
  // 2. Scattering: S_q, XRD
  // 3. Angular: BAD, PAD, DAD, CN, RD
  // 4. Trajectory: MSD, VACF, VDOS
  std::map<std::string, int> priority = {
      {"g_r", 0},  {"G_r", 1},   {"J_r", 2},  {"S_q", 10}, {"XRD", 11},
      {"BAD", 20}, {"PAD", 21},  {"DAD", 22}, {"CN", 23},  {"RD", 24},
      {"MSD", 30}, {"VACF", 31}, {"VDOS", 32}};

  std::sort(names.begin(), names.end(),
            [&](const std::string &a, const std::string &b) {
              int p1 = priority.contains(a) ? priority.at(a) : 100;
              int p2 = priority.contains(b) ? priority.at(b) : 100;
              if (p1 != p2)
                return p1 < p2;
              return a < b; // Fallback to alphabetical
            });

  available_plot_keys_ = names; // Store keys corresponding to indices

  // Build MenuItem list: {text: name, enabled: true}
  auto menu_model = std::make_shared<slint::VectorModel<MenuItem>>();
  for (const auto &name : names) {
    MenuItem item;
    const correlation::analysis::Histogram *hist = backend_.getHistogram(name);
    std::string display_text =
        (hist && !hist->title.empty()) ? hist->title : name;
    item.text = slint::SharedString(display_text);
    item.enabled = true;
    menu_model->push_back(item);
  }
  ui_.set_plot_items(menu_model);

  // Reset selection index
  ui_.set_selected_plot_index(names.empty() ? -1 : 0);
}

void AppController::handleSelectPlot(int index) {
  if (index < 0 || index >= static_cast<int>(available_plot_keys_.size()))
    return;
  const std::string &name = available_plot_keys_[index];
  const correlation::analysis::Histogram *hist = backend_.getHistogram(name);
  if (!hist)
    return;
  correlation::plotters::PlotConfig config;
  config.theme = ui_.get_is_dark()
                     ? correlation::plotters::PlotConfig::Theme::Dark
                     : correlation::plotters::PlotConfig::Theme::Light;
  std::string svg = correlation::plotters::renderHistogramAsSvg(*hist, config);

  // Load SVG directly from memory avoiding filesystem issues
  auto img = slint::private_api::load_image_from_embedded_data(
      std::span<const uint8_t>(reinterpret_cast<const uint8_t *>(svg.data()),
                               svg.size()),
      "svg");
  ui_.set_preview_plot(img);
}

void AppController::handleSavePlot() {
  int index = ui_.get_selected_plot_index();
  if (index < 0 || index >= static_cast<int>(available_plot_keys_.size()))
    return;
  const std::string &name = available_plot_keys_[index];
  const correlation::analysis::Histogram *hist = backend_.getHistogram(name);
  if (!hist)
    return;

  std::string default_path = backend_.options().output_file_base;
  if (!default_path.empty()) {
      default_path += "_" + name + ".svg";
  } else {
      default_path = name + ".svg";
  }

  current_plot_save_dialog_ = std::make_unique<pfd::save_file>(
      "Save SVG Plot", default_path,
      std::vector<std::string>{"SVG Image", "*.svg", "All Files", "*"}, pfd::opt::none);

  ui_.set_timer_running(true);
  ui_.set_text_opacity(true);
  ui_.set_analysis_status_text(
      slint::SharedString(std::string("Saving ") + name + " plot..."));
}

} // namespace correlation::app
