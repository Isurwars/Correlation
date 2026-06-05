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
#include <format>
#include <fstream>
#include <map>
#include <string>
#include <thread>
#include <vector>

namespace correlation::app {

//---------------------------------------------------------------------------//
//------------------------------- Constructors ------------------------------//
//---------------------------------------------------------------------------//

AppController::AppController(AppWindow &ui, AppBackend &backend) : ui_(ui), backend_(backend) {
  // Initialize Native File Dialog
  NFD_Init();

  // Populate the calculator groups from the factory
  populateCalculatorGroups(ui_);

  // Set export formats availability based on compile-time definitions
#ifdef CORRELATION_USE_HDF5
  ui_.set_hdf5_available(true);
#else
  ui_.set_hdf5_available(false);
#endif

#ifdef CORRELATION_USE_ARROW
  ui_.set_parquet_available(true);
#else
  ui_.set_parquet_available(false);
#endif

  // default options to UI
  handleOptionstoUI(ui_);

  // Connect the UI signals to the controller's member functions.
  // We use lambdas to capture 'this' and call the appropriate method.
  ui_.on_run_analysis([this]() { handleRunAnalysis(); });
  ui_.on_write_files([this]() { handleWriteFiles(); });
  ui_.on_browse_file([this]() { handleBrowseFile(); });
  ui_.on_validate_inputs([this]() { validateInputs(); });
  ui_.on_copy_cli_command([this]() { handleCopyCliCommand(); });

  // Handle calculator toggle: update backend options and refresh the UI model
  ui_.on_toggle_calculator(
      [this](slint::SharedString id, bool enabled) {
        backend_.setCalculatorActive(std::string(id.data()), enabled);
        updateActiveGroupFlags(ui_);
        updateCliCommand();
      });

  // Handle plot selection: generate SVG and push to UI
  ui_.on_select_plot([this](int index) { handleSelectPlot(index); });

  // Handle mouse move on preview plot
  ui_.on_mouse_move([this](float mx, float my, bool hover, float w, float h) {
    handleMouseMove(mx, my, hover, w, h);
  });

  // Handle save plot request (SVG or PDF)
  ui_.on_save_plot([this]() { handleSavePlot(); });

  // Handle pin run request
  ui_.on_pin_run([this]() { handlePinRun(); });

  // Handle clear pinned runs request
  ui_.on_clear_pinned_runs([this]() { handleClearPinnedRuns(); });

  // Handle preset load, save, delete requests
  ui_.on_load_preset([this](int index) { handleLoadPreset(index); });
  ui_.on_save_preset([this](slint::SharedString name) { handleSavePreset(std::string(name.data())); });
  ui_.on_delete_preset([this](int index) { handleDeletePreset(index); });

  // Initial load of preset list
  refreshPresetList();
}

AppController::~AppController() {
  if (analysis_thread_.joinable()) {
    analysis_thread_.join();
  }
  if (load_thread_.joinable()) {
    load_thread_.join();
  }
  // Quit Native File Dialog
  NFD_Quit();
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
  ui.set_r_bin_width(slint::SharedString(std::format("{:.2f}", opt.r_bin_width)));
  ui.set_q_max(slint::SharedString(std::format("{:.2f}", opt.q_max)));
  ui.set_q_bin_width(slint::SharedString(std::format("{:.2f}", opt.q_bin_width)));
  ui.set_r_int_max(slint::SharedString(std::format("{:.2f}", opt.r_int_max)));
  ui.set_angle_bin_width(slint::SharedString(std::format("{:.2f}", opt.angle_bin_width)));
  ui.set_dihedral_bin_width(slint::SharedString(std::format("{:.2f}", opt.dihedral_bin_width)));
  ui.set_max_ring_size(slint::SharedString(std::to_string(opt.max_ring_size)));
  ui.set_smoothing_sigma(slint::SharedString(std::format("{:.2f}", opt.smoothing_sigma)));
  ui.set_smoothing_kernel(static_cast<int>(opt.smoothing_kernel));

  ui.set_min_frame(slint::SharedString(std::to_string(opt.min_frame + 1))); // UI is 1-based
  if (opt.max_frame == -1) {
    ui.set_max_frame("End");
  } else {
    ui.set_max_frame(slint::SharedString(std::to_string(opt.max_frame)));
  }
  ui.set_time_step(slint::SharedString(std::format("{:.2f}", opt.time_step)));
  updateActiveGroupFlags(ui);
  updateCliCommand();
};

void AppController::updateActiveGroupFlags(AppWindow &ui) {
  const auto &calculators = ::correlation::calculators::CalculatorFactory::instance().getCalculators();
  const auto &opts = backend_.options();

  bool has_radial = false;
  bool has_scattering = false;
  bool has_angular = false;
  bool has_rings = false;

  for (const auto &calc : calculators) {
    const std::string &grp = calc->getGroup();
    bool enabled = true; // default on
    auto it = opts.active_calculators.find(calc->getName());
    if (it != opts.active_calculators.end()) {
      enabled = it->second;
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

  ui.set_has_radial_active(has_radial);
  ui.set_has_scattering_active(has_scattering);
  ui.set_has_angular_active(has_angular);
  ui.set_has_rings_active(has_rings);
}

ProgramOptions AppController::handleOptionsfromUI(AppWindow &ui) {
  ProgramOptions opt;
  const std::string input_path_str = ui_.get_in_file_text().data();
  std::filesystem::path full_path(input_path_str);
  std::filesystem::path output_path = full_path.parent_path() / full_path.stem();
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
  opt.angle_bin_width = safe_stof(ui_.get_angle_bin_width(), opt.angle_bin_width);
  opt.dihedral_bin_width = safe_stof(ui_.get_dihedral_bin_width(), opt.dihedral_bin_width);
  opt.max_ring_size = static_cast<size_t>(safe_stof(ui_.get_max_ring_size(), static_cast<float>(opt.max_ring_size)));

  // Collect active_calculators from the UI model
  auto groups = ui_.get_calculator_groups();
  for (size_t gi = 0; gi < groups->row_count(); ++gi) {
    auto group = groups->row_data(gi).value();
    for (size_t ci = 0; ci < group.calculators->row_count(); ++ci) {
      auto calc = group.calculators->row_data(ci).value();
      opt.active_calculators[std::string(calc.id.data())] = calc.enabled;
    }
  }

  opt.smoothing_sigma = safe_stof(ui_.get_smoothing_sigma(), opt.smoothing_sigma);
  opt.smoothing_kernel = static_cast<correlation::math::KernelType>(ui_.get_smoothing_kernel());

  // Frame Selection Logic:
  // - "Start" maps to 0
  // - "End" maps to the last frame index (or -1 for max_frame to indicate
  // 'all')
  // - Numeric values are 1-based in UI, converted to 0-based for backend.

  // Helper lambda for case-insensitive comparison
  auto to_lower = [](const std::string &s) {
    std::string data = s;
    std::transform(data.begin(), data.end(), data.begin(), [](unsigned char c) { return std::tolower(c); });
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
      slint_cutoffs->push_back({slint::SharedString(elements[i].symbol), slint::SharedString(elements[j].symbol),
                                slint::SharedString(std::format("{:.2f}", recommended[i][j]))});
    }
  }
  ui.set_bond_cutoffs(slint_cutoffs);
}

correlation::plotters::PlotConfig AppController::buildPlotConfigFromUI() {
  correlation::plotters::PlotConfig config;
  config.theme = ui_.get_is_dark() ? correlation::plotters::PlotConfig::Theme::Dark
                                   : correlation::plotters::PlotConfig::Theme::Light;

  int size_preset_val = ui_.get_export_size_preset();
  if (size_preset_val == 1) {
    config.preset_size = correlation::plotters::PlotConfig::PresetSize::SingleColumn;
  } else if (size_preset_val == 2) {
    config.preset_size = correlation::plotters::PlotConfig::PresetSize::DoubleColumn;
  } else if (size_preset_val == 3) {
    config.preset_size = correlation::plotters::PlotConfig::PresetSize::Presentation;
  } else {
    config.preset_size = correlation::plotters::PlotConfig::PresetSize::Default;
  }

  int palette_val = ui_.get_export_palette();
  if (palette_val == 1) {
    config.palette = correlation::plotters::PlotConfig::Palette::Grayscale;
  } else if (palette_val == 2) {
    config.palette = correlation::plotters::PlotConfig::Palette::Viridis;
  } else {
    config.palette = correlation::plotters::PlotConfig::Palette::OkabeIto;
  }

  config.font_scale = safe_stof(ui_.get_export_font_scale(), 1.0f);
  if (config.font_scale <= 0.0f) config.font_scale = 1.0f;

  config.line_width = safe_stof(ui_.get_export_line_width(), 3.0f);
  if (config.line_width <= 0.0f) config.line_width = 3.0f;

  config.show_grid = ui_.get_export_show_grid();
  config.show_legend = ui_.get_export_show_legend();

  return config;
}

std::vector<std::vector<double>> AppController::getBondCutoffs(AppWindow &ui) {
  auto slint_cutoffs = ui.get_bond_cutoffs();
  if (!backend_.cell())
    return {};
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
  const auto &calculators = ::correlation::calculators::CalculatorFactory::instance().getCalculators();
  const auto &opts = backend_.options();

  // Collect group names in insertion order
  std::vector<std::string> group_order;
  std::map<std::string, std::vector<const correlation::calculators::BaseCalculator *>> groups_map;
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
  if (!validateInputs()) {
    return;
  }

  ui_.set_analysis_done(false); // Reset done state
  ui_.set_analysis_running(true);
  ui_.set_progress(0.0f);
  ui_.set_analysis_status_text(slint::SharedString(AppDefaults::MSG_RUNNING_ANALYSIS));

  // Create a ProgramOptions object from the UI properties
  backend_.setOptions(handleOptionsfromUI(ui_));

  // Set the progress callback
  backend_.setProgressCallback([this](float p, const std::string &msg) { updateProgress(p, msg); });

  // run analysis in a separate thread
  if (analysis_thread_.joinable()) {
    analysis_thread_.join();
  }

  analysis_thread_ = std::thread([this]() {
    std::string err = backend_.run_analysis();

    slint::invoke_from_event_loop([this, err]() {
      ui_.set_analysis_running(false);
      if (err.empty()) {
        ui_.set_analysis_status_text(slint::SharedString(AppDefaults::MSG_ANALYSIS_ENDED));
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
  std::filesystem::path dp(default_path);
  std::string default_dir = dp.parent_path().string();
  std::string default_name = dp.filename().string();

  nfdchar_t *outPath = nullptr;
  nfdresult_t result = NFD_SaveDialogU8(&outPath, nullptr, 0,
                                        default_dir.empty() ? nullptr : default_dir.c_str(),
                                        default_name.empty() ? nullptr : default_name.c_str());

  if (result == NFD_OKAY) {
    std::string filepath(outPath);
    NFD_FreePathU8(outPath);

    std::filesystem::path p(filepath);
    if (p.has_extension()) {
      p.replace_extension("");
    }

    ProgramOptions opts = handleOptionsfromUI(ui_);
    opts.output_file_base = p.string();
    backend_.setOptions(opts);
    std::string err = backend_.write_files();
    if (err.empty()) {
      ui_.set_analysis_status_text(slint::SharedString(AppDefaults::MSG_FILES_WRITTEN));
    } else {
      ui_.set_analysis_status_text(slint::SharedString(err));
    }
  } else if (result == NFD_CANCEL) {
    ui_.set_analysis_status_text(slint::SharedString(AppDefaults::MSG_SAVE_CANCELLED));
  } else {
    std::string error_msg = "Error: ";
    error_msg += NFD_GetError();
    ui_.set_analysis_status_text(slint::SharedString(error_msg));
  }
}

void AppController::handleBrowseFile() {
  nfdfilteritem_t filterList[] = {
      {"Supported Structure Files", "arc,car,cell,cif,dat,md,outmol,poscar,contcar,vasp,xdatcar"},
      {"Materials Studio CAR", "car"},
      {"Materials Studio ARC", "arc"},
      {"CASTEP CELL", "cell"},
      {"CASTEP MD", "md"},
      {"CIF files", "cif"},
      {"ONETEP DAT", "dat"},
      {"DMol3 Outmol", "outmol"},
      {"VASP POSCAR/CONTCAR", "poscar,contcar,vasp"},
      {"VASP XDATCAR", "xdatcar"}
  };
  nfdfiltersize_t filterCount = sizeof(filterList) / sizeof(filterList[0]);

  nfdchar_t *outPath = nullptr;
  nfdresult_t result = NFD_OpenDialogU8(&outPath, filterList, filterCount, nullptr);

  if (result == NFD_OKAY) {
    std::string filepath(outPath);
    NFD_FreePathU8(outPath);

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
        message = std::string(AppDefaults::MSG_ERROR_LOADING) + std::string(e.what());
      }

      slint::invoke_from_event_loop([this, message, success]() {
        ui_.set_file_status_text(slint::SharedString(message));
        ui_.set_timer_running(false);
        ui_.set_text_opacity(false);

        if (success && backend_.cell()) {
          ui_.set_file_loaded(true);

          // Atom Counts
          auto atom_counts_map = backend_.getAtomCounts();
          auto slint_atom_counts = std::make_shared<slint::VectorModel<AtomCount>>();
          for (const auto &[symbol, count] : atom_counts_map) {
            slint_atom_counts->push_back({slint::SharedString(symbol), count});
          }
          ui_.set_atom_counts(slint_atom_counts);

          // Bond Cutoffs
          setBondCutoffs(ui_);

          // File Info
          ui_.set_num_frames(backend_.getFrameCount());
          ui_.set_total_atoms(backend_.getTotalAtomCount());
          ui_.set_removed_frames_count(static_cast<int>(backend_.getRemovedFrameCount()));
          ui_.set_time_step(slint::SharedString(std::format("{:.2f}", backend_.getRecommendedTimeStep())));

          // Update Run Analysis Card Frame Info
          ui_.set_min_frame("1");
          ui_.set_max_frame(slint::SharedString(std::to_string(backend_.getFrameCount())));
          validateInputs();
        }
      });
    });
  } else if (result == NFD_CANCEL) {
    std::string message = AppDefaults::MSG_FILE_SELECTION_CANCELLED;
    ui_.set_file_status_text(slint::SharedString(message));
    ui_.set_timer_running(false);
    ui_.set_text_opacity(false);
  } else {
    std::string error_msg = "Error opening file dialog: ";
    error_msg += NFD_GetError();
    ui_.set_file_status_text(slint::SharedString(error_msg));
    ui_.set_timer_running(false);
    ui_.set_text_opacity(false);
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
  std::map<std::string, int> priority = {{"g_r", 0},  {"G_r", 1},   {"J_r", 2},  {"S_q", 10}, {"XRD", 11},
                                         {"BAD", 20}, {"PAD", 21},  {"DAD", 22}, {"CN", 23},  {"RD", 24},
                                         {"MSD", 30}, {"VACF", 31}, {"VDOS", 32}};

  std::sort(names.begin(), names.end(), [&](const std::string &a, const std::string &b) {
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
    std::string display_text = (hist && !hist->title.empty()) ? hist->title : name;
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

  if (ui_.get_selected_plot_index() != index) {
    ui_.set_selected_plot_index(index);
  }

  correlation::plotters::PlotConfig config = buildPlotConfigFromUI();
  correlation::plotters::HoverInfo hover;
  hover.active = mouse_hover_;
  hover.mouse_x = last_mouse_x_;
  hover.mouse_y = last_mouse_y_;
  hover.widget_width = last_plot_width_;
  hover.widget_height = last_plot_height_;

  std::string svg;
  if (pinned_runs_.empty()) {
    svg = correlation::plotters::renderHistogramAsSvg(*hist, config, hover);
  } else {
    std::vector<correlation::plotters::LabeledHistogram> datasets;

    // Add the current run first
    datasets.push_back({"Current", hist});

    // Add pinned runs
    for (const auto &pr : pinned_runs_) {
      auto it = pr.histograms.find(name);
      if (it != pr.histograms.end()) {
        datasets.push_back({pr.label, &it->second});
      }
    }

    std::string key = "Total";
    const auto &partials = hist->smoothed_partials.empty() ? hist->partials : hist->smoothed_partials;
    if (!partials.empty() && partials.find(key) == partials.end()) {
      key = partials.begin()->first;
    }
    svg = correlation::plotters::renderComparisonSvg(datasets, key, config, hover);
  }

  // Load SVG directly from memory avoiding filesystem issues
  auto img = slint::private_api::load_image_from_embedded_data(
      std::span<const uint8_t>(reinterpret_cast<const uint8_t *>(svg.data()), svg.size()), "svg");
  ui_.set_preview_plot(img);
}

void AppController::handleMouseMove(float mx, float my, bool hover, float w, float h) {
  last_mouse_x_ = mx;
  last_mouse_y_ = my;
  mouse_hover_ = hover;
  last_plot_width_ = w;
  last_plot_height_ = h;

  int current_idx = ui_.get_selected_plot_index();
  if (current_idx >= 0) {
    handleSelectPlot(current_idx);
  }
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
    default_path += "_" + name;
  } else {
    default_path = name;
  }

  std::filesystem::path dp(default_path);
  std::string default_dir = dp.parent_path().string();
  std::string default_name = dp.filename().string();

  nfdfilteritem_t filterList[] = {
      {"SVG Image", "svg"},
      {"PDF Document", "pdf"}
  };
  nfdfiltersize_t filterCount = sizeof(filterList) / sizeof(filterList[0]);

  nfdchar_t *outPath = nullptr;
  nfdresult_t result = NFD_SaveDialogU8(&outPath, filterList, filterCount,
                                        default_dir.empty() ? nullptr : default_dir.c_str(),
                                        default_name.empty() ? nullptr : default_name.c_str());

  if (result == NFD_OKAY) {
    std::string filepath(outPath);
    NFD_FreePathU8(outPath);

    correlation::plotters::PlotConfig config = buildPlotConfigFromUI();

    std::string svg;
    if (pinned_runs_.empty()) {
      svg = correlation::plotters::renderHistogramAsSvg(*hist, config);
    } else {
      std::vector<correlation::plotters::LabeledHistogram> datasets;
      datasets.push_back({"Current", hist});
      for (const auto &pr : pinned_runs_) {
        auto it = pr.histograms.find(name);
        if (it != pr.histograms.end()) {
          datasets.push_back({pr.label, &it->second});
        }
      }
      std::string key = "Total";
      const auto &partials = hist->smoothed_partials.empty() ? hist->partials : hist->smoothed_partials;
      if (!partials.empty() && partials.find(key) == partials.end()) {
        key = partials.begin()->first;
      }
      svg = correlation::plotters::renderComparisonSvg(datasets, key, config);
    }

    // Check if file extension is PDF (case-insensitive)
    bool save_as_pdf = false;
    if (filepath.size() >= 4) {
      std::string ext = filepath.substr(filepath.size() - 4);
      if (ext == ".pdf" || ext == ".PDF") {
        save_as_pdf = true;
      }
    }

    if (save_as_pdf) {
      std::string temp_svg_path = filepath + ".tmp.svg";
      std::ofstream out(temp_svg_path);
      if (out.is_open()) {
        out << svg;
        out.close();

        std::string cmd = "rsvg-convert -f pdf -o \"" + filepath + "\" \"" + temp_svg_path + "\"";
        int ret = std::system(cmd.c_str());
        if (ret != 0) {
          std::string inkscape_cmd = "inkscape --export-filename=\"" + filepath + "\" \"" + temp_svg_path + "\"";
          ret = std::system(inkscape_cmd.c_str());
        }

        std::filesystem::remove(temp_svg_path);

        if (ret == 0) {
          ui_.set_analysis_status_text(slint::SharedString(name + " plot saved as PDF successfully."));
        } else {
          ui_.set_analysis_status_text(slint::SharedString("Failed to convert SVG to PDF."));
        }
      } else {
        ui_.set_analysis_status_text(slint::SharedString("Failed to write temporary SVG."));
      }
    } else {
      // Save as SVG
      std::ofstream out(filepath);
      if (out.is_open()) {
        out << svg;
        out.close();
        ui_.set_analysis_status_text(slint::SharedString(name + " plot saved successfully."));
      } else {
        ui_.set_analysis_status_text(slint::SharedString("Failed to open file for saving plot."));
      }
    }
  } else if (result == NFD_CANCEL) {
    ui_.set_analysis_status_text(slint::SharedString(AppDefaults::MSG_SAVE_CANCELLED));
  } else {
    std::string error_msg = "Error: ";
    error_msg += NFD_GetError();
    ui_.set_analysis_status_text(slint::SharedString(error_msg));
  }
}

void AppController::handlePinRun() {
  const auto &hists = backend_.getHistograms();
  if (hists.empty())
    return;

  std::string label = "Run " + std::to_string(pinned_runs_.size() + 1);
  pinned_runs_.push_back({label, hists});

  slint::invoke_from_event_loop([this]() {
    ui_.set_pinned_runs_count(static_cast<int>(pinned_runs_.size()));

    // Refresh plot if we have one selected
    int current_idx = ui_.get_selected_plot_index();
    if (current_idx >= 0) {
      handleSelectPlot(current_idx);
    }
  });
}

void AppController::handleClearPinnedRuns() {
  pinned_runs_.clear();

  slint::invoke_from_event_loop([this]() {
    ui_.set_pinned_runs_count(0);

    // Refresh plot if we have one selected
    int current_idx = ui_.get_selected_plot_index();
    if (current_idx >= 0) {
      handleSelectPlot(current_idx);
    }
  });
}

bool AppController::validateInputs() {
  bool valid = true;

  auto to_lower = [](const std::string &s) {
    std::string data = s;
    std::transform(data.begin(), data.end(), data.begin(), [](unsigned char c) { return std::tolower(c); });
    return data;
  };

  auto is_positive_float = [](const std::string &s, float &val) {
    try {
      size_t idx;
      float v = std::stof(s, &idx);
      if (idx < s.size()) {
        return false;
      }
      if (v <= 0.0f) return false;
      val = v;
      return true;
    } catch (...) {
      return false;
    }
  };

  auto is_non_negative_float = [](const std::string &s, float &val) {
    try {
      size_t idx;
      float v = std::stof(s, &idx);
      if (idx < s.size()) {
        return false;
      }
      if (v < 0.0f) return false;
      val = v;
      return true;
    } catch (...) {
      return false;
    }
  };

  auto is_positive_int = [](const std::string &s, int &val) {
    try {
      size_t idx;
      int v = std::stoi(s, &idx);
      if (idx < s.size()) {
        return false;
      }
      if (v <= 0) return false;
      val = v;
      return true;
    } catch (...) {
      return false;
    }
  };

  float r_max_val = 0.0f;
  float r_bin_val = 0.0f;
  float q_max_val = 0.0f;
  float q_bin_val = 0.0f;

  // 1. r_max
  std::string r_max_s = ui_.get_r_max().data();
  if (!is_positive_float(r_max_s, r_max_val)) {
    ui_.set_r_max_error("Must be a positive number");
    valid = false;
  } else {
    ui_.set_r_max_error("");
  }

  // 2. r_bin_width
  std::string r_bin_s = ui_.get_r_bin_width().data();
  if (!is_positive_float(r_bin_s, r_bin_val)) {
    ui_.set_r_bin_error("Must be a positive number");
    valid = false;
  } else if (r_max_val > 0.0f && r_bin_val > r_max_val) {
    ui_.set_r_bin_error("Must be ≤ r_max");
    valid = false;
  } else {
    ui_.set_r_bin_error("");
  }

  // 3. q_max
  std::string q_max_s = ui_.get_q_max().data();
  if (!is_positive_float(q_max_s, q_max_val)) {
    ui_.set_q_max_error("Must be a positive number");
    valid = false;
  } else {
    ui_.set_q_max_error("");
  }

  // 4. q_bin_width
  std::string q_bin_s = ui_.get_q_bin_width().data();
  if (!is_positive_float(q_bin_s, q_bin_val)) {
    ui_.set_q_bin_error("Must be a positive number");
    valid = false;
  } else if (q_max_val > 0.0f && q_bin_val > q_max_val) {
    ui_.set_q_bin_error("Must be ≤ q_max");
    valid = false;
  } else {
    ui_.set_q_bin_error("");
  }

  // 5. r_int_max
  float r_int_max_val = 0.0f;
  std::string r_int_max_s = ui_.get_r_int_max().data();
  if (!is_positive_float(r_int_max_s, r_int_max_val)) {
    ui_.set_r_int_max_error("Must be a positive number");
    valid = false;
  } else {
    ui_.set_r_int_max_error("");
  }

  // 6. angle_bin_width
  float angle_bin_val = 0.0f;
  std::string angle_bin_s = ui_.get_angle_bin_width().data();
  if (!is_positive_float(angle_bin_s, angle_bin_val)) {
    ui_.set_angle_bin_error("Must be a positive number");
    valid = false;
  } else if (angle_bin_val > 180.0f) {
    ui_.set_angle_bin_error("Must be ≤ 180°");
    valid = false;
  } else {
    ui_.set_angle_bin_error("");
  }

  // 7. dihedral_bin_width
  float dihedral_bin_val = 0.0f;
  std::string dihedral_bin_s = ui_.get_dihedral_bin_width().data();
  if (!is_positive_float(dihedral_bin_s, dihedral_bin_val)) {
    ui_.set_dihedral_bin_error("Must be a positive number");
    valid = false;
  } else if (dihedral_bin_val > 360.0f) {
    ui_.set_dihedral_bin_error("Must be ≤ 360°");
    valid = false;
  } else {
    ui_.set_dihedral_bin_error("");
  }

  // 8. max_ring_size
  int max_ring_val = 0;
  std::string max_ring_s = ui_.get_max_ring_size().data();
  if (!is_positive_int(max_ring_s, max_ring_val) || max_ring_val < 3) {
    ui_.set_max_ring_error("Must be an integer ≥ 3");
    valid = false;
  } else {
    ui_.set_max_ring_error("");
  }

  // 9. smoothing_sigma
  float smoothing_sigma_val = 0.0f;
  std::string smoothing_sigma_s = ui_.get_smoothing_sigma().data();
  if (!is_non_negative_float(smoothing_sigma_s, smoothing_sigma_val)) {
    ui_.set_smoothing_sigma_error("Must be a non-negative number");
    valid = false;
  } else {
    ui_.set_smoothing_sigma_error("");
  }

  // 10. time_step
  float time_step_val = 0.0f;
  std::string time_step_s = ui_.get_time_step().data();
  if (!is_positive_float(time_step_s, time_step_val)) {
    ui_.set_time_step_error("Must be a positive number");
    valid = false;
  } else {
    ui_.set_time_step_error("");
  }

  // 11. min_frame and max_frame
  int min_frame_val = -1;
  int max_frame_val = -1;
  bool min_frame_valid = true;
  bool max_frame_valid = true;

  std::string min_frame_s = ui_.get_min_frame().data();
  std::string min_frame_lower = to_lower(min_frame_s);
  int total_frames = ui_.get_num_frames();

  if (min_frame_lower == "start" || min_frame_s.empty()) {
    min_frame_val = 0;
    ui_.set_min_frame_error("");
  } else if (min_frame_lower == "end") {
    if (total_frames > 0) {
      min_frame_val = total_frames - 1;
    } else {
      min_frame_val = 0;
    }
    ui_.set_min_frame_error("");
  } else {
    int parsed_val;
    if (!is_positive_int(min_frame_s, parsed_val)) {
      ui_.set_min_frame_error("Must be positive integer, 'Start', or 'End'");
      min_frame_valid = false;
      valid = false;
    } else if (total_frames > 0 && parsed_val > total_frames) {
      ui_.set_min_frame_error(slint::SharedString(std::format("Must be ≤ total frames ({})", total_frames)));
      min_frame_valid = false;
      valid = false;
    } else {
      min_frame_val = parsed_val - 1;
      ui_.set_min_frame_error("");
    }
  }

  std::string max_frame_s = ui_.get_max_frame().data();
  std::string max_frame_lower = to_lower(max_frame_s);

  if (max_frame_lower == "end" || max_frame_s.empty()) {
    max_frame_val = total_frames > 0 ? total_frames - 1 : -1;
    ui_.set_max_frame_error("");
  } else if (max_frame_lower == "start") {
    max_frame_val = 0;
    ui_.set_max_frame_error("");
  } else {
    int parsed_val;
    if (!is_positive_int(max_frame_s, parsed_val)) {
      ui_.set_max_frame_error("Must be positive integer or 'End'");
      max_frame_valid = false;
      valid = false;
    } else if (total_frames > 0 && parsed_val > total_frames) {
      ui_.set_max_frame_error(slint::SharedString(std::format("Must be ≤ total frames ({})", total_frames)));
      max_frame_valid = false;
      valid = false;
    } else {
      max_frame_val = parsed_val - 1;
      ui_.set_max_frame_error("");
    }
  }

  // Cross-validation of min_frame and max_frame
  if (min_frame_valid && max_frame_valid && min_frame_val >= 0 && max_frame_val >= 0) {
    if (min_frame_val > max_frame_val) {
      ui_.set_min_frame_error("Start frame must be ≤ End frame");
      ui_.set_max_frame_error("End frame must be ≥ Start frame");
      valid = false;
    }
  }

  // 12. export_font_scale and export_line_width
  float font_scale_val = 0.0f;
  std::string font_scale_s = ui_.get_export_font_scale().data();
  if (!is_positive_float(font_scale_s, font_scale_val)) {
    ui_.set_export_font_scale_error("Must be a positive number");
    valid = false;
  } else {
    ui_.set_export_font_scale_error("");
  }

  float line_width_val = 0.0f;
  std::string line_width_s = ui_.get_export_line_width().data();
  if (!is_positive_float(line_width_s, line_width_val)) {
    ui_.set_export_line_width_error("Must be a positive number");
    valid = false;
  } else {
    ui_.set_export_line_width_error("");
  }

  ui_.set_has_validation_errors(!valid);
  updateCliCommand();
  return valid;
}

void AppController::updateCliCommand() {
  ProgramOptions opt = handleOptionsfromUI(ui_);
  
  std::string cmd = "correlation-cli";

  if (!opt.input_file.empty()) {
    cmd += " \"" + opt.input_file + "\"";
  } else {
    cmd += " <input_file>";
  }

  cmd += " --r-max " + std::string(ui_.get_r_max().data());
  cmd += " --r-bin " + std::string(ui_.get_r_bin_width().data());

  if (ui_.get_has_scattering_active()) {
    cmd += " --q-max " + std::string(ui_.get_q_max().data());
    cmd += " --q-bin " + std::string(ui_.get_q_bin_width().data());
    cmd += " --r-int-max " + std::string(ui_.get_r_int_max().data());
  }

  if (ui_.get_has_angular_active()) {
    cmd += " --angle-bin " + std::string(ui_.get_angle_bin_width().data());
    cmd += " --dihedral-bin " + std::string(ui_.get_dihedral_bin_width().data());
  }

  if (ui_.get_has_rings_active()) {
    cmd += " --max-ring-size " + std::string(ui_.get_max_ring_size().data());
  }

  if (ui_.get_num_frames() > 1) {
    cmd += " --time-step " + std::string(ui_.get_time_step().data());
    cmd += " --min-frame " + std::string(ui_.get_min_frame().data());
    cmd += " --max-frame " + std::string(ui_.get_max_frame().data());
  }

  if (opt.smoothing) {
    cmd += " --smoothing-sigma " + std::string(ui_.get_smoothing_sigma().data());
    std::string kernel_str = "gaussian";
    if (opt.smoothing_kernel == correlation::math::KernelType::Bump) {
      kernel_str = "bump";
    } else if (opt.smoothing_kernel == correlation::math::KernelType::Triweight) {
      kernel_str = "triweight";
    }
    cmd += " --smoothing-kernel " + kernel_str;
  } else {
    cmd += " --no-smoothing";
  }

  // Check if any groups are completely disabled in GUI to add to --disable-groups list
  std::map<std::string, std::vector<const correlation::calculators::BaseCalculator *>> groups_map;
  const auto &calculators = ::correlation::calculators::CalculatorFactory::instance().getCalculators();
  for (const auto &calc : calculators) {
    groups_map[calc->getGroup()].push_back(calc.get());
  }

  std::vector<std::string> disabled_groups;
  for (const auto &[grp_name, grp_calcs] : groups_map) {
    bool all_disabled = true;
    for (const auto *calc : grp_calcs) {
      auto it = opt.active_calculators.find(calc->getName());
      if (it != opt.active_calculators.end() && it->second) {
        all_disabled = false;
        break;
      }
    }
    if (all_disabled) {
      std::string grp_lower = grp_name;
      std::transform(grp_lower.begin(), grp_lower.end(), grp_lower.begin(), ::tolower);
      disabled_groups.push_back(grp_lower);
    }
  }

  if (!disabled_groups.empty()) {
    std::string disabled_list = "";
    for (const auto &grp : disabled_groups) {
      if (!disabled_list.empty()) {
        disabled_list += ",";
      }
      disabled_list += grp;
    }
    cmd += " --disable-groups " + disabled_list;
  }

  if (opt.use_csv) {
    cmd += " --csv";
  } else {
    cmd += " --no-csv";
  }

  if (opt.use_hdf5) {
    cmd += " --hdf5";
  } else {
    cmd += " --no-hdf5";
  }

  if (opt.use_parquet) {
    cmd += " --parquet";
  } else {
    cmd += " --no-parquet";
  }

  ui_.set_cli_command(slint::SharedString(cmd));
}

void AppController::handleCopyCliCommand() {
  std::string cmd = ui_.get_cli_command().data();
  
  // Platform-specific clipboard copy helper
#if defined(__linux__)
  FILE *pipe = popen("xclip -selection clipboard", "w");
  if (pipe) {
    fwrite(cmd.c_str(), 1, cmd.size(), pipe);
    pclose(pipe);
  }
#elif defined(_WIN32)
  if (OpenClipboard(nullptr)) {
    EmptyClipboard();
    HGLOBAL hGlob = GlobalAlloc(GMEM_MOVEABLE, cmd.size() + 1);
    if (hGlob) {
      void *pMem = GlobalLock(hGlob);
      if (pMem) {
        memcpy(pMem, cmd.c_str(), cmd.size() + 1);
        GlobalUnlock(hGlob);
        SetClipboardData(CF_TEXT, hGlob);
      }
    }
    CloseClipboard();
  }
#elif defined(__APPLE__)
  FILE *pipe = popen("pbcopy", "w");
  if (pipe) {
    fwrite(cmd.c_str(), 1, cmd.size(), pipe);
    pclose(pipe);
  }
#endif
}

void AppController::handleLoadPreset(int index) {
  if (index < 0 || index >= static_cast<int>(presets_.size()))
    return;
  const Preset &preset = presets_[index];

  // Update backend options
  ProgramOptions current_opts = backend_.options();
  // Keep input_file and output_file_base
  std::string input = current_opts.input_file;
  std::string output = current_opts.output_file_base;
  current_opts = preset.options;
  current_opts.input_file = input;
  current_opts.output_file_base = output;
  backend_.setOptions(current_opts);

  // Update UI components with the new options
  handleOptionstoUI(ui_);
  populateCalculatorGroups(ui_);
  updateActiveGroupFlags(ui_);
  validateInputs();
  updateCliCommand();

  ui_.set_analysis_status_text(slint::SharedString(std::string("Loaded preset: ") + preset.name));
}

void AppController::handleSavePreset(const std::string &name) {
  if (name.empty()) return;

  Preset preset;
  preset.name = name;
  preset.description = "Saved from Correlation UI";
  preset.options = handleOptionsfromUI(ui_);

  PresetManager::save(preset);

  ui_.set_analysis_status_text(slint::SharedString(std::string("Saved preset: ") + name));
  refreshPresetList();
}

void AppController::handleDeletePreset(int index) {
  if (index < 0 || index >= static_cast<int>(presets_.size()))
    return;

  std::string name = presets_[index].name;
  PresetManager::remove(name);

  ui_.set_analysis_status_text(slint::SharedString(std::string("Deleted preset: ") + name));
  refreshPresetList();
}

void AppController::refreshPresetList() {
  presets_ = PresetManager::loadAll();

  auto menu_model = std::make_shared<slint::VectorModel<MenuItem>>();
  for (const auto &p : presets_) {
    MenuItem item;
    item.text = slint::SharedString(p.name);
    item.enabled = true;
    menu_model->push_back(item);
  }

  ui_.set_preset_items(menu_model);
  ui_.set_selected_preset(-1);
}

} // namespace correlation::app
