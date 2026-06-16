/**
 * @file AppController.cpp
 * @brief Implementation of the application controller.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#if defined(_WIN32)
#define NOMINMAX
#include <Windows.h>
#endif

#include "app/AppBackend.hpp"
#include "app/AppController.hpp"
#include "calculators/CalculatorFactory.hpp"
#include <stdlib.h>
#include "plotters/PdfPlotter.hpp"
#include "plotters/SvgPlotter.hpp"

#include <algorithm>
#include <atomic>
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

AppController::AppController(AppWindow &window, AppBackend &backend) : window_(window), backend_(backend) {
#if defined(__linux__)
  // Prevent GTK3 from hanging for 25 seconds if xdg-desktop-portal is missing or failing
  setenv("GTK_USE_PORTAL", "0", 0);
#endif

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

  // default options to UI
  handleOptionstoUI();

  // Connect the UI signals to the controller's member functions.
  // We use lambdas to capture 'this' and call the appropriate method.
  window_.on_run_analysis([this]() { handleRunAnalysis(); });
  window_.on_cancel_analysis([this]() { backend_.cancel_analysis(); });
  window_.on_write_files([this]() { handleWriteFiles(); });
  window_.on_browse_file([this]() { handleBrowseFile(); });
  window_.on_validate_inputs([this]() { validateInputs(); });
  window_.on_copy_cli_command([this]() { handleCopyCliCommand(); });

  // Handle calculator toggle: update backend options and refresh the UI model
  window_.on_toggle_calculator([this](const slint::SharedString &calc_id, bool enabled) {
    backend_.setCalculatorActive(std::string(calc_id.data()), enabled);
    updateActiveGroupFlags();
    updateCliCommand();
  });

  // Handle plot selection: generate SVG and push to UI
  window_.on_select_plot([this](int index) { requestPlotUpdate(index, true); });

  // Handle mouse move on preview plot
  window_.on_mouse_move([this](float mouse_x, float mouse_y, bool hover, float width, float height) {
    handleMouseMove(mouse_x, mouse_y, hover, width, height);
  });

  // Handle save plot request (SVG or PDF)
  window_.on_save_plot([this]() { handleSavePlot(); });

  // Handle pin run request
  window_.on_pin_run([this]() { handlePinRun(); });

  // Handle clear pinned runs request
  window_.on_clear_pinned_runs([this]() { handleClearPinnedRuns(); });

  // Handle preset load, save, delete requests
  window_.on_load_preset([this](int index) { handleLoadPreset(index); });
  window_.on_save_preset([this](const slint::SharedString &name) { handleSavePreset(std::string(name.data())); });
  window_.on_delete_preset([this](int index) { handleDeletePreset(index); });
  window_.on_material_type_changed([this](int type) { handleMaterialTypeChanged(type); });

  // Handle plot resized callback from UI
  window_.on_plot_resized([this](float width, float height) {
    last_plot_width_ = width;
    last_plot_height_ = height;
    int current_idx = window_.get_selected_plot_index();
    if (current_idx >= 0) {
      requestPlotUpdate(current_idx, true);
    }
  });

  // Start periodic 1s timer to poll config and re-plot if anything changes
  update_timer_.start(slint::TimerMode::Repeated, std::chrono::milliseconds(1000), [this]() {
    int current_idx = window_.get_selected_plot_index();
    if (current_idx >= 0) {
      requestPlotUpdate(current_idx, false);
    }
  });

  // Initial load of preset list
  refreshPresetList();
}

AppController::~AppController() {
  if (render_thread_.joinable()) {
    render_thread_.join();
  }
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
namespace {
/**
 * @brief Safely converts a Slint SharedString to a numeric type with a default
 * fallback.
 *
 * @tparam T The numeric type to return (e.g., float, double).
 * @param str The Slint string to parse.
 * @param default_value The value to return if parsing fails.
 * @return The parsed value or default_value on error.
 */
template <typename T> T safe_parse(const slint::SharedString &str, T default_value) {
  try {
    if constexpr (std::is_same_v<T, float>) {
      return std::stof(str.data());
    } else if constexpr (std::is_same_v<T, double>) {
      return std::stod(str.data());
    } else {
      return default_value;
    }
  } catch (const std::exception &e) {
    // Optionally, log the error or update a UI status message
    return default_value;
  }
}
} // namespace

void AppController::updateProgress(float progress, const std::string &msg) {
  progress = std::clamp(progress, 0.0F, 1.0F);
  slint::invoke_from_event_loop([progress, msg, this]() {
    window_.set_progress(progress);
    if (!msg.empty()) {
      window_.set_analysis_status_text(slint::SharedString(msg));
    }
  });
}

void AppController::handleOptionstoUI() {
  ProgramOptions opt = backend_.options();
  window_.set_in_file_text(slint::SharedString(opt.input_file));
  {
    auto opts = window_.get_analysis_options();
    opts.smoothing_enabled = opt.smoothing;
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.r_max = slint::SharedString(std::format("{:.2f}", opt.r_max));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.r_bin_width = slint::SharedString(std::format("{:.2f}", opt.r_bin_width));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.q_max = slint::SharedString(std::format("{:.2f}", opt.q_max));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.q_bin_width = slint::SharedString(std::format("{:.2f}", opt.q_bin_width));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.r_int_max = slint::SharedString(std::format("{:.2f}", opt.r_int_max));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.angle_bin_width = slint::SharedString(std::format("{:.2f}", opt.angle_bin_width));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.dihedral_bin_width = slint::SharedString(std::format("{:.2f}", opt.dihedral_bin_width));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.max_ring_size = slint::SharedString(std::to_string(opt.max_ring_size));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.smoothing_sigma = slint::SharedString(std::format("{:.2f}", opt.smoothing_sigma));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.smoothing_kernel = static_cast<int>(opt.smoothing_kernel);
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.material_type = opt.material_type;
    window_.set_analysis_options(opts);
  }

  {
    auto opts = window_.get_analysis_options();
    opts.min_frame = slint::SharedString(std::to_string(opt.min_frame + 1));
    window_.set_analysis_options(opts);
  } // UI is 1-based
  if (opt.max_frame == -1) {
    {
      auto opts = window_.get_analysis_options();
      opts.max_frame = "End";
      window_.set_analysis_options(opts);
    }
  } else {
    {
      auto opts = window_.get_analysis_options();
      opts.max_frame = slint::SharedString(std::to_string(opt.max_frame));
      window_.set_analysis_options(opts);
    }
  }
  {
    auto opts = window_.get_analysis_options();
    opts.time_step = slint::SharedString(std::format("{:.2f}", opt.time_step));
    window_.set_analysis_options(opts);
  }
  updateActiveGroupFlags();
  updateCliCommand();
};

void AppController::updateActiveGroupFlags() {
  const auto &calculators = ::correlation::calculators::CalculatorFactory::instance().getCalculators();
  const auto &opts = backend_.options();

  bool has_radial = false;
  bool has_scattering = false;
  bool has_angular = false;
  bool has_rings = false;

  for (const auto &calc : calculators) {
    const std::string &grp = calc->getGroup();
    bool enabled = true; // default on
    auto calc_it = opts.active_calculators.find(calc->getName());
    if (calc_it != opts.active_calculators.end()) {
      enabled = calc_it->second;
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

  window_.set_has_radial_active(has_radial);
  window_.set_has_scattering_active(has_scattering);
  window_.set_has_angular_active(has_angular);
  window_.set_has_rings_active(has_rings);
}

ProgramOptions AppController::handleOptionsfromUI() {
  ProgramOptions opt;
  const std::string input_path_str = window_.get_in_file_text().data();
  std::filesystem::path full_path(input_path_str);
  std::filesystem::path output_path = full_path.parent_path() / full_path.stem();
  opt.input_file = input_path_str;
  opt.output_file_base = output_path.make_preferred().string();
  opt.smoothing = true;
  opt.r_max = safe_parse(window_.get_analysis_options().r_max, opt.r_max); // NOLINT(bugprone-narrowing-conversions)
  opt.r_bin_width =
      safe_parse(window_.get_analysis_options().r_bin_width, opt.r_bin_width); // NOLINT(bugprone-narrowing-conversions)
  opt.q_max = safe_parse(window_.get_analysis_options().q_max, opt.q_max);     // NOLINT(bugprone-narrowing-conversions)
  opt.q_bin_width =
      safe_parse(window_.get_analysis_options().q_bin_width, opt.q_bin_width); // NOLINT(bugprone-narrowing-conversions)
  opt.r_int_max =
      safe_parse(window_.get_analysis_options().r_int_max, opt.r_int_max); // NOLINT(bugprone-narrowing-conversions)
  opt.angle_bin_width = safe_parse(window_.get_analysis_options().angle_bin_width,
                                   opt.angle_bin_width); // NOLINT(bugprone-narrowing-conversions)
  opt.dihedral_bin_width = safe_parse(window_.get_analysis_options().dihedral_bin_width,
                                      opt.dihedral_bin_width); // NOLINT(bugprone-narrowing-conversions)
  opt.max_ring_size = static_cast<size_t>(
      safe_parse(window_.get_analysis_options().max_ring_size, static_cast<double>(opt.max_ring_size)));

  // Collect active_calculators from the UI model
  auto groups = window_.get_calculator_groups();
  for (size_t gi = 0; gi < groups->row_count(); ++gi) {
    auto group = groups->row_data(gi).value();
    for (size_t ci = 0; ci < group.calculators->row_count(); ++ci) {
      auto calc = group.calculators->row_data(ci).value();
      opt.active_calculators[std::string(calc.id.data())] = calc.enabled;
    }
  }

  opt.smoothing_sigma = safe_parse(window_.get_analysis_options().smoothing_sigma,
                                   opt.smoothing_sigma); // NOLINT(bugprone-narrowing-conversions)
  opt.smoothing_kernel = static_cast<correlation::math::KernelType>(window_.get_analysis_options().smoothing_kernel);
  opt.material_type = window_.get_analysis_options().material_type;

  // Frame Selection Logic:
  // - "Start" maps to 0
  // - "End" maps to the last frame index (or -1 for max_frame to indicate
  // 'all')
  // - Numeric values are 1-based in UI, converted to 0-based for backend.

  // Helper lambda for case-insensitive comparison
  auto to_lower = [](const std::string &str) -> std::string {
    std::string data = str;
    std::transform(data.begin(), data.end(), data.begin(), [](unsigned char chr) { return std::tolower(chr); });
    return data;
  };

  // Frame Selection
  try {
    std::string min_s = window_.get_analysis_options().min_frame.data();
    std::string min_s_lower = to_lower(min_s);

    if (min_s_lower == "start") {
      opt.min_frame = 0;
    } else if (min_s_lower == "end") {
      opt.min_frame = std::max(0, static_cast<int>(backend_.getFrameCount()) - 1);
    } else {
      opt.min_frame = std::max(0, std::stoi(min_s) - 1); // UI is 1-based
    }
  } catch (...) {
    opt.min_frame = 0;
  }

  try {
    std::string max_s = window_.get_analysis_options().max_frame.data();
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

  opt.time_step = safe_parse(window_.get_analysis_options().time_step, opt.time_step);

  // Handle Bond Cutoffs
  auto cutoffs = getBondCutoffs();
  size_t num_elements = cutoffs.size();
  opt.bond_cutoffs_sq.resize(num_elements, std::vector<double>(num_elements));
  for (size_t i = 0; i < num_elements; ++i) {
    for (size_t j = 0; j < num_elements; ++j) {
      opt.bond_cutoffs_sq[i][j] = cutoffs[i][j] * cutoffs[i][j];
    }
  }

  return opt;
};

void AppController::setBondCutoffs() {
  auto recommended = backend_.getRecommendedBondCutoffs();
  auto elements = backend_.cell()->elements();
  auto slint_cutoffs = std::make_shared<slint::VectorModel<BondCutoff>>();

  for (size_t i = 0; i < elements.size(); ++i) {
    for (size_t j = i; j < elements.size(); ++j) {
      slint_cutoffs->push_back({.element1 = slint::SharedString(elements[i].symbol),
                                .element2 = slint::SharedString(elements[j].symbol),
                                .distance = slint::SharedString(std::format("{:.2f}", recommended[i][j]))});
    }
  }
  window_.set_bond_cutoffs(slint_cutoffs);
}

correlation::plotters::PlotConfig AppController::buildPlotConfigFromUI() {
  correlation::plotters::PlotConfig config;
  config.theme = window_.get_is_dark() ? correlation::plotters::PlotConfig::Theme::Dark
                                       : correlation::plotters::PlotConfig::Theme::Light;

  int size_preset_val = window_.get_export_config().size_preset;
  if (size_preset_val == 1) {
    config.preset_size = correlation::plotters::PlotConfig::PresetSize::SingleColumn;
  } else if (size_preset_val == 2) {
    config.preset_size = correlation::plotters::PlotConfig::PresetSize::DoubleColumn;
  } else if (size_preset_val == 3) {
    config.preset_size = correlation::plotters::PlotConfig::PresetSize::Presentation;
  } else {
    config.preset_size = correlation::plotters::PlotConfig::PresetSize::Default;
    // When using the Default preset, adapt the SVG canvas to the actual
    // widget dimensions so the plot fills the preview panel without
    // letterboxing.  Fall back to the built-in 1200×900 if the widget
    // has not been laid out yet.
    if (last_plot_width_ > 1.0F && last_plot_height_ > 1.0F) {
      config.width = static_cast<double>(last_plot_width_);
      config.height = static_cast<double>(last_plot_height_);
    }
  }

  int palette_val = window_.get_export_config().palette;
  if (palette_val == 1) {
    config.palette = correlation::plotters::PlotConfig::Palette::Grayscale;
  } else if (palette_val == 2) {
    config.palette = correlation::plotters::PlotConfig::Palette::Viridis;
  } else {
    config.palette = correlation::plotters::PlotConfig::Palette::OkabeIto;
  }

  config.font_scale = safe_parse(window_.get_export_config().font_scale, 1.0F);
  if (config.font_scale <= 0.0F) {
    config.font_scale = 1.0F;
  }

  config.line_width = safe_parse(window_.get_export_config().line_width, 3.0F);
  if (config.line_width <= 0.0F) {
    config.line_width = 3.0F;
  }

  config.show_grid = window_.get_export_config().show_grid;
  config.show_legend = window_.get_export_config().show_legend;
  config.show_markers = window_.get_export_config().show_markers;
  config.fill_area = window_.get_export_config().fill_area;

  return config;
}

std::vector<std::vector<double>> AppController::getBondCutoffs() {
  auto slint_cutoffs = window_.get_bond_cutoffs();
  if (backend_.cell() == nullptr) {
    return {};
  }
  auto elements = backend_.cell()->elements();
  size_t num_elements = elements.size();
  std::vector<std::vector<double>> cutoffs(num_elements, std::vector<double>(num_elements, 0.0));

  for (size_t k = 0; k < slint_cutoffs->row_count(); ++k) {
    auto item = slint_cutoffs->row_data(k).value();
    std::string symbol1 = item.element1.data();
    std::string symbol2 = item.element2.data();
    double dist = std::stod(item.distance.data());

    int idx1 = -1;
    int idx2 = -1;
    for (size_t elem_idx = 0; elem_idx < num_elements; ++elem_idx) {
      if (elements[elem_idx].symbol == symbol1) {
        idx1 = static_cast<int>(elem_idx);
      }
      if (elements[elem_idx].symbol == symbol2) {
        idx2 = static_cast<int>(elem_idx);
      }
    }

    if (idx1 != -1 && idx2 != -1) {
      cutoffs[idx1][idx2] = dist;
      cutoffs[idx2][idx1] = dist;
    }
  }
  return cutoffs;
}

void AppController::populateCalculatorGroups() {
  const auto &calculators = ::correlation::calculators::CalculatorFactory::instance().getCalculators();
  const auto &opts = backend_.options();

  // Collect group names in insertion order
  std::vector<std::string> group_order;
  std::map<std::string, std::vector<const correlation::calculators::BaseCalculator *>> groups_map;
  for (const auto &calc : calculators) {
    const std::string &grp = calc->getGroup();
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
      auto calc_it = opts.active_calculators.find(calc->getName());
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

//---------------------------------------------------------------------------//
//---------------------------------- Methods --------------------------------//
//---------------------------------------------------------------------------//

void AppController::handleRunAnalysis() {
  if (!validateInputs()) {
    return;
  }

  window_.set_analysis_done(false); // Reset done state
  window_.set_analysis_running(true);
  window_.set_progress(0.0F);
  window_.set_analysis_status_text(slint::SharedString(AppDefaults::MSG_RUNNING_ANALYSIS));

  // Create a ProgramOptions object from the UI properties
  backend_.setOptions(handleOptionsfromUI());

  // Set the progress callback
  backend_.setProgressCallback([this](float progress, const std::string &msg) { updateProgress(progress, msg); });

  // run analysis in a separate thread
  if (analysis_thread_.joinable()) {
    analysis_thread_.join();
  }

  analysis_thread_ = std::thread([this]() {
    std::string err = backend_.run_analysis();

    slint::invoke_from_event_loop([this, err]() {
      window_.set_analysis_running(false);
      if (err.empty()) {
        if (backend_.is_cancelled()) {
          window_.set_analysis_status_text(slint::SharedString("Analysis Cancelled."));
        } else {
          window_.set_analysis_status_text(slint::SharedString(AppDefaults::MSG_ANALYSIS_ENDED));
        }
      } else {
        window_.set_analysis_status_text(slint::SharedString(err));
      }
      window_.set_analysis_done(true);
      window_.set_progress(1.0F);

      // Populate the plot dropdown and auto-preview the first histogram
      populatePlotList();
      if (!available_plot_keys_.empty()) {
        requestPlotUpdate(0, true);
      }
    });
  });

  // Detach or move is not enough, we need to keep the thread object alive.
  // We keep it as a member variable.
}

void AppController::handleWriteFiles() {
  std::string default_path = backend_.options().output_file_base;
  std::filesystem::path def_path(default_path);
  std::string default_dir = def_path.parent_path().string();
  std::string default_name = def_path.filename().string();

  std::vector<nfdfilteritem_t> filterList;
  filterList.push_back({.name = "CSV (Comma-Separated Values)", .spec = "csv"});
#ifdef CORRELATION_USE_HDF5
  filterList.push_back({.name = "HDF5 Files", .spec = "h5,hdf5"});
#endif
#ifdef CORRELATION_USE_ARROW
  filterList.push_back({.name = "Parquet Files", .spec = "parquet"});
#endif

  nfdchar_t *outPath = nullptr;
  nfdresult_t result = NFD_SaveDialogU8(&outPath, filterList.data(), filterList.size(),
                                        default_dir.empty() ? nullptr : default_dir.c_str(),
                                        default_name.empty() ? nullptr : default_name.c_str());

  if (result == NFD_OKAY) {
    std::string filepath(outPath);
    NFD_FreePathU8(outPath);

    std::filesystem::path file_path_obj(filepath);
    std::string ext = file_path_obj.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    bool use_csv = false;
    bool use_hdf5 = false;
    bool use_parquet = false;

    if (ext == ".h5" || ext == ".hdf5") {
#ifdef CORRELATION_USE_HDF5
      use_hdf5 = true;
#else
      window_.set_analysis_status_text("Error: HDF5 format is not available.");
      return;
#endif
    } else if (ext == ".parquet") {
#ifdef CORRELATION_USE_ARROW
      use_parquet = true;
#else
      window_.set_analysis_status_text("Error: Parquet format is not available.");
      return;
#endif
    } else {
      use_csv = true;
      if (ext != ".csv") {
        file_path_obj.replace_extension(".csv");
      }
    }

    if (file_path_obj.has_extension()) {
      file_path_obj.replace_extension("");
    }

    ProgramOptions opts = handleOptionsfromUI();
    opts.output_file_base = file_path_obj.string();
    opts.use_csv = use_csv;
    opts.use_hdf5 = use_hdf5;
    opts.use_parquet = use_parquet;
    backend_.setOptions(opts);

    std::string err = backend_.write_files();
    if (err.empty()) {
      window_.set_analysis_status_text(slint::SharedString(AppDefaults::MSG_FILES_WRITTEN));
    } else {
      window_.set_analysis_status_text(slint::SharedString(err));
    }
  } else if (result == NFD_CANCEL) {
    window_.set_analysis_status_text(slint::SharedString(AppDefaults::MSG_SAVE_CANCELLED));
  } else {
    std::string error_msg = "Error: ";
    error_msg += NFD_GetError();
    window_.set_analysis_status_text(slint::SharedString(error_msg));
  }
}

void AppController::handleBrowseFile() {
  std::array<nfdfilteritem_t, 10> filterList = {
      {{.name = "Supported Structure Files", .spec = "arc,car,cell,cif,dat,md,outmol,poscar,contcar,vasp,xdatcar"},
       {.name = "Materials Studio CAR", .spec = "car"},
       {.name = "Materials Studio ARC", .spec = "arc"},
       {.name = "CASTEP CELL", .spec = "cell"},
       {.name = "CASTEP MD", .spec = "md"},
       {.name = "CIF files", .spec = "cif"},
       {.name = "ONETEP DAT", .spec = "dat"},
       {.name = "DMol3 Outmol", .spec = "outmol"},
       {.name = "VASP POSCAR/CONTCAR", .spec = "poscar,contcar,vasp"},
       {.name = "VASP XDATCAR", .spec = "xdatcar"}}};
  nfdfiltersize_t filterCount = filterList.size();

  nfdchar_t *outPath = nullptr;
  nfdresult_t result = NFD_OpenDialogU8(&outPath, filterList.data(), filterCount, nullptr);

  if (result == NFD_OKAY) {
    std::string filepath(outPath);
    NFD_FreePathU8(outPath);

    window_.set_in_file_text(slint::SharedString(filepath));
    window_.set_file_status_text(slint::SharedString("Loading file..."));
    window_.set_timer_running(true);
    window_.set_text_opacity(true);
    window_.set_progress(0.0F);

    if (load_thread_.joinable()) {
      load_thread_.join();
    }

    load_thread_ = std::thread([this, filepath]() {
      backend_.setProgressCallback([this](float progress, const std::string &msg) {
        slint::invoke_from_event_loop([progress, msg, this]() {
          window_.set_progress(progress);
          if (!msg.empty()) {
            window_.set_file_status_text(slint::SharedString(msg));
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
        window_.set_file_status_text(slint::SharedString(message));
        window_.set_timer_running(false);
        window_.set_text_opacity(false);

        if (success && backend_.cell() != nullptr) {
          window_.set_file_loaded(true);

          // Atom Counts
          auto atom_counts_map = backend_.getAtomCounts();
          auto slint_atom_counts = std::make_shared<slint::VectorModel<AtomCount>>();
          for (const auto &[symbol, count] : atom_counts_map) {
            slint_atom_counts->push_back({slint::SharedString(symbol), count});
          }
          window_.set_atom_counts(slint_atom_counts);

          // Bond Cutoffs
          setBondCutoffs();

          // File Info
          window_.set_num_frames(static_cast<int>(backend_.getFrameCount()));
          window_.set_total_atoms(static_cast<int>(backend_.getTotalAtomCount()));
          window_.set_removed_frames_count(static_cast<int>(backend_.getRemovedFrameCount()));
          {
            auto opts = window_.get_analysis_options();
            opts.time_step = slint::SharedString(std::format("{:.2f}", backend_.getRecommendedTimeStep()));
            window_.set_analysis_options(opts);
          }

          // Update Run Analysis Card Frame Info
          {
            auto opts = window_.get_analysis_options();
            opts.min_frame = "1";
            window_.set_analysis_options(opts);
          }
          {
            auto opts = window_.get_analysis_options();
            opts.max_frame = slint::SharedString(std::to_string(backend_.getFrameCount()));
            window_.set_analysis_options(opts);
          }
          validateInputs();
        }
      });
    });
  } else if (result == NFD_CANCEL) {
    std::string message = AppDefaults::MSG_FILE_SELECTION_CANCELLED;
    window_.set_file_status_text(slint::SharedString(message));
    window_.set_timer_running(false);
    window_.set_text_opacity(false);
  } else {
    std::string error_msg = "Error opening file dialog: ";
    error_msg += NFD_GetError();
    window_.set_file_status_text(slint::SharedString(error_msg));
    window_.set_timer_running(false);
    window_.set_text_opacity(false);
  }
}

//---------------------------------------------------------------------------//
//---------------------------- Plot Preview Methods -------------------------//
//---------------------------------------------------------------------------//

void AppController::populatePlotList() {
  last_rendered_index_ = -1;
  auto names = backend_.getAvailableHistogramNames();

  // Custom ordering based on user request:
  // 1. Radial: g(r), G(r), J(r)
  // 2. Scattering: S_q, XRD
  // 3. Angular: BAD, PAD, DAD, CN, RD
  // 4. Trajectory: MSD, VACF, VDOS
  std::map<std::string, int> priority = {{"g_r", 0},  {"G_r", 1},   {"J_r", 2},  {"S_q", 10}, {"XRD", 11},
                                         {"BAD", 20}, {"PAD", 21},  {"DAD", 22}, {"CN", 23},  {"RD", 24},
                                         {"MSD", 30}, {"VACF", 31}, {"VDOS", 32}};

  std::ranges::sort(names, [&](const std::string &lhs, const std::string &rhs) {
    int prio_a = priority.contains(lhs) ? priority.at(lhs) : 100;
    int prio_b = priority.contains(rhs) ? priority.at(rhs) : 100;
    if (prio_a != prio_b) {
      return prio_a < prio_b;
    }
    return lhs < rhs; // Fallback to alphabetical
  });

  available_plot_keys_ = names; // Store keys corresponding to indices

  // Build MenuItem list: {text: name, enabled: true}
  auto menu_model = std::make_shared<slint::VectorModel<MenuItem>>();
  for (const auto &name : names) {
    MenuItem item;
    const correlation::analysis::Histogram *hist = backend_.getHistogram(name);
    std::string display_text = (hist != nullptr && !hist->title.empty()) ? hist->title : name;
    item.text = slint::SharedString(display_text);
    item.enabled = true;
    menu_model->push_back(item);
  }
  window_.set_plot_items(menu_model);

  // Reset selection index
  window_.set_selected_plot_index(names.empty() ? -1 : 0);
}

void AppController::handleMouseMove(float mouse_x, float mouse_y, bool hover, float width, float height) {
  // Trust hover directly from Slint since we made TouchArea visibility stable
  bool actual_hover = hover;

  // Ignore duplicate events if horizontal position, vertical position, and hover state are identical
  if (std::abs(mouse_x - last_mouse_x_) < 1e-2F && std::abs(mouse_y - last_mouse_y_) < 1e-2F &&
      actual_hover == mouse_hover_ && std::abs(width - last_plot_width_) < 1e-2F &&
      std::abs(height - last_plot_height_) < 1e-2F) {
    return;
  }

  bool hover_changed = (actual_hover != mouse_hover_);

  last_mouse_x_ = mouse_x;
  last_mouse_y_ = mouse_y;
  mouse_hover_ = actual_hover;
  last_plot_width_ = width;
  last_plot_height_ = height;

  int current_idx = window_.get_selected_plot_index();
  if (current_idx >= 0) {
    requestPlotUpdate(current_idx, hover_changed || !actual_hover);
  }
}

void AppController::requestPlotUpdate(int index, bool immediate) {
  if (index < 0 || index >= static_cast<int>(available_plot_keys_.size())) {
    return;
  }

  const std::string &name = available_plot_keys_[index];
  const correlation::analysis::Histogram *hist = backend_.getHistogram(name);
  if (hist == nullptr) {
    return;
  }

  if (immediate) {
    needs_redraw_ = true;
  }

  correlation::plotters::PlotConfig config = buildPlotConfigFromUI();
  correlation::plotters::HoverInfo hover;
  hover.active = mouse_hover_;
  hover.mouse_x = last_mouse_x_;
  hover.mouse_y = last_mouse_y_;
  hover.widget_width = last_plot_width_;
  hover.widget_height = last_plot_height_;

  if (isPlotCacheHit(index, config, hover) && !needs_redraw_) {
    return;
  }

  // If a render is already running, queue this request as pending and return
  if (is_rendering_) {
    render_pending_ = true;
    pending_plot_index_ = index;
    return;
  }

  // Start rendering
  is_rendering_ = true;
  needs_redraw_ = false;
  last_rendered_index_ = index;
  last_pinned_runs_count_ = pinned_runs_.size();
  last_config_ = config;
  last_hover_ = hover;

  // Prepare data for the thread (deep copies)
  RenderTaskData data;
  data.active_hist = *hist;
  data.config = config;
  data.hover = hover;
  data.ashcroft_weights = backend_.getAshcroftWeights();

  for (const auto &pinned_run : pinned_runs_) {
    auto hist_it = pinned_run.histograms.find(name);
    if (hist_it != pinned_run.histograms.end()) {
      data.comparison_hists.emplace_back(pinned_run.label, hist_it->second);
    }
  }

  // Join previous render thread if it exists and is done
  if (render_thread_.joinable()) {
    render_thread_.join();
  }

  // Ensure selection index in UI is correct (must be on main/UI thread)
  if (window_.get_selected_plot_index() != index) {
    window_.set_selected_plot_index(index);
  }

  executePlotRender(std::move(data));
}

bool AppController::isPlotCacheHit(int index, const correlation::plotters::PlotConfig &config,
                                   const correlation::plotters::HoverInfo &hover) const {
  return (index == last_rendered_index_ && pinned_runs_.size() == last_pinned_runs_count_ &&
          config.theme == last_config_.theme && config.preset_size == last_config_.preset_size &&
          std::abs(config.width - last_config_.width) < 1e-2 && std::abs(config.height - last_config_.height) < 1e-2 &&
          config.palette == last_config_.palette && std::abs(config.font_scale - last_config_.font_scale) < 1e-4 &&
          std::abs(config.line_width - last_config_.line_width) < 1e-4 && config.show_grid == last_config_.show_grid &&
          config.show_legend == last_config_.show_legend && config.show_markers == last_config_.show_markers &&
          config.fill_area == last_config_.fill_area && hover.active == last_hover_.active &&
          std::abs(hover.mouse_x - last_hover_.mouse_x) < 1e-2 &&
          std::abs(hover.mouse_y - last_hover_.mouse_y) < 1e-2 &&
          std::abs(hover.widget_width - last_hover_.widget_width) < 1e-2 &&
          std::abs(hover.widget_height - last_hover_.widget_height) < 1e-2);
}

void AppController::executePlotRender(RenderTaskData data) {
  render_thread_ = std::thread([this, data]() {
    std::string svg;
    if (data.comparison_hists.empty()) {
      svg =
          correlation::plotters::renderHistogramAsSvg(data.active_hist, data.config, data.hover, data.ashcroft_weights);
    } else {
      std::vector<correlation::plotters::LabeledHistogram> datasets;
      datasets.push_back({"Current", &data.active_hist});
      for (const auto &pair : data.comparison_hists) {
        datasets.push_back({pair.first, &pair.second});
      }

      std::string key = "Total";
      const auto &partials =
          data.active_hist.smoothed_partials.empty() ? data.active_hist.partials : data.active_hist.smoothed_partials;
      if (!partials.empty() && !partials.contains(key)) {
        key = partials.begin()->first;
      }
      svg = correlation::plotters::renderComparisonSvg(datasets, key, data.config, data.hover);
    }

    // Write SVG to a unique temporary file to enable asynchronous parsing by Slint's background thread.
    // (This avoids UI stutter since load_image_from_embedded_data forces synchronous parsing on the UI thread).
    // Note: /tmp is mounted as tmpfs (RAM) on Linux, so this does not cause mechanical disk thrashing.
    static std::atomic<uint64_t> file_counter{0};
    auto temp_dir = std::filesystem::temp_directory_path();
    auto temp_path = temp_dir / ("correlation_preview_" + std::to_string(file_counter++) + ".svg");

    std::string final_temp_path;
    std::ofstream out(temp_path);
    if (out) {
      out << svg;
      out.close();
      final_temp_path = temp_path.string();
    }

    // Dispatch UI updates back to the Slint event loop
    slint::invoke_from_event_loop([this, final_temp_path, svg]() {
      if (!final_temp_path.empty()) {
        auto img = slint::Image::load_from_path(slint::SharedString(final_temp_path));
        window_.set_preview_plot(img);

        // Remove the temporary file. Slint has already captured the path and deferred the load.
        std::error_code error_code;
        std::filesystem::remove(final_temp_path, error_code);
      } else {
        auto img = slint::private_api::load_image_from_embedded_data(
            std::span<const uint8_t>(reinterpret_cast<const uint8_t *>(svg.data()), svg.size()), // NOLINT
            "svg");
        window_.set_preview_plot(img);
      }

      is_rendering_ = false;

      // If another render was requested during this drawing, run it now
      if (render_pending_) {
        render_pending_ = false;
        requestPlotUpdate(pending_plot_index_, true);
      }
    });
  });
}

void AppController::handleSavePlot() {
  int index = window_.get_selected_plot_index();
  if (index < 0 || index >= static_cast<int>(available_plot_keys_.size())) {
    return;
  }
  const std::string &name = available_plot_keys_[index];
  const correlation::analysis::Histogram *hist = backend_.getHistogram(name);
  if (hist == nullptr) {
    return;
  }

  std::string default_path = backend_.options().output_file_base;
  if (!default_path.empty()) {
    default_path += "_" + name;
  } else {
    default_path = name;
  }

  std::filesystem::path def_path(default_path);
  std::string default_dir = def_path.parent_path().string();
  std::string default_name = def_path.filename().string();

  std::array<nfdfilteritem_t, 2> filterList = {
      {{.name = "SVG Image", .spec = "svg"}, {.name = "PDF Document", .spec = "pdf"}}};
  nfdfiltersize_t filterCount = filterList.size();

  nfdchar_t *outPath = nullptr;
  nfdresult_t result =
      NFD_SaveDialogU8(&outPath, filterList.data(), filterCount, default_dir.empty() ? nullptr : default_dir.c_str(),
                       default_name.empty() ? nullptr : default_name.c_str());

  if (result == NFD_OKAY) {
    std::string filepath(outPath);
    NFD_FreePathU8(outPath);

    executeSavePlot(filepath, hist, name);
  } else if (result == NFD_CANCEL) {
    window_.set_analysis_status_text(slint::SharedString(AppDefaults::MSG_SAVE_CANCELLED));
  } else {
    std::string error_msg = "Error: ";
    error_msg += NFD_GetError();
    window_.set_analysis_status_text(slint::SharedString(error_msg));
  }
}

namespace {
std::string getComparisonKey(const correlation::analysis::Histogram *hist) {
  std::string key = "Total";
  const auto &partials = hist->smoothed_partials.empty() ? hist->partials : hist->smoothed_partials;
  if (!partials.empty() && !partials.contains(key)) {
    key = partials.begin()->first;
  }
  return key;
}
} // namespace

void AppController::executeSavePlot(const std::string &filepath, const correlation::analysis::Histogram *hist,
                                    const std::string &name) {
  correlation::plotters::PlotConfig config = buildPlotConfigFromUI();
  config.use_native_text = true;

  if (config.preset_size == correlation::plotters::PlotConfig::PresetSize::Default) {
    config.width = 1200.0;
    config.height = 900.0;
  }

  std::string ext = std::filesystem::path(filepath).extension().string();
  bool save_as_pdf = (ext == ".pdf" || ext == ".PDF");

  auto build_datasets = [&]() {
    std::vector<correlation::plotters::LabeledHistogram> datasets;
    datasets.push_back({"Current", hist});
    for (const auto &pinned_run : pinned_runs_) {
      auto hist_it = pinned_run.histograms.find(name);
      if (hist_it != pinned_run.histograms.end()) {
        datasets.push_back({pinned_run.label, &hist_it->second});
      }
    }
    return datasets;
  };

  if (save_as_pdf) {
    if (pinned_runs_.empty()) {
      correlation::plotters::renderHistogramAsPdf(*hist, filepath, config);
    } else {
      correlation::plotters::renderComparisonPdf(build_datasets(), getComparisonKey(hist), filepath, config);
    }
    window_.set_analysis_status_text(slint::SharedString(name + " plot saved as PDF successfully."));
  } else {
    std::string svg;
    if (pinned_runs_.empty()) {
      svg = correlation::plotters::renderHistogramAsSvg(*hist, config, {}, backend_.getAshcroftWeights());
    } else {
      svg = correlation::plotters::renderComparisonSvg(build_datasets(), getComparisonKey(hist), config);
    }
    std::ofstream out(filepath);
    if (out.is_open()) {
      out << svg;
      out.close();
      window_.set_analysis_status_text(slint::SharedString(name + " plot saved successfully."));
    } else {
      window_.set_analysis_status_text(slint::SharedString("Failed to open file for saving plot."));
    }
  }
}

void AppController::handlePinRun() {
  const auto &hists = backend_.getHistograms();
  if (hists.empty()) {
    return;
  }

  std::string label = "Run " + std::to_string(pinned_runs_.size() + 1);
  pinned_runs_.push_back({label, hists});

  slint::invoke_from_event_loop([this]() {
    window_.set_pinned_runs_count(static_cast<int>(pinned_runs_.size()));

    // Refresh plot if we have one selected
    int current_idx = window_.get_selected_plot_index();
    if (current_idx >= 0) {
      requestPlotUpdate(current_idx, true);
    }
  });
}

void AppController::handleClearPinnedRuns() {
  pinned_runs_.clear();

  slint::invoke_from_event_loop([this]() {
    window_.set_pinned_runs_count(0);

    // Refresh plot if we have one selected
    int current_idx = window_.get_selected_plot_index();
    if (current_idx >= 0) {
      requestPlotUpdate(current_idx, true);
    }
  });
}

namespace {
std::string to_lower_str(const std::string &str) {
  std::string data = str;
  std::ranges::transform(data, data.begin(), [](const unsigned char chr) { return std::tolower(chr); });
  return data;
}

bool is_positive_float(const std::string &str, float &val) {
  try {
    size_t idx = 0;
    float parsed_val = std::stof(str, &idx);
    if (idx < str.size() || parsed_val <= 0.0F) {
      return false;
    }
    val = parsed_val;
    return true;
  } catch (...) {
    return false;
  }
}

bool is_non_negative_float(const std::string &str, float &val) {
  try {
    size_t idx = 0;
    float parsed_val = std::stof(str, &idx);
    if (idx < str.size() || parsed_val < 0.0F) {
      return false;
    }
    val = parsed_val;
    return true;
  } catch (...) {
    return false;
  }
}

bool is_positive_int(const std::string &str, int &val) {
  try {
    size_t idx = 0;
    int parsed_val = std::stoi(str, &idx);
    if (idx < str.size() || parsed_val <= 0) {
      return false;
    }
    val = parsed_val;
    return true;
  } catch (...) {
    return false;
  }
}

bool validate_min_frame_input(const std::string &frame_s, int total_frames, int &frame_val,
                              slint::SharedString &error_out) {
  std::string frame_lower = to_lower_str(frame_s);

  if (frame_lower == "start" || frame_s.empty()) {
    frame_val = 0;
    return true;
  }

  if (frame_lower == "end") {
    frame_val = total_frames > 0 ? total_frames - 1 : 0;
    return true;
  }

  int parsed_val = 0;
  if (!is_positive_int(frame_s, parsed_val)) {
    error_out = "Must be positive integer, 'Start', or 'End'";
    return false;
  }

  if (total_frames > 0 && parsed_val > total_frames) {
    error_out = slint::SharedString(std::format("Must be ≤ total frames ({})", total_frames));
    return false;
  }

  frame_val = parsed_val - 1;
  return true;
}

bool validate_max_frame_input(const std::string &frame_s, int total_frames, int &frame_val,
                              slint::SharedString &error_out) {
  std::string frame_lower = to_lower_str(frame_s);

  if (frame_lower == "end" || frame_s.empty()) {
    frame_val = total_frames > 0 ? total_frames - 1 : -1;
    return true;
  }

  if (frame_lower == "start") {
    frame_val = 0;
    return true;
  }

  int parsed_val = 0;
  if (!is_positive_int(frame_s, parsed_val)) {
    error_out = "Must be positive integer or 'End'";
    return false;
  }

  if (total_frames > 0 && parsed_val > total_frames) {
    error_out = slint::SharedString(std::format("Must be ≤ total frames ({})", total_frames));
    return false;
  }

  frame_val = parsed_val - 1;
  return true;
}

std::string getDisabledGroupsArg(const ProgramOptions &opt) {
  std::map<std::string, std::vector<const correlation::calculators::BaseCalculator *>> groups_map;
  const auto &calculators = ::correlation::calculators::CalculatorFactory::instance().getCalculators();
  for (const auto &calc : calculators) {
    groups_map[calc->getGroup()].push_back(calc.get());
  }

  std::vector<std::string> disabled_groups;
  for (const auto &[grp_name, grp_calcs] : groups_map) {
    bool all_disabled = true;
    for (const auto *calc : grp_calcs) {
      auto calc_it = opt.active_calculators.find(calc->getName());
      if (calc_it != opt.active_calculators.end() && calc_it->second) {
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

  if (disabled_groups.empty()) {
    return "";
  }

  std::string disabled_list;
  for (const auto &grp : disabled_groups) {
    if (!disabled_list.empty()) {
      disabled_list += ",";
    }
    disabled_list += grp;
  }
  return " --disable-groups " + disabled_list;
}
} // namespace

bool AppController::validateInputs() {
  bool valid = true;
  auto errs = window_.get_app_errors();

  errs.r_max_error = "";
  errs.r_bin_error = "";
  errs.q_max_error = "";
  errs.q_bin_error = "";
  errs.r_int_max_error = "";
  errs.angle_bin_error = "";
  errs.dihedral_bin_error = "";
  errs.max_ring_error = "";
  errs.smoothing_sigma_error = "";
  errs.time_step_error = "";
  errs.min_frame_error = "";
  errs.max_frame_error = "";
  errs.export_font_scale_error = "";
  errs.export_line_width_error = "";

  float r_max_val = 0.0F;
  std::string r_max_s = window_.get_analysis_options().r_max.data();
  if (!is_positive_float(r_max_s, r_max_val)) {
    errs.r_max_error = "Must be a positive number";
    valid = false;
  }

  float r_bin_val = 0.0F;
  std::string r_bin_s = window_.get_analysis_options().r_bin_width.data();
  if (!is_positive_float(r_bin_s, r_bin_val)) {
    errs.r_bin_error = "Must be a positive number";
    valid = false;
  } else if (r_max_val > 0.0F && r_bin_val > r_max_val) {
    errs.r_bin_error = "Must be ≤ r_max";
    valid = false;
  }

  float q_max_val = 0.0F;
  std::string q_max_s = window_.get_analysis_options().q_max.data();
  if (!is_positive_float(q_max_s, q_max_val)) {
    errs.q_max_error = "Must be a positive number";
    valid = false;
  }

  float q_bin_val = 0.0F;
  std::string q_bin_s = window_.get_analysis_options().q_bin_width.data();
  if (!is_positive_float(q_bin_s, q_bin_val)) {
    errs.q_bin_error = "Must be a positive number";
    valid = false;
  } else if (q_max_val > 0.0F && q_bin_val > q_max_val) {
    errs.q_bin_error = "Must be ≤ q_max";
    valid = false;
  }

  float r_int_max_val = 0.0F;
  std::string r_int_max_s = window_.get_analysis_options().r_int_max.data();
  if (!is_positive_float(r_int_max_s, r_int_max_val)) {
    errs.r_int_max_error = "Must be a positive number";
    valid = false;
  }

  float angle_bin_val = 0.0F;
  std::string angle_bin_s = window_.get_analysis_options().angle_bin_width.data();
  if (!is_positive_float(angle_bin_s, angle_bin_val)) {
    errs.angle_bin_error = "Must be a positive number";
    valid = false;
  } else if (angle_bin_val > 180.0F) {
    errs.angle_bin_error = "Must be ≤ 180°";
    valid = false;
  }

  float dihedral_bin_val = 0.0F;
  std::string dihedral_bin_s = window_.get_analysis_options().dihedral_bin_width.data();
  if (!is_positive_float(dihedral_bin_s, dihedral_bin_val)) {
    errs.dihedral_bin_error = "Must be a positive number";
    valid = false;
  } else if (dihedral_bin_val > 360.0F) {
    errs.dihedral_bin_error = "Must be ≤ 360°";
    valid = false;
  }

  int max_ring_val = 0;
  std::string max_ring_s = window_.get_analysis_options().max_ring_size.data();
  if (!is_positive_int(max_ring_s, max_ring_val) || max_ring_val < 3) {
    errs.max_ring_error = "Must be an integer ≥ 3";
    valid = false;
  }

  float smoothing_sigma_val = 0.0F;
  std::string smoothing_sigma_s = window_.get_analysis_options().smoothing_sigma.data();
  if (!is_non_negative_float(smoothing_sigma_s, smoothing_sigma_val)) {
    errs.smoothing_sigma_error = "Must be a non-negative number";
    valid = false;
  }

  float time_step_val = 0.0F;
  std::string time_step_s = window_.get_analysis_options().time_step.data();
  if (!is_positive_float(time_step_s, time_step_val)) {
    errs.time_step_error = "Must be a positive number";
    valid = false;
  }

  int min_frame_val = -1;
  int max_frame_val = -1;
  int total_frames = window_.get_num_frames();

  std::string min_frame_s = window_.get_analysis_options().min_frame.data();
  bool min_frame_valid = validate_min_frame_input(min_frame_s, total_frames, min_frame_val, errs.min_frame_error);
  if (!min_frame_valid) {
    valid = false;
  }

  std::string max_frame_s = window_.get_analysis_options().max_frame.data();
  bool max_frame_valid = validate_max_frame_input(max_frame_s, total_frames, max_frame_val, errs.max_frame_error);
  if (!max_frame_valid) {
    valid = false;
  }

  if (min_frame_valid && max_frame_valid && min_frame_val >= 0 && max_frame_val >= 0) {
    if (min_frame_val > max_frame_val) {
      errs.min_frame_error = "Start frame must be ≤ End frame";
      errs.max_frame_error = "End frame must be ≥ Start frame";
      valid = false;
    }
  }

  float font_scale_val = 0.0F;
  std::string font_scale_s = window_.get_export_config().font_scale.data();
  if (!is_positive_float(font_scale_s, font_scale_val)) {
    errs.export_font_scale_error = "Must be a positive number";
    valid = false;
  }

  float line_width_val = 0.0F;
  std::string line_width_s = window_.get_export_config().line_width.data();
  if (!is_positive_float(line_width_s, line_width_val)) {
    errs.export_line_width_error = "Must be a positive number";
    valid = false;
  }

  window_.set_app_errors(errs);
  window_.set_has_validation_errors(!valid);
  updateCliCommand();
  return valid;
}

void AppController::updateCliCommand() {
  ProgramOptions opt = handleOptionsfromUI();

  std::string cmd = "correlation-cli";

  if (!opt.input_file.empty()) {
    cmd += " \"" + opt.input_file + "\"";
  } else {
    cmd += " <input_file>";
  }

  cmd += " --r-max " + std::string(window_.get_analysis_options().r_max.data());
  cmd += " --r-bin " + std::string(window_.get_analysis_options().r_bin_width.data());

  if (window_.get_has_scattering_active()) {
    cmd += " --q-max " + std::string(window_.get_analysis_options().q_max.data());
    cmd += " --q-bin " + std::string(window_.get_analysis_options().q_bin_width.data());
    cmd += " --r-int-max " + std::string(window_.get_analysis_options().r_int_max.data());
  }

  if (window_.get_has_angular_active()) {
    cmd += " --angle-bin " + std::string(window_.get_analysis_options().angle_bin_width.data());
    cmd += " --dihedral-bin " + std::string(window_.get_analysis_options().dihedral_bin_width.data());
  }

  if (window_.get_has_rings_active()) {
    cmd += " --max-ring-size " + std::string(window_.get_analysis_options().max_ring_size.data());
  }

  if (window_.get_num_frames() > 1) {
    cmd += " --time-step " + std::string(window_.get_analysis_options().time_step.data());
    cmd += " --min-frame " + std::string(window_.get_analysis_options().min_frame.data());
    cmd += " --max-frame " + std::string(window_.get_analysis_options().max_frame.data());
  }

  if (opt.smoothing) {
    cmd += " --smoothing-sigma " + std::string(window_.get_analysis_options().smoothing_sigma.data());
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

  cmd += getDisabledGroupsArg(opt);

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

  window_.set_cli_command(slint::SharedString(cmd));
}

void AppController::handleCopyCliCommand() {
  std::string cmd = window_.get_cli_command().data();

  // Platform-specific clipboard copy helper
#if defined(__linux__)
  // NOLINTNEXTLINE(cert-env33-c)
  FILE *pipe = popen("xclip -selection clipboard", "w");
  if (pipe != nullptr) {
    (void)fwrite(cmd.c_str(), 1, cmd.size(), pipe);
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
  // NOLINTNEXTLINE(cert-env33-c)
  FILE *pipe = popen("pbcopy", "w");
  if (pipe != nullptr) {
    (void)fwrite(cmd.c_str(), 1, cmd.size(), pipe);
    pclose(pipe);
  }
#endif
}

void AppController::handleLoadPreset(int index) {
  if (index < 0 || index >= static_cast<int>(presets_.size())) {
    return;
  }

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
  handleOptionstoUI();
  populateCalculatorGroups();
  updateActiveGroupFlags();
  validateInputs();
  updateCliCommand();

  window_.set_analysis_status_text(slint::SharedString(std::string("Loaded preset: ") + preset.name));
}

void AppController::handleSavePreset(const std::string &name) {
  if (name.empty()) {
    return;
  }

  Preset preset;
  preset.name = name;
  preset.description = "Saved from Correlation UI";
  preset.options = handleOptionsfromUI();

  PresetManager::save(preset);

  window_.set_analysis_status_text(slint::SharedString(std::string("Saved preset: ") + name));
  refreshPresetList();
}

void AppController::handleDeletePreset(int index) {
  if (index < 0 || index >= static_cast<int>(presets_.size())) {
    return;
  }

  std::string name = presets_[index].name;
  PresetManager::remove(name);

  window_.set_analysis_status_text(slint::SharedString(std::string("Deleted preset: ") + name));
  refreshPresetList();
}

void AppController::refreshPresetList() {
  presets_ = PresetManager::loadAll();

  auto menu_model = std::make_shared<slint::VectorModel<MenuItem>>();
  for (const auto &preset : presets_) {
    MenuItem item;
    item.text = slint::SharedString(preset.name);
    item.enabled = true;
    menu_model->push_back(item);
  }

  window_.set_preset_items(menu_model);
  window_.set_selected_preset(-1);
}

void AppController::handleMaterialTypeChanged(int type) {
  if (type == 2) { // Crystalline
    {
      auto opts = window_.get_analysis_options();
      opts.r_bin_width = slint::SharedString(std::format("{:.3f}", AppDefaults::R_BIN_WIDTH_CRYSTAL));
      window_.set_analysis_options(opts);
    }
    {
      auto opts = window_.get_analysis_options();
      opts.q_bin_width = slint::SharedString(std::format("{:.3f}", AppDefaults::Q_BIN_WIDTH_CRYSTAL));
      window_.set_analysis_options(opts);
    }
    {
      auto opts = window_.get_analysis_options();
      opts.angle_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::ANGLE_BIN_WIDTH_CRYSTAL));
      window_.set_analysis_options(opts);
    }
    {
      auto opts = window_.get_analysis_options();
      opts.dihedral_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::ANGLE_BIN_WIDTH_CRYSTAL));
      window_.set_analysis_options(opts);
    }
    {
      auto opts = window_.get_analysis_options();
      opts.smoothing_sigma = slint::SharedString(std::format("{:.2f}", AppDefaults::SMOOTHING_SIGMA_CRYSTAL));
      window_.set_analysis_options(opts);
    }
  } else if (type == 1) { // Liquid
    {
      auto opts = window_.get_analysis_options();
      opts.r_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::R_BIN_WIDTH_LIQUID));
      window_.set_analysis_options(opts);
    }
    {
      auto opts = window_.get_analysis_options();
      opts.q_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::Q_BIN_WIDTH_LIQUID));
      window_.set_analysis_options(opts);
    }
    {
      auto opts = window_.get_analysis_options();
      opts.angle_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::ANGLE_BIN_WIDTH_LIQUID));
      window_.set_analysis_options(opts);
    }
    {
      auto opts = window_.get_analysis_options();
      opts.dihedral_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::ANGLE_BIN_WIDTH_LIQUID));
      window_.set_analysis_options(opts);
    }
    {
      auto opts = window_.get_analysis_options();
      opts.smoothing_sigma = slint::SharedString(std::format("{:.2f}", AppDefaults::SMOOTHING_SIGMA_LIQUID));
      window_.set_analysis_options(opts);
    }
  } else { // Amorphous (0)
    {
      auto opts = window_.get_analysis_options();
      opts.r_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::R_BIN_WIDTH));
      window_.set_analysis_options(opts);
    }
    {
      auto opts = window_.get_analysis_options();
      opts.q_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::Q_BIN_WIDTH));
      window_.set_analysis_options(opts);
    }
    {
      auto opts = window_.get_analysis_options();
      opts.angle_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::ANGLE_BIN_WIDTH));
      window_.set_analysis_options(opts);
    }
    {
      auto opts = window_.get_analysis_options();
      opts.dihedral_bin_width = slint::SharedString(std::format("{:.2f}", AppDefaults::ANGLE_BIN_WIDTH));
      window_.set_analysis_options(opts);
    }
    {
      auto opts = window_.get_analysis_options();
      opts.smoothing_sigma = slint::SharedString(std::format("{:.2f}", AppDefaults::SMOOTHING_SIGMA));
      window_.set_analysis_options(opts);
    }
  }
  validateInputs();
}

} // namespace correlation::app
