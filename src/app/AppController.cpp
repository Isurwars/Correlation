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

#include "app/AnalysisRunner.hpp"
#include "app/AppController.hpp"

#include "app/FileIOHandler.hpp"
#include "app/InputValidator.hpp"
#include "app/PlotController.hpp"
#include "app/PresetController.hpp"
#include "calculators/CalculatorFactory.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <map>
#include <string>
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

  analysis_runner_ = std::make_unique<AnalysisRunner>(window_, backend_, *this);
  file_io_handler_ = std::make_unique<FileIOHandler>(window_, backend_, *this);
  input_validator_ = std::make_unique<InputValidator>(window_, backend_, *this);
  plot_controller_ = std::make_unique<PlotController>(window_, backend_);
  preset_controller_ = std::make_unique<PresetController>(window_, backend_, *this);

  // default options to UI
  handleOptionstoUI();

  // Connect the UI signals to the controller's member functions.
  // We use lambdas to capture 'this' and call the appropriate method.
  window_.on_run_analysis([this]() { analysis_runner_->handleRunAnalysis(); });
  window_.on_cancel_analysis([this]() { backend_.cancel_analysis(); });
  window_.on_browse_file([this]() { file_io_handler_->handleBrowseFile(); });
  window_.on_write_files([this]() { file_io_handler_->handleWriteFiles(); });
  window_.on_validate_inputs([this]() { input_validator_->validateInputs(); });
  window_.on_copy_cli_command([this]() { input_validator_->handleCopyCliCommand(); });

  // Handle calculator toggle: update backend options and refresh the UI model
  window_.on_toggle_calculator([this](const slint::SharedString &calc_id, bool enabled) {
    backend_.setCalculatorActive(std::string(calc_id.data()), enabled);
    updateActiveGroupFlags();
    input_validator_->updateCliCommand();
  });

  // Handle plot selection: generate SVG and push to UI
  window_.on_select_plot([this](int index) { plot_controller_->requestPlotUpdate(index, true); });

  // Handle mouse move on preview plot
  window_.on_mouse_move([this](float mouse_x, float mouse_y, bool hover, float width, float height) {
    plot_controller_->handleMouseMove(mouse_x, mouse_y, hover, width, height);
  });

  // Handle save plot request (SVG or PDF)
  window_.on_save_plot([this]() { plot_controller_->handleSavePlot(); });

  // Handle pin run request
  window_.on_pin_run([this]() { plot_controller_->handlePinRun(); });

  // Handle clear pinned runs request
  window_.on_clear_pinned_runs([this]() { plot_controller_->handleClearPinnedRuns(); });

  // Handle preset load, save, delete requests
  window_.on_load_preset([this](int index) { preset_controller_->handleLoadPreset(index); });
  window_.on_save_preset(
      [this](const slint::SharedString &name) { preset_controller_->handleSavePreset(std::string(name.data())); });
  window_.on_delete_preset([this](int index) { preset_controller_->handleDeletePreset(index); });
  window_.on_material_type_changed([this](int type) { preset_controller_->handleMaterialTypeChanged(type); });

  // Handle plot resized callback from UI
  window_.on_plot_resized([this](float width, float height) {
    plot_controller_->handlePlotResized({.width = width, .height = height});
  });

  // Initial load of preset list
  preset_controller_->refreshPresetList();
}

AppController::~AppController() {
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
    opts.lef_cutoff = slint::SharedString(std::format("{:.2f}", opt.lef_cutoff));
    window_.set_analysis_options(opts);
  }
  {
    auto opts = window_.get_analysis_options();
    opts.lef_sigma = slint::SharedString(std::format("{:.2f}", opt.lef_sigma));
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
  input_validator_->updateCliCommand();
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

  opt.lef_cutoff = safe_parse(window_.get_analysis_options().lef_cutoff,
                              opt.lef_cutoff); // NOLINT(bugprone-narrowing-conversions)
  opt.lef_sigma = safe_parse(window_.get_analysis_options().lef_sigma,
                             opt.lef_sigma); // NOLINT(bugprone-narrowing-conversions)

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

} // namespace correlation::app
