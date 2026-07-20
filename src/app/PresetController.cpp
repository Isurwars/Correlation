/**
 * @file PresetController.cpp
 * @brief Implementation of PresetController.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/PresetController.hpp"
#include "AppWindow.h"
#include "app/AppController.hpp"
#include "app/InputValidator.hpp"
#include "app/PresetManager.hpp"
#include <format>

namespace correlation::app {

PresetController::PresetController(::AppWindow &window, AppBackend &backend, AppController &controller)
    : window_(window), backend_(backend), controller_(controller) {}

void PresetController::handleLoadPreset(int index) {
  if (index < 0 || static_cast<size_t>(index) >= presets_.size()) {
    return;
  }

  const Preset &preset = presets_[index];

  // Update backend options
  ProgramOptions current_opts = backend_.options();
  // Keep input_file and output_file_base
  const std::string input = current_opts.input_file;
  const std::string output = current_opts.output_file_base;
  current_opts = preset.options;
  current_opts.input_file = input;
  current_opts.output_file_base = output;
  backend_.setOptions(current_opts);

  // Update UI components with the new options
  controller_.handleOptionstoUI();
  controller_.populateCalculatorGroups();
  controller_.updateActiveGroupFlags();
  static_cast<void>(controller_.getInputValidator()->validateInputs());
  controller_.getInputValidator()->updateCliCommand();

  window_.set_analysis_status_text(slint::SharedString(std::string("Loaded preset: ") + preset.name));
}

void PresetController::handleSavePreset(const std::string &name) {
  if (name.empty()) {
    return;
  }

  Preset preset;
  preset.name = name;
  preset.description = "Saved from Correlation UI";
  preset.options = controller_.handleOptionsfromUI();

  PresetManager::save(preset);

  window_.set_analysis_status_text(slint::SharedString(std::string("Saved preset: ") + name));
  refreshPresetList();
}

void PresetController::handleDeletePreset(int index) {
  if (index < 0 || static_cast<size_t>(index) >= presets_.size()) {
    return;
  }

  const std::string name = presets_[index].name;
  PresetManager::remove(name);

  window_.set_analysis_status_text(slint::SharedString(std::string("Deleted preset: ") + name));
  refreshPresetList();
}

void PresetController::refreshPresetList() {
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

void PresetController::handleMaterialTypeChanged(int type) {
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
  } else if (type == 0) { // Amorphous (0)
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

  // Force re-validation and update the CLI equivalent
  static_cast<void>(controller_.getInputValidator()->validateInputs());
}

} // namespace correlation::app
