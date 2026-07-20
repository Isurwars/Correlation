/**
 * @file PresetController.hpp
 * @brief Preset management logic.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "app/AppBackend.hpp"
#include "app/PresetManager.hpp"

#include <string>
#include <vector>

class AppWindow;

namespace correlation::app {

class AppController;

/**
 * @class PresetController
 * @brief Manages parameter presets and material type configurations.
 */
class PresetController {
public:
  /**
   * @brief Constructs the PresetController.
   * @param[in,out] window Reference to the UI window.
   * @param[in,out] backend Reference to the application backend.
   * @param[in,out] controller Reference to the main AppController.
   */
  PresetController(::AppWindow &window, AppBackend &backend, AppController &controller);

  ~PresetController() = default;
  PresetController(const PresetController &) = delete;
  PresetController &operator=(const PresetController &) = delete;
  PresetController(PresetController &&) = delete;
  PresetController &operator=(PresetController &&) = delete;

  /**
   * @brief Loads a preset configuration by UI index.
   * @param[in] index Target index in the preset list.
   */
  void handleLoadPreset(int index);

  /**
   * @brief Saves the current UI configuration as a named preset.
   * @param[in] name Descriptive preset identifier string.
   */
  void handleSavePreset(const std::string &name);

  /**
   * @brief Deletes a saved preset by UI index.
   * @param[in] index Target index in the preset list.
   */
  void handleDeletePreset(int index);

  /**
   * @brief Refreshes the list of available presets in the UI dropdown model.
   */
  void refreshPresetList();

  /**
   * @brief Applies default parameters when material type changes.
   * @param[in] type Selected material type enum/index.
   */
  void handleMaterialTypeChanged(int type);

private:
  ::AppWindow &window_;
  AppBackend &backend_;
  AppController &controller_;

  std::vector<Preset> presets_;
};

} // namespace correlation::app
