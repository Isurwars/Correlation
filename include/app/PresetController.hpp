/**
 * @file PresetController.hpp
 * @brief Preset management logic.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "AppWindow.h"
#include "app/AppBackend.hpp"
#include "app/PresetManager.hpp"

#include <string>
#include <vector>

namespace correlation::app {

class AppController;

/**
 * @class PresetController
 * @brief Manages parameter presets and material type configurations.
 */
class PresetController {
public:
  PresetController(AppWindow &window, AppBackend &backend, AppController &controller);
  ~PresetController() = default;
  PresetController(const PresetController &) = delete;
  PresetController &operator=(const PresetController &) = delete;
  PresetController(PresetController &&) = delete;
  PresetController &operator=(PresetController &&) = delete;

  void handleLoadPreset(int index);
  void handleSavePreset(const std::string &name);
  void handleDeletePreset(int index);
  void refreshPresetList();
  void handleMaterialTypeChanged(int type);

private:
  AppWindow &window_;
  AppBackend &backend_;
  AppController &controller_;

  std::vector<Preset> presets_;
};

} // namespace correlation::app
