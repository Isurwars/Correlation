/**
 * @file SettingsManager.hpp
 * @brief Manages application window and layout settings persistence.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include <cstdint>
#include <filesystem>

namespace correlation::app {

/**
 * @struct AppSettings
 * @brief Stores layout geometry and window dimensions across app sessions.
 */
struct AppSettings {
  uint32_t window_width{1180};
  uint32_t window_height{795};
  float left_col_width{260.0F};
  float middle_col_width{260.0F};
};

/**
 * @class SettingsManager
 * @brief Handles saving and loading AppSettings to/from ~/.correlation/settings.json.
 */
class SettingsManager {
public:
  SettingsManager() = delete;

  /**
   * @brief Gets the settings file path (~/.correlation/settings.json).
   */
  [[nodiscard]] static std::filesystem::path settingsFilePath();

  /**
   * @brief Loads application settings from disk. Returns defaults if missing or invalid.
   */
  [[nodiscard]] static AppSettings load();

  /**
   * @brief Saves application settings to disk.
   */
  static bool save(const AppSettings &settings);
};

} // namespace correlation::app
