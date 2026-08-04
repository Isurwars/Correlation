// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "app/SettingsManager.hpp"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>

namespace correlation::app {

std::filesystem::path SettingsManager::settingsFilePath() {
#ifdef _WIN32
  const char *app_data = std::getenv("APPDATA");
  if (app_data != nullptr && app_data[0] != '\0') {
    return std::filesystem::path(app_data) / "Correlation" / "settings.json";
  }
  const char *user_profile = std::getenv("USERPROFILE");
  if (user_profile != nullptr && user_profile[0] != '\0') {
    return std::filesystem::path(user_profile) / ".correlation" / "settings.json";
  }
  return std::filesystem::current_path() / "settings.json";
#else
  const char *home = std::getenv("HOME");
  if (home != nullptr && home[0] != '\0') {
    return std::filesystem::path(home) / ".correlation" / "settings.json";
  }
  return std::filesystem::current_path() / "settings.json";
#endif
}

namespace {

std::string extractJsonValue(const std::string &line, const std::string &key) {
  const auto pos = line.find("\"" + key + "\"");
  if (pos == std::string::npos) {
    return "";
  }
  const auto colon = line.find(':', pos);
  if (colon == std::string::npos) {
    return "";
  }
  auto start = colon + 1;
  while (start < line.size() && (line[start] == ' ' || line[start] == '\t')) {
    start++;
  }
  auto end = line.find_first_of(",}\r\n", start);
  if (end == std::string::npos) {
    end = line.size();
  }
  return line.substr(start, end - start);
}

void parseSettingsLine(const std::string &line, AppSettings &settings) {
  if (const auto val = extractJsonValue(line, "window_width"); !val.empty()) {
    settings.window_width = std::max(800U, static_cast<uint32_t>(std::stoul(val)));
  }
  if (const auto val = extractJsonValue(line, "window_height"); !val.empty()) {
    settings.window_height = std::max(600U, static_cast<uint32_t>(std::stoul(val)));
  }
  if (const auto val = extractJsonValue(line, "left_col_width"); !val.empty()) {
    settings.left_col_width = std::clamp(std::stof(val), 200.0F, 450.0F);
  }
  if (const auto val = extractJsonValue(line, "middle_col_width"); !val.empty()) {
    settings.middle_col_width = std::clamp(std::stof(val), 200.0F, 450.0F);
  }
}

} // namespace

AppSettings SettingsManager::load() {
  AppSettings settings;
  const auto file_path = settingsFilePath();
  if (!std::filesystem::exists(file_path)) {
    return settings;
  }

  try {
    std::ifstream in_file(file_path);
    if (!in_file.is_open()) {
      return settings;
    }
    std::string line;
    while (std::getline(in_file, line)) {
      parseSettingsLine(line, settings);
    }
  } catch (const std::exception &e) {
    std::cerr << "Failed to parse settings.json: " << e.what() << "\n";
  }
  return settings;
}

bool SettingsManager::save(const AppSettings &settings) {
  try {
    const auto file_path = settingsFilePath();
    std::filesystem::create_directories(file_path.parent_path());

    std::ofstream out(file_path);
    if (!out.is_open()) {
      return false;
    }

    out << "{\n";
    out << "  \"window_width\": " << settings.window_width << ",\n";
    out << "  \"window_height\": " << settings.window_height << ",\n";
    out << "  \"left_col_width\": " << settings.left_col_width << ",\n";
    out << "  \"middle_col_width\": " << settings.middle_col_width << "\n";
    out << "}\n";
    return true;
  } catch (const std::exception &e) {
    std::cerr << "Failed to save settings.json: " << e.what() << "\n";
    return false;
  }
}

} // namespace correlation::app
