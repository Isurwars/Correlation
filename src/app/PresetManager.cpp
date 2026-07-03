// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "app/PresetManager.hpp"
#include <algorithm>
#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
#include <sstream>

namespace correlation::app {

namespace {

std::string escapeJsonString(const std::string &str) {
  std::string out;
  for (char const chr : str) {
    if (chr == '"') {
      out += "\\\"";
    } else if (chr == '\\') {
      out += "\\\\";
    } else if (chr == '\n') {
      out += "\\n";
    } else if (chr == '\r') {
      out += "\\r";
    } else if (chr == '\t') {
      out += "\\t";
    } else {
      out += chr;
    }
  }
  return out;
}

bool isTargetKeyAt(const std::string &json, size_t pos, const std::string &key) {
  size_t const key_len = key.size();
  if (pos + 1 + key_len >= json.size()) {
    return false;
  }
  if (json.compare(pos + 1, key_len, key) != 0) {
    return false;
  }
  if (json[pos + 1 + key_len] != '"') {
    return false;
  }
  // Verify it is followed by a colon (skipping whitespace)
  size_t const colon_pos = json.find_first_not_of(" \t\n\r", pos + 1 + key_len + 1);
  return (colon_pos != std::string::npos && json[colon_pos] == ':');
}

struct LevelTracker {
  int brace_level = 0;
  int bracket_level = 0;
};

void updateLevels(char chr, LevelTracker &tracker) {
  if (chr == '{') {
    tracker.brace_level++;
  } else if (chr == '}') {
    tracker.brace_level--;
  } else if (chr == '[') {
    tracker.bracket_level++;
  } else if (chr == ']') {
    tracker.bracket_level--;
  }
}

size_t findRootKey(const std::string &json, const std::string &key) {
  bool in_string = false;
  LevelTracker tracker;
  for (size_t i = 0; i < json.size(); ++i) {
    char const chr = json[i];
    if (chr == '\\' && in_string) {
      i++; // skip next char
      continue;
    }
    if (chr == '"') {
      in_string = !in_string;
      if (in_string && tracker.brace_level == 1 && tracker.bracket_level == 0) {
        if (isTargetKeyAt(json, i, key)) {
          return i;
        }
      }
      continue;
    }
    if (!in_string) {
      updateLevels(chr, tracker);
    }
  }
  return std::string::npos;
}

std::string parseStringValue(const std::string &json, const std::string &key) {
  size_t pos = findRootKey(json, key);
  if (pos == std::string::npos) {
    return "";
  }
  pos = json.find(':', pos);
  if (pos == std::string::npos) {
    return "";
  }
  size_t const start = json.find('\"', pos);
  if (start == std::string::npos) {
    return "";
  }

  // Find the end quote, skipping escaped characters
  size_t end = std::string::npos;
  for (size_t i = start + 1; i < json.size(); ++i) {
    if (json[i] == '\\') {
      i++; // skip next character as it is escaped
    } else if (json[i] == '"') {
      end = i;
      break;
    }
  }
  if (end == std::string::npos) {
    return "";
  }
  std::string val = json.substr(start + 1, end - start - 1);

  // Unescape
  std::string unescaped;
  for (size_t i = 0; i < val.size(); ++i) {
    if (val[i] == '\\' && i + 1 < val.size()) {
      char const next = val[i + 1];
      if (next == '"') {
        unescaped += '"';
        i++;
      } else if (next == '\\') {
        unescaped += '\\';
        i++;
      } else if (next == 'n') {
        unescaped += '\n';
        i++;
      } else if (next == 'r') {
        unescaped += '\r';
        i++;
      } else if (next == 't') {
        unescaped += '\t';
        i++;
      } else {
        unescaped += val[i];
      }
    } else {
      unescaped += val[i];
    }
  }
  return unescaped;
}

double parseDoubleValue(const std::string &json, const std::string &key, double fallback) {
  size_t pos = findRootKey(json, key);
  if (pos == std::string::npos) {
    return fallback;
  }
  pos = json.find(':', pos);
  if (pos == std::string::npos) {
    return fallback;
  }
  size_t const start = json.find_first_not_of(" \t\n\r:", pos + 1);
  if (start == std::string::npos) {
    return fallback;
  }
  size_t end = json.find_first_of(",}\n\r", start);
  if (end == std::string::npos) {
    end = json.size();
  }
  try {
    return std::stod(json.substr(start, end - start));
  } catch (...) {
    return fallback;
  }
}

int parseIntValue(const std::string &json, const std::string &key, int fallback) {
  size_t pos = findRootKey(json, key);
  if (pos == std::string::npos) {
    return fallback;
  }
  pos = json.find(':', pos);
  if (pos == std::string::npos) {
    return fallback;
  }
  size_t const start = json.find_first_not_of(" \t\n\r:", pos + 1);
  if (start == std::string::npos) {
    return fallback;
  }
  size_t end = json.find_first_of(",}\n\r", start);
  if (end == std::string::npos) {
    end = json.size();
  }
  try {
    return std::stoi(json.substr(start, end - start));
  } catch (...) {
    return fallback;
  }
}

bool parseBoolValue(const std::string &json, const std::string &key, bool fallback) {
  size_t pos = findRootKey(json, key);
  if (pos == std::string::npos) {
    return fallback;
  }
  pos = json.find(':', pos);
  if (pos == std::string::npos) {
    return fallback;
  }
  size_t const start = json.find_first_not_of(" \t\n\r:", pos + 1);
  if (start == std::string::npos) {
    return fallback;
  }
  std::string const val = json.substr(start, 5);
  if (val.starts_with("true")) {
    return true;
  }
  if (val.starts_with("false")) {
    return false;
  }
  return fallback;
}

std::map<std::string, bool> parseActiveCalculators(const std::string &json) {
  std::map<std::string, bool> active;
  size_t const pos = findRootKey(json, "active_calculators");
  if (pos == std::string::npos) {
    return active;
  }
  size_t const obj_start = json.find('{', pos);
  if (obj_start == std::string::npos) {
    return active;
  }
  size_t const obj_end = json.find('}', obj_start);
  if (obj_end == std::string::npos) {
    return active;
  }

  std::string const inner = json.substr(obj_start + 1, obj_end - obj_start - 1);
  size_t idx = 0;
  while (true) {
    size_t const key_start = inner.find('\"', idx);
    if (key_start == std::string::npos) {
      break;
    }
    size_t const key_end = inner.find('\"', key_start + 1);
    if (key_end == std::string::npos) {
      break;
    }
    std::string const key = inner.substr(key_start + 1, key_end - key_start - 1);

    size_t const colon = inner.find(':', key_end);
    if (colon == std::string::npos) {
      break;
    }
    size_t const val_start = inner.find_first_not_of(" \t\n\r", colon + 1);
    if (val_start == std::string::npos) {
      break;
    }

    bool val = false;
    if (inner.substr(val_start, 4) == "true") {
      val = true;
      idx = val_start + 4;
    } else {
      val = false;
      idx = val_start + 5;
    }
    active[key] = val;
  }
  return active;
}

} // namespace

std::filesystem::path PresetManager::presetsDirectory() {
  std::filesystem::path home;
#ifdef _WIN32
  const char *userprofile = std::getenv("USERPROFILE");
  if (userprofile)
    home = userprofile;
#else
  const char *home_env = std::getenv("HOME");
  if (home_env != nullptr) {
    home = home_env;
  }
#endif
  if (home.empty()) {
    home = std::filesystem::current_path();
  }
  return home / ".correlation" / "presets";
}

std::vector<Preset> PresetManager::loadAll() {
  std::vector<Preset> presets;
  std::filesystem::path const dir = presetsDirectory();
  if (!std::filesystem::exists(dir)) {
    return presets;
  }

  for (const auto &entry : std::filesystem::directory_iterator(dir)) {
    if (entry.is_regular_file() && entry.path().extension() == ".json") {
      std::ifstream in_stream(entry.path());
      if (in_stream.is_open()) {
        std::stringstream string_stream;
        string_stream << in_stream.rdbuf();
        try {
          Preset const preset = fromJson(string_stream.str());
          presets.push_back(preset);
        } catch (const std::exception &err) {
          std::cerr << "Failed to load preset " << entry.path().string() << ": " << err.what() << '\n';
        }
      }
    }
  }

  // Sort alphabetically
  std::ranges::sort(presets, [](const Preset &pre_a, const Preset &pre_b) { return pre_a.name < pre_b.name; });

  return presets;
}

void PresetManager::save(const Preset &preset) {
  std::filesystem::path const dir = presetsDirectory();
  std::filesystem::create_directories(dir);

  // Filename is safe name
  std::string filename = preset.name;
  std::ranges::replace_if(filename, [](char chr) { return !std::isalnum(chr) && chr != '-' && chr != '_'; }, '_');

  std::filesystem::path const filepath = dir / (filename + ".json");
  std::ofstream out(filepath);
  if (out.is_open()) {
    out << toJson(preset);
  }
}

void PresetManager::remove(const std::string &name) {
  std::filesystem::path const dir = presetsDirectory();
  if (!std::filesystem::exists(dir)) {
    return;
  }

  std::string filename = name;
  std::ranges::replace_if(filename, [](char chr) { return !std::isalnum(chr) && chr != '-' && chr != '_'; }, '_');

  std::filesystem::path const filepath = dir / (filename + ".json");
  if (std::filesystem::exists(filepath)) {
    std::filesystem::remove(filepath);
  }
}

std::string PresetManager::toJson(const Preset &preset) {
  std::string active_calcs_json;
  for (auto it = preset.options.active_calculators.begin(); it != preset.options.active_calculators.end(); ++it) {
    if (it != preset.options.active_calculators.begin()) {
      active_calcs_json += ", ";
    }
    active_calcs_json += std::format("\"{}\": {}", it->first, it->second ? "true" : "false");
  }

  return std::format(
      "{{\n"
      "  \"name\": \"{}\",\n"
      "  \"description\": \"{}\",\n"
      "  \"smoothing\": {},\n"
      "  \"use_csv\": {},\n"
      "  \"use_hdf5\": {},\n"
      "  \"use_parquet\": {},\n"
      "  \"r_max\": {:.6f},\n"
      "  \"r_bin_width\": {:.6f},\n"
      "  \"q_max\": {:.6f},\n"
      "  \"q_bin_width\": {:.6f},\n"
      "  \"r_int_max\": {:.6f},\n"
      "  \"angle_bin_width\": {:.6f},\n"
      "  \"dihedral_bin_width\": {:.6f},\n"
      "  \"max_ring_size\": {},\n"
      "  \"hyper_samples\": {},\n"
      "  \"smoothing_sigma\": {:.6f},\n"
      "  \"lef_cutoff\": {:.6f},\n"
      "  \"lef_sigma\": {:.6f},\n"
      "  \"smoothing_kernel\": {},\n"
      "  \"min_frame\": {},\n"
      "  \"max_frame\": {},\n"
      "  \"time_step\": {:.6f},\n"
      "  \"material_type\": {},\n"
      "  \"active_calculators\": {{{}}}\n"
      "}}",
      escapeJsonString(preset.name), escapeJsonString(preset.description), preset.options.smoothing ? "true" : "false",
      preset.options.use_csv ? "true" : "false", preset.options.use_hdf5 ? "true" : "false",
      preset.options.use_parquet ? "true" : "false", preset.options.r_max, preset.options.r_bin_width,
      preset.options.q_max, preset.options.q_bin_width, preset.options.r_int_max, preset.options.angle_bin_width,
      preset.options.dihedral_bin_width, preset.options.max_ring_size, preset.options.hyper_samples,
      preset.options.smoothing_sigma,
      preset.options.lef_cutoff, preset.options.lef_sigma,
      static_cast<int>(preset.options.smoothing_kernel), preset.options.min_frame, preset.options.max_frame,
      preset.options.time_step, preset.options.material_type, active_calcs_json);
}

Preset PresetManager::fromJson(const std::string &json) {
  Preset preset;
  preset.name = parseStringValue(json, "name");
  preset.description = parseStringValue(json, "description");

  preset.options.smoothing = parseBoolValue(json, "smoothing", true);
  preset.options.use_csv = parseBoolValue(json, "use_csv", true);
  preset.options.use_hdf5 = parseBoolValue(json, "use_hdf5", false);
  preset.options.use_parquet = parseBoolValue(json, "use_parquet", false);

  preset.options.r_max = parseDoubleValue(json, "r_max", AppDefaults::R_MAX);
  preset.options.r_bin_width = parseDoubleValue(json, "r_bin_width", AppDefaults::R_BIN_WIDTH);
  preset.options.q_max = parseDoubleValue(json, "q_max", AppDefaults::Q_MAX);
  preset.options.q_bin_width = parseDoubleValue(json, "q_bin_width", AppDefaults::Q_BIN_WIDTH);
  preset.options.r_int_max = parseDoubleValue(json, "r_int_max", AppDefaults::R_INT_MAX);
  preset.options.angle_bin_width = parseDoubleValue(json, "angle_bin_width", AppDefaults::ANGLE_BIN_WIDTH);
  preset.options.dihedral_bin_width = parseDoubleValue(json, "dihedral_bin_width", AppDefaults::ANGLE_BIN_WIDTH);

  preset.options.max_ring_size = static_cast<size_t>(parseIntValue(json, "max_ring_size", 8));
  preset.options.hyper_samples = static_cast<size_t>(parseIntValue(json, "hyper_samples", 10000));
  preset.options.smoothing_sigma = parseDoubleValue(json, "smoothing_sigma", AppDefaults::SMOOTHING_SIGMA);
  preset.options.lef_cutoff = parseDoubleValue(json, "lef_cutoff", AppDefaults::LEF_CUTOFF);
  preset.options.lef_sigma = parseDoubleValue(json, "lef_sigma", AppDefaults::LEF_SIGMA);
  preset.options.smoothing_kernel =
      static_cast<correlation::math::KernelType>(parseIntValue(json, "smoothing_kernel", 0));

  preset.options.min_frame = parseIntValue(json, "min_frame", 0);
  preset.options.max_frame = parseIntValue(json, "max_frame", -1);
  preset.options.time_step = parseDoubleValue(json, "time_step", AppDefaults::TIME_STEP);
  preset.options.material_type = parseIntValue(json, "material_type", 0);

  preset.options.active_calculators = parseActiveCalculators(json);

  return preset;
}

} // namespace correlation::app
