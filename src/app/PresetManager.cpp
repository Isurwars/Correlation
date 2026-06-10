// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "app/PresetManager.hpp"
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <format>
#include <algorithm>

namespace correlation::app {

namespace {

std::string escapeJsonString(const std::string &s) {
  std::string out;
  for (char c : s) {
    if (c == '"') out += "\\\"";
    else if (c == '\\') out += "\\\\";
    else if (c == '\n') out += "\\n";
    else if (c == '\r') out += "\\r";
    else if (c == '\t') out += "\\t";
    else out += c;
  }
  return out;
}

size_t findRootKey(const std::string &json, const std::string &key) {
  bool in_string = false;
  int brace_level = 0;
  int bracket_level = 0;
  for (size_t i = 0; i < json.size(); ++i) {
    char c = json[i];
    if (c == '\\' && in_string) {
      i++; // skip next char
      continue;
    }
    if (c == '"') {
      in_string = !in_string;
      if (in_string) {
        if (brace_level == 1 && bracket_level == 0) {
          if (i + 1 + key.size() < json.size() &&
              json.compare(i + 1, key.size(), key) == 0 &&
              json[i + 1 + key.size()] == '"') {
            // Verify it is followed by a colon (skipping whitespace)
            size_t colon_pos = json.find_first_not_of(" \t\n\r", i + 1 + key.size() + 1);
            if (colon_pos != std::string::npos && json[colon_pos] == ':') {
              return i;
            }
          }
        }
      }
      continue;
    }
    if (!in_string) {
      if (c == '{') brace_level++;
      else if (c == '}') brace_level--;
      else if (c == '[') bracket_level++;
      else if (c == ']') bracket_level--;
    }
  }
  return std::string::npos;
}

std::string parseStringValue(const std::string &json, const std::string &key) {
  size_t pos = findRootKey(json, key);
  if (pos == std::string::npos) return "";
  pos = json.find(":", pos);
  if (pos == std::string::npos) return "";
  size_t start = json.find("\"", pos);
  if (start == std::string::npos) return "";
  
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
  if (end == std::string::npos) return "";
  std::string val = json.substr(start + 1, end - start - 1);
  
  // Unescape
  std::string unescaped;
  for (size_t i = 0; i < val.size(); ++i) {
    if (val[i] == '\\' && i + 1 < val.size()) {
      char next = val[i + 1];
      if (next == '"') { unescaped += '"'; i++; }
      else if (next == '\\') { unescaped += '\\'; i++; }
      else if (next == 'n') { unescaped += '\n'; i++; }
      else if (next == 'r') { unescaped += '\r'; i++; }
      else if (next == 't') { unescaped += '\t'; i++; }
      else { unescaped += val[i]; }
    } else {
      unescaped += val[i];
    }
  }
  return unescaped;
}

double parseDoubleValue(const std::string &json, const std::string &key, double fallback) {
  size_t pos = findRootKey(json, key);
  if (pos == std::string::npos) return fallback;
  pos = json.find(":", pos);
  if (pos == std::string::npos) return fallback;
  size_t start = json.find_first_not_of(" \t\n\r:", pos + 1);
  if (start == std::string::npos) return fallback;
  size_t end = json.find_first_of(",}\n\r", start);
  if (end == std::string::npos) end = json.size();
  try {
    return std::stod(json.substr(start, end - start));
  } catch (...) {
    return fallback;
  }
}

int parseIntValue(const std::string &json, const std::string &key, int fallback) {
  size_t pos = findRootKey(json, key);
  if (pos == std::string::npos) return fallback;
  pos = json.find(":", pos);
  if (pos == std::string::npos) return fallback;
  size_t start = json.find_first_not_of(" \t\n\r:", pos + 1);
  if (start == std::string::npos) return fallback;
  size_t end = json.find_first_of(",}\n\r", start);
  if (end == std::string::npos) end = json.size();
  try {
    return std::stoi(json.substr(start, end - start));
  } catch (...) {
    return fallback;
  }
}

bool parseBoolValue(const std::string &json, const std::string &key, bool fallback) {
  size_t pos = findRootKey(json, key);
  if (pos == std::string::npos) return fallback;
  pos = json.find(":", pos);
  if (pos == std::string::npos) return fallback;
  size_t start = json.find_first_not_of(" \t\n\r:", pos + 1);
  if (start == std::string::npos) return fallback;
  std::string val = json.substr(start, 5);
  if (val.starts_with("true")) return true;
  if (val.starts_with("false")) return false;
  return fallback;
}

std::map<std::string, bool> parseActiveCalculators(const std::string &json) {
  std::map<std::string, bool> active;
  size_t pos = findRootKey(json, "active_calculators");
  if (pos == std::string::npos) return active;
  size_t obj_start = json.find("{", pos);
  if (obj_start == std::string::npos) return active;
  size_t obj_end = json.find("}", obj_start);
  if (obj_end == std::string::npos) return active;
  
  std::string inner = json.substr(obj_start + 1, obj_end - obj_start - 1);
  size_t idx = 0;
  while (true) {
    size_t key_start = inner.find("\"", idx);
    if (key_start == std::string::npos) break;
    size_t key_end = inner.find("\"", key_start + 1);
    if (key_end == std::string::npos) break;
    std::string key = inner.substr(key_start + 1, key_end - key_start - 1);
    
    size_t colon = inner.find(":", key_end);
    if (colon == std::string::npos) break;
    size_t val_start = inner.find_first_not_of(" \t\n\r", colon + 1);
    if (val_start == std::string::npos) break;
    
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
  if (userprofile) home = userprofile;
#else
  const char *home_env = std::getenv("HOME");
  if (home_env) home = home_env;
#endif
  if (home.empty()) {
    home = std::filesystem::current_path();
  }
  return home / ".correlation" / "presets";
}

std::vector<Preset> PresetManager::loadAll() {
  std::vector<Preset> presets;
  std::filesystem::path dir = presetsDirectory();
  if (!std::filesystem::exists(dir)) {
    return presets;
  }

  for (const auto &entry : std::filesystem::directory_iterator(dir)) {
    if (entry.is_regular_file() && entry.path().extension() == ".json") {
      std::ifstream in(entry.path());
      if (in.is_open()) {
        std::stringstream ss;
        ss << in.rdbuf();
        try {
          Preset p = fromJson(ss.str());
          presets.push_back(p);
        } catch (...) {
          // Skip corrupt presets
        }
      }
    }
  }

  // Sort alphabetically
  std::sort(presets.begin(), presets.end(), [](const Preset &a, const Preset &b) {
    return a.name < b.name;
  });

  return presets;
}

void PresetManager::save(const Preset &preset) {
  std::filesystem::path dir = presetsDirectory();
  std::filesystem::create_directories(dir);

  // Filename is safe name
  std::string filename = preset.name;
  std::replace_if(filename.begin(), filename.end(), [](char c) {
    return !std::isalnum(c) && c != '-' && c != '_';
  }, '_');

  std::filesystem::path filepath = dir / (filename + ".json");
  std::ofstream out(filepath);
  if (out.is_open()) {
    out << toJson(preset);
  }
}

void PresetManager::remove(const std::string &name) {
  std::filesystem::path dir = presetsDirectory();
  if (!std::filesystem::exists(dir)) return;

  std::string filename = name;
  std::replace_if(filename.begin(), filename.end(), [](char c) {
    return !std::isalnum(c) && c != '-' && c != '_';
  }, '_');

  std::filesystem::path filepath = dir / (filename + ".json");
  if (std::filesystem::exists(filepath)) {
    std::filesystem::remove(filepath);
  }
}

std::string PresetManager::toJson(const Preset &preset) {
  std::string active_calcs_json = "";
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
    "  \"smoothing_sigma\": {:.6f},\n"
    "  \"smoothing_kernel\": {},\n"
    "  \"min_frame\": {},\n"
    "  \"max_frame\": {},\n"
    "  \"time_step\": {:.6f},\n"
    "  \"material_type\": {},\n"
    "  \"active_calculators\": {{{}}}\n"
    "}}",
    escapeJsonString(preset.name),
    escapeJsonString(preset.description),
    preset.options.smoothing ? "true" : "false",
    preset.options.use_csv ? "true" : "false",
    preset.options.use_hdf5 ? "true" : "false",
    preset.options.use_parquet ? "true" : "false",
    preset.options.r_max,
    preset.options.r_bin_width,
    preset.options.q_max,
    preset.options.q_bin_width,
    preset.options.r_int_max,
    preset.options.angle_bin_width,
    preset.options.dihedral_bin_width,
    preset.options.max_ring_size,
    preset.options.smoothing_sigma,
    static_cast<int>(preset.options.smoothing_kernel),
    preset.options.min_frame,
    preset.options.max_frame,
    preset.options.time_step,
    preset.options.material_type,
    active_calcs_json
  );
}

Preset PresetManager::fromJson(const std::string &json) {
  Preset p;
  p.name = parseStringValue(json, "name");
  p.description = parseStringValue(json, "description");
  
  p.options.smoothing = parseBoolValue(json, "smoothing", true);
  p.options.use_csv = parseBoolValue(json, "use_csv", true);
  p.options.use_hdf5 = parseBoolValue(json, "use_hdf5", false);
  p.options.use_parquet = parseBoolValue(json, "use_parquet", false);
  
  p.options.r_max = parseDoubleValue(json, "r_max", AppDefaults::R_MAX);
  p.options.r_bin_width = parseDoubleValue(json, "r_bin_width", AppDefaults::R_BIN_WIDTH);
  p.options.q_max = parseDoubleValue(json, "q_max", AppDefaults::Q_MAX);
  p.options.q_bin_width = parseDoubleValue(json, "q_bin_width", AppDefaults::Q_BIN_WIDTH);
  p.options.r_int_max = parseDoubleValue(json, "r_int_max", AppDefaults::R_INT_MAX);
  p.options.angle_bin_width = parseDoubleValue(json, "angle_bin_width", AppDefaults::ANGLE_BIN_WIDTH);
  p.options.dihedral_bin_width = parseDoubleValue(json, "dihedral_bin_width", AppDefaults::ANGLE_BIN_WIDTH);
  
  p.options.max_ring_size = static_cast<size_t>(parseIntValue(json, "max_ring_size", 8));
  p.options.smoothing_sigma = parseDoubleValue(json, "smoothing_sigma", AppDefaults::SMOOTHING_SIGMA);
  p.options.smoothing_kernel = static_cast<correlation::math::KernelType>(parseIntValue(json, "smoothing_kernel", 0));
  
  p.options.min_frame = parseIntValue(json, "min_frame", 0);
  p.options.max_frame = parseIntValue(json, "max_frame", -1);
  p.options.time_step = parseDoubleValue(json, "time_step", AppDefaults::TIME_STEP);
  p.options.material_type = parseIntValue(json, "material_type", 0);
  
  p.options.active_calculators = parseActiveCalculators(json);
  
  return p;
}

} // namespace correlation::app
