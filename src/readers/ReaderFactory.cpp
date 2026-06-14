/**
 * @file ReaderFactory.cpp
 * @brief Implementation of the reader factory.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/ReaderFactory.hpp"

#include <algorithm>
#include <fstream>

namespace correlation::readers {

ReaderFactory &ReaderFactory::instance() {
  static ReaderFactory instance_val;
  return instance_val;
}

bool ReaderFactory::registerReader(std::unique_ptr<BaseReader> reader) {
  if (!reader) {
    return false;
}

  for (const auto &ext : reader->getExtensions()) {
    std::string lower_ext = ext;
    if (lower_ext.empty()) {
      continue;
}
    if (lower_ext[0] != '.') {
      lower_ext = "." + lower_ext;
}
    std::ranges::transform(lower_ext, lower_ext.begin(), ::tolower);
    extension_map_[lower_ext] = reader.get();
  }

  readers_.push_back(std::move(reader));
  return true;
}

BaseReader *ReaderFactory::getReaderForExtension(const std::string &extension, const std::string &filename) {
  std::string lower_ext = extension;
  if (lower_ext.empty()) {
    return nullptr;
}
  if (lower_ext[0] != '.') {
    lower_ext = "." + lower_ext;
}
  std::ranges::transform(lower_ext, lower_ext.begin(), ::tolower);

  if ((lower_ext == ".out" || lower_ext == ".in") && !filename.empty()) {
    std::ifstream file(filename);
    if (file.is_open()) {
      std::string line;
      // Read first 200 lines to sniff content
      for (int i = 0; i < 200 && std::getline(file, line); ++i) {
        // Convert to uppercase for matching
        std::string uline = line;
        for (auto &c : uline) {
          c = toupper(c);
}

        if (uline.contains("CELL_PARAMETERS") || uline.contains("ATOMIC_POSITIONS") ||
            uline.contains("QUANTUM ESPRESSO") || uline.contains("PWSCF") ||
            uline.contains("&CONTROL") || uline.contains("&SYSTEM")) {
          auto it = extension_map_.find(".pwo");
          if (it != extension_map_.end()) {
            return it->second;
          }
        }
        if (uline.contains("&CELL") || uline.contains("&COORD") ||
            uline.contains("&GLOBAL") || uline.contains("CP2K")) {
          auto it = extension_map_.find(".restart");
          if (it != extension_map_.end()) {
            return it->second;
          }
        }
      }
    }
  }

  auto it = extension_map_.find(lower_ext);
  if (it != extension_map_.end()) {
    return it->second;
  }
  return nullptr;
}

std::vector<std::string> ReaderFactory::getAllExtensions() const {
  std::vector<std::string> extensions;
  extensions.reserve(extension_map_.size());
  for (const auto &[ext, _] : extension_map_) {
    extensions.push_back(ext);
  }
  return extensions;
}

const std::vector<std::unique_ptr<BaseReader>> &ReaderFactory::getReaders() const { return readers_; }

} // namespace correlation::readers
