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

namespace {

std::string sniffFormatFromOutFile(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    return "";
  }

  std::string line;
  // Read first 200 lines to sniff content
  for (int i = 0; i < 200 && std::getline(file, line); ++i) {
    // Convert to uppercase for matching
    std::string uline = line;
    for (auto &chr : uline) {
      chr = static_cast<char>(std::toupper(static_cast<unsigned char>(chr)));
    }

    if (uline.contains("CELL_PARAMETERS") || uline.contains("ATOMIC_POSITIONS") || uline.contains("QUANTUM ESPRESSO") ||
        uline.contains("PWSCF") || uline.contains("&CONTROL") || uline.contains("&SYSTEM")) {
      return ".pwo";
    }
    if (uline.contains("&CELL") || uline.contains("&COORD") || uline.contains("&GLOBAL") || uline.contains("CP2K")) {
      return ".restart";
    }
  }

  return "";
}

} // namespace

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
      lower_ext.insert(0, 1, '.');
    }
    std::ranges::transform(lower_ext, lower_ext.begin(), ::tolower);
    extension_map_[lower_ext] = reader.get();
  }

  readers_.push_back(std::move(reader));
  return true;
}

BaseReader *ReaderFactory::getReaderForExtension(const ReaderExtensionQuery &query) {
  std::string lower_ext = query.extension;
  if (lower_ext.empty()) {
    return nullptr;
  }
  if (lower_ext[0] != '.') {
    lower_ext.insert(0, 1, '.');
  }
  std::ranges::transform(lower_ext, lower_ext.begin(), ::tolower);

  if ((lower_ext == ".out" || lower_ext == ".in") && !query.filename.empty()) {
    std::string sniffed = sniffFormatFromOutFile(query.filename);
    if (!sniffed.empty()) {
      auto iter = extension_map_.find(sniffed);
      if (iter != extension_map_.end()) {
        return iter->second;
      }
    }
  }

  auto iter = extension_map_.find(lower_ext);
  if (iter != extension_map_.end()) {
    return iter->second;
  }
  return nullptr;
}

std::vector<std::string> ReaderFactory::getAllExtensions() const {
  std::vector<std::string> extensions;
  extensions.reserve(extension_map_.size());
  for (const auto &[ext, unused] : extension_map_) {
    extensions.push_back(ext);
  }
  return extensions;
}

const std::vector<std::unique_ptr<BaseReader>> &ReaderFactory::getReaders() const { return readers_; }

} // namespace correlation::readers
