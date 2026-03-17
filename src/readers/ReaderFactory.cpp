// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "readers/ReaderFactory.hpp"
#include <algorithm>

namespace FileReader {

ReaderFactory &ReaderFactory::instance() {
  static ReaderFactory instance_val;
  return instance_val;
}

bool ReaderFactory::registerReader(std::unique_ptr<BaseReader> reader) {
  if (!reader) return false;

  for (const auto &ext : reader->getExtensions()) {
    std::string lower_ext = ext;
    if (lower_ext[0] != '.') lower_ext = "." + lower_ext;
    std::transform(lower_ext.begin(), lower_ext.end(), lower_ext.begin(), ::tolower);
    extension_map_[lower_ext] = reader.get();
  }

  readers_.push_back(std::move(reader));
  return true;
}

BaseReader *ReaderFactory::getReaderForExtension(const std::string &extension) {
  std::string lower_ext = extension;
  if (lower_ext[0] != '.') lower_ext = "." + lower_ext;
  std::transform(lower_ext.begin(), lower_ext.end(), lower_ext.begin(), ::tolower);

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

const std::vector<std::unique_ptr<BaseReader>> &ReaderFactory::getReaders() const {
  return readers_;
}

} // namespace FileReader
