/**
 * @file WriterFactory.cpp
 * @brief Implementation of the writer factory.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "writers/WriterFactory.hpp"
#include <algorithm>

namespace Writer {

WriterFactory &WriterFactory::instance() {
  static WriterFactory instance_val;
  return instance_val;
}

bool WriterFactory::registerWriter(std::unique_ptr<BaseWriter> writer) {
  if (!writer) return false;

  name_map_[writer->getName()] = writer.get();

  for (const auto &ext : writer->getExtensions()) {
    std::string lower_ext = ext;
    if (lower_ext[0] != '.') lower_ext = "." + lower_ext;
    std::transform(lower_ext.begin(), lower_ext.end(), lower_ext.begin(), ::tolower);
    extension_map_[lower_ext] = writer.get();
  }

  writers_.push_back(std::move(writer));
  return true;
}

BaseWriter *WriterFactory::getWriterForExtension(const std::string &extension) {
  std::string lower_ext = extension;
  if (lower_ext[0] != '.') lower_ext = "." + lower_ext;
  std::transform(lower_ext.begin(), lower_ext.end(), lower_ext.begin(), ::tolower);

  auto it = extension_map_.find(lower_ext);
  if (it != extension_map_.end()) {
    return it->second;
  }
  return nullptr;
}

const std::vector<std::unique_ptr<BaseWriter>> &WriterFactory::getWriters() const {
  return writers_;
}

BaseWriter *WriterFactory::getWriter(const std::string &name) {
  auto it = name_map_.find(name);
  if (it != name_map_.end()) {
    return it->second;
  }
  return nullptr;
}

} // namespace Writer
