/**
 * @file StructureWriterFactory.cpp
 * @brief Implementation of the structure writer factory registry.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "writers/structure/StructureWriterFactory.hpp"

#include <algorithm>
#include <cctype>

namespace correlation::writers::structure {

namespace {

std::string normalizeString(std::string_view sv) {
  std::string result(sv);
  // Strip leading dot if present
  if (!result.empty() && result.front() == '.') {
    result.erase(0, 1);
  }
  std::ranges::transform(result, result.begin(),
                         [](unsigned char chr) { return static_cast<char>(std::tolower(chr)); });
  return result;
}

} // namespace

StructureWriterFactory &StructureWriterFactory::instance() {
  static StructureWriterFactory instance_val;
  return instance_val;
}

bool StructureWriterFactory::registerWriter(std::unique_ptr<BaseStructureWriter> writer) {
  if (!writer) {
    return false;
  }

  BaseStructureWriter *raw_ptr = writer.get();

  // Register all extensions
  for (const auto &ext : writer->getExtensions()) {
    std::string const norm_ext = normalizeString(ext);
    if (!norm_ext.empty()) {
      extension_map_[norm_ext] = raw_ptr;
    }
  }

  // Register display name and lowercase format key
  std::string const norm_name = normalizeString(writer->getName());
  if (!norm_name.empty()) {
    format_map_[norm_name] = raw_ptr;
  }

  // Also register aliases from extensions as format keys
  for (const auto &ext : writer->getExtensions()) {
    std::string const norm_ext = normalizeString(ext);
    if (!norm_ext.empty()) {
      format_map_[norm_ext] = raw_ptr;
    }
  }

  writers_.push_back(std::move(writer));
  return true;
}

BaseStructureWriter *
StructureWriterFactory::getWriterForExtension(std::string_view extension) const {
  std::string const norm = normalizeString(extension);
  auto iter = extension_map_.find(norm);
  if (iter != extension_map_.end()) {
    return iter->second;
  }
  return nullptr;
}

BaseStructureWriter *StructureWriterFactory::getWriterForFormat(std::string_view format) const {
  std::string const norm = normalizeString(format);
  auto iter = format_map_.find(norm);
  if (iter != format_map_.end()) {
    return iter->second;
  }
  return nullptr;
}

std::vector<std::string> StructureWriterFactory::getSupportedFormats() const {
  std::vector<std::string> formats;
  formats.reserve(writers_.size());
  for (const auto &writer : writers_) {
    formats.push_back(writer->getName());
  }
  return formats;
}

std::vector<std::string> StructureWriterFactory::getSupportedExtensions() const {
  std::vector<std::string> extensions;
  extensions.reserve(extension_map_.size());
  for (const auto &[ext, _] : extension_map_) {
    extensions.push_back(ext);
  }
  return extensions;
}

} // namespace correlation::writers::structure
