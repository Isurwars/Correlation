/**
 * @file StructureConverter.cpp
 * @brief Implementation of the high-level atomic structure translation pipeline.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "writers/structure/StructureConverter.hpp"
#include "core/Cell.hpp"
#include "readers/ReaderFactory.hpp"
#include "writers/structure/BaseStructureWriter.hpp"
#include "writers/structure/StructureWriterFactory.hpp"

#include <filesystem>
#include <stdexcept>

namespace correlation::writers::structure {

namespace {

correlation::readers::BaseReader *resolveSourceReader(const std::string &input_file,
                                                      const std::string &from_format) {
  auto &reader_factory = correlation::readers::ReaderFactory::instance();

  if (!from_format.empty()) {
    std::string ext = from_format;
    if (ext.front() != '.') {
      ext.insert(0, 1, '.');
    }
    auto *reader = reader_factory.getReaderForExtension(ext);
    if (reader != nullptr) {
      return reader;
    }
  }

  const std::filesystem::path input_path(input_file);
  const std::string ext = input_path.extension().string();
  return reader_factory.getReaderForExtension(
      correlation::readers::ReaderExtensionQuery{.extension = ext, .filename = input_file});
}

BaseStructureWriter *resolveDestinationWriter(const std::string &output_file,
                                              const std::string &to_format) {
  auto &writer_factory = StructureWriterFactory::instance();

  if (!to_format.empty()) {
    auto *writer = writer_factory.getWriterForFormat(to_format);
    if (writer != nullptr) {
      return writer;
    }
    writer = writer_factory.getWriterForExtension(to_format);
    if (writer != nullptr) {
      return writer;
    }
  }

  const std::filesystem::path output_path(output_file);
  std::string ext = output_path.extension().string();
  if (ext.empty()) {
    ext = output_path.filename().string();
  }
  return writer_factory.getWriterForExtension(ext);
}

void validateStructure(const correlation::core::Cell &cell, const std::string &input_file) {
  if (cell.atomCount() == 0) {
    throw std::runtime_error("Structure file '" + input_file + "' contains no atoms.");
  }
  if (cell.volume() <= 1e-9) {
    throw std::runtime_error("Structure cell in '" + input_file +
                             "' has degenerate or non-periodic volume (" +
                             std::to_string(cell.volume()) + " Å³).");
  }
}

} // namespace

void convertStructure(const std::string &input_file, const std::string &output_file,
                      const StructureConvertOptions &options) {
  if (!std::filesystem::exists(input_file)) {
    throw std::runtime_error("Input file does not exist: " + input_file);
  }

  auto *reader = resolveSourceReader(input_file, options.from_format);
  if (reader == nullptr) {
    throw std::runtime_error("Unrecognized or unsupported source structure format for: " +
                             input_file);
  }

  auto *writer = resolveDestinationWriter(output_file, options.to_format);
  if (writer == nullptr) {
    throw std::runtime_error("Unrecognized or unsupported destination structure format for: " +
                             output_file);
  }

  correlation::core::Cell const cell = reader->readStructure(input_file);
  validateStructure(cell, input_file);

  StructureWriterOptions const writer_opts{
      .fractional = options.fractional,
      .sort_by_element = true,
      .title = "Converted from " + input_file,
  };
  writer->write(output_file, cell, writer_opts);
}

} // namespace correlation::writers::structure
