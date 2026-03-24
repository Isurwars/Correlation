// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include "BaseWriter.hpp"
#include <memory>
#include <map>
#include <string>
#include <vector>

namespace Writer {

/**
 * @brief Registry for all file writers, enabling automatic discovery.
 */
class WriterFactory {
public:
  static WriterFactory &instance();

  /**
   * @brief Registers a new writer.
   * @param writer Unique pointer to the writer instance.
   * @return true if registration was successful.
   */
  bool registerWriter(std::unique_ptr<BaseWriter> writer);

  /**
   * @brief Finds a writer that supports the given file extension.
   * @param extension The file extension (e.g., ".csv").
   * @return A pointer to the writer instance, or nullptr if not found.
   */
  BaseWriter *getWriterForExtension(const std::string &extension);

  /**
   * @brief Returns all registered writers.
   */
  const std::vector<std::unique_ptr<BaseWriter>> &getWriters() const;

  /**
   * @brief Gets a writer by its name.
   * @param name The name of the writer.
   * @return Pointer to the writer, or nullptr if not found.
   */
  BaseWriter *getWriter(const std::string &name);

private:
  WriterFactory() = default;
  std::vector<std::unique_ptr<BaseWriter>> writers_;
  std::map<std::string, BaseWriter *> extension_map_;
  std::map<std::string, BaseWriter *> name_map_;
};

} // namespace Writer
