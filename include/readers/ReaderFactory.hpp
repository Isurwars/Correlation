/**
 * @file ReaderFactory.hpp
 * @brief Factory for instantiating file reader objects by format.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include "BaseReader.hpp"
#include <memory>
#include <map>
#include <string>
#include <vector>

namespace FileReader {

/**
 * @brief Registry for all file readers, enabling automatic discovery.
 */
class ReaderFactory {
public:
  static ReaderFactory &instance();

  /**
   * @brief Registers a new reader.
   * @param reader Unique pointer to the reader instance.
   * @return true if registration was successful.
   */
  bool registerReader(std::unique_ptr<BaseReader> reader);

  /**
   * @brief Finds a reader that supports the given file extension.
   * @param extension The file extension (e.g., ".car").
   * @return A pointer to the reader instance, or nullptr if not found.
   */
  BaseReader *getReaderForExtension(const std::string &extension);

  /**
   * @brief Returns all registered extensions.
   */
  std::vector<std::string> getAllExtensions() const;

  /**
   * @brief Returns all registered readers.
   */
  const std::vector<std::unique_ptr<BaseReader>> &getReaders() const;

private:
  ReaderFactory() = default;
  std::vector<std::unique_ptr<BaseReader>> readers_;
  std::map<std::string, BaseReader *> extension_map_;
};

} // namespace FileReader
