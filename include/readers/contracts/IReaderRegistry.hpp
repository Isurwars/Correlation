/**
 * @file IReaderRegistry.hpp
 * @brief Abstract interface for querying file readers by format and extension.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace correlation::readers {

class BaseReader;

/**
 * @brief Parameters for querying a reader by extension and optional filename for content sniffing.
 */
struct ReaderExtensionQuery {
  std::string_view extension;
  std::string_view filename;
};

/**
 * @brief Abstract interface defining lookup capabilities for file reader discovery and instantiation.
 */
class IReaderRegistry {
public:
  IReaderRegistry() = default;
  virtual ~IReaderRegistry() = default;

  IReaderRegistry(const IReaderRegistry &) = delete;
  IReaderRegistry &operator=(const IReaderRegistry &) = delete;
  IReaderRegistry(IReaderRegistry &&) = delete;
  IReaderRegistry &operator=(IReaderRegistry &&) = delete;

  /**
   * @brief Finds a reader matching the given file extension and optional filename sniffing query.
   * @param[in] query Structure containing extension and optional file path.
   * @return Pointer to a registered BaseReader, or nullptr if unsupported.
   */
  [[nodiscard]] virtual BaseReader *
  getReaderForExtension(const ReaderExtensionQuery &query) const = 0;

  /**
   * @brief Finds a reader matching the given file extension.
   * @param[in] extension File extension (e.g., ".xyz", ".poscar").
   * @return Pointer to a registered BaseReader, or nullptr if unsupported.
   */
  [[nodiscard]] virtual BaseReader *getReaderForExtension(std::string_view extension) const = 0;

  /**
   * @brief Returns a list of all file extensions supported by the registry.
   * @return Vector of file extension strings.
   */
  [[nodiscard]] virtual std::vector<std::string> getAllExtensions() const = 0;

  /**
   * @brief Returns all registered reader instances.
   * @return Constant reference to vector of unique pointers to BaseReader.
   */
  [[nodiscard]] virtual const std::vector<std::unique_ptr<BaseReader>> &getReaders() const = 0;
};

} // namespace correlation::readers
