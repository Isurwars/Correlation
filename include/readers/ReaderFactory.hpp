/**
 * @file ReaderFactory.hpp
 * @brief Factory for instantiating file reader objects by format.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseReader.hpp"

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace correlation::readers {

/**
 * @brief Parameters for querying a reader by extension and optional filename.
 */
struct ReaderExtensionQuery {
  std::string extension;
  std::string filename;

  ReaderExtensionQuery(const char *ext) : extension(ext) {}
  ReaderExtensionQuery(std::string ext) : extension(std::move(ext)) {}
  ReaderExtensionQuery(std::string ext, std::string file)
      : extension(std::move(ext)), filename(std::move(file)) {}
};

/**
 * @brief Registry for all file readers, enabling automatic discovery.
 */
class ReaderFactory {
public:
  /** @return Singleton instance of the ReaderFactory. */
  static ReaderFactory &instance();

  /**
   * @brief Registers a new reader.
   * @param reader Unique pointer to the reader instance.
   * @return true if registration was successful.
   */
  bool registerReader(std::unique_ptr<BaseReader> reader);

  /**
   * @brief Exception-safe template helper for static auto-registration.
   * @tparam T The specific Reader class type to instantiate.
   * @param name The human-readable name of the reader (for error logging).
   * @return true if registration succeeded, false if an exception occurred.
   */
  template <typename T> static bool registerTypeSafe(const char *name) noexcept {
    try {
      return instance().registerReader(std::make_unique<T>());
    } catch (const std::exception &e) {
      std::cerr << "[FATAL] Failed to statically register reader '" << (name != nullptr ? name : "unknown")
                << "': " << e.what() << '\n';
      return false;
    } catch (...) {
      std::cerr << "[FATAL] Failed to statically register reader '" << (name != nullptr ? name : "unknown")
                << "' due to an unknown exception." << '\n';
      return false;
    }
  }

  /**
   * @brief Finds a reader that supports the given file extension.
   * @param query The reader extension query parameters.
   * @return A pointer to the reader instance, or nullptr if not found.
   */
  BaseReader *getReaderForExtension(const ReaderExtensionQuery &query);

  /**
   * @brief Returns all registered extensions.
   * @return A vector of extension strings (e.g., {".car", ".arc"}).
   */
  std::vector<std::string> getAllExtensions() const;

  /**
   * @brief Returns all registered readers.
   * @return Constant reference to the internal vector of readers.
   */
  const std::vector<std::unique_ptr<BaseReader>> &getReaders() const;

private:
  ReaderFactory() = default;
  std::vector<std::unique_ptr<BaseReader>> readers_;  ///< Storage for reader instances.
  std::map<std::string, BaseReader *> extension_map_; ///< Lookup map for file extensions.
};

} // namespace correlation::readers
