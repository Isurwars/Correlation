/**
 * @file WriterFactory.hpp
 * @brief Factory for instantiating file writer objects by format.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "BaseWriter.hpp"

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace correlation::writers {

/**
 * @brief Registry for all file writers, enabling automatic discovery.
 */
class WriterFactory {
public:
  /** @return Singleton instance of the WriterFactory. */
  static WriterFactory &instance();

  /**
   * @brief Registers a new writer.
   * @param writer Unique pointer to the writer instance.
   * @return true if registration was successful.
   */
  bool registerWriter(std::unique_ptr<BaseWriter> writer);

  /**
   * @brief Exception-safe template helper for static auto-registration.
   * @tparam T The specific Writer class type to instantiate.
   * @param name The human-readable name of the writer (for error logging).
   * @return true if registration succeeded, false if an exception occurred.
   */
  template <typename T> static bool registerTypeSafe(const char *name) noexcept {
    try {
      return instance().registerWriter(std::make_unique<T>());
    } catch (const std::exception &e) {
      std::cerr << "[FATAL] Failed to statically register writer '" << (name != nullptr ? name : "unknown")
                << "': " << e.what() << '\n';
      return false;
    } catch (...) {
      std::cerr << "[FATAL] Failed to statically register writer '" << (name != nullptr ? name : "unknown")
                << "' due to an unknown exception." << '\n';
      return false;
    }
  }

  /**
   * @brief Finds a writer that supports the given file extension.
   * @param extension The file extension (e.g., ".csv").
   * @return A pointer to the writer instance, or nullptr if not found.
   */
  BaseWriter *getWriterForExtension(const std::string &extension);

  /**
   * @brief Returns all registered writers.
   * @return Constant reference to the internal vector of writers.
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
  std::vector<std::unique_ptr<BaseWriter>> writers_;  ///< Storage for writer instances.
  std::map<std::string, BaseWriter *> extension_map_; ///< Lookup map for file extensions.
  std::map<std::string, BaseWriter *> name_map_;      ///< Lookup map for format names.
};

} // namespace correlation::writers
