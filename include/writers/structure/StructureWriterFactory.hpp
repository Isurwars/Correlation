/**
 * @file StructureWriterFactory.hpp
 * @brief Factory registry for discovering and instantiating atomic structure writers.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "writers/structure/BaseStructureWriter.hpp"

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace correlation::writers::structure {

/**
 * @brief Factory registry for all atomic structure writers.
 */
class StructureWriterFactory {
public:
  /**
   * @brief Returns the global singleton factory instance.
   * @return Reference to StructureWriterFactory.
   */
  static StructureWriterFactory &instance();

  /**
   * @brief Registers a new structure writer instance.
   * @param[in] writer Owning pointer to the writer instance.
   * @return true if successfully registered, false on collision.
   */
  bool registerWriter(std::unique_ptr<BaseStructureWriter> writer);

  /**
   * @brief Static helper to register a writer type safely without throwing during static init.
   * @tparam T Concrete writer type deriving from BaseStructureWriter.
   * @param[in] name Debug name of the writer for error diagnostics.
   * @return true if registration succeeded.
   */
  template <typename T> static bool registerTypeSafe(const char *name) noexcept {
    try {
      return instance().registerWriter(std::make_unique<T>());
    } catch (const std::exception &e) {
      std::cerr << "[FATAL] Failed to statically register structure writer '"
                << (name != nullptr ? name : "unknown") << "': " << e.what() << '\n';
      return false;
    } catch (...) {
      std::cerr << "[FATAL] Failed to statically register structure writer '"
                << (name != nullptr ? name : "unknown") << "' due to unknown exception.\n";
      return false;
    }
  }

  /**
   * @brief Finds a writer matching a file extension.
   * @param[in] extension File extension (e.g. ".poscar", "vasp", ".cell").
   * @return Non-owning pointer to writer, or nullptr if unsupported.
   */
  [[nodiscard]] BaseStructureWriter *getWriterForExtension(std::string_view extension) const;

  /**
   * @brief Finds a writer matching a format identifier.
   * @param[in] format Format name (e.g. "poscar", "vasp", "cif", "cell", "car").
   * @return Non-owning pointer to writer, or nullptr if unsupported.
   */
  [[nodiscard]] BaseStructureWriter *getWriterForFormat(std::string_view format) const;

  /**
   * @brief Returns names of all registered structure formats.
   * @return Vector of format names.
   */
  [[nodiscard]] std::vector<std::string> getSupportedFormats() const;

  /**
   * @brief Returns list of all registered file extensions.
   * @return Vector of extensions.
   */
  [[nodiscard]] std::vector<std::string> getSupportedExtensions() const;

private:
  StructureWriterFactory() = default;
  ~StructureWriterFactory() = default;
  StructureWriterFactory(const StructureWriterFactory &) = delete;
  StructureWriterFactory &operator=(const StructureWriterFactory &) = delete;
  StructureWriterFactory(StructureWriterFactory &&) = delete;
  StructureWriterFactory &operator=(StructureWriterFactory &&) = delete;

  std::vector<std::unique_ptr<BaseStructureWriter>> writers_;
  std::map<std::string, BaseStructureWriter *> extension_map_;
  std::map<std::string, BaseStructureWriter *> format_map_;
};

} // namespace correlation::writers::structure
