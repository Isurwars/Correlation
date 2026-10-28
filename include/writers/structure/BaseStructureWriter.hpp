/**
 * @file BaseStructureWriter.hpp
 * @brief Abstract base class and options for structure file format writers.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/Cell.hpp"

#include <string>
#include <vector>

namespace correlation::writers::structure {

/**
 * @brief Configuration options for serializing atomic structure files.
 */
struct StructureWriterOptions {
  bool fractional{true}; ///< Serialize coordinates in fractional (direct) space if supported.
  bool sort_by_element{
      true};           ///< Group atoms by chemical element (required by formats like POSCAR).
  std::string title{}; ///< Optional header/comment string for formats with title blocks.
};

/**
 * @brief Abstract interface for all atomic structure writers.
 */
class BaseStructureWriter {
public:
  BaseStructureWriter() = default;
  virtual ~BaseStructureWriter() = default;
  BaseStructureWriter(const BaseStructureWriter &) = delete;
  BaseStructureWriter &operator=(const BaseStructureWriter &) = delete;
  BaseStructureWriter(BaseStructureWriter &&) = delete;
  BaseStructureWriter &operator=(BaseStructureWriter &&) = delete;

  /**
   * @brief Returns the display name of the structure format.
   * @return A human-readable name (e.g. "VASP POSCAR", "CASTEP Cell").
   */
  [[nodiscard]] virtual std::string getName() const = 0;

  /**
   * @brief Returns a list of supported file extensions.
   * @return Vector of extensions including or excluding leading dot (e.g. {"vasp", "poscar"}).
   */
  [[nodiscard]] virtual std::vector<std::string> getExtensions() const = 0;

  /**
   * @brief Writes a simulation cell to a target file path.
   * @param[in] filepath The destination path for the output file.
   * @param[in] cell The atomic simulation cell to serialize.
   * @param[in] options Options governing coordinate conventions and element ordering.
   * @throws std::runtime_error If file cannot be written or cell data is invalid.
   */
  virtual void write(const std::string &filepath, const correlation::core::Cell &cell,
                     const StructureWriterOptions &options = {}) const = 0;
};

} // namespace correlation::writers::structure
