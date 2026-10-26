/**
 * @file IStructureReader.hpp
 * @brief Abstract interface for parsing static atomic structures into a Cell.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/Cell.hpp"

#include <functional>
#include <string>

namespace correlation::readers {

/**
 * @brief Progress report callback function type.
 * Invoked with normalized progress in [0.0, 1.0] and a status message.
 */
using ProgressCallback = std::function<void(float, const std::string &)>;

/**
 * @brief Segregated contract for file format parsers that read individual snapshot structures.
 */
class IStructureReader {
public:
  IStructureReader() = default;
  virtual ~IStructureReader() = default;

  IStructureReader(const IStructureReader &) = delete;
  IStructureReader &operator=(const IStructureReader &) = delete;
  IStructureReader(IStructureReader &&) = delete;
  IStructureReader &operator=(IStructureReader &&) = delete;

  /**
   * @brief Reads an atomic structure from a file.
   * @param[in] filename Path to the file to be read.
   * @param[in] progress_callback Optional callback to report progress [0.0, 1.0].
   * @return A Cell object containing the loaded structure.
   */
  [[nodiscard]] virtual correlation::core::Cell
  readStructure(const std::string &filename, ProgressCallback progress_callback = nullptr) = 0;
};

} // namespace correlation::readers
