/**
 * @file StructureConverter.hpp
 * @brief High-level atomic structure translation pipeline.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include <string>

namespace correlation::writers::structure {

/**
 * @brief Options for converting between atomic structure files.
 */
struct StructureConvertOptions {
  std::string from_format{}; ///< Optional explicit source format name or extension.
  std::string to_format{};   ///< Optional explicit destination format name or extension.
  bool fractional{true};     ///< Coordinate convention (fractional vs Cartesian).
};

/**
 * @brief Converts an atomic structure from one format to another.
 *
 * Reads the input structure, performs validation on periodicity and atom integrity,
 * and serializes it to the destination format.
 *
 * @param[in] input_file Path to the input file.
 * @param[in] output_file Path to the output file.
 * @param[in] options Translation options.
 * @throws std::runtime_error If input file is missing, format is unsupported, or writing fails.
 */
void convertStructure(const std::string &input_file, const std::string &output_file,
                      const StructureConvertOptions &options = {});

} // namespace correlation::writers::structure
