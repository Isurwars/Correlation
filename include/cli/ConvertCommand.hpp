/**
 * @file ConvertCommand.hpp
 * @brief Subcommand handler for converting between atomic structure formats.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include <string>

namespace correlation::cli {

/**
 * @brief Options for the `convert` CLI subcommand.
 */
struct ConvertOptions {
  std::string input_file;  ///< Source structure file path.
  std::string output_file; ///< Destination structure file path.
  std::string from_format; ///< Explicit source format override.
  std::string to_format;   ///< Explicit destination format override.
  bool fractional{true};   ///< Coordinate convention (true = fractional/direct, false = cartesian).
  bool quiet{false};       ///< Suppress console progress messages.
};

/**
 * @brief Executes the structure format conversion workflow.
 * @param[in] opts Parsed convert options.
 * @return 0 on success, non-zero on failure.
 */
int executeConvert(const ConvertOptions &opts);

} // namespace correlation::cli
