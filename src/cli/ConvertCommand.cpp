/**
 * @file ConvertCommand.cpp
 * @brief Implementation of the `convert` subcommand for structure format translation.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "cli/ConvertCommand.hpp"
#include "writers/structure/StructureConverter.hpp"

#include <iostream>

namespace correlation::cli {

int executeConvert(const ConvertOptions &opts) {
  try {
    correlation::writers::structure::StructureConvertOptions const convert_opts{
        .from_format = opts.from_format,
        .to_format = opts.to_format,
        .fractional = opts.fractional,
    };
    correlation::writers::structure::convertStructure(opts.input_file, opts.output_file,
                                                      convert_opts);
    if (!opts.quiet) {
      std::cout << "Successfully converted '" << opts.input_file << "' to '" << opts.output_file
                << "'\n";
    }
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
}

} // namespace correlation::cli
