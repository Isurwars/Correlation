// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "readers/OnetepDatReader.hpp"

#include <fstream>
#include <stdexcept>

namespace OnetepDatReader {

Cell read(const std::string &file_name) {
  /*
   * This is a Stub for reading a ONETEP file and return and empty cell.
   */
  std::ifstream myfile(file_name);
  // std::string   line;
  Cell tempCell;
  // std::smatch match;

  return tempCell;
}

} // namespace OnetepDatReader
