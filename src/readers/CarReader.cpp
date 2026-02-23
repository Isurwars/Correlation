// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "readers/CarReader.hpp"

#include <array>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace CarReader {

Cell read(const std::string &file_name) {
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" +
                             std::strerror(errno) + ").");
  }

  Cell tempCell;
  std::string line;

  while (std::getline(myfile, line)) {

    // Ignore empty lines or comment lines
    if (line.empty() || line[0] == '!') {
      continue;
    }

    std::stringstream line_stream(line);
    std::string first_token;
    line_stream >> first_token;

    if (first_token == "PBC") {
      std::array<double, 6> lat;
      // The token "PBC" is consumed, so we read the 6 numbers that follow.
      if (line_stream >> lat[0] >> lat[1] >> lat[2] >> lat[3] >> lat[4] >>
          lat[5]) {
        tempCell.setLatticeParameters(lat);
      }
      continue;
    }

    if (first_token == "PBC=OFF") {
      std::array<double, 6> lat = {100.0, 100.0, 100.0, 90.0, 90.0, 90.0};
      tempCell.setLatticeParameters(lat);
      continue;
    }

    // Reset the stream to parse the full line as an atom entry.
    line_stream.clear();
    line_stream.seekg(0);

    // Declare variables for all 8 columns we need to read.
    std::string u1, u5, u6, u7, element;
    double x, y, z;

    // Read exactly 8 columns to get the element.
    if (line_stream >> u1 >> x >> y >> z >> u5 >> u6 >> u7 >> element) {
      tempCell.addAtom(element, {x, y, z});
    }
  }

  return tempCell;
}

} // namespace CarReader
