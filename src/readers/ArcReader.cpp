// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "readers/ArcReader.hpp"

#include <array>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace ArcReader {

std::vector<Cell> read(const std::string &file_name) {
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" +
                             std::strerror(errno) + ").");
  }

  std::vector<Cell> frames;
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

    if (first_token == "end") {
      if (!tempCell.isEmpty()) {
        frames.push_back(std::move(tempCell));
        tempCell = Cell(); // Reset for next frame
      }
      continue;
    }

    if (first_token == "PBC") {
      std::array<double, 6> lat;
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

    if (first_token == "PBC=ON") {
      continue;
    }

    // Check if it is a single token line (Energy)
    std::string second_token;
    if (!(line_stream >> second_token)) {
      try {
        double energy = std::stod(first_token);
        tempCell.setEnergy(energy);
        continue;
      } catch (...) {
        // Not a number, move on
      }
    }

    // Attempt to parse atom
    // Reset stream to start of line
    line_stream.clear();
    line_stream.seekg(0);

    std::string u1, u5, u6, u7, element;
    double x, y, z;

    if (line_stream >> u1 >> x >> y >> z >> u5 >> u6 >> u7 >> element) {
      tempCell.addAtom(element, {x, y, z});
    }
  }

  return frames;
}

} // namespace ArcReader
