// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "readers/CellReader.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {
void toLower(std::string &s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
}
} // namespace

namespace CellReader {

Cell read(const std::string &file_name) {
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" +
                             std::strerror(errno) + ").");
  }

  Cell tempCell;
  bool in_block = false;
  bool frac_flag = false;
  std::string current_block_type;
  std::string line;
  int lattice_row_count = 0;
  std::array<double, 6> lat;

  while (std::getline(myfile, line)) {
    /* Read line by line */
    std::stringstream line_stream(line);
    std::string token;
    line_stream >> token;
    toLower(token); // For case-insensitive matching

    if (token == "%block") {
      in_block = true;
      line_stream >> current_block_type;
      toLower(current_block_type);
      lattice_row_count = 0; // Reset for a new lattice block
      continue;
    }

    if (token == "%endblock") {
      in_block = false;
      current_block_type.clear();
      continue;
    }

    if (in_block) {
      if (current_block_type == "lattice_cart") {
        double v[3][3];
        line_stream.clear();
        line_stream.seekg(0); // Reread the full line
        if (line_stream >> v[lattice_row_count][0] >> v[lattice_row_count][1] >>
            v[lattice_row_count][2]) {
          lattice_row_count++;
          if (lattice_row_count == 3) {
            tempCell =
                Cell({v[0][0], v[0][1], v[0][2]}, {v[1][0], v[1][1], v[1][2]},
                     {v[2][0], v[2][1], v[2][2]});
          }
        }
      }

      if (current_block_type == "lattice_abc") {

        line_stream.clear();
        line_stream.seekg(0);
        if (line_stream >> lat[0 + 3 * lattice_row_count] >>
            lat[1 + 3 * lattice_row_count] >> lat[2 + 3 * lattice_row_count]) {
          lattice_row_count++;
          if (lattice_row_count == 2) {
            tempCell.setLatticeParameters(lat);
          }
        }
      }

      if (current_block_type == "positions_abs") {
        std::string element;
        double x, y, z;
        line_stream.clear();
        line_stream.seekg(0);
        if (line_stream >> element >> x >> y >> z) {
          tempCell.addAtom(element, {x, y, z});
        }
      }

      if (current_block_type == "positions_frac") {
        std::string element;
        double x, y, z;

        frac_flag = true;
        line_stream.clear();
        line_stream.seekg(0);
        if (line_stream >> element >> x >> y >> z) {
          tempCell.addAtom(element, {x, y, z});
        }
      }
    }
  }

  if (frac_flag)
    tempCell.wrapPositions();
  return tempCell;
}

} // namespace CellReader
