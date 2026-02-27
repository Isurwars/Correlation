// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "readers/OutmolReader.hpp"

#include <cerrno>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace OutmolReader {

// Conversion factor from Bohr (Hartree atomic units) to Angstroms
constexpr double BOHR_TO_ANGSTROM = 0.529177210903;

std::vector<Cell> read(const std::string &file_name) {
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" +
                             std::strerror(errno) + ").");
  }

  std::vector<Cell> frames;
  std::string line;

  std::array<double, 3> h1 = {0.0, 0.0, 0.0};
  std::array<double, 3> h2 = {0.0, 0.0, 0.0};
  std::array<double, 3> h3 = {0.0, 0.0, 0.0};
  bool cell_parsed = false;

  while (std::getline(myfile, line)) {
    if (line.find("$cell vectors") != std::string::npos) {
      std::getline(myfile, line); // row 1
      std::stringstream(line) >> h1[0] >> h1[1] >> h1[2];

      std::getline(myfile, line); // row 2
      std::stringstream(line) >> h2[0] >> h2[1] >> h2[2];

      std::getline(myfile, line); // row 3
      std::stringstream(line) >> h3[0] >> h3[1] >> h3[2];

      for (int i = 0; i < 3; ++i) {
        h1[i] *= BOHR_TO_ANGSTROM;
        h2[i] *= BOHR_TO_ANGSTROM;
        h3[i] *= BOHR_TO_ANGSTROM;
      }
      cell_parsed = true;
      continue;
    }

    if (line.find("$coordinates") != std::string::npos) {
      Cell tempCell({h1[0], h1[1], h1[2]}, {h2[0], h2[1], h2[2]},
                    {h3[0], h3[1], h3[2]});
      while (std::getline(myfile, line)) {
        if (line.find("$end") != std::string::npos) {
          break;
        }
        std::stringstream ss(line);
        std::string symbol;
        double x, y, z;
        if (ss >> symbol >> x >> y >> z) {
          tempCell.addAtom(symbol, {x * BOHR_TO_ANGSTROM, y * BOHR_TO_ANGSTROM,
                                    z * BOHR_TO_ANGSTROM});
        }
      }
      if (!tempCell.isEmpty()) {
        frames.push_back(std::move(tempCell));
      }
      continue;
    }

    if (line.find("ATOMIC  COORDINATES (au)") != std::string::npos) {
      std::getline(myfile, line); // Next line is header: df x y z ...
      Cell tempCell({h1[0], h1[1], h1[2]}, {h2[0], h2[1], h2[2]},
                    {h3[0], h3[1], h3[2]});
      while (std::getline(myfile, line)) {
        // Look for the "df" prefix
        if (line.find("df") == std::string::npos) {
          break;
        }
        std::stringstream ss(line);
        std::string df, symbol;
        double x, y, z;
        if (ss >> df >> symbol >> x >> y >> z) {
          if (df == "df") {
            tempCell.addAtom(symbol,
                             {x * BOHR_TO_ANGSTROM, y * BOHR_TO_ANGSTROM,
                              z * BOHR_TO_ANGSTROM});
          }
        }
      }
      if (!tempCell.isEmpty()) {
        frames.push_back(std::move(tempCell));
      }
      continue;
    }
  }

  return frames;
}

} // namespace OutmolReader
