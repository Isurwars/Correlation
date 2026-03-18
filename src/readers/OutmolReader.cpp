// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "readers/OutmolReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <cerrno>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <memory>
#include <functional>

#include "Cell.hpp"
#include "Trajectory.hpp"
#include "PhysicalData.hpp"

namespace FileReader {

// Automatic registration
static bool registered = ReaderFactory::instance().registerReader(
    std::make_unique<OutmolReader>()
);

Cell OutmolReader::readStructure(const std::string &filename,
                                  std::function<void(float, const std::string &)>
                                      progress_callback) {
  return read(filename, progress_callback).at(0);
}

Trajectory OutmolReader::readTrajectory(const std::string &filename,
                                         std::function<void(float, const std::string &)>
                                             progress_callback) {
  return Trajectory(read(filename, progress_callback), 1.0);
}

std::vector<Cell>
OutmolReader::read(const std::string &file_name,
                   std::function<void(float, const std::string &)> progress_callback) {
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

  myfile.seekg(0, std::ios::end);
  std::streampos file_size = myfile.tellg();
  myfile.seekg(0, std::ios::beg);
  std::streampos last_progress_pos = 0;
  size_t update_interval = file_size / 100;

  while (std::getline(myfile, line)) {
    if (progress_callback) {
      std::streampos current_pos = myfile.tellg();
      if (current_pos - last_progress_pos > update_interval) {
        float p =
            static_cast<float>(current_pos) / static_cast<float>(file_size);
        progress_callback(p, "Loading Outmol file...");
        last_progress_pos = current_pos;
      }
    }
    if (line.find("$cell vectors") != std::string::npos) {
      std::getline(myfile, line); // row 1
      std::stringstream(line) >> h1[0] >> h1[1] >> h1[2];

      std::getline(myfile, line); // row 2
      std::stringstream(line) >> h2[0] >> h2[1] >> h2[2];

      std::getline(myfile, line); // row 3
      std::stringstream(line) >> h3[0] >> h3[1] >> h3[2];

      for (int i = 0; i < 3; ++i) {
        h1[i] *= constants::BOHR_TO_ANGSTROM;
        h2[i] *= constants::BOHR_TO_ANGSTROM;
        h3[i] *= constants::BOHR_TO_ANGSTROM;
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
          tempCell.addAtom(symbol, {x * constants::BOHR_TO_ANGSTROM, y * constants::BOHR_TO_ANGSTROM,
                                    z * constants::BOHR_TO_ANGSTROM});
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
                             {x * constants::BOHR_TO_ANGSTROM, y * constants::BOHR_TO_ANGSTROM,
                              z * constants::BOHR_TO_ANGSTROM});
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

} // namespace FileReader
