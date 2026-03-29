// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include "readers/OutmolReader.hpp"
#include "Cell.hpp"
#include "Trajectory.hpp"
#include "math/Constants.hpp"
#include "readers/ReaderFactory.hpp"

#include <cerrno>
#include <cstring>
#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace FileReader {

// Automatic registration
static bool registered =
    ReaderFactory::instance().registerReader(std::make_unique<OutmolReader>());

Cell OutmolReader::readStructure(
    const std::string &filename,
    std::function<void(float, const std::string &)> progress_callback) {
  return read(filename, progress_callback).at(0);
}

Trajectory OutmolReader::readTrajectory(
    const std::string &filename,
    std::function<void(float, const std::string &)> progress_callback) {
  return Trajectory(read(filename, progress_callback), 1.0);
}

std::vector<Cell> OutmolReader::read(
    const std::string &file_name,
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

  std::stringstream ss;

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
      if (!std::getline(myfile, line))
        break; // row 1
      ss.clear();
      ss.str(line);
      ss >> h1[0] >> h1[1] >> h1[2];

      if (!std::getline(myfile, line))
        break; // row 2
      ss.clear();
      ss.str(line);
      ss >> h2[0] >> h2[1] >> h2[2];

      if (!std::getline(myfile, line))
        break; // row 3
      ss.clear();
      ss.str(line);
      ss >> h3[0] >> h3[1] >> h3[2];

      for (int i = 0; i < 3; ++i) {
        h1[i] *= correlation::math::constants::bohr_to_angstrom;
        h2[i] *= correlation::math::constants::bohr_to_angstrom;
        h3[i] *= correlation::math::constants::bohr_to_angstrom;
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
        ss.clear();
        ss.str(line);
        std::string symbol;
        double x, y, z;
        if (ss >> symbol >> x >> y >> z) {
          tempCell.addAtom(
              symbol, {x * correlation::math::constants::bohr_to_angstrom,
                       y * correlation::math::constants::bohr_to_angstrom,
                       z * correlation::math::constants::bohr_to_angstrom});
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
        ss.clear();
        ss.str(line);
        std::string df, symbol;
        double x, y, z;
        if (ss >> df >> symbol >> x >> y >> z) {
          if (df == "df") {
            tempCell.addAtom(
                symbol, {x * correlation::math::constants::bohr_to_angstrom,
                         y * correlation::math::constants::bohr_to_angstrom,
                         z * correlation::math::constants::bohr_to_angstrom});
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
