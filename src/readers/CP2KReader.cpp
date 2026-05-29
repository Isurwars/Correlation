/**
 * @file CP2KReader.cpp
 * @brief Implementation of the CP2K file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "readers/CP2KReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace correlation::readers {

namespace {
bool registered = ReaderFactory::instance().registerReader(std::make_unique<CP2KReader>());
}

correlation::core::Cell CP2KReader::readStructure(const std::string &filename,
                                                  std::function<void(float, const std::string &)> progress_callback) {

  auto traj = readTrajectory(filename, progress_callback);
  if (traj.getFrameCount() == 0) {
    throw std::runtime_error("No structure found in CP2K file: " + filename);
  }
  return traj.getFrames()[traj.getFrameCount() - 1];
}

correlation::core::Trajectory
CP2KReader::readTrajectory(const std::string &filename,
                           std::function<void(float, const std::string &)> progress_callback) {

  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  std::vector<correlation::core::Cell> frames;
  std::string line;
  correlation::core::Cell current_cell;
  bool parsing_coords = false;
  bool parsing_cell = false;
  bool has_box = false;

  while (std::getline(file, line)) {
    line.erase(0, line.find_first_not_of(" \t\r\n"));
    line.erase(line.find_last_not_of(" \t\r\n") + 1);

    if (line.empty() || line[0] == '#' || line[0] == '!')
      continue;

    // To uppercase for easier matching
    std::string uline = line;
    for (auto &c : uline)
      c = toupper(c);

    if (uline.find("&CELL") == 0) {
      parsing_cell = true;
      parsing_coords = false;
    } else if (uline.find("&COORD") == 0) {
      parsing_coords = true;
      parsing_cell = false;
      if (current_cell.atomCount() > 0) {
        frames.push_back(current_cell);
        current_cell = correlation::core::Cell();
        if (has_box) {
          current_cell.setLatticeParameters(frames.back().lattice_parameters());
        }
      }
    } else if (uline.find("&END CELL") == 0 || uline.find("&END COORD") == 0) {
      parsing_cell = false;
      parsing_coords = false;
    } else if (parsing_cell) {
      std::istringstream iss(uline);
      std::string token;
      iss >> token;
      if (token == "ABC") {
        double a, b, c;
        if (iss >> a >> b >> c) {
          current_cell.setLatticeParameters({a, b, c, 90.0, 90.0, 90.0});
          has_box = true;
        }
      } else if (token == "ALPHA_BETA_GAMMA") {
        // Can parse angles if needed
      }
    } else if (parsing_coords) {
      std::istringstream iss(line); // Original line for case sensitivity in symbols
      std::string symbol;
      double x, y, z;
      if (iss >> symbol >> x >> y >> z) {
        current_cell.addAtom(symbol, correlation::math::Vector3<double>(x, y, z));
      }
    }
  }

  if (current_cell.atomCount() > 0) {
    frames.push_back(std::move(current_cell));
  } else if (!frames.empty() && frames.back().atomCount() == 0) {
    frames.pop_back();
  }

  if (frames.empty()) {
    throw std::runtime_error("No structure found in CP2K file: " + filename);
  }

  if (progress_callback)
    progress_callback(1.0f, "CP2K trajectory loaded.");

  return correlation::core::Trajectory(frames, 1.0);
}

} // namespace correlation::readers
