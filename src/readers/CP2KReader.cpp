/**
 * @file CP2KReader.cpp
 * @brief Implementation of the CP2K file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/CP2KReader.hpp"

#include "readers/ReaderFactory.hpp"
#include <math.h>

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace correlation::readers {

// Automatic registration
// NOLINTNEXTLINE(cert-err58-cpp, bugprone-throwing-static-initialization)
const bool registered = ReaderFactory::instance().registerReader(std::make_unique<CP2KReader>());

correlation::core::Cell CP2KReader::readStructure(const std::string &filename,
                                                  std::function<void(float, const std::string &)> progress_callback) {

  auto traj = readTrajectory(filename, progress_callback);
  if (traj.getFrameCount() == 0) {
    throw std::runtime_error("No structure found in CP2K file: " + filename);
  }
  return traj.getFrames()[traj.getFrameCount() - 1];
}

namespace {

struct Cp2kParserState {
  bool parsing_coords = false;
  bool parsing_cell = false;
  bool has_box = false;
};

void processCp2kLine(const std::string &line, std::vector<correlation::core::Cell> &frames,
                     correlation::core::Cell &current_cell, Cp2kParserState &state) {
  // To uppercase for easier matching
  std::string uline = line;
  for (auto &chr : uline) {
    chr = static_cast<char>(std::toupper(static_cast<unsigned char>(chr)));
  }

  if (uline.starts_with("&CELL")) {
    state.parsing_cell = true;
    state.parsing_coords = false;
  } else if (uline.starts_with("&COORD")) {
    state.parsing_coords = true;
    state.parsing_cell = false;
    if (current_cell.atomCount() > 0) {
      frames.push_back(current_cell);
      current_cell = correlation::core::Cell();
      if (state.has_box) {
        current_cell.setLatticeParameters(frames.back().lattice_parameters());
      }
    }
  } else if (uline.starts_with("&END CELL") || uline.starts_with("&END COORD")) {
    state.parsing_cell = false;
    state.parsing_coords = false;
  } else if (state.parsing_cell) {
    std::istringstream iss(uline);
    std::string token;
    iss >> token;
    if (token == "ABC") {
      real_t param_a = 0.0;
      real_t param_b = 0.0;
      real_t param_c = 0.0;
      if (iss >> param_a >> param_b >> param_c) {
        current_cell.setLatticeParameters({param_a, param_b, param_c, 90.0, 90.0, 90.0});
        state.has_box = true;
      }
    }
  } else if (state.parsing_coords) {
    std::istringstream iss(line); // Original line for case sensitivity in symbols
    std::string symbol;
    real_t pos_x = 0.0;
    real_t pos_y = 0.0;
    real_t pos_z = 0.0;
    if (iss >> symbol >> pos_x >> pos_y >> pos_z) {
      current_cell.addAtom(symbol, correlation::math::Vector3<real_t>(pos_x, pos_y, pos_z));
    }
  }
}

} // namespace

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
  Cp2kParserState parser_state;

  while (std::getline(file, line)) {
    line.erase(0, line.find_first_not_of(" \t\r\n"));
    line.erase(line.find_last_not_of(" \t\r\n") + 1);

    if (line.empty() || line[0] == '#' || line[0] == '!') {
      continue;
    }

    processCp2kLine(line, frames, current_cell, parser_state);
  }

  if (current_cell.atomCount() > 0) {
    frames.push_back(std::move(current_cell));
  } else if (!frames.empty() && frames.back().atomCount() == 0) {
    frames.pop_back();
  }

  if (frames.empty()) {
    throw std::runtime_error("No structure found in CP2K file: " + filename);
  }

  if (progress_callback) {
    progress_callback(1.0F, "CP2K trajectory loaded.");
  }

  return {frames, 1.0};
}

} // namespace correlation::readers
