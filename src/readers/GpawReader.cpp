/**
 * @file GpawReader.cpp
 * @brief Implementation of the GPAW DFT output/log file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/GpawReader.hpp"
#include "math/LinearAlgebra.hpp"
#include "math/Precision.hpp"
#include "readers/ReaderFactory.hpp"

#include <array>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace correlation::readers {

namespace {

const bool registered = ReaderFactory::registerTypeSafe<GpawReader>("GpawReader");

struct GpawTrajectoryParser {
  std::ifstream *file{nullptr};
  std::vector<correlation::core::Cell> frames;
  correlation::core::Cell current_cell;
  bool parsing_positions{false};
  bool has_lattice{false};
  correlation::math::Matrix3<real_t> current_lattice;

  explicit GpawTrajectoryParser(std::ifstream &stream) : file(&stream) {}

  void parseLatticeVectors() {
    std::array<std::array<real_t, 3>, 3> lat_vec{};
    std::string line;
    for (auto &row : lat_vec) {
      if (std::getline(*file, line)) {
        // Strip braces/brackets/commas
        for (char &chr : line) {
          if (chr == '[' || chr == ']' || chr == ',') {
            chr = ' ';
          }
        }
        std::istringstream iss(line);
        iss >> row.at(0) >> row.at(1) >> row.at(2);
      }
    }
    current_lattice = correlation::math::Matrix3<real_t>({lat_vec.at(0)[0], lat_vec.at(0)[1], lat_vec.at(0)[2]},
                                                         {lat_vec.at(1)[0], lat_vec.at(1)[1], lat_vec.at(1)[2]},
                                                         {lat_vec.at(2)[0], lat_vec.at(2)[1], lat_vec.at(2)[2]});
    has_lattice = true;
    current_cell.updateLattice(current_lattice);
  }

  void startNewPositionBlock() {
    if (current_cell.atomCount() > 0) {
      frames.push_back(std::move(current_cell));
      current_cell = correlation::core::Cell();
      if (has_lattice) {
        current_cell.updateLattice(current_lattice);
      }
    }
    parsing_positions = true;
  }

  void parseAtomLine(const std::string &line) {
    if (line.starts_with("---") || line.starts_with("===")) {
      return;
    }

    std::istringstream iss(line);
    std::string token1;
    if (!(iss >> token1)) {
      parsing_positions = false;
      return;
    }

    std::string symbol;
    real_t pos_x = 0.0;
    real_t pos_y = 0.0;
    real_t pos_z = 0.0;

    if (std::isalpha(static_cast<unsigned char>(token1[0])) != 0) {
      symbol = token1;
      if (!(iss >> pos_x >> pos_y >> pos_z)) {
        parsing_positions = false;
        return;
      }
    } else {
      if (!(iss >> symbol >> pos_x >> pos_y >> pos_z)) {
        parsing_positions = false;
        return;
      }
    }

    current_cell.addAtom(symbol, correlation::math::Vector3<real_t>(pos_x, pos_y, pos_z));
  }

  void parse() {
    std::string line;
    while (std::getline(*file, line)) {
      line.erase(0, line.find_first_not_of(" \t\r\n"));
      line.erase(line.find_last_not_of(" \t\r\n") + 1);

      if (line.empty()) {
        continue;
      }

      if (line.starts_with("Unit cell:") || line.starts_with("Lattice vectors:") || line.starts_with("Cell:")) {
        parseLatticeVectors();
      } else if (line.starts_with("Positions:") || line.starts_with("Cartesian positions:") ||
                 line.starts_with("ATOMIC_POSITIONS")) {
        startNewPositionBlock();
      } else if (parsing_positions) {
        parseAtomLine(line);
      }
    }

    if (current_cell.atomCount() > 0) {
      frames.push_back(std::move(current_cell));
    }
  }
};

} // namespace

correlation::core::Cell GpawReader::readStructure(const std::string &filename,
                                                  std::function<void(float, const std::string &)> progress_callback) {
  auto trajectory = readTrajectory(filename, std::move(progress_callback));
  if (trajectory.getFrameCount() == 0) {
    throw std::runtime_error("No atomic coordinates found in GPAW file: " + filename);
  }
  return trajectory.getFrames().back();
}

correlation::core::Trajectory
GpawReader::readTrajectory(const std::string &filename,
                           std::function<void(float, const std::string &)> progress_callback) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open GPAW file: " + filename);
  }

  if (progress_callback) {
    progress_callback(0.1F, "Parsing GPAW output...");
  }

  GpawTrajectoryParser parser(file);
  parser.parse();

  if (parser.frames.empty()) {
    throw std::runtime_error("No valid coordinate frames parsed from GPAW file: " + filename);
  }

  if (progress_callback) {
    progress_callback(1.0F, "Finished reading GPAW file");
  }

  return {parser.frames, 1.0};
}

} // namespace correlation::readers
