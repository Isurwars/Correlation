/**
 * @file QEReader.cpp
 * @brief Implementation of the Quantum ESPRESSO file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/QEReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace correlation::readers {

namespace {
bool registered = ReaderFactory::instance().registerReader(std::make_unique<QEReader>());
}

correlation::core::Cell QEReader::readStructure(const std::string &filename,
                                                std::function<void(float, const std::string &)> progress_callback) {

  auto traj = readTrajectory(filename, progress_callback);
  if (traj.getFrameCount() == 0) {
    throw std::runtime_error("No structure found in QE file: " + filename);
  }
  // Return the last frame typically, as it might be an optimization
  return traj.getFrames()[traj.getFrameCount() - 1];
}

correlation::core::Trajectory
QEReader::readTrajectory(const std::string &filename,
                         std::function<void(float, const std::string &)> progress_callback) {

  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  std::vector<correlation::core::Cell> frames;
  std::string line;
  correlation::core::Cell current_cell;
  bool parsing_atoms = false;
  bool parsing_cell = false;
  bool has_box = false;
  double a = 0;
  double b = 0;
  double c = 0;

  while (std::getline(file, line)) {
    // Basic trimming
    line.erase(0, line.find_first_not_of(" \t\r\n"));
    line.erase(line.find_last_not_of(" \t\r\n") + 1);

    if (line.empty()) {
      continue;
}

    if (line.starts_with("CELL_PARAMETERS")) {
      parsing_cell = true;
      parsing_atoms = false;
      // Note: Full lattice vector parsing is more complex, but we'll try basic diagonal
      // For a robust implementation, you parse the 3x3 matrix.
      double v[3][3];
      for (auto & i : v) {
        if (std::getline(file, line)) {
          std::istringstream iss(line);
          iss >> i[0] >> i[1] >> i[2];
        }
      }
      current_cell.updateLattice(correlation::math::Matrix3<double>(
          {v[0][0], v[0][1], v[0][2]}, {v[1][0], v[1][1], v[1][2]}, {v[2][0], v[2][1], v[2][2]}));
      has_box = true;
      parsing_cell = false;
    } else if (line.starts_with("ATOMIC_POSITIONS")) {
      parsing_atoms = true;
      parsing_cell = false;
      if (current_cell.atomCount() > 0) {
        frames.push_back(current_cell);
        current_cell = correlation::core::Cell();
        if (has_box) {
          current_cell.updateLattice(frames.back().latticeVectors());
        }
      }
    } else if (parsing_atoms) {
      std::istringstream iss(line);
      std::string symbol;
      double x;
      double y;
      double z;
      if (iss >> symbol >> x >> y >> z) {
        // Sometimes there are constraints following, ignore them
        current_cell.addAtom(symbol, correlation::math::Vector3<double>(x, y, z));
      } else {
        // Stop parsing atoms if we hit an invalid line
        parsing_atoms = false;
      }
    }
  }

  // Push the final frame if it has atoms
  if (current_cell.atomCount() > 0) {
    frames.push_back(std::move(current_cell));
  } else if (!frames.empty() && frames.back().atomCount() == 0) {
    frames.pop_back();
  }

  if (frames.empty()) {
    throw std::runtime_error("No structure found in QE file: " + filename);
  }

  if (progress_callback) {
    progress_callback(1.0F, "QE trajectory loaded.");
}

  return {frames, 1.0};
}

} // namespace correlation::readers
