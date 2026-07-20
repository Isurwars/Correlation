/**
 * @file QEReader.cpp
 * @brief Implementation of the Quantum ESPRESSO file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/QEReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <array>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace correlation::readers {

namespace {

// Automatic registration
// NOLINTNEXTLINE(cert-err58-cpp, bugprone-throwing-static-initialization)
const bool registered = ReaderFactory::instance().registerReader(std::make_unique<QEReader>());

struct QETrajectoryParser {
  std::ifstream *file = nullptr;
  std::vector<correlation::core::Cell> frames;
  correlation::core::Cell current_cell;
  bool parsing_atoms = false;
  bool has_box = false;

  QETrajectoryParser(std::ifstream &file) : file(&file) {}

  void parseCellParameters() {
    std::array<std::array<real_t, 3>, 3> lat_vec = {};
    std::string line;
    for (auto &row : lat_vec) {
      if (std::getline(*file, line)) {
        std::istringstream iss(line);
        iss >> row.at(0) >> row.at(1) >> row.at(2);
      }
    }
    current_cell.updateLattice(correlation::math::Matrix3<real_t>(
        {lat_vec.at(0)[0], lat_vec.at(0)[1], lat_vec.at(0)[2]}, {lat_vec.at(1)[0], lat_vec.at(1)[1], lat_vec.at(1)[2]},
        {lat_vec.at(2)[0], lat_vec.at(2)[1], lat_vec.at(2)[2]}));
    has_box = true;
    parsing_atoms = false;
  }

  void handleAtomicPositionsLine() {
    parsing_atoms = true;
    if (current_cell.atomCount() > 0) {
      frames.push_back(current_cell);
      current_cell = correlation::core::Cell();
      if (has_box) {
        current_cell.updateLattice(frames.back().latticeVectors());
      }
    }
  }

  void parseAtomLine(const std::string &line) {
    std::istringstream iss(line);
    std::string symbol;
    real_t pos_x = 0.0;
    real_t pos_y = 0.0;
    real_t pos_z = 0.0;
    if (iss >> symbol >> pos_x >> pos_y >> pos_z) {
      current_cell.addAtom(symbol, correlation::math::Vector3<real_t>(pos_x, pos_y, pos_z));
    } else {
      parsing_atoms = false;
    }
  }

  void parse() {
    std::string line;
    while (std::getline(*file, line)) {
      // Basic trimming
      line.erase(0, line.find_first_not_of(" \t\r\n"));
      line.erase(line.find_last_not_of(" \t\r\n") + 1);

      if (line.empty()) {
        continue;
      }

      if (line.starts_with("CELL_PARAMETERS")) {
        parseCellParameters();
      } else if (line.starts_with("ATOMIC_POSITIONS")) {
        handleAtomicPositionsLine();
      } else if (parsing_atoms) {
        parseAtomLine(line);
      }
    }

    if (current_cell.atomCount() > 0) {
      frames.push_back(std::move(current_cell));
    } else if (!frames.empty() && frames.back().atomCount() == 0) {
      frames.pop_back();
    }
  }
};

} // namespace

correlation::core::Cell QEReader::readStructure(const std::string &filename,
                                                std::function<void(float, const std::string &)> progress_callback) {
  auto traj = readTrajectory(filename, progress_callback);
  if (traj.getFrameCount() == 0) {
    throw std::runtime_error("No structure found in QE file: " + filename);
  }
  return traj.getFrames()[traj.getFrameCount() - 1];
}

correlation::core::Trajectory
QEReader::readTrajectory(const std::string &filename,
                         std::function<void(float, const std::string &)> progress_callback) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  QETrajectoryParser parser(file);
  parser.parse();

  if (parser.frames.empty()) {
    throw std::runtime_error("No structure found in QE file: " + filename);
  }

  if (progress_callback) {
    progress_callback(1.0F, "QE trajectory loaded.");
  }

  return {parser.frames, 1.0};
}

} // namespace correlation::readers
