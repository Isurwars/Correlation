/**
 * @file OrcaReader.cpp
 * @brief Implementation of the ORCA quantum chemistry output file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/OrcaReader.hpp"
#include "math/Precision.hpp"
#include "readers/ReaderFactory.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace correlation::readers {

namespace {

const bool registered = ReaderFactory::registerTypeSafe<OrcaReader>("OrcaReader");

constexpr real_t bohr_to_angstrom = static_cast<real_t>(0.52917721092);

struct OrcaTrajectoryParser {
  std::ifstream *file{nullptr};
  std::vector<correlation::core::Cell> frames;
  correlation::core::Cell current_cell;
  bool in_coord_block{false};
  bool is_bohr{false};

  explicit OrcaTrajectoryParser(std::ifstream &stream) : file(&stream) {}

  void checkHeader(const std::string &line) {
    if (line.contains("CARTESIAN COORDINATES (ANGSTROEM)") || line.contains("COORDINATES (ANGSTROEMS)") ||
        line.contains("CARTESIAN COORDINATES (ANGSTROM)")) {
      startNewFrame(false);
    } else if (line.contains("CARTESIAN COORDINATES (A.U.)")) {
      startNewFrame(true);
    }
  }

  void startNewFrame(bool bohr) {
    if (current_cell.atomCount() > 0) {
      frames.push_back(std::move(current_cell));
      current_cell = correlation::core::Cell();
    }
    in_coord_block = true;
    is_bohr = bohr;
  }

  void parseAtomLine(const std::string &line) {
    if (line.starts_with("---") || line.starts_with("***")) {
      return;
    }
    std::istringstream iss(line);
    std::string token1;
    if (!(iss >> token1)) {
      in_coord_block = false;
      return;
    }

    std::string symbol;
    real_t pos_x = 0.0;
    real_t pos_y = 0.0;
    real_t pos_z = 0.0;

    // Line could be "C 0.0 0.0 0.0" or "0 C 0.0 0.0 0.0" (index first)
    if (std::isalpha(static_cast<unsigned char>(token1[0])) != 0) {
      symbol = token1;
      if (!(iss >> pos_x >> pos_y >> pos_z)) {
        in_coord_block = false;
        return;
      }
    } else {
      if (!(iss >> symbol >> pos_x >> pos_y >> pos_z)) {
        in_coord_block = false;
        return;
      }
    }

    if (is_bohr) {
      pos_x *= bohr_to_angstrom;
      pos_y *= bohr_to_angstrom;
      pos_z *= bohr_to_angstrom;
    }

    current_cell.addAtom(symbol, correlation::math::Vector3<real_t>(pos_x, pos_y, pos_z));
  }

  void parse() {
    std::string line;
    while (std::getline(*file, line)) {
      line.erase(0, line.find_first_not_of(" \t\r\n"));
      line.erase(line.find_last_not_of(" \t\r\n") + 1);

      if (line.empty()) {
        if (in_coord_block && current_cell.atomCount() > 0) {
          in_coord_block = false;
        }
        continue;
      }

      checkHeader(line);

      if (in_coord_block && !line.contains("CARTESIAN COORDINATES") && !line.contains("COORDINATES (ANGSTROM")) {
        parseAtomLine(line);
      }
    }

    if (current_cell.atomCount() > 0) {
      frames.push_back(std::move(current_cell));
    }
  }
};

} // namespace

correlation::core::Cell OrcaReader::readStructure(const std::string &filename,
                                                  std::function<void(float, const std::string &)> progress_callback) {
  auto trajectory = readTrajectory(filename, std::move(progress_callback));
  if (trajectory.getFrameCount() == 0) {
    throw std::runtime_error("No atomic coordinates found in ORCA file: " + filename);
  }
  return trajectory.getFrames().back();
}

correlation::core::Trajectory
OrcaReader::readTrajectory(const std::string &filename,
                           std::function<void(float, const std::string &)> progress_callback) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open ORCA file: " + filename);
  }

  if (progress_callback) {
    progress_callback(0.1F, "Parsing ORCA output...");
  }

  OrcaTrajectoryParser parser(file);
  parser.parse();

  if (parser.frames.empty()) {
    throw std::runtime_error("No valid coordinate frames parsed from ORCA file: " + filename);
  }

  if (progress_callback) {
    progress_callback(1.0F, "Finished reading ORCA file");
  }

  return {parser.frames, 1.0};
}

} // namespace correlation::readers
