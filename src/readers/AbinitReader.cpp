/**
 * @file AbinitReader.cpp
 * @brief Implementation of the ABINIT input/output file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/AbinitReader.hpp"
#include "math/LinearAlgebra.hpp"
#include "math/Precision.hpp"
#include "physics/PhysicalData.hpp"
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

const bool registered = ReaderFactory::registerTypeSafe<AbinitReader>("AbinitReader");

constexpr real_t bohr_to_angstrom = static_cast<real_t>(0.52917721092);

struct AbinitTrajectoryParser {
  std::ifstream *file{nullptr};
  std::vector<correlation::core::Cell> frames;
  correlation::core::Cell current_cell;
  bool parsing_xcart{false};
  bool parsing_xangst{false};
  bool has_lattice{false};
  std::vector<std::string> symbols;
  std::array<real_t, 3> acell{1.0, 1.0, 1.0};
  std::array<std::array<real_t, 3>, 3> rprim{{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};

  explicit AbinitTrajectoryParser(std::ifstream &stream) : file(&stream) {}

  void parseAcell(const std::string &line) {
    std::string const clean = line.substr(line.find("acell") + 5);
    std::string filtered = clean;
    for (char &chr : filtered) {
      if (chr == '*' || chr == ',') {
        chr = ' ';
      }
    }
    std::istringstream iss(filtered);
    iss >> acell[0] >> acell[1] >> acell[2];
    updateCellLattice();
  }

  void parseRprim(const std::string &line) {
    std::string const clean = line.substr(line.find("rprim") + 5);
    std::istringstream iss(clean);
    for (auto &row : rprim) {
      iss >> row[0] >> row[1] >> row[2];
    }
    updateCellLattice();
  }

  void updateCellLattice() {
    correlation::math::Matrix3<real_t> const lattice(
        {rprim[0][0] * acell[0] * bohr_to_angstrom, rprim[0][1] * acell[0] * bohr_to_angstrom,
         rprim[0][2] * acell[0] * bohr_to_angstrom},
        {rprim[1][0] * acell[1] * bohr_to_angstrom, rprim[1][1] * acell[1] * bohr_to_angstrom,
         rprim[1][2] * acell[1] * bohr_to_angstrom},
        {rprim[2][0] * acell[2] * bohr_to_angstrom, rprim[2][1] * acell[2] * bohr_to_angstrom,
         rprim[2][2] * acell[2] * bohr_to_angstrom});
    current_cell.updateLattice(lattice);
    has_lattice = true;
  }

  void startNewCoordBlock(bool is_angstrom) {
    if (current_cell.atomCount() > 0) {
      frames.push_back(std::move(current_cell));
      current_cell = correlation::core::Cell();
      if (has_lattice) {
        current_cell.updateLattice(frames.back().latticeVectors());
      }
    }
    parsing_xangst = is_angstrom;
    parsing_xcart = !is_angstrom;
  }

  void parseAtomLine(const std::string &line) {
    std::istringstream iss(line);
    std::string token1;
    if (!(iss >> token1)) {
      parsing_xcart = false;
      parsing_xangst = false;
      return;
    }

    std::string symbol = "X";
    real_t pos_x = 0.0;
    real_t pos_y = 0.0;
    real_t pos_z = 0.0;

    if (std::isalpha(static_cast<unsigned char>(token1[0])) != 0) {
      symbol = token1;
      if (!(iss >> pos_x >> pos_y >> pos_z)) {
        parsing_xcart = false;
        parsing_xangst = false;
        return;
      }
    } else {
      try {
        pos_x = static_cast<real_t>(std::stod(token1));
        if (!(iss >> pos_y >> pos_z)) {
          parsing_xcart = false;
          parsing_xangst = false;
          return;
        }
      } catch (...) {
        parsing_xcart = false;
        parsing_xangst = false;
        return;
      }
    }

    if (parsing_xcart) {
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
        continue;
      }

      if (line.contains("acell")) {
        parseAcell(line);
      } else if (line.contains("rprim")) {
        parseRprim(line);
      } else if (line.contains("xangst") || line.contains("Cartesian coordinates (angstrom)")) {
        startNewCoordBlock(true);
      } else if (line.contains("xcart") || line.contains("Cartesian coordinates (bohr)")) {
        startNewCoordBlock(false);
      } else if (parsing_xcart || parsing_xangst) {
        parseAtomLine(line);
      }
    }

    if (current_cell.atomCount() > 0) {
      frames.push_back(std::move(current_cell));
    }
  }
};

} // namespace

correlation::core::Cell AbinitReader::readStructure(const std::string &filename,
                                                    std::function<void(float, const std::string &)> progress_callback) {
  auto trajectory = readTrajectory(filename, std::move(progress_callback));
  if (trajectory.getFrameCount() == 0) {
    throw std::runtime_error("No atomic coordinates found in ABINIT file: " + filename);
  }
  return trajectory.getFrames().back();
}

correlation::core::Trajectory
AbinitReader::readTrajectory(const std::string &filename,
                             std::function<void(float, const std::string &)> progress_callback) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open ABINIT file: " + filename);
  }

  if (progress_callback) {
    progress_callback(0.1F, "Parsing ABINIT output...");
  }

  AbinitTrajectoryParser parser(file);
  parser.parse();

  if (parser.frames.empty()) {
    throw std::runtime_error("No valid coordinate frames parsed from ABINIT file: " + filename);
  }

  if (progress_callback) {
    progress_callback(1.0F, "Finished reading ABINIT file");
  }

  return {parser.frames, 1.0};
}

} // namespace correlation::readers
