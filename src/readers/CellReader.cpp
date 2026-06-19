/**
 * @file CellReader.cpp
 * @brief Implementation of the CASTEP CELL file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */
#include "readers/CellReader.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/ReaderFactory.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace correlation::readers {

// Automatic registration
// NOLINTNEXTLINE(cert-err58-cpp, bugprone-throwing-static-initialization)
static const bool registered = ReaderFactory::instance().registerReader(std::make_unique<CellReader>());

correlation::core::Cell
CellReader::readStructure(const std::string &filename,
                          std::function<void(float, const std::string &)> /*progress_callback*/) {
  return read(filename);
}

correlation::core::Trajectory
CellReader::readTrajectory(const std::string & /*filename*/,
                           std::function<void(float, const std::string &)> /*progress_callback*/) {
  throw std::runtime_error("Cell files are structures, use readStructure.");
}

namespace {
void toLower(std::string &str) {
  std::transform(str.begin(), str.end(), str.begin(), [](unsigned char chr) { return std::tolower(chr); });
}

void parseLatticeCart(std::stringstream &line_stream, int &lattice_row_count,
                      std::array<std::array<double, 3>, 3> &lat_vec, correlation::core::Cell &tempCell) {
  if (lattice_row_count >= 3) {
    return;
  }
  line_stream.clear();
  line_stream.seekg(0); // Reread the full line
  if (line_stream >> lat_vec.at(lattice_row_count).at(0) >> lat_vec.at(lattice_row_count).at(1) >>
      lat_vec.at(lattice_row_count).at(2)) {
    lattice_row_count++;
    if (lattice_row_count == 3) {
      tempCell = correlation::core::Cell({lat_vec.at(0).at(0), lat_vec.at(0).at(1), lat_vec.at(0).at(2)},
                                         {lat_vec.at(1).at(0), lat_vec.at(1).at(1), lat_vec.at(1).at(2)},
                                         {lat_vec.at(2).at(0), lat_vec.at(2).at(1), lat_vec.at(2).at(2)});
    }
  }
}

void parseLatticeAbc(std::stringstream &line_stream, int &lattice_row_count, std::array<double, 6> &lat,
                     correlation::core::Cell &tempCell) {
  if (lattice_row_count >= 2) {
    return;
  }
  line_stream.clear();
  line_stream.seekg(0);
  if (line_stream >> lat.at(0 + 3 * lattice_row_count) >> lat.at(1 + 3 * lattice_row_count) >>
      lat.at(2 + 3 * lattice_row_count)) {
    lattice_row_count++;
    if (lattice_row_count == 2) {
      tempCell.setLatticeParameters(lat);
    }
  }
}

void parsePositionsAbs(std::stringstream &line_stream, correlation::core::Cell &tempCell) {
  std::string element;
  double coord_x = 0.0;
  double coord_y = 0.0;
  double coord_z = 0.0;
  line_stream.clear();
  line_stream.seekg(0);
  if (line_stream >> element >> coord_x >> coord_y >> coord_z) {
    tempCell.addAtom(element, {coord_x, coord_y, coord_z});
  }
}

void parsePositionsFrac(std::stringstream &line_stream, correlation::core::Cell &tempCell, bool &frac_flag) {
  std::string element;
  double coord_x = 0.0;
  double coord_y = 0.0;
  double coord_z = 0.0;
  line_stream.clear();
  line_stream.seekg(0);
  if (line_stream >> element >> coord_x >> coord_y >> coord_z) {
    frac_flag = true;
    const auto &lattice_vectors = tempCell.latticeVectors();
    correlation::math::Vector3<double> const pos = {
        coord_x * lattice_vectors[0].x() + coord_y * lattice_vectors[1].x() + coord_z * lattice_vectors[2].x(),
        coord_x * lattice_vectors[0].y() + coord_y * lattice_vectors[1].y() + coord_z * lattice_vectors[2].y(),
        coord_x * lattice_vectors[0].z() + coord_y * lattice_vectors[1].z() + coord_z * lattice_vectors[2].z()};
    tempCell.addAtom(element, pos);
  }
}
} // namespace

correlation::core::Cell CellReader::read(const std::string &file_name) {
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" + std::strerror(errno) + ").");
  }

  correlation::core::Cell tempCell;
  bool in_block = false;
  bool frac_flag = false;
  std::string current_block_type;
  std::string line;
  int lattice_row_count = 0;
  std::array<double, 6> lat = {0, 0, 0, 0, 0, 0};
  std::array<std::array<double, 3>, 3> lattice_matrix = {};

  while (std::getline(myfile, line)) {
    /* Read line by line */
    std::stringstream line_stream(line);
    std::string token;
    line_stream >> token;
    toLower(token); // For case-insensitive matching

    if (token == "%block") {
      in_block = true;
      line_stream >> current_block_type;
      toLower(current_block_type);
      lattice_row_count = 0; // Reset for a new lattice block
      continue;
    }

    if (token == "%endblock") {
      in_block = false;
      current_block_type.clear();
      continue;
    }

    if (in_block) {
      if (current_block_type == "lattice_cart") {
        parseLatticeCart(line_stream, lattice_row_count, lattice_matrix, tempCell);
      } else if (current_block_type == "lattice_abc") {
        parseLatticeAbc(line_stream, lattice_row_count, lat, tempCell);
      } else if (current_block_type == "positions_abs") {
        parsePositionsAbs(line_stream, tempCell);
      } else if (current_block_type == "positions_frac") {
        parsePositionsFrac(line_stream, tempCell, frac_flag);
      }
    }
  }

  if (frac_flag) {
    tempCell.wrapPositions();
  }
  return tempCell;
}

} // namespace correlation::readers
