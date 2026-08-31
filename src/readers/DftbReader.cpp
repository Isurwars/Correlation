/**
 * @file DftbReader.cpp
 * @brief Implementation of the DFTB+ (.gen, .dftb, .hsd) file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/DftbReader.hpp"
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

const bool registered = ReaderFactory::registerTypeSafe<DftbReader>("DftbReader");

struct GenHeaderData {
  size_t atom_count{0};
  char mode{'C'};
  std::vector<std::string> element_types;
};

struct TempAtom {
  size_t type_idx{0};
  real_t x{0.0};
  real_t y{0.0};
  real_t z{0.0};
};

GenHeaderData parse_gen_header(std::istream &stream) {
  std::string header_line;
  if (!std::getline(stream, header_line)) {
    throw std::runtime_error("Empty DFTB+ GenFormat file.");
  }

  std::istringstream h_iss(header_line);
  size_t atom_count = 0;
  char mode = 'C';
  if (!(h_iss >> atom_count >> mode)) {
    throw std::runtime_error("Invalid DFTB+ GenFormat header: " + header_line);
  }
  mode = static_cast<char>(std::toupper(static_cast<unsigned char>(mode)));

  std::string element_line;
  if (!std::getline(stream, element_line)) {
    throw std::runtime_error("Missing element symbol list in DFTB+ GenFormat file.");
  }

  std::istringstream elem_iss(element_line);
  std::vector<std::string> element_types;
  std::string elem_sym;
  while (elem_iss >> elem_sym) {
    element_types.push_back(elem_sym);
  }

  if (element_types.empty()) {
    throw std::runtime_error("No element symbols found in DFTB+ header.");
  }

  return {.atom_count = atom_count, .mode = mode, .element_types = std::move(element_types)};
}

std::vector<TempAtom> parse_raw_atoms(std::istream &stream, size_t atom_count) {
  std::vector<TempAtom> raw_atoms;
  raw_atoms.reserve(atom_count);

  for (size_t i = 0; i < atom_count; ++i) {
    std::string line;
    if (!std::getline(stream, line)) {
      throw std::runtime_error("Unexpected end of file while reading DFTB+ atom coordinates.");
    }
    std::istringstream iss(line);
    size_t idx = 0;
    size_t type_idx = 0;
    real_t pos_x = 0.0;
    real_t pos_y = 0.0;
    real_t pos_z = 0.0;
    if (!(iss >> idx >> type_idx >> pos_x >> pos_y >> pos_z)) {
      throw std::runtime_error("Malformed DFTB+ atom line: " + line);
    }
    raw_atoms.push_back({.type_idx = type_idx, .x = pos_x, .y = pos_y, .z = pos_z});
  }
  return raw_atoms;
}

correlation::math::Matrix3<real_t> parse_lattice(std::istream &stream) {
  std::string origin_line;
  std::getline(stream, origin_line); // origin (0 0 0)

  std::array<std::array<real_t, 3>, 3> lat_vec{};
  for (auto &row : lat_vec) {
    std::string lat_line;
    if (std::getline(stream, lat_line)) {
      std::istringstream lat_iss(lat_line);
      lat_iss >> row[0] >> row[1] >> row[2];
    }
  }
  return {
      correlation::math::Vector3<real_t>(lat_vec[0][0], lat_vec[0][1], lat_vec[0][2]),
      correlation::math::Vector3<real_t>(lat_vec[1][0], lat_vec[1][1], lat_vec[1][2]),
      correlation::math::Vector3<real_t>(lat_vec[2][0], lat_vec[2][1], lat_vec[2][2])};
}

void populate_cell_atoms(correlation::core::Cell &cell, const std::vector<TempAtom> &raw_atoms,
                         const std::vector<std::string> &element_types, char mode,
                         const correlation::math::Matrix3<real_t> &lattice) {
  for (const auto &raw : raw_atoms) {
    size_t const elem_idx = (raw.type_idx > 0 && raw.type_idx <= element_types.size()) ? (raw.type_idx - 1) : 0;
    const std::string &symbol = element_types[elem_idx];
    if (mode == 'F') {
      correlation::math::Vector3<real_t> const cart = lattice * correlation::math::Vector3<real_t>(raw.x, raw.y, raw.z);
      cell.addAtom(symbol, cart);
    } else {
      cell.addAtom(symbol, correlation::math::Vector3<real_t>(raw.x, raw.y, raw.z));
    }
  }
}

correlation::core::Cell parse_gen_file(std::istream &stream) {
  GenHeaderData const header = parse_gen_header(stream);
  std::vector<TempAtom> const raw_atoms = parse_raw_atoms(stream, header.atom_count);

  correlation::core::Cell cell;
  if (header.mode == 'S' || header.mode == 'F') {
    correlation::math::Matrix3<real_t> const lattice = parse_lattice(stream);
    cell.updateLattice(lattice);
    populate_cell_atoms(cell, raw_atoms, header.element_types, header.mode, lattice);
  } else {
    populate_cell_atoms(cell, raw_atoms, header.element_types, header.mode, correlation::math::Matrix3<real_t>());
  }
  return cell;
}

} // namespace

correlation::core::Cell DftbReader::readStructure(const std::string &filename,
                                                  std::function<void(float, const std::string &)> progress_callback) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open DFTB+ file: " + filename);
  }

  if (progress_callback) {
    progress_callback(0.1F, "Parsing DFTB+ GenFormat structure...");
  }

  correlation::core::Cell cell = parse_gen_file(file);

  if (cell.atomCount() == 0) {
    throw std::runtime_error("No atoms parsed from DFTB+ file: " + filename);
  }

  if (progress_callback) {
    progress_callback(1.0F, "Finished reading DFTB+ file");
  }

  return cell;
}

correlation::core::Trajectory
DftbReader::readTrajectory(const std::string &filename,
                           std::function<void(float, const std::string &)> progress_callback) {
  correlation::core::Cell const cell = readStructure(filename, std::move(progress_callback));
  correlation::core::Trajectory trajectory;
  trajectory.addFrame(cell);
  return trajectory;
}

} // namespace correlation::readers
