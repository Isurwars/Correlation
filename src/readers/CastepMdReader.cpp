/**
 * @file CastepMdReader.cpp
 * @brief Implementation of the CASTEP MD trajectory reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/CastepMdReader.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "math/Constants.hpp"
#include "math/LinearAlgebra.hpp"
#include "readers/ReaderFactory.hpp"

#include <cerrno>
#include <cstring>
#include <fstream>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace correlation::readers {

// Automatic registration
const bool REGISTERED = ReaderFactory::registerTypeSafe<CastepMdReader>("CastepMdReader");

correlation::core::Cell
CastepMdReader::readStructure(const std::string &filename,
                              std::function<void(float, const std::string &)> progress_callback) {
  return CastepMdReader::read(filename, progress_callback).at(0);
}

correlation::core::Trajectory
CastepMdReader::readTrajectory(const std::string &filename,
                               std::function<void(float, const std::string &)> progress_callback) {
  return {CastepMdReader::read(filename, progress_callback), 1.0};
}

std::vector<correlation::core::Cell>
CastepMdReader::read(const std::string &file_name,
                     const std::function<void(float, const std::string &)> &progress_callback) {
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" + std::strerror(errno) +
                             ").");
  }

  std::vector<correlation::core::Cell> frames;
  std::string line;

  // Parse frames
  myfile.seekg(0, std::ios::end);
  std::streampos const file_size = myfile.tellg();
  myfile.clear();
  myfile.seekg(0, std::ios::beg);
  std::streampos last_progress_pos = 0;
  size_t const update_interval = file_size / 100;

  // Skip the header
  bool in_header = true;
  while (in_header && std::getline(myfile, line)) {
    if (line.empty() || line.find_first_not_of(' ') == std::string::npos) {
      in_header = false;
    }
  }

  // Temporary storage for building the current cell
  correlation::core::Cell temp_cell;
  real_t current_energy = 0.0;
  bool cell_has_atoms = false;

  while (std::getline(myfile, line)) {
    if (progress_callback) {
      updateProgress(myfile.tellg(), file_size, last_progress_pos, update_interval,
                     progress_callback);
    }
    // If line represents energies/time
    if (line.contains("<-- E")) {
      parseEnergyLine(line, current_energy, temp_cell, cell_has_atoms, frames);
      continue;
    }

    if (line.contains("<-- h") && !line.contains("<-- hv")) {
      parseLatticeLine(myfile, line, temp_cell);
      continue;
    }

    if (line.contains("<-- R")) {
      parseAtomLine(line, current_energy, temp_cell, cell_has_atoms);
      continue;
    }
  }

  // Push the final frame if it has atoms
  if (cell_has_atoms && !temp_cell.isEmpty()) {
    frames.push_back(std::move(temp_cell));
  }

  return frames;
}

void CastepMdReader::updateProgress(
    std::streampos current_pos, std::streampos file_size, std::streampos &last_progress_pos,
    size_t update_interval,
    const std::function<void(float, const std::string &)> &progress_callback) {
  if (std::cmp_greater(current_pos - last_progress_pos, update_interval)) {
    const float progress = static_cast<float>(current_pos) / static_cast<float>(file_size);
    progress_callback(progress, "Loading CASTEP MD file...");
    last_progress_pos = current_pos;
  }
}

void CastepMdReader::parseEnergyLine(const std::string &line, real_t &current_energy,
                                     correlation::core::Cell &temp_cell, bool &cell_has_atoms,
                                     std::vector<correlation::core::Cell> &frames) {
  // For CASTEP MD, the frame usually starts with time (single float on a
  // line without tag) We skip it and read E. If we have atoms in temp_cell,
  // it means this is a new frame starting
  if (cell_has_atoms) {
    // Save lattice parameters and energy
    std::array<real_t, 6> const last_lattice = temp_cell.latticeParameters();
    real_t const last_energy = temp_cell.getEnergy();

    frames.push_back(std::move(temp_cell));

    // Re-initialize for next frame but keep lattice & energy
    temp_cell = correlation::core::Cell(last_lattice);
    temp_cell.setEnergy(last_energy);
    cell_has_atoms = false;
  }

  // We can extract energy if needed...
  std::stringstream line_stream(line);
  if (line_stream >> current_energy) {
    temp_cell.setEnergy(current_energy); // Update energy
  }
}

void CastepMdReader::parseLatticeLine(std::ifstream &myfile, const std::string &line,
                                      correlation::core::Cell &temp_cell) {
  // Lattice vectors h are given row by row in Bohr
  // The first <-- h is row 1
  std::array<real_t, 3> lattice_vector_1{};
  std::array<real_t, 3> lattice_vector_2{};
  std::array<real_t, 3> lattice_vector_3{};
  std::stringstream line_stream(line);
  line_stream >> lattice_vector_1[0] >> lattice_vector_1[1] >> lattice_vector_1[2];

  std::string next_line;
  if (std::getline(myfile, next_line)) { // row 2
    std::stringstream ss2(next_line);
    ss2 >> lattice_vector_2[0] >> lattice_vector_2[1] >> lattice_vector_2[2];
  }

  if (std::getline(myfile, next_line)) { // row 3
    std::stringstream ss3(next_line);
    ss3 >> lattice_vector_3[0] >> lattice_vector_3[1] >> lattice_vector_3[2];
  }

  // Convert Bohr to Angstroms
  lattice_vector_1[0] *= correlation::math::bohr_to_angstrom;
  lattice_vector_1[1] *= correlation::math::bohr_to_angstrom;
  lattice_vector_1[2] *= correlation::math::bohr_to_angstrom;

  lattice_vector_2[0] *= correlation::math::bohr_to_angstrom;
  lattice_vector_2[1] *= correlation::math::bohr_to_angstrom;
  lattice_vector_2[2] *= correlation::math::bohr_to_angstrom;

  lattice_vector_3[0] *= correlation::math::bohr_to_angstrom;
  lattice_vector_3[1] *= correlation::math::bohr_to_angstrom;
  lattice_vector_3[2] *= correlation::math::bohr_to_angstrom;

  temp_cell =
      correlation::core::Cell({lattice_vector_1[0], lattice_vector_1[1], lattice_vector_1[2]},
                              {lattice_vector_2[0], lattice_vector_2[1], lattice_vector_2[2]},
                              {lattice_vector_3[0], lattice_vector_3[1], lattice_vector_3[2]});
}

void CastepMdReader::parseAtomLine(const std::string &line, real_t current_energy,
                                   correlation::core::Cell &temp_cell, bool &cell_has_atoms) {
  // Positions R are given in Bohr
  // Format: Symbol ID x y z <-- R
  std::stringstream line_stream(line);
  std::string symbol;
  int atom_id = 0;
  real_t coord_x = 0.0;
  real_t coord_y = 0.0;
  real_t coord_z = 0.0;
  if (line_stream >> symbol >> atom_id >> coord_x >> coord_y >> coord_z) {
    temp_cell.addAtom(
        symbol, correlation::math::Vector3<real_t>(coord_x * correlation::math::bohr_to_angstrom,
                                                   coord_y * correlation::math::bohr_to_angstrom,
                                                   coord_z * correlation::math::bohr_to_angstrom));
    temp_cell.setEnergy(current_energy); // Assign energy once per atom or frame
    cell_has_atoms = true;
  }
}

} // namespace correlation::readers
