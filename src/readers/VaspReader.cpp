/**
 * @file VaspReader.cpp
 * @brief Implementation of the VASP POSCAR/CONTCAR file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */
#include "readers/VaspReader.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/ReaderFactory.hpp"

#include <cerrno>
#include <cmath>
#include <cstring>
#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace correlation::readers {

// Automatic registration
static bool registered =
    ReaderFactory::instance().registerReader(std::make_unique<VaspReader>());

correlation::core::Cell VaspReader::readStructure(
    const std::string &filename,
    std::function<void(float, const std::string &)> progress_callback) {
  return read(filename);
}

correlation::core::Trajectory VaspReader::readTrajectory(
    const std::string &filename,
    std::function<void(float, const std::string &)> progress_callback) {
  throw std::runtime_error(
      "POSCAR/CONTCAR files are single structures, use readStructure.");
}

correlation::core::Cell VaspReader::read(const std::string &file_name) {
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" +
                             std::strerror(errno) + ").");
  }

  std::string line;

  // Line 1: Comment
  if (!std::getline(myfile, line)) {
    throw std::runtime_error("POSCAR: unexpected end of file (comment line).");
  }

  // Line 2: Scaling factor
  if (!std::getline(myfile, line)) {
    throw std::runtime_error(
        "POSCAR: unexpected end of file (scaling factor).");
  }
  double scaling_factor = std::stod(line);

  // Lines 3-5: Lattice vectors
  double v[3][3];
  for (int i = 0; i < 3; ++i) {
    if (!std::getline(myfile, line)) {
      throw std::runtime_error(
          "POSCAR: unexpected end of file (lattice vector).");
    }
    std::istringstream iss(line);
    if (!(iss >> v[i][0] >> v[i][1] >> v[i][2])) {
      throw std::runtime_error(
          "POSCAR: failed to parse lattice vector on line " +
          std::to_string(i + 3) + ".");
    }
  }

  // Apply scaling factor
  if (scaling_factor > 0.0) {
    // Positive: multiply all vectors by scaling factor
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        v[i][j] *= scaling_factor;
      }
    }
  } else if (scaling_factor < 0.0) {
    // Negative: interpret as target volume
    // V_target = |scaling_factor|
    // V_current = a . (b x c)
    double target_volume = std::abs(scaling_factor);
    double current_volume =
        std::abs(v[0][0] * (v[1][1] * v[2][2] - v[1][2] * v[2][1]) -
                 v[0][1] * (v[1][0] * v[2][2] - v[1][2] * v[2][0]) +
                 v[0][2] * (v[1][0] * v[2][1] - v[1][1] * v[2][0]));
    double scale = std::cbrt(target_volume / current_volume);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        v[i][j] *= scale;
      }
    }
  }

  correlation::core::Cell tempCell({v[0][0], v[0][1], v[0][2]},
                                   {v[1][0], v[1][1], v[1][2]},
                                   {v[2][0], v[2][1], v[2][2]});

  // Line 6: Species names (VASP 5+ format)
  // Could also be atom counts (VASP 4 format) — detect by checking if
  // the first token is a number.
  if (!std::getline(myfile, line)) {
    throw std::runtime_error(
        "POSCAR: unexpected end of file (species/counts).");
  }

  std::vector<std::string> species;
  std::vector<int> atom_counts;

  {
    std::istringstream iss(line);
    std::string token;
    bool is_species_line = false;

    // Peek at first token to determine if this is species or counts
    if (iss >> token) {
      try {
        std::stoi(token);
        // It's a number → VASP 4 format (no species line)
        // We don't have species names, use generic labels
        atom_counts.push_back(std::stoi(token));
        while (iss >> token) {
          atom_counts.push_back(std::stoi(token));
        }
        // Generate default species names
        for (size_t i = 0; i < atom_counts.size(); ++i) {
          species.push_back("Type" + std::to_string(i + 1));
        }
      } catch (...) {
        // Not a number → species names (VASP 5+)
        is_species_line = true;
        species.push_back(token);
        while (iss >> token) {
          species.push_back(token);
        }
      }
    }

    // If we read species names, next line is atom counts
    if (is_species_line) {
      if (!std::getline(myfile, line)) {
        throw std::runtime_error(
            "POSCAR: unexpected end of file (atom counts).");
      }
      std::istringstream count_iss(line);
      int count;
      while (count_iss >> count) {
        atom_counts.push_back(count);
      }
    }
  }

  if (species.size() != atom_counts.size()) {
    throw std::runtime_error("POSCAR: species count (" +
                             std::to_string(species.size()) +
                             ") does not match atom count entries (" +
                             std::to_string(atom_counts.size()) + ").");
  }

  int total_atoms = 0;
  for (int c : atom_counts) {
    total_atoms += c;
  }

  // Next line: "Selective dynamics" (optional) or coordinate type
  if (!std::getline(myfile, line)) {
    throw std::runtime_error(
        "POSCAR: unexpected end of file (coordinate type).");
  }

  // Check for Selective Dynamics
  char first_char = ' ';
  for (char c : line) {
    if (!std::isspace(static_cast<unsigned char>(c))) {
      first_char = c;
      break;
    }
  }

  if (first_char == 'S' || first_char == 's') {
    // Skip selective dynamics line, read next line for coordinate type
    if (!std::getline(myfile, line)) {
      throw std::runtime_error(
          "POSCAR: unexpected end of file (coordinate type after selective "
          "dynamics).");
    }
    // Re-read first char
    for (char c : line) {
      if (!std::isspace(static_cast<unsigned char>(c))) {
        first_char = c;
        break;
      }
    }
  }

  bool is_direct = (first_char == 'D' || first_char == 'd');
  // Cartesian: first_char == 'C' || first_char == 'c' || first_char == 'K' ||
  // first_char == 'k'

  // Read atom positions
  int species_idx = 0;
  int atoms_in_species = 0;
  for (int i = 0; i < total_atoms; ++i) {
    if (!std::getline(myfile, line)) {
      throw std::runtime_error(
          "POSCAR: unexpected end of file (atom position " +
          std::to_string(i + 1) + " of " + std::to_string(total_atoms) + ").");
    }

    // Determine which species this atom belongs to
    while (species_idx < static_cast<int>(atom_counts.size()) &&
           atoms_in_species >= atom_counts[species_idx]) {
      atoms_in_species = 0;
      species_idx++;
    }

    std::istringstream iss(line);
    double x, y, z;
    if (!(iss >> x >> y >> z)) {
      throw std::runtime_error("POSCAR: failed to parse atom coordinates on "
                               "atom " +
                               std::to_string(i + 1) + ".");
    }

    tempCell.addAtom(species[species_idx], {x, y, z});
    atoms_in_species++;
  }

  if (is_direct) {
    tempCell.wrapPositions();
  }

  return tempCell;
}

} // namespace correlation::readers
