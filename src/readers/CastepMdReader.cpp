// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "readers/CastepMdReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <cerrno>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <memory>
#include <functional>

#include "Cell.hpp"
#include "Trajectory.hpp"

namespace FileReader {

// Conversion factor from Bohr (Hartree atomic units) to Angstroms
constexpr double BOHR_TO_ANGSTROM = 0.529177210903;

// Automatic registration
static bool registered = ReaderFactory::instance().registerReader(
    std::make_unique<CastepMdReader>()
);

Cell CastepMdReader::readStructure(const std::string &filename,
                                    std::function<void(float, const std::string &)>
                                        progress_callback) {
  return CastepMdReader::read(filename, progress_callback).at(0);
}

Trajectory CastepMdReader::readTrajectory(const std::string &filename,
                                           std::function<void(float, const std::string &)>
                                               progress_callback) {
  return Trajectory(CastepMdReader::read(filename, progress_callback), 1.0);
}

std::vector<Cell>
CastepMdReader::read(const std::string &file_name,
                     std::function<void(float, const std::string &)> progress_callback) {
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" +
                             std::strerror(errno) + ").");
  }

  std::vector<Cell> frames;
  std::string line;

  // Skip the header until we reach the blank line
  bool in_header = true;
  while (in_header && std::getline(myfile, line)) {
    if (line.empty() || line.find_first_not_of(' ') == std::string::npos) {
      in_header = false;
    }
  }

  // Temporary storage for building the current cell
  Cell tempCell;
  double current_energy = 0.0;
  bool cell_has_atoms = false;

  // Parse frames
  myfile.seekg(0, std::ios::end);
  std::streampos file_size = myfile.tellg();
  myfile.clear();
  myfile.seekg(0, std::ios::beg);
  std::streampos last_progress_pos = 0;
  size_t update_interval = file_size / 100;

  // re-skip header since we seekg to beg
  in_header = true;
  while (in_header && std::getline(myfile, line)) {
    if (line.empty() || line.find_first_not_of(' ') == std::string::npos) {
      in_header = false;
    }
  }

  while (std::getline(myfile, line)) {
    if (progress_callback) {
      std::streampos current_pos = myfile.tellg();
      if (current_pos - last_progress_pos > update_interval) {
        float p =
            static_cast<float>(current_pos) / static_cast<float>(file_size);
        progress_callback(p, "Loading CASTEP MD file...");
        last_progress_pos = current_pos;
      }
    }
    // If line represents energies/time
    if (line.find("<-- E") != std::string::npos) {
      // For CASTEP MD, the frame usually starts with time (single float on a
      // line without tag) We skip it and read E. If we have atoms in tempCell,
      // it means this is a new frame starting
      if (cell_has_atoms) {
        frames.push_back(std::move(tempCell));
        tempCell = Cell();
        cell_has_atoms = false;
      }

      // We can extract energy if needed...
      std::stringstream ss(line);
      if (ss >> current_energy) {
        // successfully read energy
      }
      continue;
    }

    if (line.find("<-- h") != std::string::npos &&
        line.find("<-- hv") == std::string::npos) {
      // Lattice vectors h are given row by row in Bohr
      // The first <-- h is row 1
      std::array<double, 3> h1, h2, h3;
      std::stringstream(line) >> h1[0] >> h1[1] >> h1[2];

      std::getline(myfile, line); // row 2
      std::stringstream(line) >> h2[0] >> h2[1] >> h2[2];

      std::getline(myfile, line); // row 3
      std::stringstream(line) >> h3[0] >> h3[1] >> h3[2];

      // Convert Bohr to Angstroms
      for (int i = 0; i < 3; ++i) {
        h1[i] *= BOHR_TO_ANGSTROM;
        h2[i] *= BOHR_TO_ANGSTROM;
        h3[i] *= BOHR_TO_ANGSTROM;
      }

      tempCell = Cell({h1[0], h1[1], h1[2]}, {h2[0], h2[1], h2[2]},
                      {h3[0], h3[1], h3[2]});
      continue;
    }

    if (line.find("<-- R") != std::string::npos) {
      // Positions R are given in Bohr
      // Format: Symbol ID x y z <-- R
      std::stringstream ss(line);
      std::string symbol, tag;
      int id;
      double x, y, z;
      if (ss >> symbol >> id >> x >> y >> z) {
        tempCell.addAtom(symbol, linalg::Vector3<double>(x * BOHR_TO_ANGSTROM,
                                                        y * BOHR_TO_ANGSTROM,
                                                        z * BOHR_TO_ANGSTROM));
        tempCell.setEnergy(
            current_energy); // Assign energy once per atom or frame
        cell_has_atoms = true;
      }
      continue;
    }
  }

  // Push the final frame if it has atoms
  if (cell_has_atoms && !tempCell.isEmpty()) {
    frames.push_back(std::move(tempCell));
  }

  return frames;
}

} // namespace FileReader
