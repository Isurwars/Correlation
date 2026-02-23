// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "readers/LammpsDumpReader.hpp"

#include <cerrno>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace LammpsDumpReader {

Cell read(const std::string &file_name) {
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" +
                             std::strerror(errno) + ").");
  }

  Cell tempCell;
  std::string line;
  int num_atoms = 0;

  // Read the LAMMPS dump file line by line.
  // We look for "ITEM: TIMESTEP" to identify frames.
  // Note: This reader currently only keeps the LAST frame found in the file,
  // effectively treating it as a single structure reader if used via
  // readStructure.
  while (std::getline(myfile, line)) {
    // Find the beginning of a new timestep block
    if (line.find("ITEM: TIMESTEP") != std::string::npos) {
      // A new timestep is starting, so we clear the old cell data
      // This ensures we only keep the data for the *last* timestep
      tempCell = Cell();

      // The next section will be the header for this timestep
      // Read NUMBER OF ATOMS
      std::getline(myfile, line); // Timestep number (we ignore it)
      std::getline(myfile, line); // "ITEM: NUMBER OF ATOMS"
      std::getline(myfile, line); // The actual number
      num_atoms = std::stoi(line);

      // Read BOX BOUNDS
      std::getline(myfile, line); // "ITEM: BOX BOUNDS..."
      double xlo, xhi, ylo, yhi, zlo, zhi;
      std::getline(myfile, line);
      std::stringstream(line) >> xlo >> xhi;
      std::getline(myfile, line);
      std::stringstream(line) >> ylo >> yhi;
      std::getline(myfile, line);
      std::stringstream(line) >> zlo >> zhi;

      // Create an orthorhombic cell from the box bounds
      tempCell = Cell({xhi - xlo, 0.0, 0.0}, {0.0, yhi - ylo, 0.0},
                      {0.0, 0.0, zhi - zlo});

      // Read ATOMS header
      std::getline(myfile, line); // "ITEM: ATOMS id type x y z"

      // Now, read the specified number of atom lines
      for (int i = 0; i < num_atoms; ++i) {
        std::getline(myfile, line);
        std::stringstream atom_stream(line);
        int id, type;
        double x, y, z;

        // We assume the format "id type x y z"
        if (atom_stream >> id >> type >> x >> y >> z) {
          // Use the numeric type as the element string
          std::string element = std::to_string(type);
          tempCell.addAtom(element, {x, y, z});
        }
      }
    }
  }

  return tempCell;
}

} // namespace LammpsDumpReader
