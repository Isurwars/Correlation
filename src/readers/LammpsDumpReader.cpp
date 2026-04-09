/**
 * @file LammpsDumpReader.cpp
 * @brief Implementation of the LAMMPS dump file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */
#include "readers/LammpsDumpReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <cerrno>
#include <cstring>
#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Cell.hpp"
#include "Trajectory.hpp"

namespace correlation::readers {

// Automatic registration
static bool registered = ReaderFactory::instance().registerReader(
    std::make_unique<LammpsDumpReader>());

Cell LammpsDumpReader::readStructure(
    const std::string &filename,
    std::function<void(float, const std::string &)> progress_callback) {
  auto frames = read(filename, progress_callback);
  if (frames.empty()) {
    throw std::runtime_error("No frames found in LAMMPS dump file: " +
                             filename);
  }
  return frames.front();
}

Trajectory LammpsDumpReader::readTrajectory(
    const std::string &filename,
    std::function<void(float, const std::string &)> progress_callback) {
  return Trajectory(read(filename, progress_callback), 1.0);
}

// ---------------------------------------------------------------------------
// Core parser
// ---------------------------------------------------------------------------

std::vector<Cell> LammpsDumpReader::read(
    const std::string &file_name,
    std::function<void(float, const std::string &)> progress_callback) {

  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" +
                             std::strerror(errno) + ").");
  }

  // Determine file size for progress reporting.
  myfile.seekg(0, std::ios::end);
  const std::streampos file_size = myfile.tellg();
  myfile.seekg(0, std::ios::beg);
  std::streampos last_progress_pos = 0;
  const size_t update_interval =
      (file_size > 100) ? static_cast<size_t>(file_size) / 100 : 1;

  std::vector<Cell> frames;
  std::string line;

  while (std::getline(myfile, line)) {
    // -----------------------------------------------------------------------
    // Progress callback
    // -----------------------------------------------------------------------
    if (progress_callback) {
      std::streampos current_pos = myfile.tellg();
      if (current_pos - last_progress_pos >
          static_cast<std::streamoff>(update_interval)) {
        float p =
            static_cast<float>(current_pos) / static_cast<float>(file_size);
        progress_callback(p, "Loading LAMMPS dump file...");
        last_progress_pos = current_pos;
      }
    }

    // -----------------------------------------------------------------------
    // Each frame begins with "ITEM: TIMESTEP"
    // -----------------------------------------------------------------------
    if (line.find("ITEM: TIMESTEP") == std::string::npos) {
      continue;
    }

    // --- Timestep value (ignored, but consumed) ---
    std::getline(myfile, line);

    // --- NUMBER OF ATOMS ---
    std::getline(myfile, line); // "ITEM: NUMBER OF ATOMS"
    std::getline(myfile, line); // number
    int num_atoms = 0;
    try {
      num_atoms = std::stoi(line);
    } catch (...) {
      throw std::runtime_error("Failed to parse atom count in LAMMPS dump: " +
                               file_name);
    }

    // --- BOX BOUNDS ---
    // Format: "ITEM: BOX BOUNDS {pp|ss|ff} {pp|ss|ff} {pp|ss|ff}"
    // May optionally have "xy xz yz" for triclinic boxes.
    std::getline(myfile, line); // "ITEM: BOX BOUNDS ..."
    const bool triclinic = (line.find("xy") != std::string::npos);

    double xlo, xhi, ylo, yhi, zlo, zhi;
    double xy = 0.0, xz = 0.0, yz = 0.0;

    if (triclinic) {
      std::getline(myfile, line);
      std::stringstream(line) >> xlo >> xhi >> xy;
      std::getline(myfile, line);
      std::stringstream(line) >> ylo >> yhi >> xz;
      std::getline(myfile, line);
      std::stringstream(line) >> zlo >> zhi >> yz;
    } else {
      std::getline(myfile, line);
      std::stringstream(line) >> xlo >> xhi;
      std::getline(myfile, line);
      std::stringstream(line) >> ylo >> yhi;
      std::getline(myfile, line);
      std::stringstream(line) >> zlo >> zhi;
    }

    // Build the Cell from the box definition.
    // For orthorhombic boxes this is straightforward. For triclinic we pass
    // the full lattice vectors following the LAMMPS convention.
    Cell frame;
    if (triclinic) {
      // LAMMPS triclinic lattice vectors:
      //   a = (lx,  0,  0)
      //   b = (xy, ly,  0)
      //   c = (xz, yz, lz)
      const double lx = xhi - xlo;
      const double ly = yhi - ylo;
      const double lz = zhi - zlo;
      frame = Cell({lx, 0.0, 0.0}, {xy, ly, 0.0}, {xz, yz, lz});
    } else {
      frame = Cell({xhi - xlo, 0.0, 0.0}, {0.0, yhi - ylo, 0.0},
                   {0.0, 0.0, zhi - zlo});
    }

    // --- ATOMS header — discover column layout ---
    // e.g. "ITEM: ATOMS id type x y z"
    //      "ITEM: ATOMS id element x y z"
    //      "ITEM: ATOMS id type xs ys zs"
    std::getline(myfile, line); // "ITEM: ATOMS ..."

    // Parse column names from the header.
    std::istringstream header_ss(line);
    std::string token;
    header_ss >> token >> token; // consume "ITEM:" and "ATOMS"

    std::vector<std::string> col_names;
    while (header_ss >> token) {
      col_names.push_back(token);
    }

    // Locate the columns we care about. -1 means not found.
    int col_id = -1, col_type = -1, col_element = -1;
    int col_x = -1, col_y = -1, col_z = -1;
    bool scaled_coords = false;

    for (int i = 0; i < static_cast<int>(col_names.size()); ++i) {
      const auto &cn = col_names[i];
      if (cn == "id")
        col_id = i;
      else if (cn == "type")
        col_type = i;
      else if (cn == "element")
        col_element = i;
      else if (cn == "x")
        col_x = i;
      else if (cn == "y")
        col_y = i;
      else if (cn == "z")
        col_z = i;
      else if (cn == "xs" || cn == "xsu") {
        col_x = i;
        scaled_coords = true;
      } else if (cn == "ys" || cn == "ysu") {
        col_y = i;
        scaled_coords = true;
      } else if (cn == "zs" || cn == "zsu") {
        col_z = i;
        scaled_coords = true;
      }
    }

    // We need at minimum x/y/z columns and some form of element identity.
    if (col_x == -1 || col_y == -1 || col_z == -1) {
      // Skip this frame — cannot parse positions.
      for (int i = 0; i < num_atoms; ++i)
        std::getline(myfile, line);
      continue;
    }

    // Pre-cache the lattice matrix for scaled-coordinate conversion.
    const auto &lv = frame.latticeVectors();

    // --- Atom lines ---
    for (int i = 0; i < num_atoms; ++i) {
      if (!std::getline(myfile, line))
        break;

      std::istringstream atom_ss(line);
      std::vector<std::string> fields;
      while (atom_ss >> token) {
        fields.push_back(token);
      }

      const int num_fields = static_cast<int>(fields.size());
      if (col_x >= num_fields || col_y >= num_fields || col_z >= num_fields)
        continue;

      double x = std::stod(fields[col_x]);
      double y = std::stod(fields[col_y]);
      double z = std::stod(fields[col_z]);

      // Determine element symbol.
      std::string element_symbol;
      if (col_element >= 0 && col_element < num_fields) {
        // Prefer explicit element name (e.g. "Si", "O").
        element_symbol = fields[col_element];
      } else if (col_type >= 0 && col_type < num_fields) {
        // Fall back to numeric type label ("1", "2", ...).
        element_symbol = fields[col_type];
      } else {
        // Last resort: use atom id.
        element_symbol = (col_id >= 0 && col_id < num_fields)
                             ? fields[col_id]
                             : std::to_string(i + 1);
      }

      // Convert scaled (fractional) coordinates to Cartesian if needed.
      correlation::math::Vector3<double> pos;
      if (scaled_coords) {
        // Cartesian = x*a + y*b + z*c  (column vectors)
        pos = {x * lv[0][0] + y * lv[0][1] + z * lv[0][2],
               x * lv[1][0] + y * lv[1][1] + z * lv[1][2],
               x * lv[2][0] + y * lv[2][1] + z * lv[2][2]};
      } else {
        pos = {x, y, z};
      }

      frame.addAtom(element_symbol, pos);
    }

    if (!frame.isEmpty()) {
      frames.push_back(std::move(frame));
    }
  }

  return frames;
}

} // namespace correlation::readers
