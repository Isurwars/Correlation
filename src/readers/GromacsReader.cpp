/**
 * @file GromacsReader.cpp
 * @brief Implementation of the GROMACS (.gro) file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "readers/GromacsReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace correlation::readers {

namespace {
bool registered =
    ReaderFactory::instance().registerReader(std::make_unique<GromacsReader>());
}

correlation::core::Cell GromacsReader::readStructure(
    const std::string &filename,
    std::function<void(float, const std::string &)> progress_callback) {

  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  if (progress_callback)
    progress_callback(0.1f, "Reading GROMACS structure...");

  std::string line;

  // Line 1: Title/Comment
  if (!std::getline(file, line))
    throw std::runtime_error("Invalid GROMACS file: empty");

  // Line 2: Number of atoms
  if (!std::getline(file, line))
    throw std::runtime_error("Invalid GROMACS file: missing atom count");
  int num_atoms = 0;
  try {
    num_atoms = std::stoi(line);
  } catch (...) {
    throw std::runtime_error("Invalid GROMACS file: non-numeric atom count");
  }

  correlation::core::Cell cell;

  // Lines 3 to num_atoms + 2: Atoms
  // Format: 5pos 5res 5atom 5id 8x 8y 8z 8vx 8vy 8vz (last 3 optional)
  for (int i = 0; i < num_atoms; ++i) {
    if (!std::getline(file, line))
      throw std::runtime_error("Invalid GROMACS file: unexpected EOF");

    if (line.length() < 44)
      continue; // Basic coordinate check

    std::string res_name = line.substr(5, 5);
    std::string atom_name = line.substr(10, 5);
    // Trim whitespace
    atom_name.erase(0, atom_name.find_first_not_of(" "));
    atom_name.erase(atom_name.find_last_not_of(" ") + 1);

    // GROMACS atom names often start with the element (e.g., "OW", "HW1", "CL")
    // We need to extract the element symbol.
    // Simple heuristic: first character, or first two if it matches an element.
    // For now, let's just take the first alphabetic characters.
    std::string symbol = "";
    for (char c : atom_name) {
      if (std::isalpha(c))
        symbol += c;
      else
        break;
    }
    // Common GROMACS fixes
    if (symbol == "OW" || symbol == "HW")
      symbol = symbol.substr(0, 1);

    double x = std::stod(line.substr(20, 8)) * 10.0; // nm to A
    double y = std::stod(line.substr(28, 8)) * 10.0;
    double z = std::stod(line.substr(36, 8)) * 10.0;

    cell.addAtom(symbol, correlation::math::Vector3<double>(x, y, z));

    if (progress_callback && i % 1000 == 0) {
      progress_callback(0.1f + 0.8f * (static_cast<float>(i) / num_atoms),
                        "Parsing atoms...");
    }
  }

  // Last line: Box vectors (v1x v2y v3z v1y v1z v2x v2z v3x v3y)
  // Most common: 3 values (orthogonal box)
  if (std::getline(file, line)) {
    std::stringstream ss(line);
    double bx, by, bz;
    if (ss >> bx >> by >> bz) {
      // GROMACS uses nm; convert to Angstroms. Assume orthogonal box.
      cell.setLatticeParameters(
          {bx * 10.0, by * 10.0, bz * 10.0, 90.0, 90.0, 90.0});
    }
  }

  if (progress_callback)
    progress_callback(1.0f, "GROMACS file loaded.");
  return cell;
}

correlation::core::Trajectory GromacsReader::readTrajectory(
    const std::string &filename,
    std::function<void(float, const std::string &)> progress_callback) {
  // For .gro files, we usually only have one frame.
  // If there are more, they just repeat the structure.
  // For now, return a single frame trajectory.
  std::vector<correlation::core::Cell> frames;
  frames.push_back(readStructure(filename, progress_callback));
  return correlation::core::Trajectory(frames, 1.0);
}

} // namespace correlation::readers
