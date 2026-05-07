/**
 * @file XdatcarReader.cpp
 * @brief Implementation of the VASP XDATCAR trajectory file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */
#include "readers/XdatcarReader.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/ReaderFactory.hpp"

#include <algorithm>
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
    ReaderFactory::instance().registerReader(std::make_unique<XdatcarReader>());

correlation::core::Cell XdatcarReader::readStructure(
    const std::string &filename,
    std::function<void(float, const std::string &)> progress_callback) {
  auto frames = read(filename, progress_callback);
  if (frames.empty()) {
    throw std::runtime_error("XDATCAR: no frames found in file.");
  }
  return std::move(frames[0]);
}

correlation::core::Trajectory XdatcarReader::readTrajectory(
    const std::string &filename,
    std::function<void(float, const std::string &)> progress_callback) {
  return correlation::core::Trajectory(read(filename, progress_callback), 1.0);
}

std::vector<correlation::core::Cell> XdatcarReader::read(
    const std::string &file_name,
    std::function<void(float, const std::string &)> progress_callback) {
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" +
                             std::strerror(errno) + ").");
  }

  // Get file size for progress reporting
  myfile.seekg(0, std::ios::end);
  std::streampos file_size = myfile.tellg();
  myfile.seekg(0, std::ios::beg);
  std::streampos last_progress_pos = 0;
  size_t update_interval = file_size > 0 ? file_size / 100 : 1;

  std::string line;

  // Line 1: Comment
  if (!std::getline(myfile, line)) {
    throw std::runtime_error("XDATCAR: unexpected end of file (comment).");
  }

  // Line 2: Scaling factor
  if (!std::getline(myfile, line)) {
    throw std::runtime_error(
        "XDATCAR: unexpected end of file (scaling factor).");
  }
  double scaling_factor = std::stod(line);

  // Lines 3-5: Lattice vectors
  double v[3][3];
  for (int i = 0; i < 3; ++i) {
    if (!std::getline(myfile, line)) {
      throw std::runtime_error(
          "XDATCAR: unexpected end of file (lattice vector).");
    }
    std::istringstream iss(line);
    if (!(iss >> v[i][0] >> v[i][1] >> v[i][2])) {
      throw std::runtime_error(
          "XDATCAR: failed to parse lattice vector on line " +
          std::to_string(i + 3) + ".");
    }
  }

  // Apply scaling factor
  if (scaling_factor > 0.0) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        v[i][j] *= scaling_factor;
      }
    }
  } else if (scaling_factor < 0.0) {
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

  // Line 6: Species names
  if (!std::getline(myfile, line)) {
    throw std::runtime_error(
        "XDATCAR: unexpected end of file (species names).");
  }
  std::vector<std::string> species;
  {
    std::istringstream iss(line);
    std::string token;
    while (iss >> token) {
      species.push_back(token);
    }
  }

  // Line 7: Atom counts
  if (!std::getline(myfile, line)) {
    throw std::runtime_error("XDATCAR: unexpected end of file (atom counts).");
  }
  std::vector<int> atom_counts;
  {
    std::istringstream iss(line);
    int count;
    while (iss >> count) {
      atom_counts.push_back(count);
    }
  }

  if (species.size() != atom_counts.size()) {
    throw std::runtime_error(
        "XDATCAR: species count does not match atom count entries.");
  }

  int total_atoms = 0;
  for (int c : atom_counts) {
    total_atoms += c;
  }

  // Build species assignment list (which species each atom index belongs to)
  std::vector<std::string> atom_species;
  atom_species.reserve(total_atoms);
  for (size_t s = 0; s < species.size(); ++s) {
    for (int a = 0; a < atom_counts[s]; ++a) {
      atom_species.push_back(species[s]);
    }
  }

  // Read frames: each starts with "Direct configuration= N"
  std::vector<correlation::core::Cell> frames;

  while (std::getline(myfile, line)) {
    // Progress reporting
    if (progress_callback) {
      std::streampos current_pos = myfile.tellg();
      if (current_pos - last_progress_pos >
          static_cast<std::streamoff>(update_interval)) {
        float p =
            static_cast<float>(current_pos) / static_cast<float>(file_size);
        progress_callback(p, "Loading XDATCAR file...");
        last_progress_pos = current_pos;
      }
    }

    // Skip empty lines
    if (line.empty()) {
      continue;
    }

    // Check for frame separator: line starts with "Direct" or "direct"
    std::string trimmed = line;
    // Trim leading whitespace
    size_t start = trimmed.find_first_not_of(" \t");
    if (start == std::string::npos) {
      continue;
    }
    trimmed = trimmed.substr(start);

    // Detect "Direct" keyword (case-insensitive check on first 6 chars)
    if (trimmed.size() >= 6) {
      std::string prefix = trimmed.substr(0, 6);
      std::transform(prefix.begin(), prefix.end(), prefix.begin(), ::tolower);
      if (prefix == "direct") {
        // Start of a new frame — read total_atoms positions
        correlation::core::Cell tempCell({v[0][0], v[0][1], v[0][2]},
                                         {v[1][0], v[1][1], v[1][2]},
                                         {v[2][0], v[2][1], v[2][2]});

        bool valid_frame = true;
        for (int i = 0; i < total_atoms; ++i) {
          if (!std::getline(myfile, line)) {
            valid_frame = false;
            break;
          }
          std::istringstream iss(line);
          double x, y, z;
          if (!(iss >> x >> y >> z)) {
            valid_frame = false;
            break;
          }
          tempCell.addAtom(atom_species[i], {x, y, z});
        }

        if (valid_frame) {
          tempCell.wrapPositions(); // XDATCAR always uses Direct coordinates
          frames.push_back(std::move(tempCell));
        }
      }
    }
  }

  if (frames.empty()) {
    throw std::runtime_error("XDATCAR: no valid frames found in file.");
  }

  return frames;
}

} // namespace correlation::readers
