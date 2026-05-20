/**
 * @file XYZReader.cpp
 * @brief Implementation of the Extended XYZ file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "readers/XYZReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace correlation::readers {

namespace {
bool registered =
    ReaderFactory::instance().registerReader(std::make_unique<XYZReader>());
}

// ---------------------------------------------------------------------------
// Lattice parser  (Extended XYZ convention)
// ---------------------------------------------------------------------------
std::optional<std::array<double, 9>>
XYZReader::parseLattice(const std::string &comment) {
  // Search for  Lattice="v1x v1y v1z v2x v2y v2z v3x v3y v3z"
  const std::string key = "Lattice=\"";
  auto pos = comment.find(key);
  if (pos == std::string::npos)
    return std::nullopt;

  auto start = pos + key.size();
  auto end = comment.find('"', start);
  if (end == std::string::npos)
    return std::nullopt;

  std::string values = comment.substr(start, end - start);
  std::istringstream iss(values);
  std::array<double, 9> lattice{};
  for (int i = 0; i < 9; ++i) {
    if (!(iss >> lattice[i]))
      return std::nullopt;
  }
  return lattice;
}

// ---------------------------------------------------------------------------
// readStructure  – returns the last frame
// ---------------------------------------------------------------------------
correlation::core::Cell XYZReader::readStructure(
    const std::string &filename,
    std::function<void(float, const std::string &)> progress_callback) {

  auto traj = readTrajectory(filename, progress_callback);
  if (traj.getFrameCount() == 0) {
    throw std::runtime_error("No structure found in XYZ file: " + filename);
  }
  return traj.getFrames()[traj.getFrameCount() - 1];
}

// ---------------------------------------------------------------------------
// readTrajectory  – reads all frames from a (possibly multi-frame) XYZ file
// ---------------------------------------------------------------------------
correlation::core::Trajectory XYZReader::readTrajectory(
    const std::string &filename,
    std::function<void(float, const std::string &)> progress_callback) {

  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  if (progress_callback)
    progress_callback(0.0f, "Reading XYZ file...");

  // Determine file size for progress reporting.
  file.seekg(0, std::ios::end);
  const auto file_size = file.tellg();
  file.seekg(0, std::ios::beg);

  std::vector<correlation::core::Cell> frames;
  std::string line;

  while (std::getline(file, line)) {
    // --- Line 1: atom count ---
    // Skip blank lines between frames.
    if (line.find_first_not_of(" \t\r\n") == std::string::npos)
      continue;

    int num_atoms = 0;
    try {
      num_atoms = std::stoi(line);
    } catch (...) {
      throw std::runtime_error(
          "Invalid XYZ file: expected atom count, got: " + line);
    }

    if (num_atoms <= 0) {
      throw std::runtime_error(
          "Invalid XYZ file: non-positive atom count: " + line);
    }

    // --- Line 2: comment / Extended XYZ header ---
    std::string comment;
    if (!std::getline(file, comment)) {
      throw std::runtime_error(
          "Invalid XYZ file: missing comment line after atom count");
    }

    correlation::core::Cell cell;

    // Try to parse Extended XYZ lattice from comment line.
    auto lattice = parseLattice(comment);
    if (lattice) {
      const auto &L = *lattice;
      // L is row-major: a = (L[0], L[1], L[2]), b = (L[3], L[4], L[5]),
      //                 c = (L[6], L[7], L[8])
      correlation::math::Vector3<double> a(L[0], L[1], L[2]);
      correlation::math::Vector3<double> b(L[3], L[4], L[5]);
      correlation::math::Vector3<double> c(L[6], L[7], L[8]);
      cell = correlation::core::Cell(a, b, c);
    }

    // --- Lines 3..N+2: atom data ---
    for (int i = 0; i < num_atoms; ++i) {
      if (!std::getline(file, line)) {
        throw std::runtime_error(
            "Invalid XYZ file: unexpected EOF while reading atom " +
            std::to_string(i + 1));
      }

      std::istringstream iss(line);
      std::string symbol;
      double x, y, z;
      if (!(iss >> symbol >> x >> y >> z)) {
        throw std::runtime_error(
            "Invalid XYZ file: malformed atom line: " + line);
      }

      cell.addAtom(symbol, correlation::math::Vector3<double>(x, y, z));
    }

    frames.push_back(std::move(cell));

    // Report progress.
    if (progress_callback && file_size > 0) {
      float progress = static_cast<float>(file.tellg()) /
                        static_cast<float>(file_size);
      progress_callback(progress, "Reading XYZ frames...");
    }
  }

  if (frames.empty()) {
    throw std::runtime_error("No frames found in XYZ file: " + filename);
  }

  if (progress_callback)
    progress_callback(1.0f, "XYZ file loaded.");

  return correlation::core::Trajectory(frames, 1.0);
}

} // namespace correlation::readers
