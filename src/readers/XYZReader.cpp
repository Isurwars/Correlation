/**
 * @file XYZReader.cpp
 * @brief Implementation of the Extended XYZ file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#include "readers/XYZReader.hpp"
#include "readers/ReaderFactory.hpp"
#include "core/MappedFile.hpp"

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
// parseXYZFrame – parses a single frame from memory
// ---------------------------------------------------------------------------
correlation::core::Cell XYZReader::parseXYZFrame(const char *data, size_t size) {
  std::string content(data, size);
  std::istringstream stream(content);
  std::string line;

  // --- Line 1: atom count ---
  while (std::getline(stream, line)) {
    if (line.find_first_not_of(" \t\r\n") != std::string::npos) {
      break;
    }
  }

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
  if (!std::getline(stream, comment)) {
    throw std::runtime_error(
        "Invalid XYZ file: missing comment line after atom count");
  }

  correlation::core::Cell cell;

  // Try to parse Extended XYZ lattice from comment line.
  auto lattice = parseLattice(comment);
  if (lattice) {
    const auto &L = *lattice;
    correlation::math::Vector3<double> a(L[0], L[1], L[2]);
    correlation::math::Vector3<double> b(L[3], L[4], L[5]);
    correlation::math::Vector3<double> c(L[6], L[7], L[8]);
    cell = correlation::core::Cell(a, b, c);
  }

  // --- Lines 3..N+2: atom data ---
  for (int i = 0; i < num_atoms; ++i) {
    if (!std::getline(stream, line)) {
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

  return cell;
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
  return traj.getFrame(traj.getFrameCount() - 1);
}

// ---------------------------------------------------------------------------
// readTrajectory  – reads all frames from a (possibly multi-frame) XYZ file
// ---------------------------------------------------------------------------
correlation::core::Trajectory XYZReader::readTrajectory(
    const std::string &filename,
    std::function<void(float, const std::string &)> progress_callback) {

  if (progress_callback)
    progress_callback(0.0f, "Reading XYZ file...");

  auto mapped_file = std::make_shared<correlation::core::MappedFile>(filename);
  const char *data = mapped_file->data();
  const size_t total_size = mapped_file->size();

  std::vector<size_t> frame_offsets;
  size_t offset = 0;

  while (offset < total_size) {
    // Skip blank lines
    while (offset < total_size) {
      size_t line_end = offset;
      while (line_end < total_size && data[line_end] != '\n' && data[line_end] != '\r') {
        line_end++;
      }
      
      bool is_blank = true;
      for (size_t i = offset; i < line_end; ++i) {
        if (data[i] != ' ' && data[i] != '\t') {
          is_blank = false;
          break;
        }
      }
      
      if (!is_blank) {
        break;
      }
      
      if (line_end < total_size && data[line_end] == '\r') {
        line_end++;
      }
      if (line_end < total_size && data[line_end] == '\n') {
        line_end++;
      }
      offset = line_end;
    }
    
    if (offset >= total_size) {
      break;
    }
    
    size_t frame_start = offset;
    
    // Read the atom count line
    size_t line_end = offset;
    while (line_end < total_size && data[line_end] != '\n' && data[line_end] != '\r') {
      line_end++;
    }
    std::string atom_count_str(data + offset, line_end - offset);
    
    int num_atoms = 0;
    try {
      num_atoms = std::stoi(atom_count_str);
    } catch (...) {
      throw std::runtime_error("Invalid XYZ file: expected atom count, got: " + atom_count_str);
    }
    
    if (num_atoms <= 0) {
      throw std::runtime_error("Invalid XYZ file: non-positive atom count: " + atom_count_str);
    }
    
    if (line_end < total_size && data[line_end] == '\r') {
      line_end++;
    }
    if (line_end < total_size && data[line_end] == '\n') {
      line_end++;
    }
    offset = line_end;
    
    // Read the comment line
    if (offset >= total_size) {
      throw std::runtime_error("Invalid XYZ file: missing comment line after atom count");
    }
    line_end = offset;
    while (line_end < total_size && data[line_end] != '\n' && data[line_end] != '\r') {
      line_end++;
    }
    if (line_end < total_size && data[line_end] == '\r') {
      line_end++;
    }
    if (line_end < total_size && data[line_end] == '\n') {
      line_end++;
    }
    offset = line_end;
    
    // Read coordinate lines
    for (int i = 0; i < num_atoms; ++i) {
      if (offset >= total_size) {
        throw std::runtime_error("Invalid XYZ file: unexpected EOF while reading atom " + std::to_string(i + 1));
      }
      line_end = offset;
      while (line_end < total_size && data[line_end] != '\n' && data[line_end] != '\r') {
        line_end++;
      }
      if (line_end < total_size && data[line_end] == '\r') {
        line_end++;
      }
      if (line_end < total_size && data[line_end] == '\n') {
        line_end++;
      }
      offset = line_end;
    }
    
    frame_offsets.push_back(frame_start);

    // Report progress.
    if (progress_callback && total_size > 0) {
      float progress = static_cast<float>(offset) / static_cast<float>(total_size);
      progress_callback(progress, "Reading XYZ frames...");
    }
  }

  if (frame_offsets.empty()) {
    throw std::runtime_error("No frames found in XYZ file: " + filename);
  }

  frame_offsets.push_back(offset);

  if (progress_callback)
    progress_callback(1.0f, "XYZ file loaded.");

  auto parser = [](const char *d, size_t s) {
    return parseXYZFrame(d, s);
  };

  return correlation::core::Trajectory(mapped_file, std::move(frame_offsets), parser, 1.0);
}

} // namespace correlation::readers
