/**
 * @file GromacsReader.cpp
 * @brief Implementation of the GROMACS (.gro) file reader with memory-mapped
 *        lazy loading for multi-frame trajectories.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/GromacsReader.hpp"

#include "core/MappedFile.hpp"
#include "readers/ReaderFactory.hpp"
#include <math.h>

#include <cstring>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace correlation::readers {

namespace {
// Automatic registration
const bool registered = ReaderFactory::registerTypeSafe<GromacsReader>("GromacsReader");

// ---------------------------------------------------------------------------
// Helper: find end of line
// ---------------------------------------------------------------------------
inline size_t findLineEnd(const char *data, size_t total, size_t pos) {
  while (pos < total && data[pos] != '\n' && data[pos] != '\r') {
    ++pos;
  }
  return pos;
}

// ---------------------------------------------------------------------------
// Helper: skip past line ending
// ---------------------------------------------------------------------------
inline size_t skipLineEnding(const char *data, size_t total, size_t pos) {
  if (pos < total && data[pos] == '\r') {
    ++pos;
  }
  if (pos < total && data[pos] == '\n') {
    ++pos;
  }
  return pos;
}

// ---------------------------------------------------------------------------
// Helper: extract a line as std::string from [pos, lineEnd)
// ---------------------------------------------------------------------------
inline std::string extractLine(const char *data, size_t pos, size_t lineEnd) {
  return std::string(data + pos, lineEnd - pos);
}

size_t skipBlankLines(const char *data, size_t total_size, size_t offset) {
  while (offset < total_size) {
    size_t const line_end = findLineEnd(data, total_size, offset);
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
    offset = skipLineEnding(data, total_size, line_end);
  }
  return offset;
}

size_t scanNextFrame(const char *data, size_t total_size, size_t offset) {
  // Line 1: title — skip
  size_t line_end = findLineEnd(data, total_size, offset);
  offset = skipLineEnding(data, total_size, line_end);
  if (offset >= total_size) {
    return 0;
  }

  // Line 2: atom count
  line_end = findLineEnd(data, total_size, offset);
  std::string const count_str(data + offset, line_end - offset);
  offset = skipLineEnding(data, total_size, line_end);

  int num_atoms = 0;
  try {
    num_atoms = std::stoi(count_str);
  } catch (...) {
    return 0; // Not a valid frame
  }

  if (num_atoms <= 0) {
    return 0;
  }

  // Sanity check: the claimed atom count must be plausible given remaining data
  const size_t remaining_bytes = (offset < total_size) ? (total_size - offset) : 0;
  if (std::cmp_greater(num_atoms, remaining_bytes) || num_atoms > 100'000'000) {
    return 0;
  }

  // Skip N atom lines
  for (int i = 0; i < num_atoms; ++i) {
    if (offset >= total_size) {
      return 0; // Incomplete frame
    }
    line_end = findLineEnd(data, total_size, offset);
    offset = skipLineEnding(data, total_size, line_end);
  }

  // Skip box vector line
  if (offset < total_size) {
    line_end = findLineEnd(data, total_size, offset);
    offset = skipLineEnding(data, total_size, line_end);
  }

  return offset;
}

} // namespace

// ---------------------------------------------------------------------------
// parseGroFrame — parses a single .gro frame from memory
// ---------------------------------------------------------------------------
correlation::core::Cell GromacsReader::parseGroFrame(const char *data, size_t size) {
  size_t offset = 0;
  size_t lineEnd = 0;

  auto nextLine = [&]() -> std::string {
    lineEnd = findLineEnd(data, size, offset);
    std::string line = extractLine(data, offset, lineEnd);
    offset = skipLineEnding(data, size, lineEnd);
    return line;
  };

  // Line 1: Title/Comment
  nextLine(); // skip title

  // Line 2: Number of atoms
  std::string line = nextLine();
  int num_atoms = 0;
  try {
    num_atoms = std::stoi(line);
  } catch (...) {
    throw std::runtime_error("Invalid GROMACS file: non-numeric atom count");
  }

  // Guard against malformed files with absurd atom counts.  Each GRO atom line
  // is at least 44 characters, so we can estimate a reasonable upper bound from
  // the remaining data.  Additionally cap at 100 million as an absolute limit.
  if (num_atoms <= 0) {
    throw std::runtime_error("Invalid GROMACS file: non-positive atom count: " + std::to_string(num_atoms));
  }
  constexpr int kMaxAtomCount = 100'000'000;
  if (num_atoms > kMaxAtomCount) {
    throw std::runtime_error("Invalid GROMACS file: atom count exceeds limit: " + std::to_string(num_atoms));
  }
  const size_t remaining = (offset < size) ? (size - offset) : 0;
  // Each atom line needs at least ~20 bytes (even short lines). If the claimed
  // count cannot possibly fit in the remaining data, reject early.
  if (std::cmp_greater(num_atoms, remaining)) {
    throw std::runtime_error("Invalid GROMACS file: atom count (" + std::to_string(num_atoms) +
                             ") exceeds remaining data (" + std::to_string(remaining) + " bytes)");
  }

  correlation::core::Cell cell;

  // Lines 3 to num_atoms + 2: Atoms
  // Format: 5pos 5res 5atom 5id 8x 8y 8z 8vx 8vy 8vz (last 3 optional)
  for (int i = 0; i < num_atoms; ++i) {
    line = nextLine();
    if (line.length() < 44) {
      continue; // Basic coordinate check
    }

    std::string atom_name = line.substr(10, 5);
    // Trim whitespace
    atom_name.erase(0, atom_name.find_first_not_of(' '));
    atom_name.erase(atom_name.find_last_not_of(' ') + 1);

    // Extract element symbol from atom name
    std::string symbol;
    for (char const chr : atom_name) {
      if (std::isalpha(chr) != 0) {
        symbol += chr;
      } else {
        break;
      }
    }
    // Common GROMACS fixes
    if (symbol == "OW" || symbol == "HW") {
      symbol = symbol.substr(0, 1);
    }

    auto const pos_x = static_cast<real_t>(std::stod(line.substr(20, 8))) * 10.0; // nm to A
    auto const pos_y = static_cast<real_t>(std::stod(line.substr(28, 8))) * 10.0;
    auto const pos_z = static_cast<real_t>(std::stod(line.substr(36, 8))) * 10.0;

    cell.addAtom(symbol, correlation::math::Vector3<real_t>(pos_x, pos_y, pos_z));
  }

  // Last line: Box vectors (v1x v2y v3z v1y v1z v2x v2z v3x v3y)
  // Most common: 3 values (orthogonal box)
  if (offset < size) {
    line = nextLine();
    std::stringstream str_stream(line);
    real_t box_x = 0.0;
    real_t box_y = 0.0;
    real_t box_z = 0.0;
    if (str_stream >> box_x >> box_y >> box_z) {
      cell.setLatticeParameters({static_cast<real_t>(box_x * 10.0), static_cast<real_t>(box_y * 10.0),
                                 static_cast<real_t>(box_z * 10.0), static_cast<real_t>(90.0),
                                 static_cast<real_t>(90.0), static_cast<real_t>(90.0)});
    }
  }

  return cell;
}

// ---------------------------------------------------------------------------
// readStructure — returns the last frame
// ---------------------------------------------------------------------------
correlation::core::Cell
GromacsReader::readStructure(const std::string &filename,
                             std::function<void(float, const std::string &)> progress_callback) {

  auto traj = readTrajectory(filename, progress_callback);
  if (traj.getFrameCount() == 0) {
    throw std::runtime_error("No structure found in GROMACS file: " + filename);
  }
  return traj.getFrame(traj.getFrameCount() - 1);
}

// ---------------------------------------------------------------------------
// readTrajectory — memory-mapped lazy loading for multi-frame .gro
// ---------------------------------------------------------------------------
correlation::core::Trajectory
GromacsReader::readTrajectory(const std::string &filename,
                              std::function<void(float, const std::string &)> progress_callback) {

  if (progress_callback) {
    progress_callback(0.0F, "Reading GROMACS file...");
  }

  auto mapped_file = std::make_shared<correlation::core::MappedFile>(filename);
  const char *data = mapped_file->data();
  const size_t total_size = mapped_file->size();

  std::vector<size_t> frame_offsets;
  size_t offset = 0;

  while (offset < total_size) {
    // Skip any blank lines between frames
    offset = skipBlankLines(data, total_size, offset);
    if (offset >= total_size) {
      break;
    }

    size_t const frame_start = offset;
    size_t const frame_end = scanNextFrame(data, total_size, offset);
    if (frame_end == 0) {
      break; // Stop scanning on invalid/incomplete frame
    }
    offset = frame_end;

    frame_offsets.push_back(frame_start);

    // Report progress
    if (progress_callback && total_size > 0) {
      float const progress = static_cast<float>(offset) / static_cast<float>(total_size);
      progress_callback(progress, "Scanning GROMACS frames...");
    }
  }

  if (frame_offsets.empty()) {
    throw std::runtime_error("No frames found in GROMACS file: " + filename);
  }

  // Sentinel offset for the last frame's end
  frame_offsets.push_back(total_size);

  if (progress_callback) {
    progress_callback(1.0F, "GROMACS file loaded.");
  }

  auto parser = [](const char *data, size_t size) { return parseGroFrame(data, size); };

  return {mapped_file, std::move(frame_offsets), parser, 1.0};
}

} // namespace correlation::readers
