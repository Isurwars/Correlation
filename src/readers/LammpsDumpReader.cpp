/**
 * @file LammpsDumpReader.cpp
 * @brief Implementation of the LAMMPS dump file reader with memory-mapped lazy
 *        loading.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */
#include "readers/LammpsDumpReader.hpp"
#include "core/MappedFile.hpp"
#include "readers/ReaderFactory.hpp"

#include <cstring>
#include <sstream>
#include <stdexcept>

namespace correlation::readers {

// Automatic registration
static bool registered = ReaderFactory::instance().registerReader(std::make_unique<LammpsDumpReader>()); // NOLINT(cert-err58-cpp, bugprone-throwing-static-initialization)

// ---------------------------------------------------------------------------
// Helper: advance past the current line ending (\r\n or \n)
// ---------------------------------------------------------------------------
static inline size_t skipLineEnding(const char *data, size_t total, size_t pos) {
  if (pos < total && data[pos] == '\r')
    ++pos;
  if (pos < total && data[pos] == '\n')
    ++pos;
  return pos;
}

// ---------------------------------------------------------------------------
// Helper: find end of the current line (position of \n or \r)
// ---------------------------------------------------------------------------
static inline size_t findLineEnd(const char *data, size_t total, size_t pos) {
  while (pos < total && data[pos] != '\n' && data[pos] != '\r')
    ++pos;
  return pos;
}

// ---------------------------------------------------------------------------
// Helper: extract a line as std::string from [pos, lineEnd)
// ---------------------------------------------------------------------------
static inline std::string extractLine(const char *data, size_t pos, size_t lineEnd) {
  return std::string(data + pos, lineEnd - pos);
}

// ---------------------------------------------------------------------------
// parseDumpFrame — parses a single frame from a memory region
// ---------------------------------------------------------------------------
correlation::core::Cell LammpsDumpReader::parseDumpFrame(const char *data, size_t size) {
  size_t offset = 0;
  size_t lineEnd = 0;

  auto nextLine = [&]() -> std::string {
    lineEnd = findLineEnd(data, size, offset);
    std::string line = extractLine(data, offset, lineEnd);
    offset = skipLineEnding(data, size, lineEnd);
    return line;
  };

  // --- ITEM: TIMESTEP ---
  std::string line = nextLine(); // "ITEM: TIMESTEP"
  nextLine();                    // timestep value (ignored)

  // --- NUMBER OF ATOMS ---
  nextLine(); // "ITEM: NUMBER OF ATOMS"
  line = nextLine();
  int num_atoms = 0;
  try {
    num_atoms = std::stoi(line);
  } catch (...) {
    throw std::runtime_error("Failed to parse atom count in LAMMPS dump frame");
  }
  if (num_atoms <= 0) {
    throw std::runtime_error("Invalid LAMMPS dump frame: non-positive atom count: " + std::to_string(num_atoms));
  }

  // --- BOX BOUNDS ---
  line = nextLine(); // "ITEM: BOX BOUNDS ..."
  const bool triclinic = (line.find("xy") != std::string::npos);

  double xlo = 0.0, xhi = 0.0, ylo = 0.0, yhi = 0.0, zlo = 0.0, zhi = 0.0;
  double xy = 0.0, xz = 0.0, yz = 0.0;

  if (triclinic) {
    line = nextLine();
    std::stringstream(line) >> xlo >> xhi >> xy;
    line = nextLine();
    std::stringstream(line) >> ylo >> yhi >> xz;
    line = nextLine();
    std::stringstream(line) >> zlo >> zhi >> yz;
  } else {
    line = nextLine();
    std::stringstream(line) >> xlo >> xhi;
    line = nextLine();
    std::stringstream(line) >> ylo >> yhi;
    line = nextLine();
    std::stringstream(line) >> zlo >> zhi;
  }

  // Build the Cell from the box definition.
  correlation::core::Cell frame;
  if (triclinic) {
    const double lx = xhi - xlo;
    const double ly = yhi - ylo;
    const double lz = zhi - zlo;
    frame = correlation::core::Cell({lx, 0.0, 0.0}, {xy, ly, 0.0}, {xz, yz, lz});
  } else {
    frame = correlation::core::Cell({xhi - xlo, 0.0, 0.0}, {0.0, yhi - ylo, 0.0}, {0.0, 0.0, zhi - zlo});
  }

  // --- ATOMS header — discover column layout ---
  line = nextLine(); // "ITEM: ATOMS ..."

  std::istringstream header_ss(line);
  std::string token;
  header_ss >> token >> token; // consume "ITEM:" and "ATOMS"

  std::vector<std::string> col_names;
  while (header_ss >> token) {
    col_names.push_back(token);
  }

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

  if (col_x == -1 || col_y == -1 || col_z == -1) {
    return frame; // Cannot parse positions — return empty frame
  }

  const auto &lv = frame.latticeVectors();

  // --- Atom lines ---
  for (int i = 0; i < num_atoms; ++i) {
    if (offset >= size)
      break;
    line = nextLine();

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
      element_symbol = fields[col_element];
    } else if (col_type >= 0 && col_type < num_fields) {
      element_symbol = fields[col_type];
    } else {
      element_symbol = (col_id >= 0 && col_id < num_fields) ? fields[col_id] : std::to_string(i + 1);
    }

    // Convert scaled (fractional) coordinates to Cartesian if needed.
    correlation::math::Vector3<double> pos;
    if (scaled_coords) {
      // pos = x*a + y*b + z*c  (correct column-vector combination of lattice rows)
      // i.e. pos[k] = x*lv[0][k] + y*lv[1][k] + z*lv[2][k]
      pos = {x * lv[0][0] + y * lv[1][0] + z * lv[2][0], x * lv[0][1] + y * lv[1][1] + z * lv[2][1],
             x * lv[0][2] + y * lv[1][2] + z * lv[2][2]};
    } else {
      pos = {x, y, z};
    }

    frame.addAtom(element_symbol, pos);
  }

  return frame;
}

// ---------------------------------------------------------------------------
// readStructure — returns the first frame
// ---------------------------------------------------------------------------
correlation::core::Cell
LammpsDumpReader::readStructure(const std::string &filename,
                                std::function<void(float, const std::string &)> progress_callback) {
  auto traj = readTrajectory(filename, progress_callback);
  if (traj.getFrameCount() == 0) {
    throw std::runtime_error("No frames found in LAMMPS dump file: " + filename);
  }
  return traj.getFrame(0);
}

// ---------------------------------------------------------------------------
// readTrajectory — memory-mapped lazy loading
// ---------------------------------------------------------------------------
correlation::core::Trajectory
LammpsDumpReader::readTrajectory(const std::string &filename,
                                 std::function<void(float, const std::string &)> progress_callback) {

  if (progress_callback)
    progress_callback(0.0f, "Reading LAMMPS dump file...");

  auto mapped_file = std::make_shared<correlation::core::MappedFile>(filename);
  const char *data = mapped_file->data();
  const size_t total_size = mapped_file->size();

  // First pass: scan for "ITEM: TIMESTEP" to find frame byte offsets.
  std::vector<size_t> frame_offsets;

  // We search for the pattern "ITEM: TIMESTEP" at the start of a line.
  const char *needle = "ITEM: TIMESTEP";
  const size_t needle_len = std::strlen(needle);

  size_t pos = 0;
  while (pos < total_size) {
    // Check if current position starts with the needle
    bool at_line_start = (pos == 0) || (pos > 0 && (data[pos - 1] == '\n'));
    if (at_line_start && pos + needle_len <= total_size && std::memcmp(data + pos, needle, needle_len) == 0) {
      frame_offsets.push_back(pos);

      // Report progress.
      if (progress_callback && total_size > 0) {
        float progress = static_cast<float>(pos) / static_cast<float>(total_size);
        progress_callback(progress, "Scanning LAMMPS dump frames...");
      }
    }
    // Advance to the next line
    pos = findLineEnd(data, total_size, pos);
    pos = skipLineEnding(data, total_size, pos);
  }

  if (frame_offsets.empty()) {
    throw std::runtime_error("No frames found in LAMMPS dump file: " + filename);
  }

  // Add sentinel offset for the last frame's end
  frame_offsets.push_back(total_size);

  if (progress_callback)
    progress_callback(1.0f, "LAMMPS dump file loaded.");

  auto parser = [](const char *d, size_t s) { return parseDumpFrame(d, s); };

  return correlation::core::Trajectory(mapped_file, std::move(frame_offsets), parser, 1.0);
}

} // namespace correlation::readers
