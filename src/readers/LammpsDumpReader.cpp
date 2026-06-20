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

namespace {

// Automatic registration
// NOLINTNEXTLINE(cert-err58-cpp, bugprone-throwing-static-initialization)
const bool registered = ReaderFactory::instance().registerReader(std::make_unique<LammpsDumpReader>());

// ---------------------------------------------------------------------------
// Helper: advance past the current line ending (\r\n or \n)
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
// Helper: find end of the current line (position of \n or \r)
// ---------------------------------------------------------------------------
inline size_t findLineEnd(const char *data, size_t total, size_t pos) {
  while (pos < total && data[pos] != '\n' && data[pos] != '\r') {
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

// Struct to encapsulate parsing a single LAMMPS dump frame.
struct LammpsFrameParser {
  const char *data = nullptr;
  size_t size = 0;
  size_t offset = 0;
  size_t lineEnd = 0;

  std::string nextLine() {
    lineEnd = findLineEnd(data, size, offset);
    std::string line = extractLine(data, offset, lineEnd);
    offset = skipLineEnding(data, size, lineEnd);
    return line;
  }

  int parseAtomCount() {
    nextLine(); // "ITEM: NUMBER OF ATOMS"
    std::string line = nextLine();
    int num_atoms = 0;
    try {
      num_atoms = std::stoi(line);
    } catch (...) {
      throw std::runtime_error("Failed to parse atom count in LAMMPS dump frame");
    }
    if (num_atoms <= 0) {
      throw std::runtime_error("Invalid LAMMPS dump frame: non-positive atom count: " + std::to_string(num_atoms));
    }
    return num_atoms;
  }

  correlation::core::Cell parseBoxBounds() {
    std::string line = nextLine(); // "ITEM: BOX BOUNDS ..."
    const bool triclinic = (line.contains("xy"));

    double xlo = 0.0;
    double xhi = 0.0;
    double ylo = 0.0;
    double yhi = 0.0;
    double zlo = 0.0;
    double zhi = 0.0;
    double x_y = 0.0;
    double x_z = 0.0;
    double y_z = 0.0;

    if (triclinic) {
      std::stringstream(nextLine()) >> xlo >> xhi >> x_y;
      std::stringstream(nextLine()) >> ylo >> yhi >> x_z;
      std::stringstream(nextLine()) >> zlo >> zhi >> y_z;
      const double l_x = xhi - xlo;
      const double l_y = yhi - ylo;
      const double l_z = zhi - zlo;
      return correlation::core::Cell({l_x, 0.0, 0.0}, {x_y, l_y, 0.0}, {x_z, y_z, l_z});
    }
    std::stringstream(nextLine()) >> xlo >> xhi;
    std::stringstream(nextLine()) >> ylo >> yhi;
    std::stringstream(nextLine()) >> zlo >> zhi;
    return correlation::core::Cell({xhi - xlo, 0.0, 0.0}, {0.0, yhi - ylo, 0.0}, {0.0, 0.0, zhi - zlo});
  }

  struct ColumnLayout {
    int col_id = -1;
    int col_type = -1;
    int col_element = -1;
    int col_x = -1;
    int col_y = -1;
    int col_z = -1;
    bool scaled_coords = false;
  };

  ColumnLayout parseAtomsHeader() {
    std::string line = nextLine(); // "ITEM: ATOMS ..."

    std::istringstream header_ss(line);
    std::string token;
    header_ss >> token >> token; // consume "ITEM:" and "ATOMS"

    std::vector<std::string> col_names;
    while (header_ss >> token) {
      col_names.push_back(token);
    }

    ColumnLayout layout;
    for (int i = 0; i < static_cast<int>(col_names.size()); ++i) {
      const std::string &col_name = col_names[i];
      if (col_name == "id") {
        layout.col_id = i;
      } else if (col_name == "type") {
        layout.col_type = i;
      } else if (col_name == "element") {
        layout.col_element = i;
      } else if (col_name == "x") {
        layout.col_x = i;
      } else if (col_name == "y") {
        layout.col_y = i;
      } else if (col_name == "z") {
        layout.col_z = i;
      } else if (col_name == "xs" || col_name == "xsu") {
        layout.col_x = i;
        layout.scaled_coords = true;
      } else if (col_name == "ys" || col_name == "ysu") {
        layout.col_y = i;
        layout.scaled_coords = true;
      } else if (col_name == "zs" || col_name == "zsu") {
        layout.col_z = i;
        layout.scaled_coords = true;
      }
    }
    return layout;
  }

  void parseAtomLines(correlation::core::Cell &frame, int num_atoms, const ColumnLayout &layout) {
    if (layout.col_x == -1 || layout.col_y == -1 || layout.col_z == -1) {
      return; // Cannot parse positions — return empty frame
    }

    const auto &lattice_vector = frame.latticeVectors();
    std::string token;

    for (int i = 0; i < num_atoms; ++i) {
      if (offset >= size) {
        break;
      }
      std::string line = nextLine();

      std::istringstream atom_ss(line);
      std::vector<std::string> fields;
      while (atom_ss >> token) {
        fields.push_back(token);
      }

      const int num_fields = static_cast<int>(fields.size());
      if (layout.col_x >= num_fields || layout.col_y >= num_fields || layout.col_z >= num_fields) {
        continue;
      }

      double frac_x = std::stod(fields[layout.col_x]);
      double frac_y = std::stod(fields[layout.col_y]);
      double frac_z = std::stod(fields[layout.col_z]);

      // Determine element symbol.
      std::string element_symbol;
      if (layout.col_element >= 0 && layout.col_element < num_fields) {
        element_symbol = fields[layout.col_element];
      } else if (layout.col_type >= 0 && layout.col_type < num_fields) {
        element_symbol = fields[layout.col_type];
      } else {
        element_symbol =
            (layout.col_id >= 0 && layout.col_id < num_fields) ? fields[layout.col_id] : std::to_string(i + 1);
      }

      // Convert scaled (fractional) coordinates to Cartesian if needed.
      correlation::math::Vector3<double> pos;
      if (layout.scaled_coords) {
        pos = {frac_x * lattice_vector[0][0] + frac_y * lattice_vector[1][0] + frac_z * lattice_vector[2][0],
               frac_x * lattice_vector[0][1] + frac_y * lattice_vector[1][1] + frac_z * lattice_vector[2][1],
               frac_x * lattice_vector[0][2] + frac_y * lattice_vector[1][2] + frac_z * lattice_vector[2][2]};
      } else {
        pos = {frac_x, frac_y, frac_z};
      }

      frame.addAtom(element_symbol, pos);
    }
  }
};

} // namespace

// ---------------------------------------------------------------------------
// parseDumpFrame — parses a single frame from a memory region
// ---------------------------------------------------------------------------
correlation::core::Cell LammpsDumpReader::parseDumpFrame(const char *data, size_t size) {
  LammpsFrameParser parser{data, size, 0, 0};

  // --- ITEM: TIMESTEP ---
  parser.nextLine(); // "ITEM: TIMESTEP"
  parser.nextLine(); // timestep value (ignored)

  // --- NUMBER OF ATOMS ---
  int num_atoms = parser.parseAtomCount();

  // --- BOX BOUNDS & Build the Cell ---
  correlation::core::Cell frame = parser.parseBoxBounds();

  // --- ATOMS header — discover column layout ---
  auto layout = parser.parseAtomsHeader();

  // --- Atom lines ---
  parser.parseAtomLines(frame, num_atoms, layout);

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

  if (progress_callback) {
    progress_callback(0.0F, "Reading LAMMPS dump file...");
  }

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

  if (progress_callback) {
    progress_callback(1.0F, "LAMMPS dump file loaded.");
  }

  auto parser = [](const char *data, size_t size) { return parseDumpFrame(data, size); };

  return correlation::core::Trajectory(mapped_file, std::move(frame_offsets), parser, 1.0);
}

} // namespace correlation::readers
