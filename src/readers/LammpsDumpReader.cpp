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
#include "readers/NaturalSort.hpp"
#include "readers/ReaderFactory.hpp"

#include <charconv>
#include <cstring>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace correlation::readers {

namespace {

// Automatic registration
const bool REGISTERED = ReaderFactory::registerTypeSafe<LammpsDumpReader>("LammpsDumpReader");

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
inline std::string extractLine(const char *data, size_t pos, size_t line_end) {
  return {data + pos, line_end - pos};
}

// Struct to encapsulate parsing a single LAMMPS dump frame.
struct LammpsFrameParser {
  const char *data = nullptr;
  size_t size = 0;
  size_t offset = 0;
  size_t line_end = 0;

  std::string nextLine() {
    line_end = findLineEnd(data, size, offset);
    const std::string line = extractLine(data, offset, line_end);
    offset = skipLineEnding(data, size, line_end);
    return line;
  }

  int parseAtomCount() {
    nextLine(); // "ITEM: NUMBER OF ATOMS"
    const std::string line = nextLine();
    int num_atoms = 0;
    try {
      num_atoms = std::stoi(line);
    } catch (const std::exception &) {
      throw std::runtime_error("Failed to parse atom count in LAMMPS dump frame");
    }
    if (num_atoms <= 0) {
      throw std::runtime_error("Invalid LAMMPS dump frame: non-positive atom count: " +
                               std::to_string(num_atoms));
    }
    return num_atoms;
  }

  correlation::core::Cell parseBoxBounds() {
    const std::string line = nextLine(); // "ITEM: BOX BOUNDS ..."
    const bool triclinic = line.contains("xy");

    real_t xlo = 0.0;
    real_t xhi = 0.0;
    real_t ylo = 0.0;
    real_t yhi = 0.0;
    real_t zlo = 0.0;
    real_t zhi = 0.0;
    real_t x_y = 0.0;
    real_t x_z = 0.0;
    real_t y_z = 0.0;

    if (triclinic) {
      std::stringstream(nextLine()) >> xlo >> xhi >> x_y;
      std::stringstream(nextLine()) >> ylo >> yhi >> x_z;
      std::stringstream(nextLine()) >> zlo >> zhi >> y_z;
      const real_t l_x = xhi - xlo;
      const real_t l_y = yhi - ylo;
      const real_t l_z = zhi - zlo;
      return correlation::core::Cell({l_x, 0.0, 0.0}, {x_y, l_y, 0.0}, {x_z, y_z, l_z});
    }
    std::stringstream(nextLine()) >> xlo >> xhi;
    std::stringstream(nextLine()) >> ylo >> yhi;
    std::stringstream(nextLine()) >> zlo >> zhi;
    return correlation::core::Cell({xhi - xlo, 0.0, 0.0}, {0.0, yhi - ylo, 0.0},
                                   {0.0, 0.0, zhi - zlo});
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
    const std::string line = nextLine(); // "ITEM: ATOMS ..."

    std::istringstream header_ss(line);
    std::string token;
    header_ss >> token >> token; // consume "ITEM:" and "ATOMS"

    std::vector<std::string> col_names;
    while (header_ss >> token) {
      col_names.push_back(token);
    }

    ColumnLayout layout;
    for (size_t i = 0; i < col_names.size(); ++i) {
      const std::string &col_name = col_names[i];
      if (col_name == "id") {
        layout.col_id = static_cast<int>(i);
      } else if (col_name == "type") {
        layout.col_type = static_cast<int>(i);
      } else if (col_name == "element") {
        layout.col_element = static_cast<int>(i);
      } else if (col_name == "x") {
        layout.col_x = static_cast<int>(i);
      } else if (col_name == "y") {
        layout.col_y = static_cast<int>(i);
      } else if (col_name == "z") {
        layout.col_z = static_cast<int>(i);
      } else if (col_name == "xs" || col_name == "xsu") {
        layout.col_x = static_cast<int>(i);
        layout.scaled_coords = true;
      } else if (col_name == "ys" || col_name == "ysu") {
        layout.col_y = static_cast<int>(i);
        layout.scaled_coords = true;
      } else if (col_name == "zs" || col_name == "zsu") {
        layout.col_z = static_cast<int>(i);
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
      const std::string line = nextLine();

      std::istringstream atom_ss(line);
      std::vector<std::string> fields;
      while (atom_ss >> token) {
        fields.push_back(token);
      }

      const int num_fields = static_cast<int>(fields.size());
      if (layout.col_x >= num_fields || layout.col_y >= num_fields || layout.col_z >= num_fields) {
        continue;
      }

      real_t frac_x = 0.0;
      real_t frac_y = 0.0;
      real_t frac_z = 0.0;
      try {
        frac_x = static_cast<real_t>(std::stod(fields[layout.col_x]));
        frac_y = static_cast<real_t>(std::stod(fields[layout.col_y]));
        frac_z = static_cast<real_t>(std::stod(fields[layout.col_z]));
      } catch (const std::exception &) {
        continue;
      }

      // Determine element symbol.
      std::string element_symbol;
      if (layout.col_element >= 0 && layout.col_element < num_fields) {
        element_symbol = fields[layout.col_element];
      } else if (layout.col_type >= 0 && layout.col_type < num_fields) {
        element_symbol = fields[layout.col_type];
      } else {
        element_symbol = (layout.col_id >= 0 && layout.col_id < num_fields) ? fields[layout.col_id]
                                                                            : std::to_string(i + 1);
      }

      // Convert scaled (fractional) coordinates to Cartesian if needed.
      correlation::math::Vector3<real_t> pos;
      if (layout.scaled_coords) {
        pos = lattice_vector * correlation::math::Vector3<real_t>{frac_x, frac_y, frac_z};
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
  LammpsFrameParser parser{.data = data, .size = size, .offset = 0, .line_end = 0};

  // --- ITEM: TIMESTEP ---
  parser.nextLine(); // "ITEM: TIMESTEP"
  parser.nextLine(); // timestep value (ignored)

  // --- NUMBER OF ATOMS ---
  const int num_atoms = parser.parseAtomCount();

  // --- BOX BOUNDS & Build the Cell ---
  correlation::core::Cell frame = parser.parseBoxBounds();

  // --- ATOMS header — discover column layout ---
  const auto layout = parser.parseAtomsHeader();

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
  const auto traj = readTrajectory(filename, progress_callback);
  if (traj.getFrameCount() == 0) {
    throw std::runtime_error("No frames found in LAMMPS dump file: " + filename);
  }
  return traj.getFrame(0);
}

namespace {

// ---------------------------------------------------------------------------
// Helper: extract integer timestep value from frame starting at frame_start
// ---------------------------------------------------------------------------
inline std::optional<int64_t> parseTimestepValue(const char *data, size_t total_size,
                                                 size_t frame_start) {
  const size_t line1_end = findLineEnd(data, total_size, frame_start);
  const size_t val_start = skipLineEnding(data, total_size, line1_end);
  if (val_start >= total_size) {
    return std::nullopt;
  }
  const size_t val_end = findLineEnd(data, total_size, val_start);
  std::string_view val_str(data + val_start, val_end - val_start);
  while (!val_str.empty() && std::isspace(static_cast<unsigned char>(val_str.front())) != 0) {
    val_str.remove_prefix(1);
  }
  while (!val_str.empty() && std::isspace(static_cast<unsigned char>(val_str.back())) != 0) {
    val_str.remove_suffix(1);
  }
  int64_t step = 0;
  const auto [ptr, ec] = std::from_chars(val_str.data(), val_str.data() + val_str.size(), step);
  if (ec == std::errc()) {
    return step;
  }
  return std::nullopt;
}

// ---------------------------------------------------------------------------
// Helper: find byte offsets of all frame starts in a memory-mapped file
// ---------------------------------------------------------------------------
inline std::vector<size_t> findFrameStarts(const char *data, size_t total_size) {
  std::vector<size_t> frame_offsets;
  const char *needle = "ITEM: TIMESTEP";
  const size_t needle_len = std::strlen(needle);

  size_t pos = 0;
  while (pos < total_size) {
    const bool at_line_start = (pos == 0) || (pos > 0 && (data[pos - 1] == '\n'));
    if (at_line_start && pos + needle_len <= total_size &&
        std::memcmp(data + pos, needle, needle_len) == 0) {
      frame_offsets.push_back(pos);
    }
    pos = findLineEnd(data, total_size, pos);
    pos = skipLineEnding(data, total_size, pos);
  }
  return frame_offsets;
}

// ---------------------------------------------------------------------------
// Helper: process frame slices in a single mapped file, checking timesteps
// ---------------------------------------------------------------------------
inline void processFileFrames(const std::shared_ptr<correlation::core::MappedFile> &mapped,
                              uint32_t file_idx, const std::string &filename,
                              std::vector<correlation::core::Trajectory::MappedFrameSlice> &slices,
                              std::optional<int64_t> &last_timestep) {
  const char *data = mapped->data();
  const size_t total_size = mapped->size();
  std::vector<size_t> frame_offsets = findFrameStarts(data, total_size);

  if (frame_offsets.empty()) {
    throw std::runtime_error("No frames found in LAMMPS dump file: " + filename);
  }

  frame_offsets.push_back(total_size);

  for (size_t i = 0; i + 1 < frame_offsets.size(); ++i) {
    const size_t frame_start = frame_offsets[i];
    const size_t frame_end = frame_offsets[i + 1];
    const auto step_opt = parseTimestepValue(data, total_size, frame_start);

    if (step_opt.has_value()) {
      const int64_t step = *step_opt;
      if (last_timestep.has_value()) {
        if (step < *last_timestep) {
          throw std::runtime_error("LAMMPS frames out of chronological order: timestep " +
                                   std::to_string(step) + " preceded by " +
                                   std::to_string(*last_timestep) + " in file: " + filename);
        }
        if (step == *last_timestep) {
          // Deduplicate consecutive identical timesteps (restart checkpoint boundary)
          continue;
        }
      }
      last_timestep = step;
    }

    slices.push_back({
        .file_index = file_idx,
        .offset = frame_start,
        .size = frame_end - frame_start,
    });
  }
}

} // namespace

// ---------------------------------------------------------------------------
// readTrajectory — multi-file memory-mapped lazy loading
// ---------------------------------------------------------------------------
correlation::core::Trajectory LammpsDumpReader::readTrajectory(
    const std::vector<std::string> &filenames,
    std::function<void(float, const std::string &)> progress_callback) {

  if (filenames.empty()) {
    throw std::runtime_error("No files provided for LAMMPS dump trajectory reading.");
  }

  std::vector<std::string> sorted_files = filenames;
  naturalSort(sorted_files);

  if (progress_callback) {
    progress_callback(0.0F, "Opening LAMMPS dump files...");
  }

  std::vector<std::shared_ptr<correlation::core::MappedFile>> mapped_files;
  mapped_files.reserve(sorted_files.size());
  std::vector<correlation::core::Trajectory::MappedFrameSlice> frame_slices;
  std::optional<int64_t> last_timestep;

  for (size_t f_idx = 0; f_idx < sorted_files.size(); ++f_idx) {
    const auto &file_path = sorted_files[f_idx];
    const auto mapped_file = std::make_shared<correlation::core::MappedFile>(file_path);
    mapped_files.push_back(mapped_file);

    processFileFrames(mapped_file, static_cast<uint32_t>(f_idx), file_path, frame_slices,
                      last_timestep);

    if (progress_callback && !sorted_files.empty()) {
      const float progress =
          static_cast<float>(f_idx + 1) / static_cast<float>(sorted_files.size());
      progress_callback(progress, "Processed " + file_path);
    }
  }

  if (frame_slices.empty()) {
    throw std::runtime_error("No frames found across LAMMPS dump files.");
  }

  if (progress_callback) {
    progress_callback(1.0F, "LAMMPS dump files loaded.");
  }

  const auto parser = [](const char *buf, size_t size) { return parseDumpFrame(buf, size); };
  return {std::move(mapped_files), std::move(frame_slices), parser, 1.0};
}

// ---------------------------------------------------------------------------
// readTrajectory — single file or wildcard pattern
// ---------------------------------------------------------------------------
correlation::core::Trajectory LammpsDumpReader::readTrajectory(
    const std::string &filename,
    std::function<void(float, const std::string &)> progress_callback) {
  if (filename.contains('*') || filename.contains('?')) {
    const auto files = expandWildcards(filename);
    return readTrajectory(files, std::move(progress_callback));
  }
  return readTrajectory(std::vector<std::string>{filename}, std::move(progress_callback));
}

} // namespace correlation::readers
