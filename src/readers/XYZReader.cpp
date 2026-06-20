/**
 * @file XYZReader.cpp
 * @brief Implementation of the Extended XYZ file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/XYZReader.hpp"
#include "core/MappedFile.hpp"
#include "readers/ReaderFactory.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace correlation::readers {

namespace {
// Automatic registration
// NOLINTNEXTLINE(cert-err58-cpp, bugprone-throwing-static-initialization)
const bool registered = ReaderFactory::instance().registerReader(std::make_unique<XYZReader>());

struct XYZParser {
  const char *data = nullptr;
  size_t total_size = 0;
  size_t offset = 0;

  XYZParser(const char *data, size_t total_size)
      : data(data), total_size(total_size) {}

  void skipBlankLines() {
    while (offset < total_size) {
      size_t line_end = offset;
      while (line_end < total_size && data[line_end] != '\n' && data[line_end] != '\r') {
        line_end++;
      }

      bool is_blank = true;
      for (size_t char_idx = offset; char_idx < line_end; ++char_idx) {
        if (data[char_idx] != ' ' && data[char_idx] != '\t') {
          is_blank = false;
          break;
        }
      }

      if (!is_blank) {
        break;
      }

      offset = skipLineEnding(line_end);
    }
  }

  size_t skipLineEnding(size_t pos) const {
    if (pos < total_size && data[pos] == '\r') {
      pos++;
    }
    if (pos < total_size && data[pos] == '\n') {
      pos++;
    }
    return pos;
  }

  int parseAtomCount() {
    size_t line_end = offset;
    while (line_end < total_size && data[line_end] != '\n' && data[line_end] != '\r') {
      line_end++;
    }
    std::string const atom_count_str(data + offset, line_end - offset);
    offset = skipLineEnding(line_end);

    int num_atoms = 0;
    try {
      num_atoms = std::stoi(atom_count_str);
    } catch (...) {
      throw std::runtime_error("Invalid XYZ file: expected atom count, got: " + atom_count_str);
    }

    if (num_atoms <= 0) {
      throw std::runtime_error("Invalid XYZ file: non-positive atom count: " + atom_count_str);
    }

    return num_atoms;
  }

  void skipCommentLine() {
    if (offset >= total_size) {
      throw std::runtime_error("Invalid XYZ file: missing comment line after atom count");
    }
    size_t line_end = offset;
    while (line_end < total_size && data[line_end] != '\n' && data[line_end] != '\r') {
      line_end++;
    }
    offset = skipLineEnding(line_end);
  }

  void skipCoordinateLines(int num_atoms) {
    for (int atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
      if (offset >= total_size) {
        throw std::runtime_error("Invalid XYZ file: unexpected EOF while reading atom " + std::to_string(atom_idx + 1));
      }
      size_t line_end = offset;
      while (line_end < total_size && data[line_end] != '\n' && data[line_end] != '\r') {
        line_end++;
      }
      offset = skipLineEnding(line_end);
    }
  }

  std::vector<size_t> findFrameOffsets(const std::function<void(float, const std::string &)> &progress_callback) {
    std::vector<size_t> frame_offsets;

    while (offset < total_size) {
      skipBlankLines();
      if (offset >= total_size) {
        break;
      }

      size_t const frame_start = offset;
      int const num_atoms = parseAtomCount();
      skipCommentLine();
      skipCoordinateLines(num_atoms);

      frame_offsets.push_back(frame_start);

      // Report progress.
      if (progress_callback && total_size > 0) {
        float const progress = static_cast<float>(offset) / static_cast<float>(total_size);
        progress_callback(progress, "Reading XYZ frames...");
      }
    }

    return frame_offsets;
  }
};

} // namespace

// ---------------------------------------------------------------------------
// parseXYZFrame – parses a single frame from memory
// ---------------------------------------------------------------------------
correlation::core::Cell XYZReader::parseXYZFrame(const char *data, size_t size) {
  std::string const content(data, size);
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
    throw std::runtime_error("Invalid XYZ file: expected atom count, got: " + line);
  }

  if (num_atoms <= 0) {
    throw std::runtime_error("Invalid XYZ file: non-positive atom count: " + line);
  }

  // --- Line 2: comment / Extended XYZ header ---
  std::string comment;
  if (!std::getline(stream, comment)) {
    throw std::runtime_error("Invalid XYZ file: missing comment line after atom count");
  }

  correlation::core::Cell cell;

  // Try to parse Extended XYZ comment line
  auto comm_data = parseCommentLine(comment);

  if (comm_data.lattice) {
    const auto &lattice_vecs = *comm_data.lattice;
    correlation::math::Vector3<double> const param_a(lattice_vecs[0], lattice_vecs[1], lattice_vecs[2]);
    correlation::math::Vector3<double> const param_b(lattice_vecs[3], lattice_vecs[4], lattice_vecs[5]);
    correlation::math::Vector3<double> const param_c(lattice_vecs[6], lattice_vecs[7], lattice_vecs[8]);
    cell = correlation::core::Cell(param_a, param_b, param_c);
  }

  if (comm_data.energy) {
    cell.setEnergy(*comm_data.energy);
  }

  // --- Lines 3..N+2: atom data ---
  for (int i = 0; i < num_atoms; ++i) {
    if (!std::getline(stream, line)) {
      throw std::runtime_error("Invalid XYZ file: unexpected EOF while reading atom " + std::to_string(i + 1));
    }

    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::string token;
    while (iss >> token) {
      tokens.push_back(token);
    }

    int const max_idx =
        (std::max)({comm_data.species_col, comm_data.pos_x_col, comm_data.pos_y_col, comm_data.pos_z_col});

    if (std::cmp_less_equal(tokens.size(), max_idx)) {
      throw std::runtime_error("Invalid XYZ file: malformed atom line: " + line);
    }

    std::string const symbol = tokens[comm_data.species_col];
    try {
      double const pos_x = std::stod(tokens[comm_data.pos_x_col]);
      double const pos_y = std::stod(tokens[comm_data.pos_y_col]);
      double const pos_z = std::stod(tokens[comm_data.pos_z_col]);
      cell.addAtom(symbol, correlation::math::Vector3<double>(pos_x, pos_y, pos_z));
    } catch (...) {
      throw std::runtime_error("Invalid XYZ file: invalid coordinates: " + line);
    }
  }

  return cell;
}

void XYZReader::parseLattice(const std::string &comment, CommentData &data) {
  const std::string lat_key = "Lattice=\"";
  auto pos = comment.find(lat_key);
  if (pos != std::string::npos) {
    auto start = pos + lat_key.size();
    auto end = comment.find('"', start);
    if (end != std::string::npos) {
      std::string const values = comment.substr(start, end - start);
      std::istringstream iss(values);
      std::array<double, 9> lattice{};
      bool flag = true;
      for (int lat_idx = 0; lat_idx < 9; ++lat_idx) {
        if (!(iss >> lattice.at(lat_idx))) {
          flag = false;
          break;
        }
      }
      if (flag) {
        data.lattice = lattice;
      }
    }
  }
}

void XYZReader::parseEnergy(const std::string &comment, CommentData &data) {
  const std::vector<std::string> energy_keys = {"energy=", "dft_energy=", "Energy="};
  for (const auto &key : energy_keys) {
    auto pos = comment.find(key);
    if (pos != std::string::npos) {
      auto start = pos + key.size();
      if (start < comment.size()) {
        size_t end = std::string::npos;
        if (comment.at(start) == '"') {
          start++;
          end = comment.find('"', start);
        } else {
          end = comment.find_first_of(" \t\r\n", start);
        }
        std::string const val_str = comment.substr(start, end == std::string::npos ? std::string::npos : end - start);
        try {
          data.energy = std::stod(val_str);
          break;
        } catch (const std::exception &err) {
          std::cerr << "Warning: Failed to parse energy value '" << val_str << "' in comment line: " << err.what()
                    << '\n';
        }
      }
    }
  }
}

void XYZReader::parseProperties(const std::string &comment, CommentData &data) {
  const std::string prop_key = "Properties=";
  auto pos = comment.find(prop_key);
  if (pos == std::string::npos) {
    return;
  }
  auto start = pos + prop_key.size();
  if (start >= comment.size()) {
    return;
  }
  size_t end = 0;
  if (comment.at(start) == '"') {
    start++;
    end = comment.find('"', start);
  } else {
    end = comment.find_first_of(" \t\r\n", start);
  }

  std::string const props = comment.substr(start, end == std::string::npos ? std::string::npos : end - start);

  // format is name:type:cols:name:type:cols...
  std::vector<std::string> parts;
  std::istringstream p_ss(props);
  std::string p_part;
  while (std::getline(p_ss, p_part, ':')) {
    parts.push_back(p_part);
  }

  parsePropertiesParts(parts, data);
}

void XYZReader::parsePropertiesParts(const std::vector<std::string> &parts, CommentData &data) {
  int col_index = 0;
  data.species_col = -1;
  data.pos_x_col = -1;

  for (size_t part_idx = 0; part_idx + 2 < parts.size(); part_idx += 3) {
    const std::string &name = parts.at(part_idx);
    int cols = 1;
    try {
      cols = std::stoi(parts.at(part_idx + 2));
    } catch (const std::exception &err) {
      std::cerr << "Warning: Failed to parse column count '" << parts.at(part_idx + 2)
                << "' in Properties header: " << err.what() << '\n';
    }

    if (name == "species" || name == "type") {
      data.species_col = col_index;
    } else if (name == "pos") {
      data.pos_x_col = col_index;
      data.pos_y_col = col_index + 1;
      data.pos_z_col = col_index + 2;
    }
    col_index += cols;
  }

  if (data.species_col == -1) {
    data.species_col = 0;
  }
  if (data.pos_x_col == -1) {
    data.pos_x_col = 1;
    data.pos_y_col = 2;
    data.pos_z_col = 3;
  }
}

// ---------------------------------------------------------------------------
// Comment parser (Extended XYZ convention)
// ---------------------------------------------------------------------------
XYZReader::CommentData XYZReader::parseCommentLine(const std::string &comment) {
  CommentData data;
  parseLattice(comment, data);
  parseEnergy(comment, data);
  parseProperties(comment, data);
  return data;
}

// ---------------------------------------------------------------------------
// readStructure  – returns the last frame
// ---------------------------------------------------------------------------
correlation::core::Cell XYZReader::readStructure(const std::string &filename,
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
correlation::core::Trajectory
XYZReader::readTrajectory(const std::string &filename,
                          std::function<void(float, const std::string &)> progress_callback) {

  if (progress_callback) {
    progress_callback(0.0F, "Reading XYZ file...");
  }

  auto mapped_file = std::make_shared<correlation::core::MappedFile>(filename);
  const char *data = mapped_file->data();
  const size_t total_size = mapped_file->size();

  XYZParser parser(data, total_size);
  std::vector<size_t> frame_offsets = parser.findFrameOffsets(progress_callback);

  if (frame_offsets.empty()) {
    throw std::runtime_error("No frames found in XYZ file: " + filename);
  }

  // Add sentinel
  frame_offsets.push_back(parser.offset);

  if (progress_callback) {
    progress_callback(1.0F, "XYZ file loaded.");
  }

  auto parser_func = [](const char *data_begin, size_t data_size) { return parseXYZFrame(data_begin, data_size); };

  return {mapped_file, std::move(frame_offsets), parser_func, 1.0};
}

} // namespace correlation::readers
