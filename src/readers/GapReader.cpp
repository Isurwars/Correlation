/**
 * @file GapReader.cpp
 * @brief Implementation of GAP (Gaussian Approximation Potentials / QUIP) Extended XYZ trajectory reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */
#include "readers/GapReader.hpp"
#include "core/MappedFile.hpp"
#include "readers/ReaderFactory.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace correlation::readers {

namespace {
// Automatic registration
const bool registered = ReaderFactory::registerTypeSafe<GapReader>("GAP Reader");

struct GapParser {
  const char *data = nullptr;
  size_t total_size = 0;
  size_t offset = 0;

  GapParser(const char *data, size_t total_size) : data(data), total_size(total_size) {}

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

  [[nodiscard]] size_t skipLineEnding(size_t pos) const {
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
      throw std::runtime_error("Invalid GAP file: malformed atom count at offset " + std::to_string(offset));
    }
    if (num_atoms <= 0) {
      throw std::runtime_error("Invalid GAP file: non-positive atom count at offset " + std::to_string(offset));
    }
    return num_atoms;
  }

  void skipCommentLine() {
    size_t line_end = offset;
    while (line_end < total_size && data[line_end] != '\n' && data[line_end] != '\r') {
      line_end++;
    }
    offset = skipLineEnding(line_end);
  }

  void skipCoordinateLines(int num_atoms) {
    for (int atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
      if (offset >= total_size) {
        throw std::runtime_error("Invalid GAP file: unexpected EOF while reading atom lines");
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

    while (true) {
      skipBlankLines();
      if (offset >= total_size) {
        break;
      }

      size_t const frame_start = offset;
      int const num_atoms = parseAtomCount();
      skipCommentLine();
      skipCoordinateLines(num_atoms);

      frame_offsets.push_back(frame_start);

      if (progress_callback && total_size > 0 && (frame_offsets.size() % 100 == 0)) {
        float const progress = static_cast<float>(offset) / static_cast<float>(total_size);
        progress_callback(progress, "Indexing GAP frames...");
      }
    }

    return frame_offsets;
  }
};

void parseAtomLine(const std::string &line, const GapReader::CommentData &comm_data, correlation::core::Cell &cell) {
  std::istringstream iss(line);
  std::vector<std::string> tokens;
  std::string token;
  while (iss >> token) {
    tokens.push_back(token);
  }

  int const max_idx =
      (std::max)({comm_data.species_col, comm_data.pos_x_col, comm_data.pos_y_col, comm_data.pos_z_col});

  if (std::cmp_less_equal(tokens.size(), max_idx)) {
    throw std::runtime_error("Invalid GAP file: malformed atom line: " + line);
  }

  std::string const symbol = tokens[comm_data.species_col];
  try {
    const auto pos_x = static_cast<real_t>(std::stod(tokens[comm_data.pos_x_col]));
    const auto pos_y = static_cast<real_t>(std::stod(tokens[comm_data.pos_y_col]));
    const auto pos_z = static_cast<real_t>(std::stod(tokens[comm_data.pos_z_col]));
    auto &atom = cell.addAtom(symbol, correlation::math::Vector3<real_t>(pos_x, pos_y, pos_z));

    if (comm_data.force_x_col >= 0 && comm_data.force_y_col >= 0 && comm_data.force_z_col >= 0) {
      int const max_force_idx = (std::max)({comm_data.force_x_col, comm_data.force_y_col, comm_data.force_z_col});
      if (std::cmp_greater(tokens.size(), max_force_idx)) {
        const auto force_x = static_cast<real_t>(std::stod(tokens[comm_data.force_x_col]));
        const auto force_y = static_cast<real_t>(std::stod(tokens[comm_data.force_y_col]));
        const auto force_z = static_cast<real_t>(std::stod(tokens[comm_data.force_z_col]));
        atom.setVelocity(correlation::math::Vector3<real_t>(force_x, force_y, force_z));
      }
    }
  } catch (const std::exception &err) {
    throw std::runtime_error("Invalid GAP file: failed to parse coordinates/forces from atom line: " + line +
                             " (" + err.what() + ")");
  }
}

} // namespace

correlation::core::Cell GapReader::parseGapFrame(const char *data, size_t size) {
  std::string const frame_str(data, size);
  std::istringstream stream(frame_str);
  std::string line;

  while (std::getline(stream, line)) {
    if (!line.empty() && line.find_first_not_of(" \t\r\n") != std::string::npos) {
      break;
    }
  }

  if (line.empty()) {
    throw std::runtime_error("Invalid GAP file: empty frame or missing atom count");
  }

  int num_atoms = 0;
  try {
    num_atoms = std::stoi(line);
  } catch (...) {
    throw std::runtime_error("Invalid GAP file: malformed atom count: " + line);
  }

  if (num_atoms <= 0) {
    throw std::runtime_error("Invalid GAP file: non-positive atom count: " + line);
  }

  std::string comment;
  if (!std::getline(stream, comment)) {
    throw std::runtime_error("Invalid GAP file: missing comment line after atom count");
  }

  correlation::core::Cell cell;
  auto comm_data = parseCommentLine(comment);

  if (comm_data.lattice) {
    const auto &lat = *comm_data.lattice;
    correlation::math::Vector3<real_t> const param_a(lat[0], lat[1], lat[2]);
    correlation::math::Vector3<real_t> const param_b(lat[3], lat[4], lat[5]);
    correlation::math::Vector3<real_t> const param_c(lat[6], lat[7], lat[8]);
    cell = correlation::core::Cell(param_a, param_b, param_c);
  }

  if (comm_data.energy) {
    cell.setEnergy(*comm_data.energy);
  }

  cell.reserveAtoms(static_cast<size_t>(num_atoms));

  for (int i = 0; i < num_atoms; ++i) {
    if (!std::getline(stream, line)) {
      throw std::runtime_error("Invalid GAP file: unexpected EOF while reading atom " + std::to_string(i + 1));
    }
    parseAtomLine(line, comm_data, cell);
  }

  return cell;
}

void GapReader::parseLattice(const std::string &comment, CommentData &data) {
  const std::string lat_key = "Lattice=\"";
  auto pos = comment.find(lat_key);
  if (pos == std::string::npos) {
    return;
  }
  auto start = pos + lat_key.size();
  auto end = comment.find('"', start);
  if (end == std::string::npos) {
    return;
  }
  std::string const values = comment.substr(start, end - start);
  std::istringstream iss(values);
  std::array<real_t, 9> lattice{};
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

void GapReader::parseEnergy(const std::string &comment, CommentData &data) {
  const std::vector<std::string> energy_keys = {
      "gap_energy=", "GAP_energy=", "energy=", "free_energy=", "dft_energy=", "Energy="};

  for (const auto &key : energy_keys) {
    auto pos = comment.find(key);
    if (pos == std::string::npos) {
      continue;
    }
    auto start = pos + key.size();
    if (start >= comment.size()) {
      continue;
    }
    size_t end = std::string::npos;
    if (comment.at(start) == '"') {
      start++;
      end = comment.find('"', start);
    } else {
      end = comment.find_first_of(" \t\r\n", start);
    }
    std::string const val_str = comment.substr(start, end == std::string::npos ? std::string::npos : end - start);
    try {
      data.energy = static_cast<real_t>(std::stod(val_str));
      break;
    } catch (const std::exception &err) {
      std::cerr << "Warning: Failed to parse GAP energy value '" << val_str << "': " << err.what() << '\n';
    }
  }
}

void GapReader::parseProperties(const std::string &comment, CommentData &data) {
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

  std::vector<std::string> parts;
  std::istringstream p_ss(props);
  std::string p_part;
  while (std::getline(p_ss, p_part, ':')) {
    parts.push_back(p_part);
  }

  parsePropertiesParts(parts, data);
}

void GapReader::parsePropertiesParts(const std::vector<std::string> &parts, CommentData &data) {
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
                << "' in GAP Properties header: " << err.what() << '\n';
    }

    if (name == "species" || name == "type") {
      data.species_col = col_index;
    } else if (name == "pos") {
      data.pos_x_col = col_index;
      data.pos_y_col = col_index + 1;
      data.pos_z_col = col_index + 2;
    } else if (name == "gap_forces" || name == "GAP_forces" || name == "gap_force" || name == "GAP_force" ||
               name == "force" || name == "forces") {
      data.force_x_col = col_index;
      data.force_y_col = col_index + 1;
      data.force_z_col = col_index + 2;
    }
    col_index += cols;
  }

  if (data.species_col == -1 || data.pos_x_col == -1) {
    data.species_col = 0;
    data.pos_x_col = 1;
    data.pos_y_col = 2;
    data.pos_z_col = 3;
    data.force_x_col = 4;
    data.force_y_col = 5;
    data.force_z_col = 6;
  }
}

GapReader::CommentData GapReader::parseCommentLine(const std::string &comment) {
  CommentData data;
  parseLattice(comment, data);
  parseEnergy(comment, data);
  parseProperties(comment, data);
  return data;
}

correlation::core::Cell GapReader::readStructure(const std::string &filename,
                                                 std::function<void(float, const std::string &)> progress_callback) {
  auto traj = readTrajectory(filename, std::move(progress_callback));
  if (traj.getFrameCount() == 0) {
    throw std::runtime_error("No structure found in GAP file: " + filename);
  }
  return traj.getFrame(traj.getFrameCount() - 1);
}

correlation::core::Trajectory
GapReader::readTrajectory(const std::string &filename,
                          std::function<void(float, const std::string &)> progress_callback) {
  if (progress_callback) {
    progress_callback(0.0F, "Reading GAP file...");
  }

  auto mapped_file = std::make_shared<correlation::core::MappedFile>(filename);
  const char *data = mapped_file->data();
  const size_t total_size = mapped_file->size();

  GapParser parser(data, total_size);
  std::vector<size_t> frame_offsets = parser.findFrameOffsets(progress_callback);

  if (frame_offsets.empty()) {
    throw std::runtime_error("No frames found in GAP file: " + filename);
  }

  frame_offsets.push_back(parser.offset);

  if (progress_callback) {
    progress_callback(1.0F, "GAP file loaded.");
  }

  auto parser_func = [](const char *data_begin, size_t data_size) { return parseGapFrame(data_begin, data_size); };

  return {mapped_file, std::move(frame_offsets), parser_func, 1.0};
}

} // namespace correlation::readers
