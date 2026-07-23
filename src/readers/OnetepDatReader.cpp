/**
 * @file OnetepDatReader.cpp
 * @brief Implementation of the ONETEP .dat file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */
#include "readers/OnetepDatReader.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/ReaderFactory.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <functional>
#include <sstream>
#include <stdexcept>

namespace correlation::readers {

namespace {

// Automatic registration
const bool registered = ReaderFactory::registerTypeSafe<OnetepDatReader>("OnetepDatReader");

void toLower(std::string &str) {
  std::transform(str.begin(), str.end(), str.begin(), [](unsigned char chr) { return std::tolower(chr); });
}

struct OnetepDatParser {
  correlation::core::Cell tempCell;
  bool in_block = false;
  bool frac_flag = false;
  std::string current_block_type;
  int lattice_row_count = 0;
  std::array<real_t, 6> lat = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  std::array<std::array<real_t, 3>, 3> v = {};

  static void cleanLine(std::string &line) {
    // Strip comments
    size_t comment_pos = line.find_first_of("#!");
    if (comment_pos != std::string::npos) {
      line = line.substr(0, comment_pos);
    }
    // Trim
    line.erase(0, line.find_first_not_of(" \t\r\n"));
    line.erase(line.find_last_not_of(" \t\r\n") + 1);
  }

  void handleBlockStart(std::stringstream &line_stream) {
    in_block = true;
    line_stream >> current_block_type;
    toLower(current_block_type);
    lattice_row_count = 0;
  }

  void handleBlockEnd() {
    in_block = false;
    current_block_type.clear();
  }

  void parseLatticeCart(const std::string &line, std::stringstream &line_stream) {
    line_stream.clear();
    line_stream.seekg(0);
    std::string first_val;
    if (!(line_stream >> first_val)) {
      return;
    }
    std::string lower_first = first_val;
    toLower(lower_first);
    if (lower_first == "angstrom" || lower_first == "bohr") {
      return;
    }
    if (lattice_row_count < 3) {
      std::stringstream str_stream(line);
      if (str_stream >> v.at(lattice_row_count)[0] >> v.at(lattice_row_count)[1] >> v.at(lattice_row_count)[2]) {
        lattice_row_count++;
        if (lattice_row_count == 3) {
          tempCell = correlation::core::Cell({v[0][0], v[0][1], v[0][2]}, {v[1][0], v[1][1], v[1][2]},
                                             {v[2][0], v[2][1], v[2][2]});
        }
      }
    }
  }

  void parseLatticeAbc(const std::string &line, std::stringstream &line_stream) {
    line_stream.clear();
    line_stream.seekg(0);
    std::string first_val;
    if (!(line_stream >> first_val)) {
      return;
    }
    std::string lower_first = first_val;
    toLower(lower_first);
    if (lower_first == "angstrom" || lower_first == "bohr") {
      return;
    }
    if (lattice_row_count < 2) {
      std::stringstream str_stream(line);
      if (str_stream >> lat.at(0 + 3 * lattice_row_count) >> lat.at(1 + 3 * lattice_row_count) >>
          lat.at(2 + 3 * lattice_row_count)) {
        lattice_row_count++;
        if (lattice_row_count == 2) {
          tempCell.setLatticeParameters(lat);
        }
      }
    }
  }

  void parsePositionsAbs(const std::string &line, std::stringstream &line_stream) {
    std::string element;
    real_t pos_x = 0.0;
    real_t pos_y = 0.0;
    real_t pos_z = 0.0;
    line_stream.clear();
    line_stream.seekg(0);
    std::string first_val;
    if (!(line_stream >> first_val)) {
      return;
    }
    std::string lower_first = first_val;
    toLower(lower_first);
    if (lower_first == "angstrom" || lower_first == "bohr") {
      return;
    }
    std::stringstream str_stream(line);
    if (str_stream >> element >> pos_x >> pos_y >> pos_z) {
      tempCell.addAtom(element, {pos_x, pos_y, pos_z});
    }
  }

  void parsePositionsFrac(const std::string &line) {
    std::string element;
    real_t frac_x = 0.0;
    real_t frac_y = 0.0;
    real_t frac_z = 0.0;
    frac_flag = true;
    std::istringstream str_stream(line);
    if (str_stream >> element >> frac_x >> frac_y >> frac_z) {
      const auto &lattice_vectors = tempCell.latticeVectors();
      correlation::math::Vector3<real_t> pos = {
          frac_x * lattice_vectors[0][0] + frac_y * lattice_vectors[1][0] + frac_z * lattice_vectors[2][0],
          frac_x * lattice_vectors[0][1] + frac_y * lattice_vectors[1][1] + frac_z * lattice_vectors[2][1],
          frac_x * lattice_vectors[0][2] + frac_y * lattice_vectors[1][2] + frac_z * lattice_vectors[2][2]};
      tempCell.addAtom(element, pos);
    }
  }

  void processBlockLine(const std::string &line, std::stringstream &line_stream) {
    if (current_block_type == "lattice_cart") {
      parseLatticeCart(line, line_stream);
    } else if (current_block_type == "lattice_abc") {
      parseLatticeAbc(line, line_stream);
    } else if (current_block_type == "positions_abs") {
      parsePositionsAbs(line, line_stream);
    } else if (current_block_type == "positions_frac") {
      parsePositionsFrac(line);
    }
  }

  void parse(std::ifstream &myfile) {
    std::string line;
    while (std::getline(myfile, line)) {
      cleanLine(line);
      if (line.empty()) {
        continue;
      }

      std::stringstream line_stream(line);
      std::string token;
      line_stream >> token;
      toLower(token);

      if (token == "%block") {
        handleBlockStart(line_stream);
        continue;
      }

      if (token == "%endblock") {
        handleBlockEnd();
        continue;
      }

      if (in_block) {
        processBlockLine(line, line_stream);
      }
    }
  }
};

} // namespace

correlation::core::Cell
OnetepDatReader::readStructure(const std::string &filename,
                               std::function<void(float, const std::string &)> /*progress_callback*/) {
  return read(filename);
}

correlation::core::Trajectory
OnetepDatReader::readTrajectory(const std::string & /*filename*/,
                                std::function<void(float, const std::string &)> /*progress_callback*/) {
  throw std::runtime_error("ONETEP DAT files are structures, use readStructure.");
}

correlation::core::Cell OnetepDatReader::read(const std::string &file_name) {
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" + std::strerror(errno) + ").");
  }

  OnetepDatParser parser;
  parser.parse(myfile);

  if (parser.frac_flag) {
    parser.tempCell.wrapPositions();
  }
  return parser.tempCell;
}

} // namespace correlation::readers
