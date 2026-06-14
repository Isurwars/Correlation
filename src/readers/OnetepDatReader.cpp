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
#include <memory>
#include <sstream>
#include <stdexcept>

namespace correlation::readers {

// Automatic registration
static bool registered = ReaderFactory::instance().registerReader(std::make_unique<OnetepDatReader>());

correlation::core::Cell
OnetepDatReader::readStructure(const std::string &filename,
                               std::function<void(float, const std::string &)> progress_callback) {
  return read(filename);
}

correlation::core::Trajectory
OnetepDatReader::readTrajectory(const std::string &filename,
                                std::function<void(float, const std::string &)> progress_callback) {
  throw std::runtime_error("ONETEP DAT files are structures, use readStructure.");
}

namespace {
void toLower(std::string &s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
}
} // namespace

correlation::core::Cell OnetepDatReader::read(const std::string &file_name) {
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" + std::strerror(errno) + ").");
  }

  correlation::core::Cell tempCell;
  bool in_block = false;
  bool frac_flag = false;
  std::string current_block_type;
  std::string line;
  int lattice_row_count = 0;
  std::array<double, 6> lat = {0, 0, 0, 0, 0, 0};
  double v[3][3] = {{0.0}};

  while (std::getline(myfile, line)) {
    // Strip comments
    size_t comment_pos = line.find_first_of("#!");
    if (comment_pos != std::string::npos) {
      line = line.substr(0, comment_pos);
    }
    // Trim
    line.erase(0, line.find_first_not_of(" \t\r\n"));
    line.erase(line.find_last_not_of(" \t\r\n") + 1);

    if (line.empty())
      continue;

    std::stringstream line_stream(line);
    std::string token;
    line_stream >> token;
    toLower(token);

    if (token == "%block") {
      in_block = true;
      line_stream >> current_block_type;
      toLower(current_block_type);
      lattice_row_count = 0;
      continue;
    }

    if (token == "%endblock") {
      in_block = false;
      current_block_type.clear();
      continue;
    }

    if (in_block) {
      if (current_block_type == "lattice_cart") {
        line_stream.clear();
        line_stream.seekg(0);
        std::string first_val;
        if (line_stream >> first_val) {
          std::string lower_first = first_val;
          toLower(lower_first);
          if (lower_first == "angstrom" || lower_first == "bohr") {
            continue;
          }
          if (lattice_row_count < 3) {
            std::stringstream ss(line);
            if (ss >> v[lattice_row_count][0] >> v[lattice_row_count][1] >> v[lattice_row_count][2]) {
              lattice_row_count++;
              if (lattice_row_count == 3) {
                tempCell = correlation::core::Cell({v[0][0], v[0][1], v[0][2]}, {v[1][0], v[1][1], v[1][2]},
                                                   {v[2][0], v[2][1], v[2][2]});
              }
            }
          }
        }
      }

      if (current_block_type == "lattice_abc") {
        line_stream.clear();
        line_stream.seekg(0);
        std::string first_val;
        if (line_stream >> first_val) {
          std::string lower_first = first_val;
          toLower(lower_first);
          if (lower_first == "angstrom" || lower_first == "bohr") {
            continue;
          }
          if (lattice_row_count < 2) {
            std::stringstream ss(line);
            if (ss >> lat[0 + 3 * lattice_row_count] >> lat[1 + 3 * lattice_row_count] >>
                lat[2 + 3 * lattice_row_count]) {
              lattice_row_count++;
              if (lattice_row_count == 2) {
                tempCell.setLatticeParameters(lat);
              }
            }
          }
        }
      }

      if (current_block_type == "positions_abs") {
        std::string element;
        double x, y, z;
        line_stream.clear();
        line_stream.seekg(0);
        std::string first_val;
        if (line_stream >> first_val) {
          std::string lower_first = first_val;
          toLower(lower_first);
          if (lower_first == "angstrom" || lower_first == "bohr") {
            continue;
          }
          std::stringstream ss(line);
          if (ss >> element >> x >> y >> z) {
            tempCell.addAtom(element, {x, y, z});
          }
        }
      }

      if (current_block_type == "positions_frac") {
        std::string element;
        double x, y, z;
        frac_flag = true;
        line_stream.clear();
        line_stream.seekg(0);
        std::istringstream ss(line);
        if (ss >> element >> x >> y >> z) {
          const auto &lv = tempCell.latticeVectors();
          correlation::math::Vector3<double> pos = {x * lv[0][0] + y * lv[1][0] + z * lv[2][0],
                                                    x * lv[0][1] + y * lv[1][1] + z * lv[2][1],
                                                    x * lv[0][2] + y * lv[1][2] + z * lv[2][2]};
          tempCell.addAtom(element, pos);
        }
      }
    }
  }

  if (frac_flag)
    tempCell.wrapPositions();
  return tempCell;
}

} // namespace correlation::readers
