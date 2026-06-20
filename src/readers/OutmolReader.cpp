/**
 * @file OutmolReader.cpp
 * @brief Implementation of the CASTEP .outmol file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "readers/OutmolReader.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "math/Constants.hpp"
#include "readers/ReaderFactory.hpp"

#include <cerrno>
#include <cstring>
#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace correlation::readers {

namespace {

// Automatic registration
// NOLINTNEXTLINE(cert-err58-cpp, bugprone-throwing-static-initialization)
const bool registered = ReaderFactory::instance().registerReader(std::make_unique<OutmolReader>());

struct OutmolParser {
  std::ifstream *myfile = nullptr;
  std::function<void(float, const std::string &)> progress_callback;
  std::vector<correlation::core::Cell> frames;

  std::array<double, 3> h_one = {0.0, 0.0, 0.0};
  std::array<double, 3> h_two = {0.0, 0.0, 0.0};
  std::array<double, 3> h_three = {0.0, 0.0, 0.0};
  bool cell_parsed = false;

  std::streampos file_size = 0;
  std::streampos last_progress_pos = 0;
  size_t update_interval = 0;
  std::stringstream str_stream;

  OutmolParser(std::ifstream &file, std::function<void(float, const std::string &)> progress)
      : myfile(&file), progress_callback(std::move(progress)), update_interval(file_size / 100) {
    myfile->seekg(0, std::ios::end);
    file_size = myfile->tellg();
    myfile->seekg(0, std::ios::beg);
  }

  void reportProgress() {
    if (progress_callback) {
      std::streampos const current_pos = myfile->tellg();
      if (std::cmp_greater(current_pos - last_progress_pos, update_interval)) {
        float const pos_fraction = static_cast<float>(current_pos) / static_cast<float>(file_size);
        progress_callback(pos_fraction, "Loading Outmol file...");
        last_progress_pos = current_pos;
      }
    }
  }

  bool parseCellVectors() {
    std::string line;
    if (!std::getline(*myfile, line)) {
      return false; // row 1
    }
    str_stream.clear();
    str_stream.str(line);
    str_stream >> h_one[0] >> h_one[1] >> h_one[2];

    if (!std::getline(*myfile, line)) {
      return false; // row 2
    }
    str_stream.clear();
    str_stream.str(line);
    str_stream >> h_two[0] >> h_two[1] >> h_two[2];

    if (!std::getline(*myfile, line)) {
      return false; // row 3
    }
    str_stream.clear();
    str_stream.str(line);
    str_stream >> h_three[0] >> h_three[1] >> h_three[2];

    for (int i = 0; i < 3; ++i) {
      h_one.at(i) *= correlation::math::bohr_to_angstrom;
      h_two.at(i) *= correlation::math::bohr_to_angstrom;
      h_three.at(i) *= correlation::math::bohr_to_angstrom;
    }
    cell_parsed = true;
    return true;
  }

  void parseCoordinates() {
    correlation::core::Cell tempCell({h_one[0], h_one[1], h_one[2]}, {h_two[0], h_two[1], h_two[2]},
                                     {h_three[0], h_three[1], h_three[2]});
    std::string line;
    while (std::getline(*myfile, line)) {
      if (line.contains("$end")) {
        break;
      }
      str_stream.clear();
      str_stream.str(line);
      std::string symbol;
      double pos_x = 0.0;
      double pos_y = 0.0;
      double pos_z = 0.0;
      if (str_stream >> symbol >> pos_x >> pos_y >> pos_z) {
        tempCell.addAtom(symbol,
                         {pos_x * correlation::math::bohr_to_angstrom, pos_y * correlation::math::bohr_to_angstrom,
                          pos_z * correlation::math::bohr_to_angstrom});
      }
    }
    if (!tempCell.isEmpty()) {
      frames.push_back(std::move(tempCell));
    }
  }

  void parseAtomicCoordinates() {
    std::string line;
    std::getline(*myfile, line); // Next line is header: df x y z ...
    correlation::core::Cell tempCell({h_one[0], h_one[1], h_one[2]}, {h_two[0], h_two[1], h_two[2]},
                                     {h_three[0], h_three[1], h_three[2]});
    while (std::getline(*myfile, line)) {
      // Look for the "df" prefix
      if (!line.contains("df")) {
        break;
      }
      str_stream.clear();
      str_stream.str(line);
      std::string data_block;
      std::string symbol;
      double pos_x = 0.0;
      double pos_y = 0.0;
      double pos_z = 0.0;
      if (str_stream >> data_block >> symbol >> pos_x >> pos_y >> pos_z) {
        if (data_block == "df") {
          tempCell.addAtom(symbol,
                           {pos_x * correlation::math::bohr_to_angstrom, pos_y * correlation::math::bohr_to_angstrom,
                            pos_z * correlation::math::bohr_to_angstrom});
        }
      }
    }
    if (!tempCell.isEmpty()) {
      frames.push_back(std::move(tempCell));
    }
  }

  void parse() {
    std::string line;
    while (std::getline(*myfile, line)) {
      reportProgress();

      if (line.contains("$cell vectors")) {
        if (!parseCellVectors()) {
          break;
        }
        continue;
      }

      if (line.contains("$coordinates")) {
        parseCoordinates();
        continue;
      }

      if (line.contains("ATOMIC  COORDINATES (au)")) {
        parseAtomicCoordinates();
        continue;
      }
    }
  }
};

} // namespace

correlation::core::Cell OutmolReader::readStructure(const std::string &filename,
                                                    std::function<void(float, const std::string &)> progress_callback) {
  return read(filename, progress_callback).at(0);
}

correlation::core::Trajectory
OutmolReader::readTrajectory(const std::string &filename,
                             std::function<void(float, const std::string &)> progress_callback) {
  return {read(filename, progress_callback), 1.0};
}

std::vector<correlation::core::Cell>
OutmolReader::read(const std::string &file_name, std::function<void(float, const std::string &)> progress_callback) {
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" + std::strerror(errno) + ").");
  }

  OutmolParser parser(myfile, std::move(progress_callback));
  parser.parse();

  return parser.frames;
}

} // namespace correlation::readers
