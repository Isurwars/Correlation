/**
 * @file CarReader.cpp
 * @brief Implementation of the CAR file format reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */
#include "readers/CarReader.hpp"
#include "core/Cell.hpp"
#include "core/Trajectory.hpp"
#include "readers/ReaderFactory.hpp"

#include <array>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace correlation::readers {

// Automatic registration
// NOLINTNEXTLINE(cert-err58-cpp, bugprone-throwing-static-initialization)
static const bool registered = ReaderFactory::instance().registerReader(std::make_unique<CarReader>());

correlation::core::Cell
CarReader::readStructure(const std::string &filename,
                         std::function<void(float, const std::string &)> /*progress_callback*/) {
  return read(filename);
}

correlation::core::Trajectory
CarReader::readTrajectory(const std::string & /*filename*/,
                          std::function<void(float, const std::string &)> /*progress_callback*/) {
  throw std::runtime_error("CAR files are structures, use readStructure.");
}

correlation::core::Cell CarReader::read(const std::string &file_name) {
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" + std::strerror(errno) + ").");
  }

  correlation::core::Cell tempCell;
  std::string line;

  while (std::getline(myfile, line)) {

    // Ignore empty lines or comment lines
    if (line.empty() || line[0] == '!') {
      continue;
    }

    std::stringstream line_stream(line);
    std::string first_token;
    line_stream >> first_token;

    if (first_token == "PBC") {
      std::array<real_t, 6> lattice_params{};
      // The token "PBC" is consumed, so we read the 6 numbers that follow.
      if (line_stream >> lattice_params[0] >> lattice_params[1] >> lattice_params[2] >> lattice_params[3] >>
          lattice_params[4] >> lattice_params[5]) {
        tempCell.setLatticeParameters(lattice_params);
      }
      continue;
    }

    if (first_token == "PBC=OFF") {
      const std::array<real_t, 6> lattice_params = {100.0, 100.0, 100.0, 90.0, 90.0, 90.0};
      tempCell.setLatticeParameters(lattice_params);
      continue;
    }

    // Reset the stream to parse the full line as an atom entry.
    line_stream.clear();
    line_stream.seekg(0);

    // Declare variables for all 8 columns we need to read.
    std::string dummy_token_1;
    std::string dummy_token_5;
    std::string dummy_token_6;
    std::string dummy_token_7;
    std::string element;
    real_t coord_x = 0.0;
    real_t coord_y = 0.0;
    real_t coord_z = 0.0;

    // Read exactly 8 columns to get the element.
    if (line_stream >> dummy_token_1 >> coord_x >> coord_y >> coord_z >> dummy_token_5 >> dummy_token_6 >>
        dummy_token_7 >> element) {
      tempCell.addAtom(element, {coord_x, coord_y, coord_z});
    }
  }

  return tempCell;
}

} // namespace correlation::readers
