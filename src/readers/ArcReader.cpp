/**
 * @file ArcReader.cpp
 * @brief Implementation of the ARC file format reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */
#include "readers/ArcReader.hpp"
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
#include <utility>

namespace correlation::readers {

// Automatic registration
// NOLINTNEXTLINE(cert-err58-cpp, bugprone-throwing-static-initialization)
static const bool registered = ReaderFactory::instance().registerReader(std::make_unique<ArcReader>());

correlation::core::Cell
ArcReader::readStructure(const std::string & /*filename*/,
                         std::function<void(float, const std::string &)> /*progress_callback*/) {
  throw std::runtime_error("ARC files are trajectories, use readTrajectory.");
}

correlation::core::Trajectory
ArcReader::readTrajectory(const std::string &filename,
                          std::function<void(float, const std::string &)> progress_callback) {
  return {read(filename, progress_callback), 1.0};
}

std::vector<correlation::core::Cell>
ArcReader::read(const std::string &file_name,
                const std::function<void(float, const std::string &)> &progress_callback) {
  std::ifstream myfile(file_name);
  if (!myfile.is_open()) {
    throw std::runtime_error("Unable to read file: " + file_name + " (" + std::strerror(errno) + ").");
  }

  std::vector<correlation::core::Cell> frames;
  correlation::core::Cell tempCell;
  std::string line;

  myfile.seekg(0, std::ios::end);
  std::streampos const file_size = myfile.tellg();
  myfile.seekg(0, std::ios::beg);
  std::streampos last_progress_pos = 0;
  size_t const update_interval = file_size / 100;

  while (std::getline(myfile, line)) {
    if (progress_callback) {
      updateProgress(myfile.tellg(), file_size, last_progress_pos, update_interval, progress_callback);
    }
    parseLine(line, tempCell, frames);
  }

  return frames;
}

void ArcReader::updateProgress(std::streampos current_pos, std::streampos file_size, std::streampos &last_progress_pos,
                               size_t update_interval,
                               const std::function<void(float, const std::string &)> &progress_callback) {
  if (std::cmp_greater(current_pos - last_progress_pos, update_interval)) {
    const float progress = static_cast<float>(current_pos) / static_cast<float>(file_size);
    progress_callback(progress, "Loading ARC file...");
    last_progress_pos = current_pos;
  }
}

void ArcReader::parseLine(const std::string &line, correlation::core::Cell &tempCell,
                          std::vector<correlation::core::Cell> &frames) {
  // Ignore empty lines or comment lines
  if (line.empty() || line[0] == '!') {
    return;
  }

  std::stringstream line_stream(line);
  std::string first_token;
  line_stream >> first_token;

  if (first_token == "end") {
    if (!tempCell.isEmpty()) {
      frames.push_back(std::move(tempCell));
      tempCell = correlation::core::Cell(); // Reset for next frame
    }
    return;
  }

  if (first_token == "PBC") {
    std::array<real_t, 6> lattice_params{};
    if (line_stream >> lattice_params[0] >> lattice_params[1] >> lattice_params[2] >> lattice_params[3] >>
        lattice_params[4] >> lattice_params[5]) {
      tempCell.setLatticeParameters(lattice_params);
    }
    return;
  }

  if (first_token == "PBC=OFF") {
    const std::array<real_t, 6> lattice_params = {100.0, 100.0, 100.0, 90.0, 90.0, 90.0};
    tempCell.setLatticeParameters(lattice_params);
    return;
  }

  if (first_token == "PBC=ON") {
    return;
  }

  // Check if it is a single token line (Energy)
  std::string second_token;
  if (!(line_stream >> second_token)) {
    real_t energy = 0.0;
    std::istringstream parse_stream(first_token);
    if (parse_stream >> energy) {
      tempCell.setEnergy(energy);
      return;
    }
  }

  // Attempt to parse atom
  // Reset stream to start of line
  line_stream.clear();
  line_stream.seekg(0);

  std::string dummy_token_1;
  std::string dummy_token_5;
  std::string dummy_token_6;
  std::string dummy_token_7;
  std::string element;
  real_t coord_x = 0.0;
  real_t coord_y = 0.0;
  real_t coord_z = 0.0;

  if (line_stream >> dummy_token_1 >> coord_x >> coord_y >> coord_z >> dummy_token_5 >> dummy_token_6 >>
      dummy_token_7 >> element) {
    tempCell.addAtom(element, {coord_x, coord_y, coord_z});
  }
}

} // namespace correlation::readers
