/**
 * @file OnetepDatReader.cpp
 * @brief Implementation of the ONETEP .dat file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */
#include "readers/OnetepDatReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <fstream>
#include <stdexcept>
#include <memory>
#include <functional>

#include "Cell.hpp"
#include "Trajectory.hpp"

namespace correlation::readers {

// Automatic registration
static bool registered = ReaderFactory::instance().registerReader(
    std::make_unique<OnetepDatReader>()
);

Cell OnetepDatReader::readStructure(const std::string &filename,
                                     std::function<void(float, const std::string &)>
                                         progress_callback) {
  return read(filename);
}

Trajectory OnetepDatReader::readTrajectory(const std::string &filename,
                                            std::function<void(float, const std::string &)>
                                                progress_callback) {
  throw std::runtime_error("ONETEP DAT files are structures, use readStructure.");
}

Cell OnetepDatReader::read(const std::string &file_name) {
  /*
   * This is a Stub for reading a ONETEP file and return and empty cell.
   */
  std::ifstream myfile(file_name);
  // std::string   line;
  Cell tempCell;
  // std::smatch match;

  return tempCell;
}

} // namespace OnetepDatReader
