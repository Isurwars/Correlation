// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "FileReader.hpp"

#include <algorithm>
#include <filesystem>
#include <stdexcept>

#include "readers/ArcReader.hpp"
#include "readers/CarReader.hpp"
#include "readers/CastepMdReader.hpp"
#include "readers/CellReader.hpp"
#include "readers/CifReader.hpp"
#include "readers/LammpsDumpReader.hpp"
#include "readers/OnetepDatReader.hpp"
#include "readers/OutmolReader.hpp"

namespace FileReader {

FileType determineFileType(const std::string &filename) {
  std::string ext = std::filesystem::path(filename).extension().string();
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

  if (ext == ".car")
    return FileType::Car;
  if (ext == ".cell")
    return FileType::Cell;
  if (ext == ".cif")
    return FileType::Cif;
  if (ext == ".arc")
    return FileType::Arc;
  if (ext == ".dump" || ext == ".lammpstrj")
    return FileType::LammpsDump;
  if (ext == ".dat")
    return FileType::OnetepDat;
  if (ext == ".md")
    return FileType::CastepMd;
  if (ext == ".outmol")
    return FileType::Outmol;

  throw std::runtime_error("Unsupported file extension: " + ext);
}

Cell readStructure(const std::string &filename, FileType type) {
  switch (type) {
  case FileType::Car:
    return CarReader::read(filename);
  case FileType::Cif:
    return CifReader::read(filename);
  case FileType::Cell:
    return CellReader::read(filename);
  case FileType::OnetepDat:
    return OnetepDatReader::read(filename);
  case FileType::LammpsDump:
    return LammpsDumpReader::read(filename);
  case FileType::Arc:
    throw std::runtime_error("ARC files are trajectories, use readTrajectory.");
  case FileType::CastepMd:
    throw std::runtime_error(
        "CASTEP .md files are trajectories, use readTrajectory.");
  default:
    throw std::invalid_argument("Unknown file type specified.");
  }
}

Trajectory readTrajectory(const std::string &filename, FileType type) {
  switch (type) {
  case FileType::Arc: {
    std::vector<Cell> frames = ArcReader::read(filename);
    return Trajectory(frames, 1.0); // Default time_step 1.0 for now
  }
  case FileType::CastepMd: {
    std::vector<Cell> frames = CastepMdReader::read(filename);
    return Trajectory(frames, 1.0); // Default time_step 1.0 for now
  }
  case FileType::Outmol: {
    std::vector<Cell> frames = OutmolReader::read(filename);
    return Trajectory(frames, 1.0); // Default time_step 1.0 for now
  }
  default:
    throw std::runtime_error("Unsupported trajectory format.");
  }
}

} // namespace FileReader
