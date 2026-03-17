// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include "FileReader.hpp"
#include "readers/ReaderFactory.hpp"

#include <algorithm>
#include <filesystem>
#include <stdexcept>

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

  return FileType::Unknown;
}

Cell readStructure(
    const std::string &filename, FileType type,
    std::function<void(float, const std::string &)> progress_callback) {
  
  std::string ext = std::filesystem::path(filename).extension().string();
  auto reader = ReaderFactory::instance().getReaderForExtension(ext);
  
  if (reader) {
    return reader->readStructure(filename, progress_callback);
  }

  throw std::runtime_error("No reader found for extension: " + ext);
}

Trajectory readTrajectory(
    const std::string &filename, FileType type,
    std::function<void(float, const std::string &)> progress_callback) {
  
  std::string ext = std::filesystem::path(filename).extension().string();
  auto reader = ReaderFactory::instance().getReaderForExtension(ext);
  
  if (reader) {
    if (reader->isTrajectory()) {
      return reader->readTrajectory(filename, progress_callback);
    } else {
      // If it's not a trajectory reader but we asked for a trajectory,
      // maybe we just want it as a single frame trajectory?
      // BaseReader doesn't guarantee readTrajectory works for structure-only files.
      // But some might.
      try {
        Cell c = reader->readStructure(filename, progress_callback);
        std::vector<Cell> frames;
        frames.push_back(std::move(c));
        return Trajectory(frames, 1.0);
      } catch (...) {
        throw std::runtime_error("Reader for " + ext + " does not support trajectory reading.");
      }
    }
  }

  throw std::runtime_error("No reader found for extension: " + ext);
}

} // namespace FileReader
