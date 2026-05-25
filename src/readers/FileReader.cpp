/**
 * @file FileReader.cpp
 * @brief Implementation of the unified file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */
#include "readers/FileReader.hpp"
#include "readers/ReaderFactory.hpp"
#include "core/MappedFile.hpp"

#include <algorithm>
#include <filesystem>
#include <stdexcept>

namespace correlation::readers {

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
  if (ext == ".poscar" || ext == ".contcar" || ext == ".vasp")
    return FileType::Vasp;
  if (ext == ".xdatcar")
    return FileType::Xdatcar;
  if (ext == ".gro")
    return FileType::Gromacs;
  if (ext == ".pdb" || ext == ".ent")
    return FileType::Pdb;
  if (ext == ".xyz" || ext == ".exyz")
    return FileType::Xyz;

  // Check basename for extensionless VASP files (POSCAR, CONTCAR, XDATCAR)
  if (ext.empty() || ext == ".") {
    std::string basename = std::filesystem::path(filename).filename().string();
    std::transform(basename.begin(), basename.end(), basename.begin(), ::tolower);
    if (basename == "poscar" || basename == "contcar")
      return FileType::Vasp;
    if (basename == "xdatcar")
      return FileType::Xdatcar;
  }

  return FileType::Unknown;
}

/**
 * @brief Finds a reader for the given filename.
 *
 * Tries the file extension first via ReaderFactory.  When the extension is
 * empty (extensionless VASP files such as POSCAR, CONTCAR, XDATCAR) falls
 * back to a case-insensitive basename check.
 *
 * @param filename Path to the file.
 * @return Pointer to the matching reader, or nullptr if none found.
 */
static BaseReader *findReaderForFile(const std::string &filename) {
  std::string ext = std::filesystem::path(filename).extension().string();

  if (!ext.empty()) {
    auto *reader = ReaderFactory::instance().getReaderForExtension(ext);
    if (reader)
      return reader;
  }

  // Extensionless files: try basename (handles POSCAR, CONTCAR, XDATCAR)
  std::string basename =
      std::filesystem::path(filename).filename().string();
  std::transform(basename.begin(), basename.end(), basename.begin(),
                 ::tolower);

  if (basename == "poscar" || basename == "contcar")
    return ReaderFactory::instance().getReaderForExtension(".poscar");
  if (basename == "xdatcar")
    return ReaderFactory::instance().getReaderForExtension(".xdatcar");

  return nullptr;
}

correlation::core::Cell readStructure(
    const std::string &filename, FileType type,
    std::function<void(float, const std::string &)> progress_callback) {

  auto *reader = findReaderForFile(filename);

  if (reader) {
    return reader->readStructure(filename, progress_callback);
  }

  throw std::runtime_error("No reader found for file: " + filename);
}

correlation::core::Trajectory readTrajectory(
    const std::string &filename, FileType type,
    std::function<void(float, const std::string &)> progress_callback) {

  auto *reader = findReaderForFile(filename);

  if (reader) {
    // Enforce 4 GiB trajectory file size limit.
    auto file_size = std::filesystem::file_size(filename);
    if (static_cast<std::uint64_t>(file_size) >
        correlation::core::kMaxTrajectoryBytes) {
      throw std::runtime_error(
          "Trajectory file exceeds the 4 GiB memory limit: " + filename);
    }

    if (reader->isTrajectory()) {
      return reader->readTrajectory(filename, progress_callback);
    } else {
      // If it's not a trajectory reader but we asked for a trajectory,
      // wrap a single structure in a one-frame trajectory.
      try {
        correlation::core::Cell c =
            reader->readStructure(filename, progress_callback);
        std::vector<correlation::core::Cell> frames;
        frames.push_back(std::move(c));
        return correlation::core::Trajectory(frames, 1.0);
      } catch (...) {
        throw std::runtime_error("Reader for \"" + filename +
                                 "\" does not support trajectory reading.");
      }
    }
  }

  throw std::runtime_error("No reader found for file: " + filename);
}

} // namespace correlation::readers
