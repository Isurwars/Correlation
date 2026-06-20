/**
 * @file FileReader.cpp
 * @brief Implementation of the unified file reader.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */
#include "readers/FileReader.hpp"
#include "core/MappedFile.hpp"
#include "readers/ReaderFactory.hpp"

#include <algorithm>
#include <filesystem>
#include <stdexcept>
#include <utility>

namespace correlation::readers {

FileType determineFileType(const std::string &filename) {
  std::string ext = std::filesystem::path(filename).extension().string();
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

  if (ext == ".car") {
    return FileType::Car;
  }
  if (ext == ".cell") {
    return FileType::Cell;
  }
  if (ext == ".cif") {
    return FileType::Cif;
  }
  if (ext == ".arc") {
    return FileType::Arc;
  }
  if (ext == ".dump" || ext == ".lammpstrj") {
    return FileType::LammpsDump;
  }
  if (ext == ".dat") {
    return FileType::OnetepDat;
  }
  if (ext == ".md") {
    return FileType::CastepMd;
  }
  if (ext == ".outmol") {
    return FileType::Outmol;
  }
  if (ext == ".poscar" || ext == ".contcar" || ext == ".vasp") {
    return FileType::Vasp;
  }
  if (ext == ".xdatcar") {
    return FileType::Xdatcar;
  }
  if (ext == ".gro") {
    return FileType::Gromacs;
  }
  if (ext == ".pdb" || ext == ".ent") {
    return FileType::Pdb;
  }
  if (ext == ".xyz" || ext == ".exyz") {
    return FileType::Xyz;
  }

  // Check basename for extensionless VASP files (POSCAR, CONTCAR, XDATCAR)
  if (ext.empty() || ext == ".") {
    std::string basename = std::filesystem::path(filename).filename().string();
    std::transform(basename.begin(), basename.end(), basename.begin(), ::tolower);
    if (basename == "poscar" || basename == "contcar") {
      return FileType::Vasp;
    }
    if (basename == "xdatcar") {
      return FileType::Xdatcar;
    }
  }

  return FileType::Unknown;
}

namespace {

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
BaseReader *findReaderForFile(const std::string &filename) {
  std::string const ext = std::filesystem::path(filename).extension().string();

  if (!ext.empty()) {
    auto *reader = ReaderFactory::instance().getReaderForExtension({ext, filename});
    if (reader != nullptr) {
      return reader;
    }
  }

  // Extensionless files: try basename (handles POSCAR, CONTCAR, XDATCAR)
  std::string basename = std::filesystem::path(filename).filename().string();
  std::transform(basename.begin(), basename.end(), basename.begin(), ::tolower);

  if (basename == "poscar" || basename == "contcar") {
    return ReaderFactory::instance().getReaderForExtension(".poscar");
  }
  if (basename == "xdatcar") {
    return ReaderFactory::instance().getReaderForExtension(".xdatcar");
  }

  return nullptr;
}

BaseReader *findReaderForType(FileType type) {
  switch (type) {
  case FileType::Car:
    return ReaderFactory::instance().getReaderForExtension(".car");
  case FileType::Cell:
    return ReaderFactory::instance().getReaderForExtension(".cell");
  case FileType::Cif:
    return ReaderFactory::instance().getReaderForExtension(".cif");
  case FileType::Arc:
    return ReaderFactory::instance().getReaderForExtension(".arc");
  case FileType::LammpsDump:
    return ReaderFactory::instance().getReaderForExtension(".dump");
  case FileType::OnetepDat:
    return ReaderFactory::instance().getReaderForExtension(".dat");
  case FileType::CastepMd:
    return ReaderFactory::instance().getReaderForExtension(".md");
  case FileType::Outmol:
    return ReaderFactory::instance().getReaderForExtension(".outmol");
  case FileType::Vasp:
    return ReaderFactory::instance().getReaderForExtension(".poscar");
  case FileType::Xdatcar:
    return ReaderFactory::instance().getReaderForExtension(".xdatcar");
  case FileType::Gromacs:
    return ReaderFactory::instance().getReaderForExtension(".gro");
  case FileType::Pdb:
    return ReaderFactory::instance().getReaderForExtension(".pdb");
  case FileType::Xyz:
    return ReaderFactory::instance().getReaderForExtension(".xyz");
  default:
    return nullptr;
  }
}

} // namespace

correlation::core::Cell readStructure(const std::string &filename, FileType type,
                                      std::function<void(float, const std::string &)> progress_callback) {

  BaseReader *reader = nullptr;
  if (type != FileType::Unknown) {
    reader = findReaderForType(type);
  }
  if (reader == nullptr) {
    reader = findReaderForFile(filename);
  }

  if (reader != nullptr) {
    return reader->readStructure(filename, std::move(progress_callback));
  }

  throw std::runtime_error("No reader found for file: " + filename);
}

correlation::core::Trajectory readTrajectory(const std::string &filename, FileType type,
                                             const std::function<void(float, const std::string &)> &progress_callback) {

  BaseReader *reader = nullptr;
  if (type != FileType::Unknown) {
    reader = findReaderForType(type);
  }
  if (reader == nullptr) {
    reader = findReaderForFile(filename);
  }

  if (reader != nullptr) {
    // Enforce 4 GiB trajectory file size limit.
    auto file_size = std::filesystem::file_size(filename);
    if (static_cast<std::uint64_t>(file_size) > correlation::core::kMaxTrajectoryBytes) {
      throw std::runtime_error("Trajectory file exceeds the 4 GiB memory limit: " + filename);
    }

    if (reader->isTrajectory()) {
      return reader->readTrajectory(filename, progress_callback);
    }
    // wrap a single structure in a one-frame trajectory.
    try {
      correlation::core::Cell cell = reader->readStructure(filename, progress_callback);
      std::vector<correlation::core::Cell> frames;
      frames.push_back(std::move(cell));
      return {frames, 1.0};
    } catch (...) {
      throw std::runtime_error("Reader for \"" + filename + "\" does not support trajectory reading.");
    }
  }

  throw std::runtime_error("No reader found for file: " + filename);
}

} // namespace correlation::readers
