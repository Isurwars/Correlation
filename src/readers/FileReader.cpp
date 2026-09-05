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

namespace {

FileType getExtensionlessVaspType(const std::filesystem::path &filepath) {
  std::string basename = filepath.filename().string();
  std::ranges::transform(basename, basename.begin(), ::tolower);
  if (basename == "poscar" || basename == "contcar") {
    return FileType::Vasp;
  }
  if (basename == "xdatcar") {
    return FileType::Xdatcar;
  }
  return FileType::Unknown;
}

const std::unordered_map<std::string, FileType> &getExtensionTypeMap() {
  static const std::unordered_map<std::string, FileType> k_extension_map = {
      {".car", FileType::Car},           {".cell", FileType::Cell},        {".cif", FileType::Cif},
      {".arc", FileType::Arc},           {".dump", FileType::LammpsDump},  {".lammpstrj", FileType::LammpsDump},
      {".dat", FileType::OnetepDat},     {".md", FileType::CastepMd},      {".outmol", FileType::Outmol},
      {".poscar", FileType::Vasp},       {".contcar", FileType::Vasp},     {".vasp", FileType::Vasp},
      {".xdatcar", FileType::Xdatcar},   {".gro", FileType::Gromacs},      {".pdb", FileType::Pdb},
      {".ent", FileType::Pdb},           {".xyz", FileType::Xyz},          {".exyz", FileType::Xyz},
      {".mace", FileType::Mace},         {".extxyz", FileType::Mace},      {".chgnet", FileType::Chgnet},
      {".gap", FileType::Gap},           {".quip", FileType::Gap},
  };
  return k_extension_map;
}

} // namespace

FileType determineFileType(const std::string &filename) {
  const std::filesystem::path filepath(filename);
  std::string ext = filepath.extension().string();
  std::ranges::transform(ext, ext.begin(), ::tolower);

  const auto &map = getExtensionTypeMap();
  const auto iter = map.find(ext);
  if (iter != map.end()) {
    return iter->second;
  }

  // Check basename for extensionless VASP files (POSCAR, CONTCAR, XDATCAR)
  if (ext.empty() || ext == ".") {
    return getExtensionlessVaspType(filepath);
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
    auto *reader = ReaderFactory::instance().getReaderForExtension({.extension = ext, .filename = filename});
    if (reader != nullptr) {
      return reader;
    }
  }

  // Extensionless files: try basename (handles POSCAR, CONTCAR, XDATCAR)
  std::string basename = std::filesystem::path(filename).filename().string();
  std::ranges::transform(basename, basename.begin(), ::tolower);

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
  case FileType::Mace:
    return ReaderFactory::instance().getReaderForExtension(".mace");
  case FileType::Chgnet:
    return ReaderFactory::instance().getReaderForExtension(".chgnet");
  case FileType::Gap:
    return ReaderFactory::instance().getReaderForExtension(".gap");
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
