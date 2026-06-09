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
#include <fstream>
#include <sstream>
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
 * @brief Sniffs the content of a file to determine the correct reader.
 * @param filename Path to the file.
 * @return Pointer to the matching reader, or nullptr if none found.
 */
static BaseReader *sniffReaderFromContent(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open())
    return nullptr;

  std::string line;
  std::vector<std::string> lines;
  for (int i = 0; i < 200 && std::getline(file, line); ++i) {
    lines.push_back(line);
  }

  if (lines.empty())
    return nullptr;

  // 1. LAMMPS Dump: look for ITEM: TIMESTEP or ITEM: NUMBER OF ATOMS
  for (const auto &l : lines) {
    if (l.find("ITEM: TIMESTEP") != std::string::npos || l.find("ITEM: NUMBER OF ATOMS") != std::string::npos) {
      return ReaderFactory::instance().getReaderForExtension(".dump");
    }
  }

  // 2. PDB: look for HEADER or CRYST1 or ATOM
  for (const auto &l : lines) {
    if (l.find("HEADER") == 0 || l.find("CRYST1") == 0 || l.find("ATOM  ") == 0) {
      return ReaderFactory::instance().getReaderForExtension(".pdb");
    }
  }

  // 3. Quantum Espresso
  for (const auto &l : lines) {
    std::string uline = l;
    std::transform(uline.begin(), uline.end(), uline.begin(), ::toupper);
    if (uline.find("CELL_PARAMETERS") != std::string::npos || uline.find("ATOMIC_POSITIONS") != std::string::npos ||
        uline.find("QUANTUM ESPRESSO") != std::string::npos || uline.find("PWSCF") != std::string::npos ||
        uline.find("&CONTROL") != std::string::npos || uline.find("&SYSTEM") != std::string::npos) {
      return ReaderFactory::instance().getReaderForExtension(".pwo");
    }
  }

  // 4. CP2K
  for (const auto &l : lines) {
    std::string uline = l;
    std::transform(uline.begin(), uline.end(), uline.begin(), ::toupper);
    if (uline.find("&CELL") != std::string::npos || uline.find("&COORD") != std::string::npos ||
        uline.find("&GLOBAL") != std::string::npos || uline.find("CP2K") != std::string::npos) {
      return ReaderFactory::instance().getReaderForExtension(".restart");
    }
  }

  // 5. VASP (POSCAR/CONTCAR or XDATCAR)
  for (const auto &l : lines) {
    std::string uline = l;
    std::transform(uline.begin(), uline.end(), uline.begin(), ::toupper);
    if (uline.find("DIRECT CONFIGURATION") != std::string::npos) {
      return ReaderFactory::instance().getReaderForExtension(".xdatcar");
    }
  }
  bool has_direct_or_cartesian = false;
  for (const auto &l : lines) {
    std::string uline = l;
    std::transform(uline.begin(), uline.end(), uline.begin(), ::toupper);
    if (uline.find("DIRECT") != std::string::npos || uline.find("CARTESIAN") != std::string::npos) {
      has_direct_or_cartesian = true;
      break;
    }
  }
  if (has_direct_or_cartesian) {
    return ReaderFactory::instance().getReaderForExtension(".poscar");
  }

  // 6. XYZ or Gromacs
  if (lines.size() >= 3) {
    try {
      size_t count1 = std::stoul(lines[0]);
      std::istringstream iss(lines[2]);
      std::string sym;
      double x, y, z;
      if (iss >> sym >> x >> y >> z) {
        return ReaderFactory::instance().getReaderForExtension(".xyz");
      }
    } catch (...) {}
  }

  return nullptr;
}

/**
 * @brief Finds a reader for the given filename.
 *
 * Tries the file extension first via ReaderFactory.  When the extension is
 * empty (extensionless VASP files such as POSCAR, CONTCAR, XDATCAR) falls
 * back to a case-insensitive basename check. If both fail, performs content
 * sniffing as a fallback.
 *
 * @param filename Path to the file.
 * @return Pointer to the matching reader, or nullptr if none found.
 */
static BaseReader *findReaderForFile(const std::string &filename) {
  std::string ext = std::filesystem::path(filename).extension().string();

  if (!ext.empty()) {
    auto *reader = ReaderFactory::instance().getReaderForExtension(ext, filename);
    if (reader)
      return reader;
  }

  // Extensionless files: try basename (handles POSCAR, CONTCAR, XDATCAR)
  std::string basename = std::filesystem::path(filename).filename().string();
  std::transform(basename.begin(), basename.end(), basename.begin(), ::tolower);

  if (basename == "poscar" || basename == "contcar")
    return ReaderFactory::instance().getReaderForExtension(".poscar");
  if (basename == "xdatcar")
    return ReaderFactory::instance().getReaderForExtension(".xdatcar");

  // Fallback to content sniffing
  return sniffReaderFromContent(filename);
}

static BaseReader *findReaderForType(FileType type) {
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

correlation::core::Cell readStructure(const std::string &filename, FileType type,
                                      std::function<void(float, const std::string &)> progress_callback) {

  BaseReader *reader = nullptr;
  if (type != FileType::Unknown) {
    reader = findReaderForType(type);
  }
  if (!reader) {
    reader = findReaderForFile(filename);
  }

  if (reader) {
    return reader->readStructure(filename, progress_callback);
  }

  throw std::runtime_error("No reader found for file: " + filename);
}

correlation::core::Trajectory readTrajectory(const std::string &filename, FileType type,
                                             std::function<void(float, const std::string &)> progress_callback) {

  BaseReader *reader = nullptr;
  if (type != FileType::Unknown) {
    reader = findReaderForType(type);
  }
  if (!reader) {
    reader = findReaderForFile(filename);
  }

  if (reader) {
    // Enforce 4 GiB trajectory file size limit.
    auto file_size = std::filesystem::file_size(filename);
    if (static_cast<std::uint64_t>(file_size) > correlation::core::kMaxTrajectoryBytes) {
      throw std::runtime_error("Trajectory file exceeds the 4 GiB memory limit: " + filename);
    }

    if (reader->isTrajectory()) {
      return reader->readTrajectory(filename, progress_callback);
    } else {
      // If it's not a trajectory reader but we asked for a trajectory,
      // wrap a single structure in a one-frame trajectory.
      try {
        correlation::core::Cell c = reader->readStructure(filename, progress_callback);
        std::vector<correlation::core::Cell> frames;
        frames.push_back(std::move(c));
        return correlation::core::Trajectory(frames, 1.0);
      } catch (...) {
        throw std::runtime_error("Reader for \"" + filename + "\" does not support trajectory reading.");
      }
    }
  }

  throw std::runtime_error("No reader found for file: " + filename);
}

} // namespace correlation::readers
