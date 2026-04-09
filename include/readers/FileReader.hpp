/**
 * @file FileReader.hpp
 * @brief Unified file reading interface for structure and trajectory files.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <functional>
#include <string>

#include "Cell.hpp"
#include "Trajectory.hpp"

// The correlation::readers namespace encapsulates all functionality related to
// reading structure and trajectory files. It delegates to specialized readers.
namespace correlation::readers {

/**
 * @brief A type-safe enum to specify the format of a structure or trajectory
 * file.
 */
enum class FileType {
  Car,        ///< Accelrys .car file format
  Cell,       ///< CASTEP .cell physical structure format
  Cif,        ///< Crystallographic Information File (.cif)
  OnetepDat,  ///< ONETEP input/output .dat format
  Arc,        ///< Accelrys .arc trajectory format
  LammpsDump, ///< LAMMPS atomic dump trajectory format
  CastepMd,   ///< CASTEP molecular dynamics .md trajectory format
  Outmol,     ///< DMol3 .outmol format
  Unknown     ///< Unrecognized or unsupported file format
};

/**
 * @brief Reads an atomic structure from a file.
 *
 * @param filename The path to the structure file.
 * @param type The format of the file.
 * @param progress_callback Optional callback invoked with progress (0.0–1.0)
 * and a status message.
 * @return A Cell object containing the parsed atomic structure.
 * @throws std::runtime_error if the file cannot be opened or is malformed.
 */
Cell readStructure(const std::string &filename, FileType type,
                   std::function<void(float, const std::string &)>
                       progress_callback = nullptr);

/**
 * @brief Reads a trajectory from a file.
 *
 * @param filename The path to the trajectory file.
 * @param type The format of the file.
 * @param progress_callback Optional callback invoked with progress (0.0–1.0)
 * and a status message.
 * @return A Trajectory object.
 * @throws std::runtime_error if the file cannot be opened or is malformed.
 */
Trajectory readTrajectory(const std::string &filename, FileType type,
                          std::function<void(float, const std::string &)>
                              progress_callback = nullptr);

/**
 * @brief Determines the FileType from a given filename extension.
 * @param filename The full path of the file.
 * @return The corresponding FileType enum value.
 * @throws std::runtime_error if the file has no extension or an unsupported
 * one.
 */
FileType determineFileType(const std::string &filename);

} // namespace correlation::readers
