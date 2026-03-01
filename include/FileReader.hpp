// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <string>

#include "Cell.hpp"
#include "Trajectory.hpp"

// The FileReader namespace encapsulates all functionality related to reading
// structure and trajectory files. It delegates to specialized readers.
namespace FileReader {

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
 * @return A Cell object containing the parsed atomic structure.
 * @throws std::runtime_error if the file cannot be opened or is malformed.
 */
Cell readStructure(const std::string &filename, FileType type);

/**
 * @brief Reads a trajectory from a file.
 *
 * @param filename The path to the trajectory file.
 * @param type The format of the file.
 * @return A Trajectory object.
 * @throws std::runtime_error if the file cannot be opened or is malformed.
 */
Trajectory readTrajectory(const std::string &filename, FileType type);

/**
 * @brief Determines the FileType from a given filename extension.
 * @param filename The full path of the file.
 * @return The corresponding FileType enum value.
 * @throws std::runtime_error if the file has no extension or an unsupported
 * one.
 */
FileType determineFileType(const std::string &filename);

} // namespace FileReader
