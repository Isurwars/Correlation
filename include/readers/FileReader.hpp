/**
 * @file FileReader.hpp
 * @brief Unified file reading interface for structure and trajectory files.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "core/Cell.hpp"
#include "core/Trajectory.hpp"

#include "readers/BaseReader.hpp"

#include <cstdint>
#include <expected>
#include <string>

// The correlation::readers namespace encapsulates all functionality related to
// reading structure and trajectory files. It delegates to specialized readers.
namespace correlation::readers {

/**
 * @brief A type-safe enum to specify the format of a structure or trajectory
 * file.
 */
enum class FileType : std::uint8_t {
  Car,        ///< Accelrys .car file format
  Cell,       ///< CASTEP .cell physical structure format
  Cif,        ///< Crystallographic Information File (.cif)
  OnetepDat,  ///< ONETEP input/output .dat format
  Arc,        ///< Accelrys .arc trajectory format
  LammpsDump, ///< LAMMPS atomic dump trajectory format
  CastepMd,   ///< CASTEP molecular dynamics .md trajectory format
  Outmol,     ///< DMol3 .outmol format
  Vasp,       ///< VASP POSCAR/CONTCAR single-structure format
  Xdatcar,    ///< VASP XDATCAR trajectory format
  Gromacs,    ///< GROMACS .gro format
  Pdb,        ///< PDB .pdb format
  Xyz,        ///< Extended XYZ .xyz format
  Mace,       ///< MACE Extended XYZ format (.mace, .extxyz)
  Chgnet,     ///< CHGNet Extended XYZ format (.chgnet)
  Gap,        ///< GAP Extended XYZ format (.gap, .quip)
  Nequip,     ///< NequIP Extended XYZ format (.nequip)
  Unknown,    ///< Unrecognized or unsupported file format
};

/**
 * @brief Reads an atomic structure from a file.
 *
 * @param[in] filename The path to the structure file.
 * @param[in] type The format of the file.
 * @param[in] progress_callback Optional callback invoked with progress (0.0–1.0)
 * and a status message.
 * @return A correlation::core::Cell object containing the parsed atomic
 * structure.
 * @throws std::runtime_error if the file cannot be opened or is malformed.
 */
correlation::core::Cell readStructure(const std::string &filename, FileType type,
                                      const ProgressCallback &progress_callback = nullptr);

/**
 * @brief Reads a trajectory from a file.
 *
 * @param[in] filename The path to the trajectory file.
 * @param[in] type The format of the file.
 * @param[in] progress_callback Optional callback invoked with progress (0.0–1.0)
 * and a status message.
 * @return A correlation::core::Trajectory object.
 * @throws std::runtime_error if the file cannot be opened or is malformed.
 */
correlation::core::Trajectory readTrajectory(const std::string &filename, FileType type,
                                             const ProgressCallback &progress_callback = nullptr);

/**
 * @brief Reads a structure from the given file, returning Cell or diagnostic error string.
 *
 * @param[in] filename Path to the file containing the structure.
 * @param[in] type File format type used to decode the file.
 * @param[in] progress_callback Optional callback invoked during reading.
 * @return Expected containing the parsed Cell on success, or diagnostic string on failure.
 */
[[nodiscard]] std::expected<correlation::core::Cell, std::string>
readStructureExpected(const std::string &filename, FileType type = FileType::Unknown,
                      const ProgressCallback &progress_callback = nullptr);

/**
 * @brief Reads a trajectory from the given file, returning Trajectory or diagnostic error string.
 *
 * @param[in] filename Path to the file containing the trajectory.
 * @param[in] type File format type used to decode the file.
 * @param[in] progress_callback Optional callback invoked during reading.
 * @return Expected containing the parsed Trajectory on success, or diagnostic string on failure.
 */
[[nodiscard]] std::expected<correlation::core::Trajectory, std::string>
readTrajectoryExpected(const std::string &filename, FileType type = FileType::Unknown,
                       const ProgressCallback &progress_callback = nullptr);

/**
 * @brief Determines the FileType from a given filename extension.
 * @param[in] filename The full path of the file.
 * @return The corresponding FileType enum value.
 * @throws std::runtime_error if the file has no extension or an unsupported
 * one.
 */
FileType determineFileType(const std::string &filename);

} // namespace correlation::readers
