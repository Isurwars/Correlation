// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#include <string>

#include "Cell.hpp"
#include "Trajectory.hpp"

// The FileIO namespace encapsulates all functionality related to reading and
// writing structure files.
namespace FileIO {

// A type-safe enum to specify the format of a structure file.
// This is more robust than using raw string extensions.
enum class FileType { Car, Cell, Cif, OnetepDat, Arc, Unknown };

// --- Main Public Interface ---

/**
 * @brief Reads an atomic structure from a file.
 *
 * This function serves as a single entry point for reading all supported file
 * formats. It dispatches to the appropriate specialized reader based on the
 * provided FileType.
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
 * This function reads a series of frames from a file (e.g., ARC) and returns
 * a Trajectory object.
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

} // namespace FileIO
