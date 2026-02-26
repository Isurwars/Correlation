// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#pragma once

#include <string>
#include <vector>

#include "../Cell.hpp"

namespace CastepMdReader {

/**
 * @brief Reads a sequence of frames from a CASTEP .md trajectory file.
 * @param file_name The path to the .md file.
 * @return A vector of Cell objects representing the trajectory frames.
 * @throws std::runtime_error if the file cannot be read or parsed correctly.
 */
std::vector<Cell> read(const std::string &file_name);

} // namespace CastepMdReader
