// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#pragma once

#include <functional>
#include <string>
#include <vector>

#include "../Cell.hpp"

namespace CastepMdReader {

/**
 * @brief Reads a sequence of frames from a CASTEP .md trajectory file.
 * @param file_name The path to the .md file.
 * @param progress_callback Optional callback to report loading progress [0, 1].
 * @return A vector of Cell objects representing the trajectory frames.
 * @throws std::runtime_error if the file cannot be opened or parsed.
 */
std::vector<Cell> read(const std::string &file_name,
                       std::function<void(float, const std::string &)>
                           progress_callback = nullptr);

} // namespace CastepMdReader
