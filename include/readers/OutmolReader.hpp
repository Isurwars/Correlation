// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#pragma once

#include <string>
#include <vector>

#include "Cell.hpp"

namespace OutmolReader {
std::vector<Cell> read(const std::string &file_name);
}
