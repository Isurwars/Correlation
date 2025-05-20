#ifndef INCLUDE_WRITE_FILES_HPP_
#define INCLUDE_WRITE_FILES_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include <string>

#include "Cell.hpp"      // Cell actions

//-------------------------------------------------------------------------//
//------------------------------- Methods ---------------------------------//
//-------------------------------------------------------------------------//

// Write CSV (Cell, out_file_name, Smoothing)
void WriteCSV(Cell, std::string, bool = false);

#endif // INCLUDE_WRITE_FILES_HPP_
