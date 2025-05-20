#ifndef INCLUDE_READFILES_HPP_
#define INCLUDE_READFILES_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include <string>

#include "../include/Cell.hpp"
//---------------------------------------------------------------------------//
//-------------------------------- Methods ----------------------------------//
//---------------------------------------------------------------------------//
std::vector<std::vector<double>> readBond(std::string, Cell);
Cell readCar(std::string);
Cell readCell(std::string);
Cell readOnetepDat(std::string);
Cell readCif(std::string);
Cell readLammpsDump(std::string);

#endif // INCLUDE_READFILES_HPP_
