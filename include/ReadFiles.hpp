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
void toLower(std::string &);
static std::string cleanToken(std::string);
static void replaceAll(std::string &, const std::string &, const std::string &);
std::vector<std::vector<double>> readBond(const std::string &, const Cell &);
Cell readCar(const std::string &);
Cell readCell(const std::string &);
Cell readCif(const std::string &);
Cell readLammpsDump(const std::string &);
Cell readOnetepDat(std::string);

#endif // INCLUDE_READFILES_HPP_
