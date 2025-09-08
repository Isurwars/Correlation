#ifndef INCLUDE_FILEIO_HPP_
#define INCLUDE_FILEIO_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#include <string>

#include "Cell.hpp"

namespace FileIO {
// --- Helper Struct for Symmetry Operations ---
struct SymmetryOp {
  linalg::Matrix3<double> rotation{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  linalg::Vector3<double> translation{0, 0, 0};

  // Applies the operation: new_pos = rotation * old_pos + translation
  linalg::Vector3<double> apply(const linalg::Vector3<double> &pos) const {
    return rotation * pos + translation;
  }
};
// --- Helper Functions ---
void toLower(std::string &);
std::string cleanCifValue(std::string);
void parseSymmetryComponent(std::string, int, SymmetryOp &);
SymmetryOp parseSymmetryString(const std::string &);
std::vector<std::string> tokenizeCifLine(const std::string &);
// --- Main Reader Functions ---
Cell readCar(const std::string &);
Cell readCell(const std::string &);
Cell readCif(const std::string &);
Cell readLammpsDump(const std::string &);
Cell readOnetepDat(const std::string &);

} // namespace FileIO

#endif // INCLUDE_FILEIO_HPP_
