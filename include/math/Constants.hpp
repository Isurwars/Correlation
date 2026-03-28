// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#if __cplusplus >= 202002L
#include <numbers>
#endif

/**
 * @namespace correlation::math::constants
 * @brief Defines universal mathematical and physical constants.
 */
namespace correlation::math::constants {

// Use the standard library's high-precision pi from C++20 if available
// or a high-precision constant.
#if __cplusplus >= 202002L
constexpr double pi = std::numbers::pi;
#else
constexpr double pi =
    3.1415926535897932384626433832795028841971693993751058209749445923078164062;
#endif

constexpr double rad2deg = 180.0 / pi;    // RADIANS -> DEGREES
constexpr double deg2rad = pi / 180.0;    // DEGREES -> RADIANS
constexpr double THz_to_cmInv = 33.35641; // THz -> cm^-1
constexpr double THz_to_meV = 4.135667;   // THz -> meV
constexpr double BOHR_TO_ANGSTROM = 0.529177210903; // Bohr -> Angstroms

} // namespace correlation::math::constants
