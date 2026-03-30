// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

#if __cplusplus >= 202002L
#include <numbers>
#endif

/**
 * @namespace correlation::math
 * @brief Defines universal mathematical and physical constants.
 */
namespace correlation::math {

// Use the standard library's high-precision pi from C++20 if available
// or a high-precision constant.
#if __cplusplus >= 202002L
constexpr double pi = std::numbers::pi;
#else
constexpr double pi =
    3.1415926535897932384626433832795028841971693993751058209749445923078164062;
#endif

constexpr double rad_to_deg = 180.0 / pi;           // RADIANS -> DEGREES
constexpr double deg_to_rad = pi / 180.0;           // DEGREES -> RADIANS
constexpr double thz_to_cminv = 33.35641;           // THz -> cm^-1
constexpr double thz_to_mev = 4.135667;             // THz -> meV
constexpr double bohr_to_angstrom = 0.529177210903; // Bohr -> Angstroms
constexpr double angstrom_to_bohr = 1.0 / bohr_to_angstrom;
constexpr double two_pi = 2.0 * pi;
constexpr double four_pi = 4.0 * pi;
constexpr double kb_ev_per_k = 8.617333262e-5;      // eV / K
constexpr double hbar_ev_ps = 6.582119569e-4;       // eV·ps

} // namespace correlation::math
