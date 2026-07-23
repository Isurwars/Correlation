/**
 * @file Constants.hpp
 * @brief Mathematical and physical constants for the Correlation tool.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/Precision.hpp"

#if __cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)
#include <numbers>
#endif

/**
 * @namespace correlation::math
 * @brief Defines universal mathematical and physical constants.
 */
namespace correlation::math {

/** @brief High-precision mathematical pi constant. */
#if __cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)
constexpr real_t pi = static_cast<real_t>(std::numbers::pi);
#else
constexpr real_t pi = static_cast<real_t>(3.1415926535897932384626433832795028841971693993751058209749445923078164062);
#endif

constexpr real_t rad_to_deg = static_cast<real_t>(180.0) / pi;           ///< Conversion factor: Radians to Degrees
constexpr real_t deg_to_rad = pi / static_cast<real_t>(180.0);           ///< Conversion factor: Degrees to Radians
constexpr real_t thz_to_cminv = static_cast<real_t>(33.35641);           ///< Conversion factor: THz to cm^-1
constexpr real_t thz_to_mev = static_cast<real_t>(4.135667);             ///< Conversion factor: THz to meV
constexpr real_t bohr_to_angstrom = static_cast<real_t>(0.529177210903); ///< Conversion factor: Bohr to Angstroms
constexpr real_t angstrom_to_bohr =
    static_cast<real_t>(1.0) / bohr_to_angstrom;                    ///< Conversion factor: Angstroms to Bohr
constexpr real_t two_pi = static_cast<real_t>(2.0) * pi;            ///< 2 * pi
constexpr real_t four_pi = static_cast<real_t>(4.0) * pi;           ///< 4 * pi
constexpr real_t kb_ev_per_k = static_cast<real_t>(8.617333262e-5); ///< Boltzmann constant in eV/K
constexpr real_t hbar_ev_ps = static_cast<real_t>(6.582119569e-4);  ///< Reduced Planck constant in eV*ps

} // namespace correlation::math
