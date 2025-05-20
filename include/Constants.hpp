#ifndef INCLUDE_CONSTANTS_HPP_
#define INCLUDE_CONSTANTS_HPP_
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright (c) 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
#include <string>
#include <vector>

/*-----------------------------------------------------------------------------
 * Constants definitions
 *-----------------------------------------------------------------------------
 */

namespace constants {
constexpr double pi{3.141592653589793238463};      // PI value
constexpr double rad2deg{57.29577951308232087679}; // RADIANS -> DEGREES
constexpr double deg2rad{0.017453292519943295769}; // DEGREES -> RADIANS
} // namespace constants

/*-----------------------------------------------------------------------------
 * Covalent Radii definitions in Angstroms.
 * Based on the following sources:
 * Molecular single-bond Covalent Radii for elements 1-118.
 * Pyykko, P. and Atsumi, M., DOI: 10.1002/chem.200800987
 * Colvant Radii revisited, Cordero B, et al.
 * DOI: 10.1039/b801115j
 *-----------------------------------------------------------------------------
 */

double covalentRadii(const std::string&);

/*-----------------------------------------------------------------------------
 * Atomic form factors parameters are adimensional.
 * The atomic form factor parameters were taken from the International Tables
 * for Crystallography: http://it.iucr.org/CB/ch6o1v0001/
 *-----------------------------------------------------------------------------
 */
std::vector<double> atomicFormFactorParameters(const std::string&);

#endif // INCLUDE_CONSTANTS_HPP_
