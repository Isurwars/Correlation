#ifndef SRC_CONSTANTS_H_
#define SRC_CONSTANTS_H_

/*
 * Copyright [2021] <@isurwars>
 */
#include <string>

/*
 * Constants definitions
 */

namespace constants {
constexpr double pi { 3.141592653589793238463 };      // PI value
constexpr double rad2deg { 57.29577951308232087679 };  // RADIANS -> DEGREES
constexpr double deg2rad { 0.017453292519943295769 };  // DEGREES -> RADIANS
}

/*
 * Covalent Radii definitions in Angstroms.
 * Based on the following sources:
 * Molecular single-bond Covalent Radii for elements 1-118.
 * Pyykko, P. and Atsumi, M., DOI: 10.1002/chem.200800987
 * Colvant Radii revisited, Cordero B, et al.
 * DOI: 10.1039/b801115j
 */

double
Covalent_Radii(
  std::string);


#endif  // SRC_CONSTANTS_H_
