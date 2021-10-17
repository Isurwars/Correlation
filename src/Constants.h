#ifndef SRC_CONSTANTS_H_
#define SRC_CONSTANTS_H_
/* ---------------------------------------------------------------------
 * Correlation: An Analysis Tool for Liquids and for Amorphous Solids
 * Copyright (c) 2013-2021 Isaías Rodríguez <isurwars@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the MIT License version as published in:
 * https://github.com/Isurwars/Correlation/blob/main/LICENSE
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * ----------------------------------------------------------------------
 */
#include <string>

/*
 * Constants definitions
 */

namespace constants {
  constexpr double pi { 3.141592653589793238463 };      // PI value
  constexpr double rad2deg { 57.29577951308232087679 }; // RADIANS -> DEGREES
  constexpr double deg2rad { 0.017453292519943295769 }; // DEGREES -> RADIANS
}

/*
 * Covalent Radii definitions in Angstroms.
 * Based on the following sources:
 * Molecular single-bond Covalent Radii for elements 1-118.
 * Pyykko, P. and Atsumi, M., DOI: 10.1002/chem.200800987
 * Colvant Radii revisited, Cordero B, et al.
 * DOI: 10.1039/b801115j
 */

double Covalent_Radii(std::string);

#endif  // SRC_CONSTANTS_H_
