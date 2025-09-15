// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2025 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: MIT
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE

#pragma once

/**
 * @namespace constants
 * @brief Defines universal mathematical and physical constants.
 */
namespace constants {
// Use the standard library's high-precision pi from C++20
constexpr double pi =
    3.1415926535897932384626433832795028841971693993751058209749445923078164062;
constexpr double rad2deg = 180.0 / pi; // RADIANS -> DEGREES
constexpr double deg2rad = pi / 180.0; // DEGREES -> RADIANS
} // namespace constants
