/**
 * @file ElementColors.hpp
 * @brief CPK color palette and visualization radii for atomic 3D rendering.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "physics/PhysicalData.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <string_view>

namespace correlation::render3d {

/**
 * @struct RgbaColor
 * @brief 8-bit RGBA color representation for vertex colors and framebuffers.
 */
struct RgbaColor {
  uint8_t r{255};
  uint8_t g{255};
  uint8_t b{255};
  uint8_t a{255};

  [[nodiscard]] constexpr bool operator==(const RgbaColor &) const noexcept = default;

  [[nodiscard]] constexpr uint32_t toU32() const noexcept {
    return (static_cast<uint32_t>(r) << 24) | (static_cast<uint32_t>(g) << 16) |
           (static_cast<uint32_t>(b) << 8) | static_cast<uint32_t>(a);
  }
};

/**
 * @brief Retrieves standard CPK color for a given atomic number.
 * @param atomic_number Atomic number (Z).
 * @return RGBA color representation.
 */
[[nodiscard]] constexpr RgbaColor getElementColor(int atomic_number) noexcept {
  switch (atomic_number) {
  case 1:
    return RgbaColor{245, 245, 245, 255}; // H
  case 2:
    return RgbaColor{217, 255, 255, 255}; // He
  case 3:
    return RgbaColor{204, 128, 255, 255}; // Li
  case 4:
    return RgbaColor{194, 255, 0, 255}; // Be
  case 5:
    return RgbaColor{255, 181, 181, 255}; // B
  case 6:
    return RgbaColor{64, 64, 64, 255}; // C
  case 7:
    return RgbaColor{48, 80, 248, 255}; // N
  case 8:
    return RgbaColor{255, 13, 13, 255}; // O
  case 9:
    return RgbaColor{144, 224, 80, 255}; // F
  case 10:
    return RgbaColor{179, 227, 245, 255}; // Ne
  case 11:
    return RgbaColor{171, 92, 242, 255}; // Na
  case 12:
    return RgbaColor{138, 255, 0, 255}; // Mg
  case 13:
    return RgbaColor{191, 166, 166, 255}; // Al
  case 14:
    return RgbaColor{240, 200, 160, 255}; // Si
  case 15:
    return RgbaColor{255, 128, 0, 255}; // P
  case 16:
    return RgbaColor{255, 255, 48, 255}; // S
  case 17:
    return RgbaColor{31, 240, 31, 255}; // Cl
  case 18:
    return RgbaColor{128, 209, 227, 255}; // Ar
  case 19:
    return RgbaColor{143, 64, 212, 255}; // K
  case 20:
    return RgbaColor{61, 255, 0, 255}; // Ca
  case 21:
    return RgbaColor{230, 230, 230, 255}; // Sc
  case 22:
    return RgbaColor{191, 194, 199, 255}; // Ti
  case 23:
    return RgbaColor{166, 166, 171, 255}; // V
  case 24:
    return RgbaColor{138, 153, 199, 255}; // Cr
  case 25:
    return RgbaColor{156, 122, 199, 255}; // Mn
  case 26:
    return RgbaColor{224, 102, 51, 255}; // Fe
  case 27:
    return RgbaColor{240, 144, 160, 255}; // Co
  case 28:
    return RgbaColor{80, 208, 80, 255}; // Ni
  case 29:
    return RgbaColor{200, 115, 51, 255}; // Cu
  case 30:
    return RgbaColor{125, 128, 176, 255}; // Zn
  case 31:
    return RgbaColor{194, 143, 143, 255}; // Ga
  case 32:
    return RgbaColor{102, 143, 143, 255}; // Ge
  case 33:
    return RgbaColor{189, 128, 227, 255}; // As
  case 34:
    return RgbaColor{255, 161, 0, 255}; // Se
  case 35:
    return RgbaColor{166, 41, 41, 255}; // Br
  case 36:
    return RgbaColor{92, 184, 209, 255}; // Kr
  case 47:
    return RgbaColor{192, 192, 192, 255}; // Ag
  case 53:
    return RgbaColor{148, 0, 148, 255}; // I
  case 54:
    return RgbaColor{66, 158, 176, 255}; // Xe
  case 78:
    return RgbaColor{208, 208, 224, 255}; // Pt
  case 79:
    return RgbaColor{255, 209, 35, 255}; // Au
  case 80:
    return RgbaColor{184, 184, 208, 255}; // Hg
  case 82:
    return RgbaColor{87, 89, 97, 255}; // Pb
  case 92:
    return RgbaColor{0, 143, 255, 255}; // U
  default:
    return RgbaColor{255, 20, 147, 255}; // Deep Pink fallback
  }
}

namespace detail {

struct SymbolZMapping {
  std::string_view symbol;
  int z;
};

inline constexpr std::array<SymbolZMapping, 44> symbol_z_table = {
    {{"H", 1},   {"He", 2},  {"Li", 3},  {"Be", 4},  {"B", 5},   {"C", 6},   {"N", 7},   {"O", 8},
     {"F", 9},   {"Ne", 10}, {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15},  {"S", 16},
     {"Cl", 17}, {"Ar", 18}, {"K", 19},  {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23},  {"Cr", 24},
     {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}, {"Ga", 31}, {"Ge", 32},
     {"As", 33}, {"Se", 34}, {"Br", 35}, {"Kr", 36}, {"Ag", 47}, {"I", 53},  {"Xe", 54}, {"Pt", 78},
     {"Au", 79}, {"Hg", 80}, {"Pb", 82}, {"U", 92}}};

} // namespace detail

/**
 * @brief Retrieves standard CPK color for a chemical symbol.
 * @param symbol Element symbol (e.g. "H", "Fe", "Si").
 * @return RGBA color representation.
 */
[[nodiscard]] constexpr RgbaColor getElementColor(std::string_view symbol) noexcept {
  const auto it = std::ranges::find_if(
      detail::symbol_z_table,
      [symbol](const detail::SymbolZMapping &entry) noexcept { return entry.symbol == symbol; });
  if (it != detail::symbol_z_table.end()) {
    return getElementColor(it->z);
  }
  return RgbaColor{255, 20, 147, 255};
}

/**
 * @brief Retrieves visualization sphere radius scaled from covalent radius.
 * @param symbol Element symbol.
 * @param scale Scaling factor (default 0.4 for ball-and-stick).
 * @return Sphere radius in Angstroms.
 */
[[nodiscard]] inline float getElementDisplayRadius(std::string_view symbol, float scale = 0.4f) {
  try {
    const auto r = correlation::physics::getCovalentRadius(symbol);
    return static_cast<float>(r) * scale;
  } catch (...) {
    return 0.5f * scale;
  }
}

} // namespace correlation::render3d
