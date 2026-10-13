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
#include <cstddef>
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
 * @brief Full 118-element CPK lookup table based on official ASE/Jmol standards.
 *
 * Index 0 holds the fallback color (Deep Pink {255, 20, 147, 255}) for unassigned or invalid Z.
 * Indices 1 to 118 map directly to elements H (Z=1) through Og (Z=118) in O(1) time.
 */
inline constexpr std::array<RgbaColor, 119> kElementColors = {
    {/*   0: Xx */ RgbaColor{255, 20, 147, 255},
     /*   1: H  */ RgbaColor{255, 255, 255, 255},
     /*   2: He */ RgbaColor{217, 255, 255, 255},
     /*   3: Li */ RgbaColor{204, 128, 255, 255},
     /*   4: Be */ RgbaColor{194, 255, 0, 255},
     /*   5: B  */ RgbaColor{255, 181, 181, 255},
     /*   6: C  */ RgbaColor{144, 144, 144, 255},
     /*   7: N  */ RgbaColor{48, 80, 248, 255},
     /*   8: O  */ RgbaColor{255, 13, 13, 255},
     /*   9: F  */ RgbaColor{144, 224, 80, 255},
     /*  10: Ne */ RgbaColor{179, 227, 245, 255},
     /*  11: Na */ RgbaColor{171, 92, 242, 255},
     /*  12: Mg */ RgbaColor{138, 255, 0, 255},
     /*  13: Al */ RgbaColor{191, 166, 166, 255},
     /*  14: Si */ RgbaColor{240, 200, 160, 255},
     /*  15: P  */ RgbaColor{255, 128, 0, 255},
     /*  16: S  */ RgbaColor{255, 255, 48, 255},
     /*  17: Cl */ RgbaColor{31, 240, 31, 255},
     /*  18: Ar */ RgbaColor{128, 209, 227, 255},
     /*  19: K  */ RgbaColor{143, 64, 212, 255},
     /*  20: Ca */ RgbaColor{61, 255, 0, 255},
     /*  21: Sc */ RgbaColor{230, 230, 230, 255},
     /*  22: Ti */ RgbaColor{191, 194, 199, 255},
     /*  23: V  */ RgbaColor{166, 166, 171, 255},
     /*  24: Cr */ RgbaColor{138, 153, 199, 255},
     /*  25: Mn */ RgbaColor{156, 122, 199, 255},
     /*  26: Fe */ RgbaColor{224, 102, 51, 255},
     /*  27: Co */ RgbaColor{240, 144, 160, 255},
     /*  28: Ni */ RgbaColor{80, 208, 80, 255},
     /*  29: Cu */ RgbaColor{200, 128, 51, 255},
     /*  30: Zn */ RgbaColor{125, 128, 176, 255},
     /*  31: Ga */ RgbaColor{194, 143, 143, 255},
     /*  32: Ge */ RgbaColor{102, 143, 143, 255},
     /*  33: As */ RgbaColor{189, 128, 227, 255},
     /*  34: Se */ RgbaColor{255, 161, 0, 255},
     /*  35: Br */ RgbaColor{166, 41, 41, 255},
     /*  36: Kr */ RgbaColor{92, 184, 209, 255},
     /*  37: Rb */ RgbaColor{112, 46, 176, 255},
     /*  38: Sr */ RgbaColor{0, 255, 0, 255},
     /*  39: Y  */ RgbaColor{148, 255, 255, 255},
     /*  40: Zr */ RgbaColor{148, 224, 224, 255},
     /*  41: Nb */ RgbaColor{115, 194, 201, 255},
     /*  42: Mo */ RgbaColor{84, 181, 181, 255},
     /*  43: Tc */ RgbaColor{59, 158, 158, 255},
     /*  44: Ru */ RgbaColor{36, 143, 143, 255},
     /*  45: Rh */ RgbaColor{10, 125, 140, 255},
     /*  46: Pd */ RgbaColor{0, 105, 133, 255},
     /*  47: Ag */ RgbaColor{192, 192, 192, 255},
     /*  48: Cd */ RgbaColor{255, 217, 143, 255},
     /*  49: In */ RgbaColor{166, 117, 115, 255},
     /*  50: Sn */ RgbaColor{102, 128, 128, 255},
     /*  51: Sb */ RgbaColor{158, 99, 181, 255},
     /*  52: Te */ RgbaColor{212, 122, 0, 255},
     /*  53: I  */ RgbaColor{148, 0, 148, 255},
     /*  54: Xe */ RgbaColor{66, 158, 176, 255},
     /*  55: Cs */ RgbaColor{87, 23, 143, 255},
     /*  56: Ba */ RgbaColor{0, 201, 0, 255},
     /*  57: La */ RgbaColor{112, 212, 255, 255},
     /*  58: Ce */ RgbaColor{255, 255, 199, 255},
     /*  59: Pr */ RgbaColor{217, 255, 199, 255},
     /*  60: Nd */ RgbaColor{199, 255, 199, 255},
     /*  61: Pm */ RgbaColor{163, 255, 199, 255},
     /*  62: Sm */ RgbaColor{143, 255, 199, 255},
     /*  63: Eu */ RgbaColor{97, 255, 199, 255},
     /*  64: Gd */ RgbaColor{69, 255, 199, 255},
     /*  65: Tb */ RgbaColor{48, 255, 199, 255},
     /*  66: Dy */ RgbaColor{31, 255, 199, 255},
     /*  67: Ho */ RgbaColor{0, 255, 156, 255},
     /*  68: Er */ RgbaColor{0, 230, 117, 255},
     /*  69: Tm */ RgbaColor{0, 212, 82, 255},
     /*  70: Yb */ RgbaColor{0, 191, 56, 255},
     /*  71: Lu */ RgbaColor{0, 171, 36, 255},
     /*  72: Hf */ RgbaColor{77, 194, 255, 255},
     /*  73: Ta */ RgbaColor{77, 166, 255, 255},
     /*  74: W  */ RgbaColor{33, 148, 214, 255},
     /*  75: Re */ RgbaColor{38, 125, 171, 255},
     /*  76: Os */ RgbaColor{38, 102, 150, 255},
     /*  77: Ir */ RgbaColor{23, 84, 135, 255},
     /*  78: Pt */ RgbaColor{208, 208, 224, 255},
     /*  79: Au */ RgbaColor{255, 209, 35, 255},
     /*  80: Hg */ RgbaColor{184, 184, 208, 255},
     /*  81: Tl */ RgbaColor{166, 84, 77, 255},
     /*  82: Pb */ RgbaColor{87, 89, 97, 255},
     /*  83: Bi */ RgbaColor{158, 79, 181, 255},
     /*  84: Po */ RgbaColor{171, 92, 0, 255},
     /*  85: At */ RgbaColor{117, 79, 69, 255},
     /*  86: Rn */ RgbaColor{66, 130, 150, 255},
     /*  87: Fr */ RgbaColor{66, 0, 102, 255},
     /*  88: Ra */ RgbaColor{0, 125, 0, 255},
     /*  89: Ac */ RgbaColor{112, 171, 250, 255},
     /*  90: Th */ RgbaColor{0, 186, 255, 255},
     /*  91: Pa */ RgbaColor{0, 161, 255, 255},
     /*  92: U  */ RgbaColor{0, 143, 255, 255},
     /*  93: Np */ RgbaColor{0, 128, 255, 255},
     /*  94: Pu */ RgbaColor{0, 107, 255, 255},
     /*  95: Am */ RgbaColor{84, 92, 242, 255},
     /*  96: Cm */ RgbaColor{120, 92, 227, 255},
     /*  97: Bk */ RgbaColor{138, 79, 227, 255},
     /*  98: Cf */ RgbaColor{161, 54, 212, 255},
     /*  99: Es */ RgbaColor{179, 31, 212, 255},
     /* 100: Fm */ RgbaColor{179, 31, 186, 255},
     /* 101: Md */ RgbaColor{179, 13, 166, 255},
     /* 102: No */ RgbaColor{189, 13, 135, 255},
     /* 103: Lr */ RgbaColor{199, 0, 102, 255},
     /* 104: Rf */ RgbaColor{204, 0, 89, 255},
     /* 105: Db */ RgbaColor{209, 0, 79, 255},
     /* 106: Sg */ RgbaColor{217, 0, 69, 255},
     /* 107: Bh */ RgbaColor{224, 0, 56, 255},
     /* 108: Hs */ RgbaColor{230, 0, 46, 255},
     /* 109: Mt */ RgbaColor{235, 0, 38, 255},
     /* 110: Ds */ RgbaColor{232, 0, 41, 255},
     /* 111: Rg */ RgbaColor{235, 0, 36, 255},
     /* 112: Cn */ RgbaColor{237, 0, 31, 255},
     /* 113: Nh */ RgbaColor{239, 0, 26, 255},
     /* 114: Fl */ RgbaColor{242, 0, 20, 255},
     /* 115: Mc */ RgbaColor{244, 0, 15, 255},
     /* 116: Lv */ RgbaColor{247, 0, 10, 255},
     /* 117: Ts */ RgbaColor{250, 0, 5, 255},
     /* 118: Og */ RgbaColor{252, 0, 0, 255}}};

namespace detail {

/**
 * @struct SymbolZMapping
 * @brief Mapping between chemical symbol and atomic number (Z).
 */
struct SymbolZMapping {
  std::string_view symbol;
  int z;
};

/**
 * @brief Full chemical symbol to atomic number mapping for all 118 elements.
 */
inline constexpr std::array<SymbolZMapping, 118> symbol_z_table = {
    {{"H", 1},    {"He", 2},   {"Li", 3},   {"Be", 4},   {"B", 5},    {"C", 6},    {"N", 7},
     {"O", 8},    {"F", 9},    {"Ne", 10},  {"Na", 11},  {"Mg", 12},  {"Al", 13},  {"Si", 14},
     {"P", 15},   {"S", 16},   {"Cl", 17},  {"Ar", 18},  {"K", 19},   {"Ca", 20},  {"Sc", 21},
     {"Ti", 22},  {"V", 23},   {"Cr", 24},  {"Mn", 25},  {"Fe", 26},  {"Co", 27},  {"Ni", 28},
     {"Cu", 29},  {"Zn", 30},  {"Ga", 31},  {"Ge", 32},  {"As", 33},  {"Se", 34},  {"Br", 35},
     {"Kr", 36},  {"Rb", 37},  {"Sr", 38},  {"Y", 39},   {"Zr", 40},  {"Nb", 41},  {"Mo", 42},
     {"Tc", 43},  {"Ru", 44},  {"Rh", 45},  {"Pd", 46},  {"Ag", 47},  {"Cd", 48},  {"In", 49},
     {"Sn", 50},  {"Sb", 51},  {"Te", 52},  {"I", 53},   {"Xe", 54},  {"Cs", 55},  {"Ba", 56},
     {"La", 57},  {"Ce", 58},  {"Pr", 59},  {"Nd", 60},  {"Pm", 61},  {"Sm", 62},  {"Eu", 63},
     {"Gd", 64},  {"Tb", 65},  {"Dy", 66},  {"Ho", 67},  {"Er", 68},  {"Tm", 69},  {"Yb", 70},
     {"Lu", 71},  {"Hf", 72},  {"Ta", 73},  {"W", 74},   {"Re", 75},  {"Os", 76},  {"Ir", 77},
     {"Pt", 78},  {"Au", 79},  {"Hg", 80},  {"Tl", 81},  {"Pb", 82},  {"Bi", 83},  {"Po", 84},
     {"At", 85},  {"Rn", 86},  {"Fr", 87},  {"Ra", 88},  {"Ac", 89},  {"Th", 90},  {"Pa", 91},
     {"U", 92},   {"Np", 93},  {"Pu", 94},  {"Am", 95},  {"Cm", 96},  {"Bk", 97},  {"Cf", 98},
     {"Es", 99},  {"Fm", 100}, {"Md", 101}, {"No", 102}, {"Lr", 103}, {"Rf", 104}, {"Db", 105},
     {"Sg", 106}, {"Bh", 107}, {"Hs", 108}, {"Mt", 109}, {"Ds", 110}, {"Rg", 111}, {"Cn", 112},
     {"Nh", 113}, {"Fl", 114}, {"Mc", 115}, {"Lv", 116}, {"Ts", 117}, {"Og", 118}}};

} // namespace detail

/**
 * @brief Retrieves the atomic number (Z) for a chemical symbol.
 * @param symbol Element symbol (e.g. "H", "Fe", "Si").
 * @return Atomic number Z in [1, 118], or 0 if unknown.
 */
[[nodiscard]] constexpr int getAtomicNumber(std::string_view symbol) noexcept {
  const auto it = std::ranges::find_if(
      detail::symbol_z_table,
      [symbol](const detail::SymbolZMapping &entry) noexcept { return entry.symbol == symbol; });
  if (it != detail::symbol_z_table.end()) {
    return it->z;
  }
  return 0;
}

/**
 * @brief Retrieves standard CPK color for a given atomic number in O(1) time.
 * @param atomic_number Atomic number (Z).
 * @return RGBA color representation.
 */
[[nodiscard]] constexpr RgbaColor getElementColor(int atomic_number) noexcept {
  if (atomic_number >= 1 && atomic_number <= 118) {
    return kElementColors[static_cast<size_t>(atomic_number)];
  }
  return kElementColors[0];
}

/**
 * @brief Retrieves standard CPK color for a chemical symbol.
 * @param symbol Element symbol (e.g. "H", "Fe", "Si").
 * @return RGBA color representation.
 */
[[nodiscard]] constexpr RgbaColor getElementColor(std::string_view symbol) noexcept {
  return getElementColor(getAtomicNumber(symbol));
}

/**
 * @brief Retrieves visualization sphere radius scaled from covalent radius.
 * @param symbol Element symbol.
 * @param scale Scaling factor (default 0.4 for ball-and-stick).
 * @return Sphere radius in Angstroms.
 */
[[nodiscard]] inline float getElementDisplayRadius(std::string_view symbol, float scale = 0.4F) {
  try {
    const auto r = correlation::physics::getCovalentRadius(symbol);
    return static_cast<float>(r) * scale;
  } catch (...) {
    return 0.5F * scale;
  }
}

} // namespace correlation::render3d
