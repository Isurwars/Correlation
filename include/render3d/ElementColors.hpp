/**
 * @file ElementColors.hpp
 * @brief CPK color palette and visualization radii for atomic 3D rendering.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "physics/PhysicalData.hpp"
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
    return (static_cast<uint32_t>(r) << 24U) | (static_cast<uint32_t>(g) << 16U) |
           (static_cast<uint32_t>(b) << 8U) | static_cast<uint32_t>(a);
  }
};

/**
 * @brief Full 118-element CPK lookup table based on official ASE/Jmol standards.
 *
 * Index 0 holds the fallback color (Deep Pink {255, 20, 147, 255}) for unassigned or invalid Z.
 * Indices 1 to 118 map directly to elements H (Z=1) through Og (Z=118) in O(1) time.
 */
inline constexpr std::array<RgbaColor, 119> ELEMENT_COLORS = {
    {
        /*   0: Xx */ RgbaColor{.r = 255, .g = 20, .b = 147, .a = 255},
        /*   1: H  */ RgbaColor{.r = 255, .g = 255, .b = 255, .a = 255},
        /*   2: He */ RgbaColor{.r = 217, .g = 255, .b = 255, .a = 255},
        /*   3: Li */ RgbaColor{.r = 204, .g = 128, .b = 255, .a = 255},
        /*   4: Be */ RgbaColor{.r = 194, .g = 255, .b = 0, .a = 255},
        /*   5: B  */ RgbaColor{.r = 255, .g = 181, .b = 181, .a = 255},
        /*   6: C  */ RgbaColor{.r = 144, .g = 144, .b = 144, .a = 255},
        /*   7: N  */ RgbaColor{.r = 48, .g = 80, .b = 248, .a = 255},
        /*   8: O  */ RgbaColor{.r = 255, .g = 13, .b = 13, .a = 255},
        /*   9: F  */ RgbaColor{.r = 144, .g = 224, .b = 80, .a = 255},
        /*  10: Ne */ RgbaColor{.r = 179, .g = 227, .b = 245, .a = 255},
        /*  11: Na */ RgbaColor{.r = 171, .g = 92, .b = 242, .a = 255},
        /*  12: Mg */ RgbaColor{.r = 138, .g = 255, .b = 0, .a = 255},
        /*  13: Al */ RgbaColor{.r = 191, .g = 166, .b = 166, .a = 255},
        /*  14: Si */ RgbaColor{.r = 240, .g = 200, .b = 160, .a = 255},
        /*  15: P  */ RgbaColor{.r = 255, .g = 128, .b = 0, .a = 255},
        /*  16: S  */ RgbaColor{.r = 255, .g = 255, .b = 48, .a = 255},
        /*  17: Cl */ RgbaColor{.r = 31, .g = 240, .b = 31, .a = 255},
        /*  18: Ar */ RgbaColor{.r = 128, .g = 209, .b = 227, .a = 255},
        /*  19: K  */ RgbaColor{.r = 143, .g = 64, .b = 212, .a = 255},
        /*  20: Ca */ RgbaColor{.r = 61, .g = 255, .b = 0, .a = 255},
        /*  21: Sc */ RgbaColor{.r = 230, .g = 230, .b = 230, .a = 255},
        /*  22: Ti */ RgbaColor{.r = 191, .g = 194, .b = 199, .a = 255},
        /*  23: V  */ RgbaColor{.r = 166, .g = 166, .b = 171, .a = 255},
        /*  24: Cr */ RgbaColor{.r = 138, .g = 153, .b = 199, .a = 255},
        /*  25: Mn */ RgbaColor{.r = 156, .g = 122, .b = 199, .a = 255},
        /*  26: Fe */ RgbaColor{.r = 224, .g = 102, .b = 51, .a = 255},
        /*  27: Co */ RgbaColor{.r = 240, .g = 144, .b = 160, .a = 255},
        /*  28: Ni */ RgbaColor{.r = 80, .g = 208, .b = 80, .a = 255},
        /*  29: Cu */ RgbaColor{.r = 200, .g = 128, .b = 51, .a = 255},
        /*  30: Zn */ RgbaColor{.r = 125, .g = 128, .b = 176, .a = 255},
        /*  31: Ga */ RgbaColor{.r = 194, .g = 143, .b = 143, .a = 255},
        /*  32: Ge */ RgbaColor{.r = 102, .g = 143, .b = 143, .a = 255},
        /*  33: As */ RgbaColor{.r = 189, .g = 128, .b = 227, .a = 255},
        /*  34: Se */ RgbaColor{.r = 255, .g = 161, .b = 0, .a = 255},
        /*  35: Br */ RgbaColor{.r = 166, .g = 41, .b = 41, .a = 255},
        /*  36: Kr */ RgbaColor{.r = 92, .g = 184, .b = 209, .a = 255},
        /*  37: Rb */ RgbaColor{.r = 112, .g = 46, .b = 176, .a = 255},
        /*  38: Sr */ RgbaColor{.r = 0, .g = 255, .b = 0, .a = 255},
        /*  39: Y  */ RgbaColor{.r = 148, .g = 255, .b = 255, .a = 255},
        /*  40: Zr */ RgbaColor{.r = 148, .g = 224, .b = 224, .a = 255},
        /*  41: Nb */ RgbaColor{.r = 115, .g = 194, .b = 201, .a = 255},
        /*  42: Mo */ RgbaColor{.r = 84, .g = 181, .b = 181, .a = 255},
        /*  43: Tc */ RgbaColor{.r = 59, .g = 158, .b = 158, .a = 255},
        /*  44: Ru */ RgbaColor{.r = 36, .g = 143, .b = 143, .a = 255},
        /*  45: Rh */ RgbaColor{.r = 10, .g = 125, .b = 140, .a = 255},
        /*  46: Pd */ RgbaColor{.r = 0, .g = 105, .b = 133, .a = 255},
        /*  47: Ag */ RgbaColor{.r = 192, .g = 192, .b = 192, .a = 255},
        /*  48: Cd */ RgbaColor{.r = 255, .g = 217, .b = 143, .a = 255},
        /*  49: In */ RgbaColor{.r = 166, .g = 117, .b = 115, .a = 255},
        /*  50: Sn */ RgbaColor{.r = 102, .g = 128, .b = 128, .a = 255},
        /*  51: Sb */ RgbaColor{.r = 158, .g = 99, .b = 181, .a = 255},
        /*  52: Te */ RgbaColor{.r = 212, .g = 122, .b = 0, .a = 255},
        /*  53: I  */ RgbaColor{.r = 148, .g = 0, .b = 148, .a = 255},
        /*  54: Xe */ RgbaColor{.r = 66, .g = 158, .b = 176, .a = 255},
        /*  55: Cs */ RgbaColor{.r = 87, .g = 23, .b = 143, .a = 255},
        /*  56: Ba */ RgbaColor{.r = 0, .g = 201, .b = 0, .a = 255},
        /*  57: La */ RgbaColor{.r = 112, .g = 212, .b = 255, .a = 255},
        /*  58: Ce */ RgbaColor{.r = 255, .g = 255, .b = 199, .a = 255},
        /*  59: Pr */ RgbaColor{.r = 217, .g = 255, .b = 199, .a = 255},
        /*  60: Nd */ RgbaColor{.r = 199, .g = 255, .b = 199, .a = 255},
        /*  61: Pm */ RgbaColor{.r = 163, .g = 255, .b = 199, .a = 255},
        /*  62: Sm */ RgbaColor{.r = 143, .g = 255, .b = 199, .a = 255},
        /*  63: Eu */ RgbaColor{.r = 97, .g = 255, .b = 199, .a = 255},
        /*  64: Gd */ RgbaColor{.r = 69, .g = 255, .b = 199, .a = 255},
        /*  65: Tb */ RgbaColor{.r = 48, .g = 255, .b = 199, .a = 255},
        /*  66: Dy */ RgbaColor{.r = 31, .g = 255, .b = 199, .a = 255},
        /*  67: Ho */ RgbaColor{.r = 0, .g = 255, .b = 156, .a = 255},
        /*  68: Er */ RgbaColor{.r = 0, .g = 230, .b = 117, .a = 255},
        /*  69: Tm */ RgbaColor{.r = 0, .g = 212, .b = 82, .a = 255},
        /*  70: Yb */ RgbaColor{.r = 0, .g = 191, .b = 56, .a = 255},
        /*  71: Lu */ RgbaColor{.r = 0, .g = 171, .b = 36, .a = 255},
        /*  72: Hf */ RgbaColor{.r = 77, .g = 194, .b = 255, .a = 255},
        /*  73: Ta */ RgbaColor{.r = 77, .g = 166, .b = 255, .a = 255},
        /*  74: W  */ RgbaColor{.r = 33, .g = 148, .b = 214, .a = 255},
        /*  75: Re */ RgbaColor{.r = 38, .g = 125, .b = 171, .a = 255},
        /*  76: Os */ RgbaColor{.r = 38, .g = 102, .b = 150, .a = 255},
        /*  77: Ir */ RgbaColor{.r = 23, .g = 84, .b = 135, .a = 255},
        /*  78: Pt */ RgbaColor{.r = 208, .g = 208, .b = 224, .a = 255},
        /*  79: Au */ RgbaColor{.r = 255, .g = 209, .b = 35, .a = 255},
        /*  80: Hg */ RgbaColor{.r = 184, .g = 184, .b = 208, .a = 255},
        /*  81: Tl */ RgbaColor{.r = 166, .g = 84, .b = 77, .a = 255},
        /*  82: Pb */ RgbaColor{.r = 87, .g = 89, .b = 97, .a = 255},
        /*  83: Bi */ RgbaColor{.r = 158, .g = 79, .b = 181, .a = 255},
        /*  84: Po */ RgbaColor{.r = 171, .g = 92, .b = 0, .a = 255},
        /*  85: At */ RgbaColor{.r = 117, .g = 79, .b = 69, .a = 255},
        /*  86: Rn */ RgbaColor{.r = 66, .g = 130, .b = 150, .a = 255},
        /*  87: Fr */ RgbaColor{.r = 66, .g = 0, .b = 102, .a = 255},
        /*  88: Ra */ RgbaColor{.r = 0, .g = 125, .b = 0, .a = 255},
        /*  89: Ac */ RgbaColor{.r = 112, .g = 171, .b = 250, .a = 255},
        /*  90: Th */ RgbaColor{.r = 0, .g = 186, .b = 255, .a = 255},
        /*  91: Pa */ RgbaColor{.r = 0, .g = 161, .b = 255, .a = 255},
        /*  92: U  */ RgbaColor{.r = 0, .g = 143, .b = 255, .a = 255},
        /*  93: Np */ RgbaColor{.r = 0, .g = 128, .b = 255, .a = 255},
        /*  94: Pu */ RgbaColor{.r = 0, .g = 107, .b = 255, .a = 255},
        /*  95: Am */ RgbaColor{.r = 84, .g = 92, .b = 242, .a = 255},
        /*  96: Cm */ RgbaColor{.r = 120, .g = 92, .b = 227, .a = 255},
        /*  97: Bk */ RgbaColor{.r = 138, .g = 79, .b = 227, .a = 255},
        /*  98: Cf */ RgbaColor{.r = 161, .g = 54, .b = 212, .a = 255},
        /*  99: Es */ RgbaColor{.r = 179, .g = 31, .b = 212, .a = 255},
        /* 100: Fm */ RgbaColor{.r = 179, .g = 31, .b = 186, .a = 255},
        /* 101: Md */ RgbaColor{.r = 179, .g = 13, .b = 166, .a = 255},
        /* 102: No */ RgbaColor{.r = 189, .g = 13, .b = 135, .a = 255},
        /* 103: Lr */ RgbaColor{.r = 199, .g = 0, .b = 102, .a = 255},
        /* 104: Rf */ RgbaColor{.r = 204, .g = 0, .b = 89, .a = 255},
        /* 105: Db */ RgbaColor{.r = 209, .g = 0, .b = 79, .a = 255},
        /* 106: Sg */ RgbaColor{.r = 217, .g = 0, .b = 69, .a = 255},
        /* 107: Bh */ RgbaColor{.r = 224, .g = 0, .b = 56, .a = 255},
        /* 108: Hs */ RgbaColor{.r = 230, .g = 0, .b = 46, .a = 255},
        /* 109: Mt */ RgbaColor{.r = 235, .g = 0, .b = 38, .a = 255},
        /* 110: Ds */ RgbaColor{.r = 232, .g = 0, .b = 41, .a = 255},
        /* 111: Rg */ RgbaColor{.r = 235, .g = 0, .b = 36, .a = 255},
        /* 112: Cn */ RgbaColor{.r = 237, .g = 0, .b = 31, .a = 255},
        /* 113: Nh */ RgbaColor{.r = 239, .g = 0, .b = 26, .a = 255},
        /* 114: Fl */ RgbaColor{.r = 242, .g = 0, .b = 20, .a = 255},
        /* 115: Mc */ RgbaColor{.r = 244, .g = 0, .b = 15, .a = 255},
        /* 116: Lv */ RgbaColor{.r = 247, .g = 0, .b = 10, .a = 255},
        /* 117: Ts */ RgbaColor{.r = 250, .g = 0, .b = 5, .a = 255},
        /* 118: Og */ RgbaColor{.r = 252, .g = 0, .b = 0, .a = 255},
    },
};

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
inline constexpr std::array<SymbolZMapping, 118> SYMBOL_Z_TABLE = {
    {
        {.symbol = "H", .z = 1},    {.symbol = "He", .z = 2},   {.symbol = "Li", .z = 3},
        {.symbol = "Be", .z = 4},   {.symbol = "B", .z = 5},    {.symbol = "C", .z = 6},
        {.symbol = "N", .z = 7},    {.symbol = "O", .z = 8},    {.symbol = "F", .z = 9},
        {.symbol = "Ne", .z = 10},  {.symbol = "Na", .z = 11},  {.symbol = "Mg", .z = 12},
        {.symbol = "Al", .z = 13},  {.symbol = "Si", .z = 14},  {.symbol = "P", .z = 15},
        {.symbol = "S", .z = 16},   {.symbol = "Cl", .z = 17},  {.symbol = "Ar", .z = 18},
        {.symbol = "K", .z = 19},   {.symbol = "Ca", .z = 20},  {.symbol = "Sc", .z = 21},
        {.symbol = "Ti", .z = 22},  {.symbol = "V", .z = 23},   {.symbol = "Cr", .z = 24},
        {.symbol = "Mn", .z = 25},  {.symbol = "Fe", .z = 26},  {.symbol = "Co", .z = 27},
        {.symbol = "Ni", .z = 28},  {.symbol = "Cu", .z = 29},  {.symbol = "Zn", .z = 30},
        {.symbol = "Ga", .z = 31},  {.symbol = "Ge", .z = 32},  {.symbol = "As", .z = 33},
        {.symbol = "Se", .z = 34},  {.symbol = "Br", .z = 35},  {.symbol = "Kr", .z = 36},
        {.symbol = "Rb", .z = 37},  {.symbol = "Sr", .z = 38},  {.symbol = "Y", .z = 39},
        {.symbol = "Zr", .z = 40},  {.symbol = "Nb", .z = 41},  {.symbol = "Mo", .z = 42},
        {.symbol = "Tc", .z = 43},  {.symbol = "Ru", .z = 44},  {.symbol = "Rh", .z = 45},
        {.symbol = "Pd", .z = 46},  {.symbol = "Ag", .z = 47},  {.symbol = "Cd", .z = 48},
        {.symbol = "In", .z = 49},  {.symbol = "Sn", .z = 50},  {.symbol = "Sb", .z = 51},
        {.symbol = "Te", .z = 52},  {.symbol = "I", .z = 53},   {.symbol = "Xe", .z = 54},
        {.symbol = "Cs", .z = 55},  {.symbol = "Ba", .z = 56},  {.symbol = "La", .z = 57},
        {.symbol = "Ce", .z = 58},  {.symbol = "Pr", .z = 59},  {.symbol = "Nd", .z = 60},
        {.symbol = "Pm", .z = 61},  {.symbol = "Sm", .z = 62},  {.symbol = "Eu", .z = 63},
        {.symbol = "Gd", .z = 64},  {.symbol = "Tb", .z = 65},  {.symbol = "Dy", .z = 66},
        {.symbol = "Ho", .z = 67},  {.symbol = "Er", .z = 68},  {.symbol = "Tm", .z = 69},
        {.symbol = "Yb", .z = 70},  {.symbol = "Lu", .z = 71},  {.symbol = "Hf", .z = 72},
        {.symbol = "Ta", .z = 73},  {.symbol = "W", .z = 74},   {.symbol = "Re", .z = 75},
        {.symbol = "Os", .z = 76},  {.symbol = "Ir", .z = 77},  {.symbol = "Pt", .z = 78},
        {.symbol = "Au", .z = 79},  {.symbol = "Hg", .z = 80},  {.symbol = "Tl", .z = 81},
        {.symbol = "Pb", .z = 82},  {.symbol = "Bi", .z = 83},  {.symbol = "Po", .z = 84},
        {.symbol = "At", .z = 85},  {.symbol = "Rn", .z = 86},  {.symbol = "Fr", .z = 87},
        {.symbol = "Ra", .z = 88},  {.symbol = "Ac", .z = 89},  {.symbol = "Th", .z = 90},
        {.symbol = "Pa", .z = 91},  {.symbol = "U", .z = 92},   {.symbol = "Np", .z = 93},
        {.symbol = "Pu", .z = 94},  {.symbol = "Am", .z = 95},  {.symbol = "Cm", .z = 96},
        {.symbol = "Bk", .z = 97},  {.symbol = "Cf", .z = 98},  {.symbol = "Es", .z = 99},
        {.symbol = "Fm", .z = 100}, {.symbol = "Md", .z = 101}, {.symbol = "No", .z = 102},
        {.symbol = "Lr", .z = 103}, {.symbol = "Rf", .z = 104}, {.symbol = "Db", .z = 105},
        {.symbol = "Sg", .z = 106}, {.symbol = "Bh", .z = 107}, {.symbol = "Hs", .z = 108},
        {.symbol = "Mt", .z = 109}, {.symbol = "Ds", .z = 110}, {.symbol = "Rg", .z = 111},
        {.symbol = "Cn", .z = 112}, {.symbol = "Nh", .z = 113}, {.symbol = "Fl", .z = 114},
        {.symbol = "Mc", .z = 115}, {.symbol = "Lv", .z = 116}, {.symbol = "Ts", .z = 117},
        {.symbol = "Og", .z = 118},
    },
};

/**
 * @brief High-contrast categorical palette for numeric LAMMPS atom types ("1", "2", "3", ...).
 */
inline constexpr std::array<RgbaColor, 10> NUMERIC_TYPE_COLORS = {
    {
        {.r = 56, .g = 189, .b = 248},  // Type 1: Sky Blue (#38BDF8)
        {.r = 251, .g = 146, .b = 60},  // Type 2: Warm Orange (#FB923C)
        {.r = 52, .g = 211, .b = 153},  // Type 3: Emerald Green (#34D399)
        {.r = 250, .g = 204, .b = 21},  // Type 4: Gold Yellow (#FACC15)
        {.r = 167, .g = 139, .b = 250}, // Type 5: Soft Violet (#A78BFA)
        {.r = 248, .g = 113, .b = 113}, // Type 6: Coral Red (#F87171)
        {.r = 45, .g = 212, .b = 191},  // Type 7: Teal (#2DD4BF)
        {.r = 244, .g = 114, .b = 182}, // Type 8: Rose Pink (#F472B6)
        {.r = 163, .g = 230, .b = 53},  // Type 9: Lime (#A3E635)
        {.r = 96, .g = 165, .b = 250},  // Type 10: Cobalt Blue (#60A5FA)
    },
};

} // namespace detail

/**
 * @brief Retrieves the atomic number (Z) for a chemical symbol.
 * @param symbol Element symbol (e.g. "H", "Fe", "Si").
 * @return Atomic number Z in [1, 118], or 0 if unknown.
 */
[[nodiscard]] constexpr int getAtomicNumber(std::string_view symbol) noexcept {
  for (const auto &entry : detail::SYMBOL_Z_TABLE) {
    if (entry.symbol == symbol) {
      return entry.z;
    }
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
    return ELEMENT_COLORS.at(static_cast<size_t>(atomic_number));
  }
  return ELEMENT_COLORS[0];
}

/**
 * @brief Retrieves standard CPK color for a chemical symbol, or a categorical color for numeric
 *        LAMMPS atom types ("1", "2", ...).
 * @param symbol Element symbol (e.g. "H", "Fe", "Si", or "1").
 * @return RGBA color representation.
 */
[[nodiscard]] constexpr RgbaColor getElementColor(std::string_view symbol) noexcept {
  const int z = getAtomicNumber(symbol);
  if (z != 0) {
    return getElementColor(z);
  }
  if (!symbol.empty()) {
    size_t numeric_val = 0;
    bool all_digits = true;
    for (const char chr : symbol) {
      if (chr < '0' || chr > '9') {
        all_digits = false;
        break;
      }
      numeric_val = numeric_val * 10 + static_cast<size_t>(chr - '0');
    }
    if (all_digits && numeric_val >= 1) {
      return detail::NUMERIC_TYPE_COLORS.at((numeric_val - 1) % detail::NUMERIC_TYPE_COLORS.size());
    }
  }
  return ELEMENT_COLORS[0];
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
    return 1.25F * scale;
  }
}

} // namespace correlation::render3d
