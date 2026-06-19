/**
 * @file PathFont.hpp
 * @brief High-fidelity Roboto vector font for SVG path rendering.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace correlation::plotters {

/**
 * @struct Glyph
 * @brief Represents a single vector character glyph.
 *
 * Contains the spacing information (left/right bounds) and a collection
 * of strokes defined as a series of (x,y) point pairs.
 */
struct Glyph {
  double left{0.0};                                            ///< Left side bearing.
  double right{0.0};                                           ///< Right side bearing.
  std::vector<std::vector<std::pair<double, double>>> strokes; ///< Vector paths for the character.
};

/**
 * @brief High-fidelity Roboto vector font implementation.
 */
class Roboto {
public:
  /** @return Singleton instance of the Roboto font. */
  static Roboto &instance();

  /**
   * @brief Renders a string of text as SVG path data.
   *
   * @param text The UTF-8 string to render.
   * @param start_x The starting x-coordinate.
   * @param start_y The starting y-coordinate (baseline).
   * @param font_size The font size (height in pixels).
   * @param anchor Text anchor: "start", "middle", or "end".
   * @return A string containing SVG path instructions (M/L/Z commands).
   */
  // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
  std::string render(const std::string &text, double start_x, double start_y, double font_size, const std::string &anchor = "start");

private:
  std::map<uint32_t, Glyph> glyphs_; ///< Local storage for character glyph data.

  Roboto();

  // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
  void add(uint32_t code_point, double left_bearing, double right_bearing,
           std::vector<std::vector<std::pair<double, double>>> stroke_paths);

  static std::vector<uint32_t> utf8ToUnicode(const std::string &utf8_string);

  void registerBasic_1();
  void registerBasic_2();
  void registerExtended();
};

} // namespace correlation::plotters
