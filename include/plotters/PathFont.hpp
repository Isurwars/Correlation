/**
 * @file PathFont.hpp
 * @brief High-fidelity Roboto vector font for SVG path rendering.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#pragma once

#include "math/Precision.hpp"
#include <cstdint>
#include <map>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace correlation::plotters {

using Point2D = std::pair<real_t, real_t>;

/**
 * @struct Glyph
 * @brief Represents a single vector character glyph.
 *
 * Contains the spacing information (left/right bounds) and non-owning views
 * over contiguous point coordinates and stroke slice lengths.
 */
struct Glyph {
  real_t left{static_cast<real_t>(0.0)};  ///< Left side bearing.
  real_t right{static_cast<real_t>(0.0)}; ///< Right side bearing.
  std::span<const Point2D> points;        ///< Contiguous point coordinates.
  std::span<const size_t> stroke_lengths; ///< Number of points per stroke.
};

/**
 * @struct TextRenderParameters
 * @brief Parameters for vector text rendering.
 */
struct TextRenderParameters {
  std::string text;
  real_t start_x{static_cast<real_t>(0.0)};
  real_t start_y{static_cast<real_t>(0.0)};
  real_t font_size{static_cast<real_t>(12.0)};
  std::string anchor{"start"};
};

/**
 * @struct GlyphParameters
 * @brief Parameters for registering a vector character glyph.
 */
struct GlyphParameters {
  uint32_t code_point{0};
  real_t left_bearing{static_cast<real_t>(0.0)};
  real_t right_bearing{static_cast<real_t>(0.0)};
  std::span<const Point2D> points;
  std::span<const size_t> stroke_lengths;
};

/**
 * @brief High-fidelity Roboto vector font implementation.
 */
class Roboto {
public:
  /** @return Singleton instance of the Roboto font. */
  static Roboto &instance();

  /**
   * @brief Renders a string of text as SVG path data using parameters struct.
   *
   * @param params Rendering options including text, coordinates, font size, and text anchor.
   * @return A string containing SVG path instructions (M/L/Z commands).
   */
  std::string render(TextRenderParameters const &params);

private:
  std::map<uint32_t, Glyph> glyphs_; ///< Local storage for character glyph data.

  Roboto();

  void add(GlyphParameters const &params);
  void add(uint32_t code_point, real_t left_bearing, real_t right_bearing,
           std::span<const Point2D> points, std::span<const size_t> stroke_lengths);

  static std::vector<uint32_t> utf8ToUnicode(const std::string &utf8_string);

  void registerBasic1();
  void registerBasic2();
  void registerExtended();
};

} // namespace correlation::plotters
