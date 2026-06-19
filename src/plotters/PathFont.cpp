/**
 * @file PathFont.cpp
 * @brief High-fidelity Roboto vector font implementation details.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "plotters/PathFont.hpp"
#include <format>
#include <utility>

namespace correlation::plotters {

Roboto &Roboto::instance() {
  static Roboto inst;
  return inst;
}

Roboto::Roboto() {
  registerBasic_1();
  registerBasic_2();
  registerExtended();
}

void Roboto::add(const uint32_t code_point, const double left_bearing, const double right_bearing,
                 std::vector<std::vector<std::pair<double, double>>> stroke_paths) {
  glyphs_[code_point] = Glyph{
    .left = left_bearing,
    .right = right_bearing,
    .strokes = std::move(stroke_paths)
  };
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
std::string Roboto::render(const std::string &text, const double start_x, const double start_y, const double font_size, const std::string &anchor) {
  const double scale = font_size; // Glyphs are normalized to EM square of size 1.0
  const std::vector<uint32_t> codes = utf8ToUnicode(text);

  double total_width = 0.0;
  for (const uint32_t code_point : codes) {
    if (glyphs_.contains(code_point)) {
      total_width += (glyphs_.at(code_point).right - glyphs_.at(code_point).left) * scale;
    }
  }

  double cur_x = start_x;
  if (anchor == "middle") {
    cur_x -= total_width / 2.0;
  } else if (anchor == "end") {
    cur_x -= total_width;
  }

  std::string path_data;
  for (const uint32_t code_point : codes) {
    if (glyphs_.contains(code_point)) {
      const auto &glyph = glyphs_.at(code_point);
      const double offset_x = cur_x; // since left is 0.0, offset_x is just cur_x
      for (const auto &stroke : glyph.strokes) {
        for (size_t idx = 0; idx < stroke.size(); ++idx) {
          const double point_x = offset_x + stroke[idx].first * scale;
          const double point_y = start_y - stroke[idx].second * scale; // Flip Y direction (TTF Y is up, SVG Y is down)
          if (idx == 0) {
            path_data += std::format("M {:.2f} {:.2f} ", point_x, point_y);
          } else {
            path_data += std::format("L {:.2f} {:.2f} ", point_x, point_y);
          }
        }
        if (!stroke.empty()) {
          path_data += "Z ";
        }
      }
      cur_x += glyph.right * scale; // Advance by the width
    }
  }

  return path_data;
}

std::vector<uint32_t> Roboto::utf8ToUnicode(const std::string &utf8_string) {
  std::vector<uint32_t> unicode_points;
  for (size_t index = 0; index < utf8_string.length();) {
    const unsigned char byte_val = utf8_string[index];
    if (byte_val < 0x80) {
      unicode_points.push_back(byte_val);
      index++;
    } else if ((byte_val & 0xE0) == 0xC0) {
      if (index + 1 < utf8_string.length()) {
        unicode_points.push_back(((byte_val & 0x1F) << 6) | (utf8_string[index + 1] & 0x3F));
        index += 2;
      } else {
        index++;
      }
    } else if ((byte_val & 0xF0) == 0xE0) {
      if (index + 2 < utf8_string.length()) {
        const uint32_t code_point = ((byte_val & 0x0F) << 12) | ((utf8_string[index + 1] & 0x3F) << 6) | (utf8_string[index + 2] & 0x3F);
        unicode_points.push_back(code_point);
        index += 3;
      } else {
        index++;
      }
    } else {
      index++;
    }
  }
  return unicode_points;
}

} // namespace correlation::plotters
