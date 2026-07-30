/**
 * @file PathFontTests.cpp
 * @brief Unit tests for Roboto vector font path renderer.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "plotters/PathFont.hpp"

#include <gtest/gtest.h>
#include <string>

namespace correlation::plotters::testing {

TEST(PathFontTests, RenderEmptyStringReturnsEmptyPath) {
  TextRenderParameters params;
  params.text = "";
  params.start_x = 0.0;
  params.start_y = 0.0;
  params.font_size = 12.0;
  params.anchor = "start";

  std::string svg_path = Roboto::instance().render(params);
  EXPECT_TRUE(svg_path.empty());
}

TEST(PathFontTests, RenderBasicAsciiStringProducesValidPathCommands) {
  TextRenderParameters params;
  params.text = "123.45 RDF g(r)";
  params.start_x = 10.0;
  params.start_y = 20.0;
  params.font_size = 14.0;
  params.anchor = "start";

  std::string svg_path = Roboto::instance().render(params);
  EXPECT_FALSE(svg_path.empty());
  // Path data should contain SVG move ('M') or line ('L') commands
  EXPECT_NE(svg_path.find('M'), std::string::npos);
}

TEST(PathFontTests, TextAnchorsModifyStartingOffset) {
  TextRenderParameters params_start;
  params_start.text = "Test";
  params_start.font_size = 12.0;
  params_start.anchor = "start";

  TextRenderParameters params_middle = params_start;
  params_middle.anchor = "middle";

  TextRenderParameters params_end = params_start;
  params_end.anchor = "end";

  std::string path_start = Roboto::instance().render(params_start);
  std::string path_middle = Roboto::instance().render(params_middle);
  std::string path_end = Roboto::instance().render(params_end);

  EXPECT_NE(path_start, path_middle);
  EXPECT_NE(path_start, path_end);
  EXPECT_NE(path_middle, path_end);
}

TEST(PathFontTests, UnknownGlyphFallbackDoesNotCrash) {
  TextRenderParameters params;
  params.text = "Special glyph: \xCE\xB1 \xCE\xB2"; // UTF-8 alpha, beta
  params.font_size = 10.0;

  EXPECT_NO_THROW({
    std::string path = Roboto::instance().render(params);
    (void)path;
  });
}

} // namespace correlation::plotters::testing
