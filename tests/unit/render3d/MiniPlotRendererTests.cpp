/**
 * @file MiniPlotRendererTests.cpp
 * @brief Unit tests for MiniPlotRenderer offscreen CPU curve rendering.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "render3d/MiniPlotRenderer.hpp"
#include <gtest/gtest.h>

namespace correlation::render3d {

using namespace correlation::analysis;
using correlation::real_t;

namespace {

class MiniPlotRendererTests : public ::testing::Test {
protected:
  [[nodiscard]] static MetricCurveData createMockData() {
    MetricCurveData data;
    data.bins = {1.0, 2.0, 3.0, 4.0, 5.0};
    data.mean_curve = {static_cast<real_t>(0.2), 1.5, static_cast<real_t>(3.8),
                       static_cast<real_t>(1.2), static_cast<real_t>(0.4)};
    data.frame_curves.resize(3);
    data.frame_curves[0] = {static_cast<real_t>(0.1), static_cast<real_t>(1.1), 3.5,
                            static_cast<real_t>(1.3), static_cast<real_t>(0.3)};
    data.frame_curves[1] = {static_cast<real_t>(0.3), static_cast<real_t>(1.9),
                            static_cast<real_t>(4.1), 1.0, 0.5};
    data.frame_curves[2] = {static_cast<real_t>(0.2), static_cast<real_t>(1.4),
                            static_cast<real_t>(3.2), 1.5, static_cast<real_t>(0.4)};
    data.name = "g(r)";
    data.title = "Radial Distribution Function";
    data.x_unit = "Å";
    data.y_unit = "";
    return data;
  }
};

} // namespace

TEST_F(MiniPlotRendererTests, IdentifiesPeakAccurately) {
  const std::vector<real_t> bins = {1.0, 2.0, 3.0, 4.0, 5.0};
  const std::vector<real_t> values = {static_cast<real_t>(0.1), 1.5, static_cast<real_t>(4.2),
                                      static_cast<real_t>(1.8), static_cast<real_t>(0.3)};

  const auto peak = MiniPlotRenderer::findPeak(bins, values);
  EXPECT_TRUE(peak.valid);
  EXPECT_NEAR(peak.x, static_cast<real_t>(3.0), static_cast<real_t>(1e-4));
  EXPECT_NEAR(peak.y, static_cast<real_t>(4.2), static_cast<real_t>(1e-4));
}

TEST_F(MiniPlotRendererTests, RendersValidBufferDimensions) {
  const auto data = createMockData();
  const auto buf = MiniPlotRenderer::render(data, 1, 260, 130);

  EXPECT_EQ(buf.width, 260U);
  EXPECT_EQ(buf.height, 130U);
  EXPECT_EQ(buf.pixels.size(), 260U * 130U);

  // Background should not be entirely black or transparent
  bool has_content = false;
  for (const auto &pixel : buf.pixels) {
    if (pixel.r != 0 || pixel.g != 0 || pixel.b != 0) {
      has_content = true;
      break;
    }
  }
  EXPECT_TRUE(has_content);
}

TEST_F(MiniPlotRendererTests, HandlesEmptyDataGracefully) {
  const MetricCurveData empty_data;
  const auto buf = MiniPlotRenderer::render(empty_data, std::nullopt, 200, 100);

  EXPECT_EQ(buf.width, 200U);
  EXPECT_EQ(buf.height, 100U);
  EXPECT_EQ(buf.pixels.size(), 200U * 100U);
}

TEST_F(MiniPlotRendererTests, RendersWithNullOptFrameIndex) {
  const auto data = createMockData();
  // Render with only mean curve
  const auto buf = MiniPlotRenderer::render(data, std::nullopt, 260, 130);
  EXPECT_EQ(buf.pixels.size(), 260U * 130U);
}

TEST_F(MiniPlotRendererTests, RendersWithManualBoundsAndLogScale) {
  const auto data = createMockData();
  const MiniPlotConfig config{
      .manual_x_min = 2.0,
      .manual_x_max = 4.0,
      .manual_y_min = 0.5,
      .manual_y_max = 5.0,
      .y_log_scale = true,
  };
  const auto buf = MiniPlotRenderer::render(data, 1, config, 260, 130);
  EXPECT_EQ(buf.width, 260U);
  EXPECT_EQ(buf.height, 130U);
  EXPECT_EQ(buf.pixels.size(), 260U * 130U);
}

} // namespace correlation::render3d
