/**
// Correlation - Liquid and Amorphous Solid Analysis Tool
// Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
// SPDX-License-Identifier: AGPL-3.0-only
// Full license: https://github.com/Isurwars/Correlation/blob/main/LICENSE
*/

#include "render3d/MiniPlotRenderer.hpp"

#include <gtest/gtest.h>

using namespace correlation::render3d;
using namespace correlation::analysis;

class MiniPlotRendererTests : public ::testing::Test {
protected:
  MetricCurveData createMockData() {
    MetricCurveData data;
    data.bins = {1.0, 2.0, 3.0, 4.0, 5.0};
    data.mean_curve = {0.2, 1.5, 3.8, 1.2, 0.4};
    data.frame_curves.resize(3);
    data.frame_curves[0] = {0.1, 1.1, 3.5, 1.3, 0.3};
    data.frame_curves[1] = {0.3, 1.9, 4.1, 1.0, 0.5};
    data.frame_curves[2] = {0.2, 1.4, 3.2, 1.5, 0.4};
    data.name = "g(r)";
    data.title = "Radial Distribution Function";
    data.x_unit = "Å";
    data.y_unit = "";
    return data;
  }
};

TEST_F(MiniPlotRendererTests, IdentifiesPeakAccurately) {
  const std::vector<real_t> bins = {1.0, 2.0, 3.0, 4.0, 5.0};
  const std::vector<real_t> values = {0.1, 1.5, 4.2, 1.8, 0.3};

  const auto peak = MiniPlotRenderer::findPeak(bins, values);
  EXPECT_TRUE(peak.valid);
  EXPECT_NEAR(peak.x, 3.0, 1e-4);
  EXPECT_NEAR(peak.y, 4.2, 1e-4);
}

TEST_F(MiniPlotRendererTests, RendersValidBufferDimensions) {
  const auto data = createMockData();
  const auto buf = MiniPlotRenderer::render(data, 1, 260, 130);

  EXPECT_EQ(buf.width, 260U);
  EXPECT_EQ(buf.height, 130U);
  EXPECT_EQ(buf.pixels.size(), 260U * 130U);

  // Background should not be entirely black or transparent
  bool has_content = false;
  for (const auto &p : buf.pixels) {
    if (p.r != 0 || p.g != 0 || p.b != 0) {
      has_content = true;
      break;
    }
  }
  EXPECT_TRUE(has_content);
}

TEST_F(MiniPlotRendererTests, HandlesEmptyDataGracefully) {
  MetricCurveData empty_data;
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

