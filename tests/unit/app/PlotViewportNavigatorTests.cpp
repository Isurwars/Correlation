/**
 * @file PlotViewportNavigatorTests.cpp
 * @brief Unit tests for PlotViewportNavigator.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/viewmodel/PlotViewportNavigator.hpp"
#include <gtest/gtest.h>

namespace correlation::app {
namespace {

correlation::analysis::Histogram createSampleHistogram() {
  correlation::analysis::Histogram hist;
  hist.bins = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  hist.partials["Total"] = {0.0, 0.5, 1.5, 3.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.0, 0.0};
  return hist;
}

TEST(PlotViewportNavigatorTests, InitialStateHasNoActiveZoom) {
  const PlotViewportNavigator nav;
  EXPECT_FALSE(nav.hasActiveZoom());
  const auto limits = nav.limits();
  EXPECT_FALSE(limits.x_min.has_value());
  EXPECT_FALSE(limits.x_max.has_value());
  EXPECT_FALSE(limits.y_min.has_value());
  EXPECT_FALSE(limits.y_max.has_value());
}

TEST(PlotViewportNavigatorTests, ComputesScalesCorrectly) {
  const auto hist = createSampleHistogram();
  const correlation::plotters::PlotConfig config;

  const auto x_scale = PlotViewportNavigator::computeXScale(hist, config);
  EXPECT_LE(x_scale.min, 0.0);
  EXPECT_GE(x_scale.max, 10.0);

  const auto y_scale = PlotViewportNavigator::computeYScale(hist, config);
  EXPECT_LE(y_scale.min, 0.0);
  EXPECT_GE(y_scale.max, 3.0);
}

TEST(PlotViewportNavigatorTests, ManualLimitsOverrideScales) {
  const auto hist = createSampleHistogram();
  correlation::plotters::PlotConfig config;
  config.manual_x_min = 2.0;
  config.manual_x_max = 6.0;
  config.manual_y_min = 0.5;
  config.manual_y_max = 2.5;

  const auto x_scale = PlotViewportNavigator::computeXScale(hist, config);
  EXPECT_DOUBLE_EQ(x_scale.min, 2.0);
  EXPECT_DOUBLE_EQ(x_scale.max, 6.0);

  const auto y_scale = PlotViewportNavigator::computeYScale(hist, config);
  EXPECT_DOUBLE_EQ(y_scale.min, 0.5);
  EXPECT_DOUBLE_EQ(y_scale.max, 2.5);
}

TEST(PlotViewportNavigatorTests, IgnoresSmallZoomRect) {
  PlotViewportNavigator nav;
  const auto hist = createSampleHistogram();
  const correlation::plotters::PlotConfig config;

  nav.applyZoomRect(100.0F, 100.0F, 102.0F, 102.0F, 800.0F, 600.0F, hist, config);
  EXPECT_FALSE(nav.hasActiveZoom());
}

TEST(PlotViewportNavigatorTests, AppliesValidZoomRect) {
  PlotViewportNavigator nav;
  const auto hist = createSampleHistogram();
  const correlation::plotters::PlotConfig config;

  nav.applyZoomRect(100.0F, 100.0F, 300.0F, 400.0F, 800.0F, 600.0F, hist, config);
  EXPECT_TRUE(nav.hasActiveZoom());

  const auto limits = nav.limits();
  ASSERT_TRUE(limits.x_min.has_value());
  ASSERT_TRUE(limits.x_max.has_value());
  ASSERT_TRUE(limits.y_min.has_value());
  ASSERT_TRUE(limits.y_max.has_value());
  EXPECT_LT(*limits.x_min, *limits.x_max);
  EXPECT_LT(*limits.y_min, *limits.y_max);
}

TEST(PlotViewportNavigatorTests, ResetZoomRestoresState) {
  PlotViewportNavigator nav;
  const auto hist = createSampleHistogram();
  const correlation::plotters::PlotConfig config;

  nav.applyZoomRect(100.0F, 100.0F, 300.0F, 400.0F, 800.0F, 600.0F, hist, config);
  EXPECT_TRUE(nav.hasActiveZoom());

  nav.resetZoom();
  EXPECT_FALSE(nav.hasActiveZoom());
  const auto limits = nav.limits();
  EXPECT_FALSE(limits.x_min.has_value());
  EXPECT_FALSE(limits.x_max.has_value());
}

TEST(PlotViewportNavigatorTests, AppliesPanWithDeltaShift) {
  PlotViewportNavigator nav;
  const auto hist = createSampleHistogram();
  const correlation::plotters::PlotConfig config;

  nav.applyPan(50.0F, -30.0F, 400.0F, 300.0F, 800.0F, 600.0F, hist, config);
  EXPECT_TRUE(nav.hasActiveZoom());

  const auto limits = nav.limits();
  EXPECT_TRUE(limits.x_min.has_value());
  EXPECT_TRUE(limits.x_max.has_value());
}

TEST(PlotViewportNavigatorTests, AppliesWheelZoomInAndOut) {
  PlotViewportNavigator nav;
  const auto hist = createSampleHistogram();
  const correlation::plotters::PlotConfig config;

  // Zoom in
  nav.applyWheelZoom(400.0F, 300.0F, 1.0F, 800.0F, 600.0F, hist, config);
  EXPECT_TRUE(nav.hasActiveZoom());
  const auto limits_in = nav.limits();
  ASSERT_TRUE(limits_in.x_max.has_value());
  ASSERT_TRUE(limits_in.x_min.has_value());
  const real_t range_in = limits_in.x_max.value() - limits_in.x_min.value();

  // Zoom out
  nav.applyWheelZoom(400.0F, 300.0F, -1.0F, 800.0F, 600.0F, hist, config);
  const auto limits_out = nav.limits();
  ASSERT_TRUE(limits_out.x_max.has_value());
  ASSERT_TRUE(limits_out.x_min.has_value());
  const real_t range_out = limits_out.x_max.value() - limits_out.x_min.value();

  EXPECT_GT(range_out, range_in);
}

TEST(PlotViewportNavigatorTests, AppliesLimitsToPlotConfig) {
  PlotViewportNavigator nav;
  nav.setLimits(ViewportLimits{
      .x_min = 1.5,
      .x_max = 8.5,
      .y_min = -0.2,
      .y_max = 4.2,
  });

  correlation::plotters::PlotConfig config;
  nav.applyToConfig(config);

  ASSERT_TRUE(config.manual_x_min.has_value());
  ASSERT_TRUE(config.manual_x_max.has_value());
  ASSERT_TRUE(config.manual_y_min.has_value());
  ASSERT_TRUE(config.manual_y_max.has_value());

  EXPECT_NEAR(*config.manual_x_min, static_cast<real_t>(1.5), 1e-5);
  EXPECT_NEAR(*config.manual_x_max, static_cast<real_t>(8.5), 1e-5);
  EXPECT_NEAR(*config.manual_y_min, static_cast<real_t>(-0.2), 1e-5);
  EXPECT_NEAR(*config.manual_y_max, static_cast<real_t>(4.2), 1e-5);
}

} // namespace
} // namespace correlation::app
