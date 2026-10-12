/**
 * @file PlotMarkerManagerTests.cpp
 * @brief Unit tests for PlotMarkerManager.
 * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)
 * @par License
 * SPDX-License-Identifier: AGPL-3.0-only
 */

#include "app/viewmodel/PlotMarkerManager.hpp"
#include <gtest/gtest.h>

namespace correlation::app {
namespace {

TEST(PlotMarkerManagerTests, InitialStateIsEmpty) {
  const PlotMarkerManager manager;
  EXPECT_EQ(manager.markersCount(), 0);
  EXPECT_TRUE(manager.markers().empty());

  correlation::plotters::PlotConfig config;
  manager.applyToConfig(config);
  EXPECT_TRUE(config.reference_lines.empty());
}

TEST(PlotMarkerManagerTests, AddsVerticalMarkerCorrectly) {
  PlotMarkerManager manager;
  manager.addMarker(static_cast<real_t>(2.5), true);

  EXPECT_EQ(manager.markersCount(), 1);
  const auto markers = manager.markers();
  ASSERT_EQ(markers.size(), 1);

  EXPECT_NEAR(markers[0].value, static_cast<real_t>(2.5), 1e-5);
  EXPECT_TRUE(markers[0].is_vertical);
  EXPECT_EQ(markers[0].color_hex, "#0284C7");
  EXPECT_EQ(markers[0].label, "X: 2.500");
}

TEST(PlotMarkerManagerTests, AddsHorizontalMarkerCorrectly) {
  PlotMarkerManager manager;
  manager.addMarker(static_cast<real_t>(1.234), false);

  EXPECT_EQ(manager.markersCount(), 1);
  const auto markers = manager.markers();
  ASSERT_EQ(markers.size(), 1);

  EXPECT_NEAR(markers[0].value, static_cast<real_t>(1.234), 1e-5);
  EXPECT_FALSE(markers[0].is_vertical);
  EXPECT_EQ(markers[0].color_hex, "#E11D48");
  EXPECT_EQ(markers[0].label, "Y: 1.234");
}

TEST(PlotMarkerManagerTests, ClearsMarkersCorrectly) {
  PlotMarkerManager manager;
  manager.addMarker(static_cast<real_t>(1.0), true);
  manager.addMarker(static_cast<real_t>(2.0), false);
  EXPECT_EQ(manager.markersCount(), 2);

  manager.clearMarkers();
  EXPECT_EQ(manager.markersCount(), 0);
  EXPECT_TRUE(manager.markers().empty());
}

TEST(PlotMarkerManagerTests, AppliesMarkersToPlotConfig) {
  PlotMarkerManager manager;
  manager.addMarker(static_cast<real_t>(3.5), true);
  manager.addMarker(static_cast<real_t>(0.5), false);

  correlation::plotters::PlotConfig config;
  manager.applyToConfig(config);

  ASSERT_EQ(config.reference_lines.size(), 2);
  EXPECT_TRUE(config.reference_lines[0].is_vertical);
  EXPECT_NEAR(config.reference_lines[0].value, static_cast<real_t>(3.5), 1e-5);
  EXPECT_FALSE(config.reference_lines[1].is_vertical);
  EXPECT_NEAR(config.reference_lines[1].value, static_cast<real_t>(0.5), 1e-5);
}

} // namespace
} // namespace correlation::app
